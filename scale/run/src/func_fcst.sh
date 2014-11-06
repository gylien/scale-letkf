#!/bin/bash
#===============================================================================
#
#  Steps of 'fcst.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=4
stepname[1]='Prepare boundary files'
stepfunc[1]='boundary'
stepname[2]='Perturb boundaries'
stepfunc[2]='pertbdy'
stepname[3]='Run ensemble forecasts'
stepfunc[3]='ensfcst'
stepname[4]='Run verification'
stepfunc[4]='verf'

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run ensemble forecasts and (optional) verifications.

Configuration files:
  config.all
  config.$myname1 (optional)

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  MEMBERS     List of forecast members ('mean' for ensemble mean)
               all:     Run all members including ensemble mean (default)
               mems:    Run all members but not including ensemble mean
               '2 4 6': Run members 2, 4, 6
  CYCLE       Number of forecast cycles run in parallel
               (default: 1)
  CYCLE_SKIP  Run forecasts every ? cycles
               (default: 1)
  IF_VERF     Run verification? [Not finished!]
               0: No (default)
               1: Yes
              * to run the verification, a shared disk storing observations
                and reference model analyses needs to be used
  IF_EFSO     Use EFSO forecast length and output interval? [Not finished!]
               0: No (default)
               1: Yes
  ISTEP       The initial step in the first cycle from which this script starts
               (default: the first step)
  FSTEP       The final step in the last cycle by which this script ends
               (default: the last step)
  TIME_LIMIT  Requested time limit (only used when using a job scheduler)
               (default: 30 minutes)
"

if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
  echo "$USAGE"
  exit 0
fi

#-------------------------------------------------------------------------------
# set parameters from command line

STIME=${1:-$STIME}; shift
ETIME=${1:-$ETIME}; shift
MEMBERS=${1:-$MEMBERS}; shift
CYCLE=${1:-$CYCLE}; shift
CYCLE_SKIP=${1:-$CYCLE_SKIP}; shift
IF_VERF=${1:-$IF_VERF}; shift
IF_EFSO=${1:-$IF_EFSO}; shift
ISTEP=${1:-$ISTEP}; shift
FSTEP=${1:-$FSTEP}; shift
TIME_LIMIT="${1:-$TIME_LIMIT}"

#-------------------------------------------------------------------------------
# if some necessary parameters are not given, print the usage help and exit

if [ -z "$STIME" ]; then
  echo "$USAGE" >&2
  exit 1
fi

#-------------------------------------------------------------------------------
# error detection

if ((MACHINE_TYPE == 10 && ONLINE_STGOUT != 0)); then
  echo "[Error] $myname: When \$MACHINE_TYPE = 10, \$ONLINE_STGOUT needs to be 0." >&2
  exit 1
fi

#... more detections...

#-------------------------------------------------------------------------------
# assign default values to and standardize the parameters

STIME=$(datetime $STIME)
ETIME=$(datetime ${ETIME:-$STIME})
if [ -z "$MEMBERS" ] || [ "$MEMBERS" = 'all' ]; then
  MEMBERS="mean $(printf "$MEMBER_FMT " $(seq $MEMBER))"
elif [ "$MEMBERS" = 'mems' ]; then
  MEMBERS=$(printf "$MEMBER_FMT " $(seq $MEMBER))
else
  tmpstr=''
  for m in $MEMBERS; do
    if [ "$m" = 'mean' ] || [ "$m" = 'sprd' ]; then
      tmpstr="$tmpstr$m "
    else
      tmpstr="$tmpstr$(printf $MEMBER_FMT $((10#$m))) "
      (($? != 0)) && exit 1
    fi
  done
  MEMBERS="$tmpstr"
fi
CYCLE=${CYCLE:-1}
CYCLE_SKIP=${CYCLE_SKIP:-1}
IF_VERF=${IF_VERF:-0}
IF_EFSO=${IF_EFSO:-0}
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-$nsteps}
TIME_LIMIT=${TIME_LIMIT:-"0:30:00"}

#-------------------------------------------------------------------------------
}

#===============================================================================

staging_list () {
#-------------------------------------------------------------------------------

if ((TMPDAT_MODE == 1 && MACHINE_TYPE != 10)); then
#-------------------
  safe_init_tmpdir $TMPDAT
  safe_init_tmpdir $TMPDAT/exec
  ln -fs $MODELDIR/scale-les $TMPDAT/exec
  ln -fs $MODELDIR/scale-les_init $TMPDAT/exec
  ln -fs $MODELDIR/scale-les_pp $TMPDAT/exec
  ln -fs $COMMON_DIR/pdbash $TMPDAT/exec
  ln -fs $DATADIR/rad $TMPDAT
  ln -fs $ANLWRF $TMPDAT/wrf

  safe_init_tmpdir $TMPDAT/conf
  ln -fs $SCRP_DIR/*.conf $TMPDAT/conf

  if ((PREP_TOPO == 1)); then
    ln -fs $DATADIR/topo_prep $TMPDAT/topo_prep
  else
    ln -fs $DATADIR/topo $TMPDAT
  fi
  if ((PREP_LANDUSE == 1)); then
    ln -fs $DATADIR/landuse_prep $TMPDAT/landuse_prep
  else
    ln -fs $DATADIR/landuse $TMPDAT
  fi
#-------------------
else
#-------------------

  cat >> $STAGING_DIR/stagein.dat << EOF
${MODELDIR}/scale-les|exec/scale-les
${MODELDIR}/scale-les_init|exec/scale-les_init
${MODELDIR}/scale-les_pp|exec/scale-les_pp
${COMMON_DIR}/pdbash|exec/pdbash
${SCRP_DIR}/scale.conf|conf/scale.conf
${SCRP_DIR}/scale_init.conf|conf/scale_init.conf
${DATADIR}/rad|rad
EOF

  time=$STIME
  etime_anlwrf=$(datetime $ETIME $((FCSTLEN+ANLWRF_INT)) s)
  while ((time <= etime_anlwrf)); do
    path="wrfout_d01_${time}"
    echo "${ANLWRF}/${path}|wrf/${path}" >> $STAGING_DIR/stagein.dat
    if ((time == etime_anlwrf)); then
      break
    fi
    time=$(datetime $time $ANLWRF_INT s)
  done

  if ((PREP_TOPO == 1)); then
    for q in $(seq $mem_np); do
      path="topo_prep/topo$(printf $SCALE_SFX $((q-1)))"
      echo "${DATADIR}/${path}|${path}" >> $STAGING_DIR/stagein.dat
    done
  else
    cat >> $STAGING_DIR/stagein.dat << EOF
${SCRP_DIR}/scale_pp_topo.conf|conf/scale_pp_topo.conf
${DATADIR}/topo/DEM50M/Products|topo/DEM50M/Products
EOF
  fi

  if ((PREP_LANDUSE == 1)); then
    for q in $(seq $mem_np); do
      path="landuse_prep/landuse$(printf $SCALE_SFX $((q-1)))"
      echo "${DATADIR}/${path}|${path}" >> $STAGING_DIR/stagein.dat
    done
  else
    cat >> $STAGING_DIR/stagein.dat << EOF
${SCRP_DIR}/scale_pp_landuse.conf|conf/scale_pp_landuse.conf
${DATADIR}/landuse/LU100M/Products|landuse/LU100M/Products
EOF
  fi

  if ((MACHINE_TYPE == 10)); then
    cat >> $STAGING_DIR/stagein.dat << EOF
${COMMON_DIR}/datetime|exec/datetime
EOF
  fi
#-------------------
fi

#-------------------------------------------------------------------------------

if ((TMPOUT_MODE == 1 && MACHINE_TYPE != 10)); then
#-------------------
  mkdir -p $(dirname $TMPOUT)
  ln -fs $OUTDIR $TMPOUT
#-------------------
else
#-------------------
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  loop=0
  while ((time <= ETIME)); do
    loop=$((loop+1))
    if ((ONLINE_STGOUT == 1)); then
      stgoutstep="stageout.loop.${loop}"
    else
      stgoutstep='stageout.out'
    fi

    for c in $(seq $CYCLE); do
      time2=$(datetime $time $((lcycles * (c-1))) s)
      if ((time2 <= ETIME)); then
        for m in $(seq $fmember); do
          mm=$(((c-1) * fmember + m))
          for q in $(seq $mem_np); do
            #-------------------
            # stage-in

            path="${time2}/anal/${name_m[$mm]}/init$(printf $SCALE_SFX $((q-1)))"
            echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/stagein.out.${mem2proc[$(((mm-1)*mem_np+q))]}

            #-------------------
            # stage-out

            if ((OUT_OPT <= 2)); then
              path="${time2}/fcst/${name_m[$mm]}/history$(printf $SCALE_SFX $((q-1)))"
              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$(((mm-1)*mem_np+q))]}
            fi
            if ((OUT_OPT <= 1)); then
              path="${time2}/fcst/${name_m[$mm]}/init_$(datetime ${time2} $FCSTLEN s)$(printf $SCALE_SFX $((q-1)))"
              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$(((mm-1)*mem_np+q))]}
            fi

            #-------------------
          done

          #-------------------

          if ((LOG_OPT <= 3)); then
            path="${time2}/log/scale/${name_m[$mm]}_LOG${SCALE_LOG_SFX}"
            echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$(((mm-1)*mem_np+1))]}
          fi
#          if ((LOG_OPT <= 1)); then
#            # perturb bdy log
#          fi

          #-------------------
        done

        #-------------------

        if ((LOG_OPT <= 2)); then
          if ((repeat_mems <= fmember)); then
            tmpidx=1                              # mm=1
          else
            tmpidx=$((((c-1)*fmember)*mem_np+1))  # mm=$(((c-1) * fmember + 1))
          fi
          path="${time2}/log/scale_topo/pp_LOG${SCALE_LOG_SFX}"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$tmpidx]}
          path="${time2}/log/scale_landuse/pp_LOG${SCALE_LOG_SFX}"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$tmpidx]}
          path="${time2}/log/scale_bdy/init_LOG${SCALE_LOG_SFX}"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2proc[$tmpidx]}
        fi

        #-------------------
      fi
    done

    time=$(datetime $time $((lcycles * CYCLE)) s)
  done
#-------------------
fi



#for c in `seq $CYCLES`; do
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig" >> $tmpstageout/out.${node_m[$mt]}
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#    fh=0
#    while [ "$fh" -le "$FCSTLEN" ]; do
#      fhhh=`printf '%03d' $fh`
#      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
#      if [ "$OUT_OPT" -le 1 ]; then
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$mt]}
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      if [ "$OUT_OPT" -le 2 ]; then
#        echo "mv|fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|verfo1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa2/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#    fh=$((fh + FCSTOUT))
#    done
#  done
#done

#-------------------------------------------------------------------------------
}

#===============================================================================

boundary_sub () {
#-------------------------------------------------------------------------------
# Run a series of scripts (topo/landuse/init) to make the boundary files.
#
# Usage: make_boundary NODEFILE PDBASH_PROC_OPT
#
#   NODEFILE          $NODEFILE in functions 'pdbash' and 'mpirunf'
#   PDBASH_PROC_OPT   $PROC_OPT in function 'pdbash'
#
# Other input variables:
#   $c       Cycle number
#   $cf      Formatted cycle number
#   $stimes  Start time of this cycle
#   $SCRP_DIR
#   $TMPRUN
#   $TMPDAT
#   $mem_np
#   $PREP_TOPO
#   $PREP_LANDUSE
#-------------------------------------------------------------------------------

if (($# < 2)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PDBASH_PROC_OPT="$1"; shift

#-------------------------------------------------------------------------------
# topo

if ((PREP_TOPO != 1)); then
  pdbash $NODEFILE $PDBASH_PROC_OPT \
    $SCRP_DIR/src/pre_scale_pp_topo.sh ${stimes[$c]} $TMPRUN/scale_pp_topo/${cf} $TMPDAT/exec $TMPDAT
  mpirunf $NODEFILE \
    $TMPRUN/scale_pp_topo/${cf} ./scale-les_pp pp.conf
  if ((LOG_OPT <= 2)); then
    pdbash $NODEFILE $PDBASH_PROC_OPT \
      $SCRP_DIR/src/post_scale_pp_topo.sh ${stimes[$c]} $TMPRUN/scale_pp_topo/${cf}
  fi
fi

#-------------------------------------------------------------------------------
# landuse

if ((PREP_LANDUSE != 1)); then
  pdbash $NODEFILE $PDBASH_PROC_OPT \
    $SCRP_DIR/src/pre_scale_pp_landuse.sh ${stimes[$c]} $TMPRUN/scale_pp_landuse/${cf} $TMPDAT/exec $TMPDAT
  mpirunf $NODEFILE \
    $TMPRUN/scale_pp_landuse/${cf} ./scale-les_pp pp.conf
  if ((LOG_OPT <= 2)); then
    pdbash $NODEFILE $PDBASH_PROC_OPT \
      $SCRP_DIR/src/post_scale_pp_landuse.sh ${stimes[$c]} $TMPRUN/scale_pp_landuse/${cf}
  fi
fi

#-------------------------------------------------------------------------------
# init

if ((PREP_TOPO == 1)); then
  local topo_base="$TMPDAT/topo_prep/topo"
else
  local topo_base="$TMPRUN/scale_pp_topo/${cf}/topo"
fi
if ((PREP_LANDUSE == 1)); then
  local landuse_base="$TMPDAT/landuse_prep/landuse"
else
  local landuse_base="$TMPRUN/scale_pp_landuse/${cf}/landuse"
fi
pdbash $NODEFILE $PDBASH_PROC_OPT \
  $SCRP_DIR/src/pre_scale_init.sh $mem_np $topo_base $landuse_base $TMPDAT/wrf/wrfout_d01 \
  ${stimes[$c]} $FCSTLEN $TMPRUN/scale_init/${cf} $TMPDAT/exec $TMPDAT
mpirunf $NODEFILE \
  $TMPRUN/scale_init/${cf} ./scale-les_init init.conf
if ((LOG_OPT <= 2)); then
  pdbash $NODEFILE $PDBASH_PROC_OPT \
    $SCRP_DIR/src/post_scale_init.sh ${stimes[$c]} $FCSTLEN $TMPRUN/scale_init/${cf}
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

boundary () {
#-------------------------------------------------------------------------------

echo
if ((PREP_BDY == 1)); then
  echo "  ... skip this step (use prepared boundary files)"
  return 0
fi

if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
#-------------------
  ipm=0
  for c in $(seq $rcycle); do
    cf=$(printf $CYCLE_FMT $c)
    ipm=$((ipm+1))
    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    cfr=$(printf $CYCLE_FMT $(((ipm-1)/fmember+1))) # try to use processes in parallel
    echo "  ${stimesfmt[$c]}: node ${node_m[$ipm]} [$(datetime_now)]"

    boundary_sub proc.${cfr}.${name_m[$ipm]} one &
    sleep $BGJOB_INT
  done
  wait
#-------------------
else # local run directory: run multiple members as needed
#-------------------
  if ((repeat_mems <= fmember)); then
    ipm=0
    for c in $(seq $rcycle); do
      cf=$(printf $CYCLE_FMT $c)
      for m in $(seq $repeat_mems); do
        ipm=$((ipm+1))
        if ((ipm > parallel_mems)); then wait; ipm=1; fi
        echo "  ${stimesfmt[$c]}: node ${node_m[$m]} [$(datetime_now)]"

        boundary_sub proc.$(printf $CYCLE_FMT 1).${name_m[$m]} alln &
        sleep $BGJOB_INT
      done
    done
    wait
  else
    ipm=0
    for c in $(seq $rcycle); do
      cf=$(printf $CYCLE_FMT $c)
      for m in $(seq $fmember); do
        mm=$(((c-1) * fmember + m))
        ipm=$((ipm+1))
        if ((ipm > parallel_mems)); then wait; ipm=1; fi
        echo "  ${stimesfmt[$c]}: node ${node_m[$mm]} [$(datetime_now)]"

        boundary_sub proc.${cf}.${name_m[$mm]} alln &
        sleep $BGJOB_INT
      done
    done
    wait
  fi
#-------------------
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

pertbdy () {
#-------------------------------------------------------------------------------

echo
if ((PERTURB_BDY == 0)); then
  echo "  ... skip this step (do not perturb boundaries)"
  return 1
fi

###### not finished yet...
echo "pertbdy..."
######

ipm=0
for c in $(seq $rcycle); do
  cf=$(printf $CYCLE_FMT $c)
  if ((PREP_BDY == 1)); then
    bdy_base="$TMPDAT/bdy_prep/bdy_${stimes[$c]}"
  else
    bdy_base="$TMPRUN/scale_init/${cf}/boundary"
  fi

  for m in $(seq $fmember); do
    mm=$(((c-1) * fmember + m))
    ipm=$((ipm+1))
    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    echo "  ${stimesfmt[$c]}, member ${name_m[$mm]}: node ${node_m[$mm]} [$(datetime_now)]"

    if ((TMPRUN_MODE <= 2)); then
      proc_opt='one'
    else
      proc_opt='alln'
    fi
#   ......
#    ( pdbash proc.${cf}.${name_m[$mm]} $proc_opt $SCRP_DIR/src/pre_scale.sh $mem_np \
#        $TMPOUT/${stimes[$c]}/anal/${name_m[$mm]}/init $bdy_base $topo_base $landuse_base \
#        ${stimes[$c]} $FCSTLEN $FCSTOUT $TMPRUN/scale/${cf}_${name_m[$mm]} $TMPDAT/exec $TMPDAT ;
#      mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \
#        ./scale-les run.conf ) &
#   ......
#   $TMPRUN/pertbdy/${cf}_${name_m[$mm]}

    sleep $BGJOB_INT
  done
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst () {
#-------------------------------------------------------------------------------

echo

ipm=0
for c in $(seq $rcycle); do
  cf=$(printf $CYCLE_FMT $c)
  if ((PREP_TOPO == 1)); then
    topo_base="$TMPDAT/topo_prep/topo"
  else
    topo_base="$TMPRUN/scale_pp_topo/${cf}/topo"
  fi
  if ((PREP_LANDUSE == 1)); then
    landuse_base="$TMPDAT/landuse_prep/landuse"
  else
    landuse_base="$TMPRUN/scale_pp_landuse/${cf}/landuse"
  fi

  for m in $(seq $fmember); do
    mm=$(((c-1) * fmember + m))
    ipm=$((ipm+1))
    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    echo "  ${stimesfmt[$c]}, member ${name_m[$mm]}: node ${node_m[$mm]} [$(datetime_now)]"

    if ((PERTURB_BDY == 1)); then
######
      echo "not finished yet..."
#      bdy_base="$TMPRUN/pertbdy/${cf}_${name_m[$mm]}/boundary"
######
    elif ((PREP_BDY == 1)); then
      bdy_base="$TMPDAT/bdy_prep/bdy_${stimes[$c]}"
    else
      bdy_base="$TMPRUN/scale_init/${cf}/boundary"
    fi
    if ((TMPRUN_MODE <= 2)); then
      proc_opt='one'
    else
      proc_opt='alln'
    fi
    ( pdbash proc.${cf}.${name_m[$mm]} $proc_opt $SCRP_DIR/src/pre_scale.sh $mem_np \
        $TMPOUT/${stimes[$c]}/anal/${name_m[$mm]}/init $bdy_base $topo_base $landuse_base \
        ${stimes[$c]} $FCSTLEN $FCSTOUT $TMPRUN/scale/${cf}_${name_m[$mm]} $TMPDAT/exec $TMPDAT ;
      mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \
        ./scale-les run.conf ;
      pdbash proc.${cf}.${name_m[$mm]} $proc_opt $SCRP_DIR/src/post_scale.sh $mem_np \
        ${stimes[$c]} ${name_m[$mm]} $FCSTLEN $TMPRUN/scale/${cf}_${name_m[$mm]} ) &

    sleep $BGJOB_INT
  done
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

verf () {
#-------------------------------------------------------------------------------

echo
if ((IF_VERF == 0)); then
  echo "  ... skip this step (do not run verification)"
  return 1
fi

echo "verf..."


#cd $TMPMPI
#cp -f $RUNDIR/datetime.sh .

#cat > fcst_31.sh << EOF
#. datetime.sh
#mem="\$1"
#cyc="\$2"
#stime="\$3"
#Syyyymmddhh=\${stime:0:10}
#mkdir -p $ltmpout/fcst/\${Syyyymmddhh}/\${mem}
#mkdir -p $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}
#mkdir -p $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}

#mkdir -p $ltmpssio/\${cyc}_\${mem}
#cd $ltmpssio/\${cyc}_\${mem}
#ln -fs $ltmpprog/ss2grd .
#ln -fs $ltmpprog/ss2grdp .
#fh=0
#while [ "\$fh" -le "$FCSTLEN" ]; do
#  fhh="\$(printf '%02d' \$fh)"
#  fhhh="\$(printf '%03d' \$fh)"
#  Fyyyymmddhh=\$(datetime \$stime \$fh h | cut -c 1-10)
#  cd $ltmpssio/\${cyc}_\${mem}
#  rm -f fort.*
#  if [ "\$fh" -eq 0 ]; then
#    ln -fs $ltmpgfs/\${cyc}_\${mem}/sig_ini fort.11
#    ln -fs $ltmpgfs/\${cyc}_\${mem}/sfc_ini fort.12
#  else
#    ln -fs $ltmpgfs/\${cyc}_\${mem}/SIG.F\${fhh} fort.11
#    ln -fs $ltmpgfs/\${cyc}_\${mem}/SFC.F\${fhh} fort.12
#  fi
#  ./ss2grd
#  mv -f fort.31 $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd
#EOF
#if [ "$OUT_OPT" -le 2 ]; then
#  cat >> fcst_31.sh << EOF
#  ./ss2grdp
#  mv -f fort.31 $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd
#EOF
#fi
#if [ "$OUT_OPT" -le 1 ]; then
#  cat >> fcst_31.sh << EOF
#  if [ "\$fh" -eq 0 ]; then
#    cp -fL $ltmpgfs/\${cyc}_\${mem}/sig_ini $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sig
#    cp -fL $ltmpgfs/\${cyc}_\${mem}/sfc_ini $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sfc
#  else
#    mv -f $ltmpgfs/\${cyc}_\${mem}/SIG.F\${fhh} $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sig
#    mv -f $ltmpgfs/\${cyc}_\${mem}/SFC.F\${fhh} $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sfc
#  fi
#EOF
#fi
#cat >> fcst_31.sh << EOF
#fh=\$((fh+$FCSTOUT))
#done
#EOF

##-------------------------------------------------------------------------------

#echo
#ppnl=$((ppn*2))
#np=1
#pcount=0
#for c in `seq $CYCLES`; do
#  cf=`printf $CYCLE_FMT $c`
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    if [ "${node_m[$mt]}" = "${node[1]}" ]; then
#      pcount=$((pcount+np))
#      if [ "$pcount" -gt "$ppnl" ]; then
#        echo "    wait..."
#        wait
#        pcount=$np
#      fi
#    fi
#    echo "  ${stimef[$c]}, member ${name_m[$mt]} on node '${node_m[$mt]}'"
#    $MPIBIN/mpiexec -host ${node_m[$mt]} bash fcst_31.sh ${name_m[$mt]} $cf "${STIME[$c]}" &
#    sleep $BGJOB_INT
#  done
#done
#echo "    wait..."
#wait

#======================###

#cd $TMPMPI
#cp -f $RUNDIR/datetime.sh .

#cat > fcst_41.sh << EOF
#. datetime.sh
#mem="\$1"
#cyc="\$2"
#stime="\$3"
#Syyyymmddhh=\${stime:0:10}

#mkdir -p $ltmpverify/\${cyc}_\${mem}
#cd $ltmpverify/\${cyc}_\${mem}
#ln -fs $ltmpprog/verify .
#fh=0
#while [ "\$fh" -le "$FCSTLEN" ]; do
#  fhhh=\$(printf '%03d' \$fh)
#  Fyyyymmddhh=\$(datetime \$stime \$fh h | cut -c 1-10)

#  if [ -s "$ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd" ] &&
#     [ -s "$ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd" ]; then
#    cd $ltmpverify/\${cyc}_\${mem}
#    rm -f fcst.grd fcstp.grd obs??.dat ana??.grd
#    ln -s $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd fcst.grd
#    ln -s $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd fcstp.grd
####### only support shared disk
#    cat $OBS/obs\${Fyyyymmddhh}/t.dat > obs01.dat
#    cat $OBS/obs\${Fyyyymmddhh}/t-1.dat >> obs01.dat
#    cat $OBS/obs\${Fyyyymmddhh}/t+1.dat >> obs01.dat
#    ln -s $ANLGRDP/\${Fyyyymmddhh}.grd ana01.grd
##    ln -s $ANLGRDP2/\${Fyyyymmddhh}.grd ana02.grd
#######
#    ./verify > /dev/null 2>&1

#    mkdir -p $ltmpout/verfo1/\${fhhh}/\${mem}
#    mkdir -p $ltmpout/verfa1/\${fhhh}/\${mem}
#    mkdir -p $ltmpout/verfa2/\${fhhh}/\${mem}
#    mv -f vrfobs01.dat $ltmpout/verfo1/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
#    mv -f vrfana01.dat $ltmpout/verfa1/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
##    mv -f vrfana02.dat $ltmpout/verfa2/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
#  fi
#fh=\$((fh+$FCSTOUT))
#done
#EOF

##-------------------------------------------------------------------------------

#echo
#ppnl=$ppn
##ppnl=$((ppn*2))
#np=1
#pcount=0
#for c in `seq $CYCLES`; do
#  cf=`printf $CYCLE_FMT $c`
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    if [ "${node_m[$mt]}" = "${node[1]}" ]; then
#      pcount=$((pcount+np))
#      if [ "$pcount" -gt "$ppnl" ]; then
#        echo "    wait..."
#        wait
#        pcount=$np
#      fi
#    fi
#    echo "  ${stimef[$c]}, member ${name_m[$mt]} on node '${node_m[$mt]}'"
#    $MPIBIN/mpiexec -host ${node_m[$mt]} bash fcst_41.sh ${name_m[$mt]} $cf "${STIME[$c]}" &
#    sleep $BGJOB_INT
#  done
#done
#echo "    wait..."
#wait

#-------------------------------------------------------------------------------
}

#===============================================================================
