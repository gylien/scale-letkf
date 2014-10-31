#!/bin/bash
#===============================================================================
#
#  Steps of 'fcst.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

function set_topo_landuse_bases {
#-------------------------------------------------------------------------------

if ((FIX_TOPO == 1)); then
  topo_base="$TMPDAT/topo_fix/topo"
else
  topo_base="$TMPRUN/scale_pp_topo/${cf}/topo"
fi
if ((FIX_LANDUSE == 1)); then
  landuse_base="$TMPDAT/landuse_fix/landuse"
else
  landuse_base="$TMPRUN/scale_pp_landuse/${cf}/landuse"
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function setting {
#-------------------------------------------------------------------------------
# define steps

nsteps=6
stepname[1]='Prepare topo files'
stepfunc[1]='topo'
stepname[2]='Prepare landuse files'
stepfunc[2]='landuse'
stepname[3]='Prepare boundary files'
stepfunc[3]='boundary'
stepname[4]='Perturb boundaries'
stepfunc[4]='pertbdy'
stepname[5]='Run ensemble forecasts'
stepfunc[5]='ensfcst'
stepname[6]='Run verification'
stepfunc[6]='verf'

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run ensemble forecasts and (optional) verifications.

Configuration files:
  config.all
  config.$myname1 (optional)

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP]

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
FSTEP=${1:-$FSTEP}

#-------------------------------------------------------------------------------
# if some necessary parameters are not given, print the usage help and exit

if [ -z "$STIME" ]; then
  echo "$USAGE" >&2
  exit 1
fi

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

#-------------------------------------------------------------------------------
}

#===============================================================================

function staging_list {
#-------------------------------------------------------------------------------

if ((INPUT_MODE == 1 && MACHINE_TYPE != 10)); then
#---------------------------------------
  mkdir -p $TMPDAT
  ln -fs $MODELDIR $TMPDAT/exec
  ln -fs $DATADIR $TMPDAT/data
  ln -fs $ANLWRF $TMPDAT/wrf
#
######...... more
#
#---------------------------------------
else
#---------------------------------------

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
    echo "${ANLWRF}/wrfout_d01_${time}|wrf/wrfout_d01_${time}" >> $STAGING_DIR/stagein.dat
    if ((time == etime_anlwrf)); then
      break
    fi
    time=$(datetime $time $ANLWRF_INT s)
  done

  if ((FIX_TOPO == 1)); then
    for q in $(seq $mem_np); do
      echo "${DATADIR}/topo_fix/topo$(printf $SCALE_SFX $((q-1)))|topo_fix/topo$(printf $SCALE_SFX $((q-1)))" >> $STAGING_DIR/stagein.dat
    done
  else
    cat >> $STAGING_DIR/stagein.dat << EOF
${SCRP_DIR}/scale_pp_topo.conf|conf/scale_pp_topo.conf
${DATADIR}/topo/DEM50M/Products|topo/DEM50M/Products
EOF
  fi

  if ((FIX_LANDUSE == 1)); then
    for q in $(seq $mem_np); do
      echo "${DATADIR}/landuse_fix/landuse$(printf $SCALE_SFX $((q-1)))|landuse_fix/landuse$(printf $SCALE_SFX $((q-1)))" >> $STAGING_DIR/stagein.dat
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
#---------------------------------------
fi

if ((OUTPUT_MODE == 1 && MACHINE_TYPE != 10)); then
#---------------------------------------
  mkdir -p $(cd $TMPOUT/.. && pwd)
  ln -fs $OUTDIR $TMPOUT
#---------------------------------------
else
#---------------------------------------
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do

    for c in $(seq $CYCLE); do
      time2=$(datetime $time $((lcycles * (c-1))) s)
      if ((time2 <= ETIME)); then
        for m in $(seq $fmember); do
          mm=$(((c-1) * fmember + m))
          for q in $(seq $mem_np); do
            opath="anal/${name_m[$mm]}/${time2}$(printf $SCALE_SFX $((q-1)))"
            echo "${OUTDIR}/${opath}|${opath}" >> $STAGING_DIR/stagein.out.${mem2proc[$(((mm-1)*mem_np+q))]}
          done
        done
      fi
    done

    time=$(datetime $time $((lcycles * CYCLE)) s)
  done
#---------------------------------------
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
###### combine topo, landuse, boundary!!!!!

function topo {
#-------------------------------------------------------------------------------

echo

if ((FIX_TOPO == 1)); then

  echo "  ... skip this step (use fixed topo files)"

else

  #-------------------
  if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  #-------------------
    ipm=0
    for c in $(seq $rcycle); do
      cf=$(printf $CYCLE_FMT $c)
      ipm=$((ipm+1))
      if ((ipm > parallel_mems)); then wait; ipm=1; fi
      cfr=$(printf $CYCLE_FMT $(((ipm-1)/fmember+1))) # try to use processes in parallel
      echo "  ${stimesfmt[$c]}: node ${node_m[$ipm]} [$(datetime_now)]"
      ( pdbash proc.${cfr}.${name_m[$ipm]} one $SCRP_DIR/src/pre_scale_pp_topo.sh \
          ${stimes[$c]} $TMPRUN/scale_pp_topo/${cf} \
          $TMPDAT/exec $TMPDAT &&
        mpirunf proc.${cfr}.${name_m[$ipm]} $TMPRUN/scale_pp_topo/${cf} \
          ./scale-les_pp pp.conf ) &
      sleep $BGJOB_INT
    done
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
          ( pdbash proc.$(printf $CYCLE_FMT 1).${name_m[$m]} alln $SCRP_DIR/src/pre_scale_pp_topo.sh \
              ${stimes[$c]} $TMPRUN/scale_pp_topo/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.$(printf $CYCLE_FMT 1).${name_m[$m]} $TMPRUN/scale_pp_topo/${cf} \
              ./scale-les_pp pp.conf ) &
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
          ( pdbash proc.${cf}.${name_m[$mm]} alln $SCRP_DIR/src/pre_scale_pp_topo.sh \
              ${stimes[$c]} $TMPRUN/scale_pp_topo/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale_pp_topo/${cf} \
              ./scale-les_pp pp.conf ) &
          sleep $BGJOB_INT
        done
      done
      wait
    fi
  #-------------------
  fi
  #-------------------

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function landuse {
#-------------------------------------------------------------------------------

echo

if ((FIX_LANDUSE == 1)); then

  echo "  ... skip this step (use fixed landuse files)"

else

  #-------------------
  if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  #-------------------
    ipm=0
    for c in $(seq $rcycle); do
      cf=$(printf $CYCLE_FMT $c)
      ipm=$((ipm+1))
      if ((ipm > parallel_mems)); then wait; ipm=1; fi
      cfr=$(printf $CYCLE_FMT $(((ipm-1)/fmember+1))) # try to use processes in parallel
      echo "  ${stimesfmt[$c]}: node ${node_m[$ipm]} [$(datetime_now)]"
      ( pdbash proc.${cfr}.${name_m[$ipm]} one $SCRP_DIR/src/pre_scale_pp_landuse.sh \
          ${stimes[$c]} $TMPRUN/scale_pp_landuse/${cf} \
          $TMPDAT/exec $TMPDAT &&
        mpirunf proc.${cfr}.${name_m[$ipm]} $TMPRUN/scale_pp_landuse/${cf} \
          ./scale-les_pp pp.conf ) &
      sleep $BGJOB_INT
    done
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
          ( pdbash proc.$(printf $CYCLE_FMT 1).${name_m[$m]} alln $SCRP_DIR/src/pre_scale_pp_landuse.sh \
              ${stimes[$c]} $TMPRUN/scale_pp_landuse/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.$(printf $CYCLE_FMT 1).${name_m[$m]} $TMPRUN/scale_pp_landuse/${cf} \
              ./scale-les_pp pp.conf ) &
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
          ( pdbash proc.${cf}.${name_m[$mm]} alln $SCRP_DIR/src/pre_scale_pp_landuse.sh \
              ${stimes[$c]} $TMPRUN/scale_pp_landuse/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale_pp_landuse/${cf} \
              ./scale-les_pp pp.conf ) &
          sleep $BGJOB_INT
        done
      done
      wait
    fi
  #-------------------
  fi
  #-------------------

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function boundary {
#-------------------------------------------------------------------------------

echo

if ((FIX_LANDUSE == 1)); then

  echo "  ... skip this step (use fixed landuse files)"

else

  #-------------------
  if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  #-------------------
    ipm=0
    for c in $(seq $rcycle); do
      cf=$(printf $CYCLE_FMT $c)
      set_topo_landuse_bases
      ipm=$((ipm+1))
      if ((ipm > parallel_mems)); then wait; ipm=1; fi
      cfr=$(printf $CYCLE_FMT $(((ipm-1)/fmember+1))) # try to use processes in parallel
      echo "  ${stimesfmt[$c]}: node ${node_m[$ipm]} [$(datetime_now)]"
      ( pdbash proc.${cf}.${name_m[$ipm]} one $SCRP_DIR/src/pre_scale_init.sh $mem_np \
          $topo_base $landuse_base $TMPDAT/wrf/wrfout_d01 \
          ${stimes[$c]} $FCSTLEN $TMPRUN/scale_init/${cf} \
          $TMPDAT/exec $TMPDAT &&
        mpirunf proc.${cf}.${name_m[$ipm]} $TMPRUN/scale_init/${cf} \
          ./scale-les_init init.conf ) &
      sleep $BGJOB_INT
    done
  #-------------------
  else # local run directory: run multiple members as needed
  #-------------------
    if ((repeat_mems <= fmember)); then
      ipm=0
      for c in $(seq $rcycle); do
        cf=$(printf $CYCLE_FMT $c)
        set_topo_landuse_bases
        for m in $(seq $repeat_mems); do
          ipm=$((ipm+1))
          if ((ipm > parallel_mems)); then wait; ipm=1; fi
          echo "  ${stimesfmt[$c]}: node ${node_m[$m]} [$(datetime_now)]"
          ( pdbash proc.$(printf $CYCLE_FMT 1).${name_m[$m]} alln $SCRP_DIR/src/pre_scale_init.sh $mem_np \
              $topo_base $landuse_base $TMPDAT/wrf/wrfout_d01 \
              ${stimes[$c]} $FCSTLEN $TMPRUN/scale_init/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.$(printf $CYCLE_FMT 1).${name_m[$m]} $TMPRUN/scale_init/${cf} \
              ./scale-les_init init.conf ) &
          sleep $BGJOB_INT
        done
      done
      wait
    else
      ipm=0
      for c in $(seq $rcycle); do
        cf=$(printf $CYCLE_FMT $c)
        set_topo_landuse_bases
        for m in $(seq $fmember); do
          mm=$(((c-1) * fmember + m))
          ipm=$((ipm+1))
          if ((ipm > parallel_mems)); then wait; ipm=1; fi
          echo "  ${stimesfmt[$c]}: node ${node_m[$mm]} [$(datetime_now)]"
          ( pdbash proc.${cf}.${name_m[$mm]} alln $SCRP_DIR/src/pre_scale_init.sh $mem_np \
              $topo_base $landuse_base $TMPDAT/wrf/wrfout_d01 \
              ${stimes[$c]} $FCSTLEN $TMPRUN/scale_init/${cf} \
              $TMPDAT/exec $TMPDAT &&
            mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale_init/${cf} \
              ./scale-les_init init.conf ) &
          sleep $BGJOB_INT
        done
      done
      wait
    fi
  #-------------------
  fi
  #-------------------

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function pertbdy {
#-------------------------------------------------------------------------------

echo

if ((PERTURB_BDY == 0)); then

  echo "  ... skip this step (do not perturb boundaries)"

else

  echo "pertbdy..."

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function ensfcst {
#-------------------------------------------------------------------------------

echo

ipm=0
for c in $(seq $rcycle); do
  cf=$(printf $CYCLE_FMT $c)
  set_topo_landuse_bases
  for m in $(seq $fmember); do
    mm=$(((c-1) * fmember + m))
    ipm=$((ipm+1))
    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    echo "  ${stimesfmt[$c]}, member ${name_m[$mm]}: node ${node_m[$mm]} [$(datetime_now)]"
    if ((TMPRUN_MODE <= 2)); then
      ( pdbash proc.${cf}.${name_m[$mm]} one $SCRP_DIR/src/pre_scale.sh $mem_np \
          $TMPOUT/anal/${name_m[$mm]}/${stimes[$c]} \
          $TMPRUN/scale_init/${cf}/boundary \
          $topo_base $landuse_base \
          ${stimes[$c]} $FCSTLEN $FCSTOUT $TMPRUN/scale/${cf}_${name_m[$mm]} \
          $TMPDAT/exec $TMPDAT &&
        mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \
          ./scale-les run.conf ) &
    else
      ( pdbash proc.${cf}.${name_m[$mm]} alln $SCRP_DIR/src/pre_scale.sh $mem_np \
          $TMPOUT/anal/${name_m[$mm]}/${stimes[$c]} \
          $TMPRUN/scale_init/${cf}/boundary \
          $topo_base $landuse_base \
          ${stimes[$c]} $FCSTLEN $FCSTOUT $TMPRUN/scale/${cf}_${name_m[$mm]} \
          $TMPDAT/exec $TMPDAT &&
        mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \
          ./scale-les run.conf ) &
    fi
    sleep $BGJOB_INT
  done
done
wait



#for c in $(seq $rcycle); do
#  cf=$(printf $CYCLE_FMT $c)
#  for m in $(seq $fmember); do
#    mm=$(((c-1) * fmember + m))

#    pdbash proc.${cf}.${name_m[$mm]} alln $SCRP_DIR/src/pre_scale.sh \
#      $TMPOUT/anal/${name_m[$mm]}/${stimes[$c]} \
#      /data1/gylien/scale/scale-les/test/case_real/ctl_4_small_domain/initial_WRF/boundary \
#      /data1/gylien/scale/scale-les/test/case_real/ctl_4_small_domain/topo/topo \
#      /data1/gylien/scale/scale-les/test/case_real/ctl_4_small_domain/landuse/landuse \
#      ${stimes[$c]} $FCSTLEN $FCSTOUT $TMPRUN/scale/${cf}_${name_m[$mm]} \
#      $TMPDAT/exec $TMPDAT &

#    sleep $BGJOB_INT
#  done
#done
#wait

#-------------------------------------------------------------------------------

#echo
#ppnl=$((ppn*mem_nodes))
#pcount=0
#for c in `seq $rcycle`; do
##  if [ "$LOG_OPT" -le 2 ]; then
#    mkdir -p $OUTDIR/log/scalefcst/${stimes[$c]}
##  fi
#  cf=$(printf $CYCLE_FMT $c)
#  for m in `seq $fmember`; do
#    mm=$(((c-1) * fmember + m))
##    np=`cat $tmpnode/machinefile.gfs.${cf}.${name_m[$mt]} | wc -l`
####    if [ "${node[${mem2proc[$(((mm-1)*mem_np+1))]}]}" = "${node[1]}" ]; then
####      pcount=$((pcount+mem_np))
####      if [ "$pcount" -gt "$ppnl" ]; then
####        echo "    wait..."
####        wait
####        pcount=$mem_np
####      fi
####    fi
#    echo "  run SCALE forecast: ${stimes[$c]}, member ${name_m[$mm]} on node ${node_m[$mm]}"
##    if [ "$LOG_OPT" -le 2 ]; then

##echo "    mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \\"
##echo "            ./scale-les run.conf > $OUTDIR/log/scalefcst/${stimes[$c]}/scale_${name_m[$mm]}.log 2>&1 &"
#    mpirunf proc.${cf}.${name_m[$mm]} $TMPRUN/scale/${cf}_${name_m[$mm]} \
#            ./scale-les run.conf > $OUTDIR/log/scalefcst/${stimes[$c]}/scale_${name_m[$mm]}.log 2>&1 &
##    else
##      $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gfs.${cf}.${name_m[$mt]} -n $mem_np \
##                      -wdir $ltmpgfs/${cf}_${name_m[$mt]} \
##                      ./global_fcst > /dev/null 2>&1 &
##    fi
#    sleep $BGJOB_INT
#  done
#done
#echo "    wait..."
#wait

#-------------------------------------------------------------------------------
}

#===============================================================================

function verf {
#-------------------------------------------------------------------------------

echo

if ((IF_VERF == 0)); then

  echo "  ... skip this step (do not run verification)"

else

  echo "verf..."

fi


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
#if [ "$FOUT_OPT" -le 2 ]; then
#  cat >> fcst_31.sh << EOF
#  ./ss2grdp
#  mv -f fort.31 $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd
#EOF
#fi
#if [ "$FOUT_OPT" -le 1 ]; then
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

function final {
#-------------------------------------------------------------------------------

echo

echo "final..."


#for c in `seq $CYCLES`; do
#  STIMEgrads=$(datetimegrads ${STIME[$c]})
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    if [ "$FOUT_OPT" -le 1 ]; then
#      mkdir -p $OUTDIR/fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}
#    fi
#    mkdir -p $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 x > \
#                     $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhx.ctl
#    if [ "$FOUT_OPT" -le 2 ]; then
#      mkdir -p $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}
#      $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 p > \
#                       $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhp.ctl
#    fi

#    fh=0
#    while [ "$fh" -le "$FCSTLEN" ]; do
#      fhhh=`printf '%03d' $fh`
#      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
#      mkdir -p $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
#      cd $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
#      ln -fs ../../../fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
#      if [ ! -s 'yyyymmddhhx.ctl' ]; then
#        $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
#                         yyyymmddhhx.ctl
#      fi
#      if [ "$FOUT_OPT" -le 2 ]; then
#        mkdir -p $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
#        cd $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
#        ln -fs ../../../fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
#        if [ ! -s 'yyyymmddhhp.ctl' ]; then
#          $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                           yyyymmddhhp.ctl
#        fi
#      fi
#    fh=$((fh + FCSTOUT))
#    done
#  done
#done

##-------------------------------------------------------------------------------
#if [ "$SHAREDISK" = '0' ]; then
##-------------------------------------------------------------------------------

#cd $TMPMPI
#mkdir -p $tmpstageout
#rm -f $tmpstageout/*

#for c in `seq $CYCLES`; do
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig" >> $tmpstageout/out.${node_m[$mt]}
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#    fh=0
#    while [ "$fh" -le "$FCSTLEN" ]; do
#      fhhh=`printf '%03d' $fh`
#      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
#      if [ "$FOUT_OPT" -le 1 ]; then
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$mt]}
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      if [ "$FOUT_OPT" -le 2 ]; then
#        echo "mv|fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|verfo1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa2/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#    fh=$((fh + FCSTOUT))
#    done
#  done
#done

##-------------------------------------------------------------------------------

#stageout $ltmpout 0  # clean stageout
#$MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
#                rm -fr $LTMP1/${tmpsubdir} $LTMP2/${tmpsubdir} &
#wait

##-------------------------------------------------------------------------------
#elif [ "$SHAREDISK" = '1' ]; then
##-------------------------------------------------------------------------------

#rm -fr $ltmprun1 $ltmprun2

##-------------------------------------------------------------------------------
#fi
##-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
}

#===============================================================================
