#!/bin/bash
#===============================================================================
#
#  Steps of 'cycle.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=6
stepname[1]='Prepare boundary files'
stepfunc[1]='boundary'
stepname[2]='Perturb boundaries'
stepfunc[2]='pertbdy'
stepname[3]='Run ensemble forecasts'
stepfunc[3]='ensfcst'
stepname[4]='Thin observations'
stepfunc[4]='obsthin'
stepname[5]='Run observation operator'
stepfunc[5]='obsope'
stepname[6]='Run LETKF'
stepfunc[6]='letkf'

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run data assimilation cycles.

Configuration files:
  config.main
  config.$myname1

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME ISTEP FSTEP TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
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
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-$nsteps}
TIME_LIMIT=${TIME_LIMIT:-"0:30:00"}

#-------------------------------------------------------------------------------
# common variables

CYCLEFLEN=$WINDOW_E   # Model forecast length in a cycle (hour)
CYCLEFOUT=$LTIMESLOT  # Model forecast output interval (hour)

BUILTIN_STAGING=$((MACHINE_TYPE != 10 && MACHINE_TYPE != 11))

if ((TMPRUN_MODE <= 2)); then
  PROC_OPT='one'
else
  PROC_OPT='alln'
fi

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
  ln -fs $OBSUTIL_DIR/obsope $TMPDAT/exec
  ln -fs $LETKF_DIR/letkf $TMPDAT/exec
  ln -fs $DATADIR/rad $TMPDAT
  ln -fs $DATADIR/land $TMPDAT
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
${OBSUTIL_DIR}/obsope|exec/obsope
${LETKF_DIR}/letkf|exec/letkf
${SCRP_DIR}/config.nml.scale|conf/config.nml.scale
${SCRP_DIR}/config.nml.scale_init|conf/config.nml.scale_init
${SCRP_DIR}/config.nml.obsope|conf/config.nml.obsope
${SCRP_DIR}/config.nml.letkf|conf/config.nml.letkf
${DATADIR}/rad|rad
${DATADIR}/land|land
EOF

  time=$STIME
  etime_anlwrf=$(datetime $ETIME $((CYCLEFLEN+ANLWRF_INT)) s)
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
${SCRP_DIR}/config.nml.scale_pp_topo|conf/config.nml.scale_pp_topo
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
${SCRP_DIR}/config.nml.scale_pp_landuse|conf/config.nml.scale_pp_landuse
${DATADIR}/landuse/LU100M/Products|landuse/LU100M/Products
EOF
  fi

  if ((MACHINE_TYPE == 10)); then
    cat >> $STAGING_DIR/stagein.dat << EOF
${COMMON_DIR}/datetime|exec/datetime
EOF
  fi

  time=$(datetime $STIME $LCYCLE s)
  while ((time <= $(datetime $ETIME $LCYCLE s))); do
    cat >> $STAGING_DIR/stagein.dat << EOF
${OBS}/obs_${time}.dat|obs/obs_${time}.dat
${OBS}/radar_${time}.dat|obs/radar_${time}.dat
EOF
    time=$(datetime $time $LCYCLE s)
  done
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
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  loop=0
  while ((time <= ETIME)); do
    loop=$((loop+1))
    if ((ONLINE_STGOUT == 1)); then
      stgoutstep="stageout.loop.${loop}"
    else
      stgoutstep='stageout.out'
    fi

    #-------------------
    # stage-in

    for m in $(seq $mmean); do
      for q in $(seq $mem_np); do

        path="${time}/anal/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+q))]}

      done
    done

    #-------------------
    # stage-out

    for m in $(seq $MEMBER); do
      for q in $(seq $mem_np); do

        if ((OUT_OPT <= 1)); then
          path="${atime}/gues/${name_m[$m]}/history$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi
        if ((OUT_OPT <= 2)); then
          path="${atime}/gues/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        if ((OUT_OPT <= 3)); then
          path="${atime}/anal/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        if ((OBSOUT_OPT <= 2)); then
          path="${atime}/obsgues/${name_m[$m]}/obsda.${name_m[$m]}.$(printf $PROCESS_FMT $((q-1))).dat"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

      done

      #-------------------

      if ((LOG_OPT <= 3)); then
        path="${time}/log/scale/${name_m[$m]}_LOG${SCALE_LOG_SFX}"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
      fi
#          if ((LOG_OPT <= 1)); then
#            # perturb bdy log
#          fi

      #-------------------
    done

    #-------------------

    for q in $(seq $mem_np); do
      #-------------------
      # mean/sprd

      if ((OUT_OPT <= 4)); then
        path="${atime}/gues/mean/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
        path="${atime}/gues/sprd/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}

        path="${atime}/anal/mean/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
        path="${atime}/anal/sprd/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
      fi

      #-------------------
      # meanf

      if ((OUT_OPT <= 1)); then
        path="${atime}/gues/meanf/history$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
      fi
      if ((OUT_OPT <= 4)); then
        path="${atime}/gues/meanf/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
      fi

      #-------------------
    done

    #-------------------

    if ((LOG_OPT <= 4)); then
      path="${atime}/log/obsope/NOUT-000000" ###### 'NOUT-000000' as a variable
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
      path="${atime}/log/letkf/NOUT-000000" ###### 'NOUT-000000' as a variable
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
    fi

    if ((LOG_OPT <= 2)); then
      path="${time}/log/scale_topo/pp_LOG${SCALE_LOG_SFX}"
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
      path="${time}/log/scale_landuse/pp_LOG${SCALE_LOG_SFX}"
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
      path="${time}/log/scale_bdy/init_LOG${SCALE_LOG_SFX}"
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
    fi

    #-------------------

    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
#-------------------
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

boundary_sub () {
#-------------------------------------------------------------------------------
# Run a series of scripts (topo/landuse/init) to make the boundary files.
#
# Usage: boundary_sub NODEFILE
#
#   NODEFILE  $NODEFILE in functions 'pdbash' and 'mpirunf'
#
# Other input variables:
#   $time
#   $SCRP_DIR
#   $TMPRUN
#   $TMPDAT
#   $mem_np
#   $PREP_TOPO
#   $PREP_LANDUSE
#   $PROC_OPT
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"

#-------------------------------------------------------------------------------
# topo

if ((PREP_TOPO != 1)); then
  pdbash $NODEFILE $PROC_OPT \
    $SCRP_DIR/src/pre_scale_pp_topo.sh $time $TMPRUN/scale_pp_topo $TMPDAT/exec $TMPDAT
  mpirunf $NODEFILE \
    $TMPRUN/scale_pp_topo ./scale-les_pp pp.conf
  if ((LOG_OPT <= 2)); then
    pdbash $NODEFILE $PROC_OPT \
      $SCRP_DIR/src/post_scale_pp_topo.sh $time $TMPRUN/scale_pp_topo
  fi
fi

#-------------------------------------------------------------------------------
# landuse

if ((PREP_LANDUSE != 1)); then
  pdbash $NODEFILE $PROC_OPT \
    $SCRP_DIR/src/pre_scale_pp_landuse.sh $time $TMPRUN/scale_pp_landuse $TMPDAT/exec $TMPDAT
  mpirunf $NODEFILE \
    $TMPRUN/scale_pp_landuse ./scale-les_pp pp.conf
  if ((LOG_OPT <= 2)); then
    pdbash $NODEFILE $PROC_OPT \
      $SCRP_DIR/src/post_scale_pp_landuse.sh $time $TMPRUN/scale_pp_landuse
  fi
fi

#-------------------------------------------------------------------------------
# init

if ((PREP_TOPO == 1)); then
  local topo_base="$TMPDAT/topo_prep/topo"
else
  local topo_base="$TMPRUN/scale_pp_topo/topo"
fi
if ((PREP_LANDUSE == 1)); then
  local landuse_base="$TMPDAT/landuse_prep/landuse"
else
  local landuse_base="$TMPRUN/scale_pp_landuse/landuse"
fi
pdbash $NODEFILE $PROC_OPT \
  $SCRP_DIR/src/pre_scale_init.sh $mem_np $topo_base $landuse_base $TMPDAT/wrf/wrfout_d01 \
  $time $CYCLEFLEN $TMPRUN/scale_init $TMPDAT/exec $TMPDAT
mpirunf $NODEFILE \
  $TMPRUN/scale_init ./scale-les_init init.conf
if ((LOG_OPT <= 2)); then
  pdbash $NODEFILE $PROC_OPT \
    $SCRP_DIR/src/post_scale_init.sh $time $TMPRUN/scale_init
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
  echo "  ${timefmt}: node ${node_m[1]} [$(datetime_now)]"

  boundary_sub proc.${name_m[1]}
#-------------------
else # local run directory: run multiple members as needed
#-------------------
#  ipm=0
  for m in $(seq $((repeat_mems <= MEMBER ? repeat_mems : MEMBER))); do
#    ipm=$((ipm+1))
#    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    echo "  ${timefmt}: node ${node_m[$m]} [$(datetime_now)]"

    boundary_sub proc.${name_m[$m]} &
    sleep $BGJOB_INT
  done
  wait
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

if ((PREP_BDY == 1)); then
  bdy_base="$TMPDAT/bdy_prep/bdy_${time}"
else
  bdy_base="$TMPRUN/scale_init/boundary"
fi

ipm=0
for m in $(seq $MEMBER); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

#   ......
#    ( pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_scale.sh $mem_np \
#        $TMPOUT/${time}/anal/${name_m[$m]}/init $bdy_base $topo_base $landuse_base \
#        ${time} $CYCLEFLEN $CYCLEFOUT $TMPRUN/scale_${name_m[$m]} $TMPDAT/exec $TMPDAT ;
#      mpirunf proc.${name_m[$m]} $TMPRUN/scale/${name_m[$m]} \
#        ./scale-les run.conf ) &
#   ......
#   $TMPRUN/pertbdy/${name_m[$m]}

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst () {
#-------------------------------------------------------------------------------

echo

if ((PREP_TOPO == 1)); then
  topo_base="$TMPDAT/topo_prep/topo"
else
  topo_base="$TMPRUN/scale_pp_topo/topo"
fi
if ((PREP_LANDUSE == 1)); then
  landuse_base="$TMPDAT/landuse_prep/landuse"
else
  landuse_base="$TMPRUN/scale_pp_landuse/landuse"
fi

ipm=0
for m in $(seq $mmean); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  if ((PERTURB_BDY == 1)); then
######
    echo "not finished yet..."
#      bdy_base="$TMPRUN/pertbdy/${name_m[$m]}/boundary"
######
  elif ((PREP_BDY == 1)); then
    bdy_base="$TMPDAT/bdy_prep/bdy_${time}"
  else
    bdy_base="$TMPRUN/scale_init/boundary"
  fi
  ( pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_scale.sh $mem_np \
      $TMPOUT/${time}/anal/${name_m[$m]}/init $bdy_base $topo_base $landuse_base \
      ${time} $CYCLEFLEN $LCYCLE $CYCLEFOUT $TMPRUN/scale/${name_m[$m]} $TMPDAT/exec $TMPDAT ;
    mpirunf proc.${name_m[$m]} $TMPRUN/scale/${name_m[$m]} \
      ./scale-les run.conf > /dev/null ;
    pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/post_scale.sh $mem_np \
      ${time} ${name_m[$m]} $CYCLEFLEN $TMPRUN/scale/${name_m[$m]} $myname1 ) &

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

obsthin () {
#-------------------------------------------------------------------------------

echo
if ((THINNING == 0)); then
  echo "  ... skip this step (do not thin observations)"
  return 1
fi

###### not finished yet...
echo "obsthin..."
######


#-------------------------------------------------------------------------------
}

#===============================================================================

obsope () {
#-------------------------------------------------------------------------------

echo

pdbash node $PROC_OPT $SCRP_DIR/src/pre_obsope_node.sh \
  $atime $TMPRUN/obsope $TMPDAT/exec $TMPDAT/obs \
  $mem_nodes $mem_np $slot_s $slot_e $slot_b

ipm=0
for m in $(seq $MEMBER); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_obsope.sh \
    $atime ${name_m[$m]} $TMPRUN/obsope &

  sleep $BGJOB_INT
done
wait

mpirunf proc $TMPRUN/obsope ./obsope obsope.conf > /dev/null

ipm=0
for m in $(seq $MEMBER); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
#  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/post_obsope.sh \
    ${atime} ${name_m[$m]} $TMPRUN/obsope &

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

letkf () {
#-------------------------------------------------------------------------------

echo

if ((PREP_TOPO == 1)); then
  topo_base="$TMPDAT/topo_prep/topo"
else
  topo_base="$TMPRUN/scale_pp_topo/topo"
fi

pdbash node $PROC_OPT $SCRP_DIR/src/pre_letkf_node.sh \
  $atime $TMPRUN/letkf $TMPDAT/exec $TMPDAT/obs \
  $mem_nodes $mem_np $slot_s $slot_e $slot_b

ipm=0
for m in $(seq $mmean); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_letkf.sh \
    $topo_base $atime ${name_m[$m]} $TMPRUN/letkf &

  sleep $BGJOB_INT
done
wait

#exit

mpirunf proc $TMPRUN/letkf ./letkf letkf.conf > /dev/null

ipm=0
for m in $(seq $mmean); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
#  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/post_letkf.sh \
    ${atime} ${name_m[$m]} $TMPRUN/letkf &

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================
