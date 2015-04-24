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
  ln -fs $DATADIR/rad $TMPDAT/rad
  ln -fs $DATADIR/land $TMPDAT/land
  ln -fs $DATADIR/topo $TMPDAT
  ln -fs $DATADIR/landuse $TMPDAT

  ln -fs $OBS $TMPDAT/obs

  safe_init_tmpdir $TMPDAT/conf
  ln -fs $SCRP_DIR/config.* $TMPDAT/conf
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

  if [ "$TOPO_FORMAT" != 'prep' ] || [ "$LANDUSE_FORMAT" != 'prep' ]; then
    echo "${SCRP_DIR}/config.nml.scale_pp|conf/config.nml.scale_pp" >> $STAGING_DIR/stagein.dat
  fi
  if [ "$TOPO_FORMAT" != 'prep' ]; then
    echo "${DATADIR}/topo/${TOPO_FORMAT}/Products|topo/${TOPO_FORMAT}/Products" >> $STAGING_DIR/stagein.dat
  fi
  if [ "$LANDUSE_FORMAT" != 'prep' ]; then
    echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products|landuse/${LANDUSE_FORMAT}/Products" >> $STAGING_DIR/stagein.dat
  fi

  time=$(datetime $STIME $LCYCLE s)
  while ((time <= $(datetime $ETIME $LCYCLE s))); do
    for iobs in $(seq $OBSNUM); do
      if [ "${OBSNAME[$iobs]}" != '' ]; then
        echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs/${OBSNAME[$iobs]}_${time}.dat" >> $STAGING_DIR/stagein.dat
      fi
    done
    time=$(datetime $time $LCYCLE s)
  done

  if ((MACHINE_TYPE == 10)); then
    echo "${COMMON_DIR}/datetime|exec/datetime" >> $STAGING_DIR/stagein.dat
  fi
#-------------------
fi

#-------------------------------------------------------------------------------

if ((TMPOUT_MODE == 1 && MACHINE_TYPE != 10)); then
#-------------------
  mkdir -p $(dirname $TMPOUT)
  ln -fs $OUTDIR $TMPOUT

  time=$STIME
  while ((time <= ETIME)); do
    #-------------------
    if [ "$TOPO_FORMAT" = 'prep' ]; then
      ln -fs ${DATA_TOPO} $TMPOUT/${time}/topo
    fi
    if [ "$LANDUSE_FORMAT" = 'prep' ]; then
      if ((LANDUSE_UPDATE == 1)); then
        ln -fs ${DATA_LANDUSE}/${time} $TMPOUT/${time}/landuse
      else
        ln -fs ${DATA_LANDUSE} $TMPOUT/${time}/landuse
      fi
    fi
    time=$(datetime $time $LCYCLE s)
    #-------------------
  done

  if ((BDY_FORMAT == 2)); then
    ln -fs $DATA_BDY_WRF $TMPOUT/bdywrf
  fi
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
    #-------------------

    # anal
    #-------------------
    if ((loop == 1)) && ((MAKEINIT != 1)); then
      for m in $(seq $mmean); do
        for q in $(seq $mem_np); do
          path="${time}/anal/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    fi

    # topo
    #-------------------
    if [ "$TOPO_FORMAT" = 'prep' ]; then
      for m in $(seq $mmean); do
        for q in $(seq $mem_np); do
          pathin="${DATA_TOPO}/topo$(printf $SCALE_SFX $((q-1)))"
          path="${time}/topo/topo$(printf $SCALE_SFX $((q-1)))"
          echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    fi

    # landuse
    #-------------------
    if [ "$LANDUSE_FORMAT" = 'prep' ]; then
      if [ -d "${DATA_LANDUSE}/${time}" ]; then
        pathin="${DATA_LANDUSE}/${time}"
      else
        pathin="${DATA_LANDUSE}"
      fi
      for m in $(seq $mmean); do
        for q in $(seq $mem_np); do
          pathin="${pathin}/landuse$(printf $SCALE_SFX $((q-1)))"
          path="${time}/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
          echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    fi

    # bdy
    #-------------------
    if ((BDY_FORMAT == 2)); then
      for m in $(seq $mmean); do
        for q in $(seq $mem_np); do
          time_dby=${time}
          etime_bdy=$(datetime ${time} $((CYCLEFLEN+BDYINT)) s)
          while ((time_dby < etime_bdy)); do
            if ((BDY_ENS == 1)); then
              pathin="$DATA_BDY_WRF/${name_m[$m]}/wrfout_${time_dby}"
              path="bdywrf/${name_m[$m]}/wrfout_${time_dby}"
            else
              pathin="$DATA_BDY_WRF/mean/wrfout_${time_dby}"
              path="bdywrf/mean/wrfout_${time_dby}"
            fi
            echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+q))]}
            time_dby=$(datetime $time_dby $BDYINT s)
          done
        done
      done
    fi

    #-------------------
    # stage-out
    #-------------------

    for m in $(seq $MEMBER); do
      #-------------------

      for q in $(seq $mem_np); do
        #-------------------

        # bdy [members]
        #-------------------
        if ((BDYOUT_OPT <= 1)) && ((BDY_ENS == 1)); then
          path="${time}/bdy/${name_m[$m]}/boundary$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        # gues [history]
        #-------------------
        if ((OUT_OPT <= 1)); then
          path="${atime}/gues/${name_m[$m]}/history$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        # gues [restart]
        #-------------------
        if ((OUT_OPT <= 2)); then
          path="${atime}/gues/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        # anal
        #-------------------
        if ((loop == 1)) && ((MAKEINIT == 1)); then
          path="${time}/anal/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi
        if ((OUT_OPT <= 3)); then
          path="${atime}/anal/${name_m[$m]}/init$(printf $SCALE_SFX $((q-1)))"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        # obsgues
        #-------------------
        if ((OBSOUT_OPT <= 2)); then
          path="${atime}/obsgues/${name_m[$m]}/obsda.${name_m[$m]}.$(printf $PROCESS_FMT $((q-1))).dat"
          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+q))]}
        fi

        #-------------------
      done

#      if ((LOG_OPT <= 1)); then
#        # perturb bdy log
#      fi

      #-------------------
    done

    #-------------------

    for q in $(seq $mem_np); do
      #-------------------

      # topo
      #-------------------
      if ((TOPOOUT_OPT <= 1)); then
        path="${time}/topo/topo$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$q]} # ues m=1 instead of m=mmean to enhance parallelization
      fi

      # landuse
      #-------------------
      if ((LANDUSEOUT_OPT <= 1)); then
        path="${time}/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$q]} # ues m=1 instead of m=mmean to enhance parallelization
      fi

      # bdy [mean]
      #-------------------
      if ((BDYOUT_OPT <= 2)) && ((BDY_ENS != 1)); then
        path="${time}/bdy/mean/boundary$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$q]} # ues m=1 instead of m=mmean to enhance parallelization
      fi

      # mean/sprd
      #-------------------
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

      # meanf
      #-------------------
      if ((OUT_OPT <= 1)); then
        path="${atime}/gues/meanf/history$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
      fi
      if ((OUT_OPT <= 4)); then
        path="${atime}/gues/meanf/init$(printf $SCALE_SFX $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+q))]}
      fi

      # log [obsope/letkf]
      #-------------------
      if ((LOG_OPT <= 4)); then
        path="${atime}/log/obsope/NOUT-$(printf $PROCESS_FMT $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
        path="${atime}/log/letkf/NOUT-$(printf $PROCESS_FMT $((q-1)))"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
      fi

      #-------------------
    done

    # log [scale_pp]
    #-------------------
    if ((LOG_OPT <= 2)); then
      path="${time}/log/scale_pp/pp_LOG${SCALE_LOG_SFX}"
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
    fi

    # log [scale_init: mean]
    #-------------------
    if ((LOG_OPT <= 2)) && ((BDY_ENS != 1)); then
      path="${time}/log/scale_init/mean_init_LOG${SCALE_LOG_SFX}"
      echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[1]}
    fi

    #-------------------

    for m in $(seq $mmean); do
      #-------------------

      # log [scale_init: members]
      #-------------------
      if ((LOG_OPT <= 2)) && ((BDY_ENS == 1)); then
        path="${time}/log/scale_init/${name_m[$m]}_init_LOG${SCALE_LOG_SFX}"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
      fi

      # log [scale]
      #-------------------
      if ((LOG_OPT <= 3)); then
        path="${time}/log/scale/${name_m[$m]}_LOG${SCALE_LOG_SFX}"
        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
      fi

      #-------------------
    done

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
# pp (topo/landuse)

if [ "$TOPO_FORMAT" != 'prep' ] || [ "$LANDUSE_FORMAT" != 'prep' ]; then
  pdbash $NODEFILE $PROC_OPT \
    $SCRP_DIR/src/pre_scale_pp.sh $time $TMPRUN/scale_pp $TMPDAT/exec $TMPDAT
  mpirunf $NODEFILE \
    $TMPRUN/scale_pp ./scale-les_pp pp.conf
  pdbash $NODEFILE $PROC_OPT \
    $SCRP_DIR/src/post_scale_pp.sh $time $TMPRUN/scale_pp
fi

#-------------------------------------------------------------------------------
# init

if ((BDY_ENS != 1)); then
  if ((BDY_FORMAT == 2)); then
    pdbash $NODEFILE $PROC_OPT \
      $SCRP_DIR/src/pre_scale_init.sh $mem_np \
      $TMPOUT/${time}/topo/topo $TMPOUT/${time}/landuse/landuse \
      $TMPOUT/bdywrf/mean/wrfout \
      $time $CYCLEFLEN mean $TMPRUN/scale_init/mean $TMPDAT/exec $TMPDAT
    mpirunf $NODEFILE \
      $TMPRUN/scale_init/mean ./scale-les_init init.conf
    pdbash $NODEFILE $PROC_OPT \
      $SCRP_DIR/src/post_scale_init.sh $time mean $TMPRUN/scale_init/mean
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

boundary () {
#-------------------------------------------------------------------------------

echo
if ((BDY_ENS == 1)); then
  echo "     -- topo/landuse"
  echo
fi

if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
#-------------------
  echo "  ${timefmt}: node ${node_m[1]} [$(datetime_now)]"

  boundary_sub proc.${name_m[1]}
#-------------------
else # local run directory: run multiple members as needed
#-------------------
  for m in $(seq $((repeat_mems <= mmean ? repeat_mems : mmean))); do
    echo "  ${timefmt}: node ${node_m[$m]} [$(datetime_now)]"

    boundary_sub proc.${name_m[$m]} &
    sleep $BGJOB_INT
  done
  wait
#-------------------
fi

#-------------------------------------------------------------------------------

if ((BDY_ENS == 1)); then
  echo
  echo "     -- boundary"
  echo

  ipm=0
  for m in $(seq $mmean); do
    ipm=$((ipm+1))
    if ((ipm > parallel_mems)); then wait; ipm=1; fi
    echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

#    if ((PERTURB_BDY == 1)); then
#      ...
#    fi

    if ((BDY_FORMAT == 2)); then
      ( pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_scale_init.sh $mem_np \
          $TMPOUT/${time}/topo/topo $TMPOUT/${time}/landuse/landuse \
          $TMPOUT/bdywrf/mean/wrfout $time $CYCLEFLEN ${name_m[$m]} \
          $TMPRUN/scale_init/${name_m[$m]} $TMPDAT/exec $TMPDAT ;
        mpirunf proc.${name_m[$m]} \
          $TMPRUN/scale_init/${name_m[$m]} ./scale-les_init init.conf ;
        pdbash proc.${name_m[$m]} $PROC_OPT \
          $SCRP_DIR/src/post_scale_init.sh $time ${name_m[$m]} \
          $TMPRUN/scale_init/${name_m[$m]} ) &
    fi

    sleep $BGJOB_INT
  done
  wait
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

ipm=0
for m in $(seq $mmean); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

#  if ((PERTURB_BDY == 1)); then
#    ...
#  fi

  if ((BDY_ENS == 1)); then
    bdy_base="$TMPOUT/${time}/bdy/${name_m[$m]}/boundary"
  else
    bdy_base="$TMPOUT/${time}/bdy/mean/boundary"
  fi
  ( pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_scale.sh $mem_np \
      $TMPOUT/${time}/anal/${name_m[$m]}/init $bdy_base \
      $TMPOUT/${time}/topo/topo $TMPOUT/${time}/landuse/landuse \
      $time $CYCLEFLEN $LCYCLE $CYCLEFOUT $TMPRUN/scale/${name_m[$m]} $TMPDAT/exec $TMPDAT ;
    mpirunf proc.${name_m[$m]} $TMPRUN/scale/${name_m[$m]} \
      ./scale-les run.conf > /dev/null ;
    pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/post_scale.sh $mem_np \
      $time ${name_m[$m]} $CYCLEFLEN $TMPRUN/scale/${name_m[$m]} $myname1 ) &

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
    $mem_np ${atime} ${name_m[$m]} $TMPRUN/obsope &

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================

letkf () {
#-------------------------------------------------------------------------------

echo

pdbash node $PROC_OPT $SCRP_DIR/src/pre_letkf_node.sh \
  $atime $TMPRUN/letkf $TMPDAT/exec $TMPDAT/obs \
  $mem_nodes $mem_np $slot_s $slot_e $slot_b

ipm=0
for m in $(seq $mmean); do
  ipm=$((ipm+1))
  if ((ipm > parallel_mems)); then wait; ipm=1; fi
  echo "  ${timefmt}, member ${name_m[$m]}: node ${node_m[$m]} [$(datetime_now)]"

  pdbash proc.${name_m[$m]} $PROC_OPT $SCRP_DIR/src/pre_letkf.sh \
    $TMPOUT/${time}/topo/topo $atime ${name_m[$m]} $TMPRUN/letkf &

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
    $mem_np ${atime} ${name_m[$m]} $TMPRUN/letkf &

  sleep $BGJOB_INT
done
wait

#-------------------------------------------------------------------------------
}

#===============================================================================
