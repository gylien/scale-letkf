#!/bin/bash
#===============================================================================
#
#  Run data assimilation cycles.
#
#  November 2014, modified from GFS-LETKF, Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle.sh [..]
#
#  Use settings:
#    config.main
#    config.cycle
#    config.nml.scale_pp
#    config.nml.scale_init
#    config.nml.scale
#    config.nml.scale_user
#    config.nml.grads_boundary
#    config.nml.obsope
#    config.nml.letkf
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='cycle'

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

echo "[$(datetime_now)] ### 1" >&2

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@" >&2

setting "$@" || exit $?

if [ "$CONF_MODE" = 'static' ]; then
  . src/func_common_static.sh || exit $?
  . src/func_${job}_static.sh || exit $?
fi

echo
print_setting || exit $?

echo "[$(datetime_now)] ### 2" >&2

#===============================================================================
# Initialize temporary directories

if ((RUN_LEVEL <= 2)) && ((ISTEP == 1)); then
  safe_init_tmpdir $TMP || exit $?
fi

echo "[$(datetime_now)] ### 3" >&2

#===============================================================================
# Determine the distibution schemes

declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

#if ((RUN_LEVEL <= 2)) && ((ISTEP == 1)); then
if ((RUN_LEVEL <= 2)); then
  safe_init_tmpdir $NODEFILE_DIR || exit $?
fi
#distribute_da_cycle - $NODEFILE_DIR || exit $? 
distribute_da_cycle - - || exit $?  # Do not use distr

echo "[$(datetime_now)] ### 4" >&2

#===============================================================================
# Determine the staging list and then stage in

if ((RUN_LEVEL <= 1)) && ((ISTEP == 1)); then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
  if [ "$CONF_MODE" = 'static' ]; then
    staging_list_static || exit $?
    if ((DISK_MODE == 3)); then
      config_file_list $TMP/config || exit $?
    else
      config_file_list || exit $?
    fi
  else
    staging_list || exit $?
  fi

  stage_in node || exit $?
fi

echo "[$(datetime_now)] ### 5" >&2

#===============================================================================
# Run data assimilation cycles

function online_stgout_bgjob () {
  local ILOOP="$1"; shift
  local ITIME="$1"
  touch lock.$ILOOP
  echo "[$(datetime_now)] ${ITIME}: Stage-out (background job)" >&2
  while [ -e "lock.$((ILOOP-1))" ]; do
    sleep 1s
  done

  stage_out node $ILOOP || exit $?

  echo "[$(datetime_now)] ${ITIME}: Stage-out (background job completed)" >&2
  rm -f lock.$ILOOP
}

#-------------------------------------------------------------------------------

cd $TMPROOT

#-------------------------------------------------------------------------------

s_flag=1
e_flag=0
time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

  timefmt="$(datetime_fmt ${time})"
  loop=$((loop+1))
  if (($(datetime $time $LCYCLE s) > ETIME)); then
    e_flag=1
  fi
  obstime $time || exit $?

#-------------------------------------------------------------------------------
# Write the header of the log file

  echo "[$(datetime_now)] ### 7" >&2

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                          SCALE-LETKF                           |"
  echo " +----------------------------------------------------------------+"
  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then
      printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
    fi
  done
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Start time:               ${timefmt}"
  echo "  Forecast length:          $CYCLEFLEN s"
  echo "  Assimilation window:      $WINDOW_S - $WINDOW_E s ($((WINDOW_E-WINDOW_S)) s)"
  echo
  echo "  Observation timeslots:"
  for is in $(seq $slot_s $slot_e); do
    if ((is == slot_b)); then
      printf "  %4d - %s [base]\n" ${is} "${timefmt_sl[$is]}"
    else
      printf "  %4d - %s\n" ${is} "${timefmt_sl[$is]}"
    fi
  done
  echo
  echo "  Nodes used:               $NNODES_APPAR"
#  for n in $(seq $NNODES_APPAR); do
#    echo "    ${node[$n]}"
#  done
  echo
  echo "  Processes per node:       $PPN_APPAR"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Ensemble size:            $MEMBER"
#  for m in $(seq $mtot); do
#    echo "      ${name_m[$m]}: ${node_m[$m]}"
#  done
  echo

#-------------------------------------------------------------------------------
# Call functions to run the job

  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

      ######
      if ((s == 1)); then
        if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared topo and landuse files)" >&2
          continue
        elif ((BDY_FORMAT == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared boundary files)" >&2
          continue
        elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (already done in the first cycle)" >&2
          continue
        fi
      fi
      if ((s == 2)); then
        if ((BDY_FORMAT == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared boundary files)" >&2
          continue
        fi
      fi
      ######

      echo "[$(datetime_now)] ${time}: ${stepname[$s]}" >&2

      enable_iter=0
      nit=1

      if ((s <= 3)); then
        conf_time=$time
      else
        conf_time=$atime
      fi

      noit='-'

      $MPIRUN -n ${PJM_MPI_PROC} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf $noit log/${stepexecname[$s]}.NOUT_${conf_time} || exit $?

    fi
  done

#-------------------------------------------------------------------------------
# Online stage out

  if ((RUN_LEVEL <= 3)); then
    if ((ONLINE_STGOUT == 1)); then
      online_stgout_bgjob $loop $time &
    fi
  fi

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo " +----------------------------------------------------------------+"
  echo " |               SCALE-LETKF successfully completed               |"
  echo " +----------------------------------------------------------------+"
  echo

#-------------------------------------------------------------------------------

  if ((DACYCLE == 1)); then
    if ((DACYCLE_RUN_FCST == 1)); then
      cp -f ${TMP}/*.png ${OUTDIR}/${STIME}/dafcst
    fi
    break
  fi

#-------------------------------------------------------------------------------

  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
  s_flag=0

#-------------------------------------------------------------------------------
done
#-------------------------------------------------------------------------------

#===============================================================================
# Stage out

if ((RUN_LEVEL <= 3)); then
  if ((ONLINE_STGOUT == 1)); then
    wait
  else
    echo "[$(datetime_now)] Finalization (stage out)" >&2

    stage_out node || exit $?
  fi

  if ((CLEAR_TMP == 1 && USE_TMPL == 1)); then
    pdbash node all $SCRP_DIR/src/stage_out_rm_stgdir_node.sh $TMPL local # || exit $?
  fi
fi

if ((RUN_LEVEL <= 1)); then
  if [ "$CONF_MODE" = 'static' ]; then
    if ((DISK_MODE == 3)); then
      config_file_save $TMP/config || exit $?
    else
      config_file_save || exit $?
    fi
  fi
fi

#===============================================================================
# Remove temporary directories

if ((RUN_LEVEL <= 3)); then
  if ((CLEAR_TMP == 1)); then
    safe_rm_tmpdir $TMP || exit $?
  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
