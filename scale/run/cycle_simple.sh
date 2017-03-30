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
#    cycle.sh [STIME ETIME ISTEP FSTEP]
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
#-------------------------------------------------------------------------------
#
#  RUN_LEVEL:
#    0: Run everything (default)
#    1: Staging list files are ready; skip generating them (not implemented yet...)
#    2: File staging has been done; skip staging
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"

#===============================================================================
# Configuration

. config.main || exit $?
. config.cycle || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_cycle.sh || exit $?
. src/func_cycle_simple.sh || exit $?

RUN_LEVEL=${RUN_LEVEL:-0}

echo "[$(datetime_now)] ### 1" >&2

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@" >&2

setting "$1" "$2" "$3" "$4" "$5" || exit $?

echo
print_setting || exit $?

echo "[$(datetime_now)] ### 2" >&2

#===============================================================================
# Initialize temporary directories

if ((RUN_LEVEL <= 1)) && ((ISTEP == 1)); then
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

#if ((RUN_LEVEL <= 1)) && ((ISTEP == 1)); then
if ((RUN_LEVEL <= 1)); then
  safe_init_tmpdir $NODEFILE_DIR || exit $?
  distribute_da_cycle machinefile $NODEFILE_DIR || exit $?
else
  distribute_da_cycle - - $NODEFILE_DIR/distr || exit $?
fi

echo "[$(datetime_now)] ### 4" >&2

#===============================================================================
# Determine the staging list and then stage in

if ((RUN_LEVEL <= 1)) && ((ISTEP == 1)); then
  echo "[$(datetime_now)] Initialization (stage-in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
  staging_list_simple || exit $?

  config_file_list || exit $?

  stage_in node || exit $?
fi

echo "[$(datetime_now)] ### 5" >&2

#===============================================================================
# Run initialization scripts on all nodes

#if ((TMPRUN_MODE <= 2)); then
#  pdbash node one $SCRP_DIR/src/init_all_node.sh cycle || exit $?
#else
#  pdbash node all $SCRP_DIR/src/init_all_node.sh cycle || exit $?
#fi

echo "[$(datetime_now)] ### 6" >&2

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

cd $TMP_DIR

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
  for n in $(seq $NNODES_APPAR); do
    echo "    ${node[$n]}"
  done
  echo
  echo "  Processes per node:       $PPN_APPAR"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Ensemble size:            $MEMBER"
  for m in $(seq $mtot); do
    echo "      ${name_m[$m]}: ${node_m[$m]}"
  done
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
      if ((s == 4)); then
        if ((OBSOPE_RUN == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (only use integrated observation operators)" >&2
          continue
        fi
      fi
      ######

      echo "[$(datetime_now)] ${time}: ${stepname[$s]}" >&2

      nodestr=proc
######      if ((IO_ARB == 1)); then
######        if ((s == 3)); then
######          nodestr='set1.proc'
######        elif ((s == 5)); then
######          nodestr='set2.proc'
######        fi
######      fi

      if ((s <= 3)); then
        conf_time=$time
      else
        conf_time=$atime
      fi

######      if ((IO_ARB == 1)); then ##                                 
######        mpirunf ${nodestr} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf log/${stepexecname[$s]}.NOUT_${conf_time} &
######      else ##
        mpirunf ${nodestr} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf log/${stepexecname[$s]}.NOUT_${conf_time}
######      fi ##

    fi
  done

######  if ((IO_ARB == 1)); then ##                                 
######    wait                   ##
######  fi                       ##

#-------------------------------------------------------------------------------
# Online stage out

  if ((RUN_LEVEL <= 1)); then
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

  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
  s_flag=0

#-------------------------------------------------------------------------------
done
#-------------------------------------------------------------------------------

#===============================================================================
# Stage out

if ((RUN_LEVEL <= 1)); then
  if ((ONLINE_STGOUT == 1)); then
    wait
  else
    echo "[$(datetime_now)] Finalization (stage-out)" >&2

    stage_out node || exit $?
  fi

  if ((CLEAR_TMP == 1 && USE_TMPL == 1)); then
    pdbash node all $SCRP_DIR/src/stage_out_rm_stgdir_node.sh $TMPL local # || exit $?
  fi
fi

#===============================================================================
# Remove temporary directories

if ((RUN_LEVEL <= 1)); then
  if ((CLEAR_TMP == 1)); then
    safe_rm_tmpdir $TMP || exit $?
  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
