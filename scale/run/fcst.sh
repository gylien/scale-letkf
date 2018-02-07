#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#
#  August  2014, modified from GFS-LETKF, Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst.sh [..]
#
#  Use settings:
#    config.main
#    config.fcst
#    config.nml.scale_pp
#    config.nml.scale_init
#    config.nml.scale
#    config.nml.scale_user
#    config.nml.grads_boundary
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='fcst'

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@" >&2

setting "$@" || exit $?

###if [ "$CONF_MODE" = 'static' ]; then
###  . src/func_${job}_static.sh || exit $?
###fi

echo
print_setting || exit $?

#===============================================================================
# Initialize temporary directories

if ((RUN_LEVEL <= 2)) && ((ISTEP == 1)); then
  safe_init_tmpdir $TMP || exit $?
fi

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
distribute_fcst "$MEMBERS" $CYCLE "$NODELIST_TYPE" $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list and then stage in

if ((RUN_LEVEL <= 1)) && ((ISTEP == 1)); then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
###  if [ "$CONF_MODE" = 'static' ]; then
###    staging_list_static || exit $?
###    if ((DISK_MODE == 3)); then
###      config_file_list $TMP/config || exit $?
###    else
###      config_file_list || exit $?
###    fi
###  else
    staging_list || exit $?
###  fi

  stage_in node || exit $?
fi

#===============================================================================
# Run cycles of forecasts

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

declare -a stimes
declare -a stimesfmt
lcycles=$((LCYCLE * CYCLE_SKIP))
s_flag=1
e_flag=0
time=$STIME
loop=0

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

  loop=$((loop+1))

  for c in $(seq $CYCLE); do
    time2=$(datetime $time $((lcycles * (c-1))) s)
    if (($(datetime $time2 $lcycles s) > ETIME)); then
      e_flag=1
    fi

    if ((time2 <= ETIME)); then
      stimes[$c]=$time2
      stimesfmt[$c]="$(datetime_fmt ${stimes[$c]})"
      rcycle=$c  # The "real" number of cycles
    else
      stimes[$c]=
      stimesfmt[$c]=
    fi
  done

#-------------------------------------------------------------------------------
# Write the header of the log file

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                        SCALE-Forecasts                         |"
  echo " +----------------------------------------------------------------+"
  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then
      printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
    fi
  done
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Number of cycles:         $rcycle"
  echo "  Forecast start time:"
  for c in $(seq $rcycle); do
    printf "    Cycle %-5s %s\n" "$c:" "${stimesfmt[$c]}"
  done
  echo
  echo "  Forecast length:          $FCSTLEN s"
  echo "  Output interval:          $FCSTOUT s"
  echo
  echo "  Nodes used:               $NNODES_APPAR"
#  if ((MTYPE == 1)); then
    for n in $(seq $NNODES_APPAR); do
      echo "    ${node[$n]}"
    done
#  fi
  echo
  echo "  Processes per node:       $PPN_APPAR"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Number of members:        $fmember"
  for c in $(seq $rcycle); do
    echo "    Cycle $c:"
    for m in $(seq $fmember); do
      mm=$(((c-1) * fmember + m))
      echo "      ${name_m[$m]}: ${node_m[$mm]}"
    done
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
      ######

      echo "[$(datetime_now)] ${time}: ${stepname[$s]}" >&2

      enable_iter=0
      if ((s == 2 && BDY_ENS == 1)); then
        enable_iter=1
      elif ((s == 3)); then
        enable_iter=1
      fi

###      if [ "$CONF_MODE" = 'static' ]; then

###        mpirunf ${nodestr} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf log/${stepexecname[$s]}.NOUT_${conf_time} || exit $?

###      else

        execpath="${stepexecdir[$s]}/${stepexecname[$s]}"
        stdout_dir="$TMPOUT/${stimes[1]}/log/fcst_$(basename ${stepexecdir[$s]})"
        if ((enable_iter == 1)); then
          for it in $(seq $nitmax); do
            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: start" >&2

            mpirunf proc $execpath ${execpath}.conf "${stdout_dir}/NOUT-${it}" "$SCRP_DIR/fcst_step.sh" $loop $it || exit $?

            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: end" >&2
          done
        else
          mpirunf proc $execpath ${execpath}.conf "${stdout_dir}/NOUT" "$SCRP_DIR/fcst_step.sh" $loop || exit $?
        fi

###      fi

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
  echo " |             SCALE-Forecasts successfully completed             |"
  echo " +----------------------------------------------------------------+"
  echo

#-------------------------------------------------------------------------------

  time=$(datetime $time $((lcycles * CYCLE)) s)
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

###if ((RUN_LEVEL <= 1)); then
###  if [ "$CONF_MODE" = 'static' ]; then
###    if ((DISK_MODE == 3)); then
###      config_file_save $TMP/config || exit $?
###    else
###      config_file_save || exit $?
###    fi
###  fi
###fi

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
