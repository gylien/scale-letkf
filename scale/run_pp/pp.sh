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
job='pp'

#===============================================================================
# Configuration

. config.main || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_pp.sh || exit $?

#-------------------------------------------------------------------------------

CYCLE=0
ISTEP=1
FSTEP=1
echo "[$(datetime_now)] Start $myname $@" >&2

setting "$@" || exit $?

if [ "$CONF_MODE" = 'static' ]; then
  . src/func_common_static.sh || exit $?
  . src/func_fcst_static.sh || exit $?
fi

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
distribute_fcst "mdet" 1 "$NODELIST_TYPE" $NODEFILE_DIR || exit $?

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
#-------------------------------------------------------------------------------
# Write the header of the log file

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                        SCALE-PP                                |"
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Nodes used:               $NNODES_APPAR"
 echo
  echo "  Processes per node:       $PPN_APPAR"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
 echo

#-------------------------------------------------------------------------------
# Call functions to run the job

      ######
      s=1
      loop=1
      rcycle=1
      echo "[$(datetime_now)] ${time}: ${stepname[$s]}" >&2

      nodestr=proc

      if [ "$CONF_MODE" = 'static' ]; then

        mkdir -p log/${stepexecname[$s]}
          mpirunf ${nodestr} ./${stepexecname[$s]} fcst_${stepexecname[$s]}_${stimes[1]}.conf log/${stepexecname[$s]}/NOUT_${stimes[1]} || exit $?

      else

        execpath="${stepexecdir[$s]}/${stepexecname[$s]}"
        stdout_dir="$TMPOUT/${stimes[1]}/log/fcst_$(basename ${stepexecdir[$s]})"
          mpirunf proc $execpath ${execpath}.conf "${stdout_dir}/NOUT" "$SCRP_DIR/pp_step.sh" $loop || exit $?

      fi

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo " +----------------------------------------------------------------+"
  echo " |             SCALE-PP        successfully completed             |"
  echo " +----------------------------------------------------------------+"
  echo

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
