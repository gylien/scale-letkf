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
#===============================================================================

cd "$(dirname "$0")"
myname='cycle.sh'
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.main || exit $?
. config.$myname1 || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_$myname1.sh || exit $?
. src/func_${myname1}_simple.sh || exit $?

echo "[$(datetime_now)] ### 1" >&2

#-------------------------------------------------------------------------------

#if [ "$STG_TYPE" = 'K_rankdir' ]; then
#  SCRP_DIR="."
#  if ((TMPDAT_MODE <= 2)); then
#    TMPDAT="../dat"
#  else
#    TMPDAT="./dat"
#  fi
#  if ((TMPRUN_MODE <= 2)); then
#    TMPRUN="../run"
#  else
#    TMPRUN="./run"
#  fi
#  if ((TMPOUT_MODE <= 2)); then
#    TMPOUT="../out"
#  else
#    TMPOUT="./out"
#  fi
#fi

echo "[$(datetime_now)] Start $myname $@" >&2

setting "$1" "$2" "$3" "$4" "$5" || exit $?

echo
print_setting || exit $?

echo "[$(datetime_now)] ### 2" >&2

#-------------------------------------------------------------------------------

if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
  safe_init_tmpdir $TMP || exit $?
fi

echo "[$(datetime_now)] ### 3" >&2

#===============================================================================
# Determine the distibution schemes

declare -a node
declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

#if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
if [ "$STG_TYPE" = 'builtin' ]; then
  safe_init_tmpdir $NODEFILE_DIR || exit $?
  distribute_da_cycle machinefile $NODEFILE_DIR || exit $?
else
  distribute_da_cycle - - $NODEFILE_DIR/distr || exit $?
fi

echo "[$(datetime_now)] ### 4" >&2

#===============================================================================
# Determine the staging list and then stage in

if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
  echo "[$(datetime_now)] Initialization (stage-in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
  staging_list_simple || exit $?

  config_file_list || exit $?

  if [ -s "$STAGING_DIR/stagein_link.1" ] || [ -s "$STAGING_DIR/stagein_link" ]; then
#    safe_init_tmpdir $TMP || exit $?
    errmsg=$(bash $SCRP_DIR/src/stage_in_ln.sh $NNODES $STAGING_DIR/stagein_link $TMP 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
      exit 1
    fi
  fi
  if [ -s "$STAGING_DIR/stagein_share.1" ] || [ -s "$STAGING_DIR/stagein_share" ]; then
#    safe_init_tmpdir $TMP || exit $?
#    pdbash node all $SCRP_DIR/src/stage_in_init_stgdir_node.sh $TMP share || exit $?
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_in_cp_node.sh $NNODES $STAGING_DIR/stagein_share $TMP share $SCP_THREAD 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
      exit 1
    fi
  fi
  if [ -s "$STAGING_DIR/stagein_local.1" ] || [ -s "$STAGING_DIR/stagein_local" ]; then
    pdbash node all $SCRP_DIR/src/stage_in_init_stgdir_node.sh $TMPL local || exit $?
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_in_cp_node.sh $NNODES $STAGING_DIR/stagein_local $TMPL local $SCP_THREAD 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
      exit 1
    fi
  fi

  if ((DISK_MODE == 1)); then
    if ((ONLINE_STGOUT == 1)); then
      time=$STIME
      loop=0
      while ((time <= ETIME)); do
        loop=$((loop+1))
        if [ -s "$STAGING_DIR/stageout_link_loop${loop}.1" ] || [ -s "$STAGING_DIR/stageout_link_loop${loop}" ]; then
          errmsg=$(bash $SCRP_DIR/src/stage_out_ln.sh $NNODES $STAGING_DIR/stageout_link_loop${loop} $TMP 2>&1)
          if [ -n "$errmsg" ]; then
            echo "$errmsg" >&2
            exit 1
          fi
        fi
        time=$(datetime $time $LCYCLE s)
      done
    else
      if [ -s "$STAGING_DIR/stageout_link.1" ] || [ -s "$STAGING_DIR/stageout_link" ]; then
        errmsg=$(bash $SCRP_DIR/src/stage_out_ln.sh $NNODES $STAGING_DIR/stageout_link $TMP 2>&1)
        if [ -n "$errmsg" ]; then
          echo "$errmsg" >&2
          exit 1
        fi
      fi
    fi
  fi
fi

echo "[$(datetime_now)] ### 5" >&2

#===============================================================================
# Run initialization scripts on all nodes

#if ((TMPRUN_MODE <= 2)); then
#  pdbash node one $SCRP_DIR/src/init_all_node.sh $myname1 || exit $?
#else
#  pdbash node all $SCRP_DIR/src/init_all_node.sh $myname1 || exit $?
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

  if [ -s "$STAGING_DIR/stageout_share_loop${ILOOP}.1" ] || [ -s "$STAGING_DIR/stageout_share_loop${ILOOP}" ]; then
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES $STAGING_DIR/stageout_share_loop${ILOOP} $TMP $SCP_THREAD 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
#      exit 1
    fi
  fi
  if [ -s "$STAGING_DIR/stageout_local_loop${ILOOP}.1" ] || [ -s "$STAGING_DIR/stageout_local_loop${ILOOP}" ]; then
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES $STAGING_DIR/stageout_local_loop${ILOOP} $TMPL $SCP_THREAD 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
#      exit 1
    fi
  fi

  echo "[$(datetime_now)] ${ITIME}: Stage-out (background job completed)" >&2
  rm -f lock.$ILOOP
}

#-------------------------------------------------------------------------------

if [ "$STG_TYPE" = 'builtin' ]; then
  cd $TMP_EXE
fi

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
      if ((IO_ARB == 1)); then
        if ((s == 3)); then
          nodestr='set1.proc'
        elif ((s == 5)); then
          nodestr='set2.proc'
        fi
      fi

#      if ((s <= 3)); then
#        stdout_dir="$TMPOUT/${time}/log/$(basename ${stepexecdir[$s]})"
#      else
#        stdout_dir="$TMPOUT/${atime}/log/$(basename ${stepexecdir[$s]})"
#      fi
      if ((s <= 3)); then
        conf_time=$time
      else
        conf_time=$atime
      fi

#echo "$stdout_dir" >&2
#echo ${stepexecdir[$s]} >&2
#echo $(rev_path ${stepexecdir[$s]}) >&2

#      NNP=$(cat ${NODEFILE_DIR}/${nodestr} | wc -l)
##      mpiexec -n $NNP -vcoordfile "${NODEFILE_DIR}/${nodestr}" -of-proc log/${stepexecname[$s]}.NOUT_${conf_time} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf || exit $?
#echo "      mpiexec -n $NNP -vcoordfile \"${NODEFILE_DIR}/${nodestr}\" -of-proc log/${stepexecname[$s]}.NOUT_${conf_time} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf || exit \$?"

      HOSTLIST=$(cat ${NODEFILE_DIR}/${nodestr})
      HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')
      $MPIRUN $HOSTLIST 1 ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf log/${stepexecname[$s]}.NOUT_${conf_time} || exit $?
#      echo "$MPIRUN $HOSTLIST 1 ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf log/${stepexecname[$s]}.NOUT_${conf_time} || exit \$?"

#      if [ "$STG_TYPE" = 'K_rankdir' ]; then

#        mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" ${stepexecdir[$s]} \
#                "$(rev_path ${stepexecdir[$s]})/${myname1}_step.sh" "$time" "$loop" || exit $?
#      else

#        if ((IO_ARB == 1)); then ##                                 
#          mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" . \
#                  "$SCRP_DIR/${myname1}_step.sh" "$time" "$loop" || exit $? &
#        else ##
#          mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" . \
#                  "$SCRP_DIR/${myname1}_step.sh" "$time" "$loop" || exit $?
#        fi ##
#      fi

    fi
  done

######  if ((IO_ARB == 1)); then ##                                 
######    wait                   ##
######  fi                       ##

#-------------------------------------------------------------------------------
# Online stage out

  if [ "$STG_TYPE" = 'builtin' ]; then
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

if [ "$STG_TYPE" = 'builtin' ]; then
  if ((ONLINE_STGOUT == 1)); then
    wait
  else
    echo "[$(datetime_now)] Finalization (stage-out)" >&2

    if [ -s "$STAGING_DIR/stageout_share.1" ] || [ -s "$STAGING_DIR/stageout_share" ]; then
      errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES $STAGING_DIR/stageout_share $TMP $SCP_THREAD 2>&1)
      if [ -n "$errmsg" ]; then
        echo "$errmsg" >&2
#        exit 1
      fi
    fi
    if [ -s "$STAGING_DIR/stageout_local.1" ] || [ -s "$STAGING_DIR/stageout_local" ]; then
      errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES $STAGING_DIR/stageout_local $TMPL $SCP_THREAD 2>&1)
      if [ -n "$errmsg" ]; then
        echo "$errmsg" >&2
#        exit 1
      fi
    fi
  fi

  if ((CLEAR_TMP == 1 && USE_TMPL == 1)); then
    pdbash node all $SCRP_DIR/src/stage_out_rm_stgdir_node.sh $TMPL local # || exit $?
  fi
#  if ((CLEAR_TMP == 1)); then
#    safe_rm_tmpdir $TMP || exit $?
#  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
