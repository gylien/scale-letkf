#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#
#  August  2014, modified from GFS-LETKF, Guo-Yuan Lien
#  October 2014, modified,                Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP]
#
#  Use settings:
#    config.all
#    config.fcst
#    scale_pp_topo.conf
#    scale_pp_landuse.conf
#    scale_init.conf
#    scale.conf
#
#===============================================================================

cd "$(dirname "$0")"
myname=$(basename "$0")
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.all
(($? != 0)) && exit $?
. config.$myname1
(($? != 0)) && exit $?

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}"

builtin_staging=$((MACHINE_TYPE != 10 && MACHINE_TYPE != 11))

#-------------------------------------------------------------------------------

mkdir -p $LOGDIR
#exec 2>> $LOGDIR/${myname1}.err
exec 2> >(tee -a $LOGDIR/${myname1}.err >&2)

echo "[$(datetime_now)] Start $myname $@" >&2

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN THREADS \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT OUT_OPT LOG_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '                      %-10s = %s\n' $vname "${!vname}" >&2
done

#-------------------------------------------------------------------------------

if ((builtin_staging)); then
  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
    safe_init_tmpdir $TMP
  fi
  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
    safe_init_tmpdir $TMPL
  fi
fi

#===============================================================================
# Determine the distibution schemes

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m

if ((builtin_staging)); then
  safe_init_tmpdir $NODEFILE_DIR
  distribute_fcst "$MEMBERS" $CYCLE machinefile $NODEFILE_DIR
else
  distribute_fcst "$MEMBERS" $CYCLE - -
fi

#===============================================================================
# Determine the staging list and then stage in

if ((builtin_staging)); then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR
  staging_list
  if ((TMPDAT_MODE >= 2 || TMPOUT_MODE >= 2)); then
    pdbash node all $SCRP_DIR/src/stage_in.sh
  fi
fi

#===============================================================================
# Run cycles of forecasts

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

#  exec > $LOGDIR/${myname1}_${stimes[1]}.log
  exec > >(tee $LOGDIR/${myname1}_${stimes[1]}.log)

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
  echo "  Nodes used:               $NNODES"
#  if ((MTYPE == 1)); then
    for n in $(seq $NNODES); do
      echo "    ${node[$n]}"
    done
#  fi
  echo
  echo "  Processes per node:       $PPN"
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
  echo "===================================================================="

#-------------------------------------------------------------------------------
# Call functions to run the job

  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

      echo "[$(datetime_now)] ${stimes[1]}: ${stepname[$s]}" >&2
      echo
      printf " %2d. %-55s\n" $s "${stepname[$s]}"

      ${stepfunc[$s]}

      echo
      echo "===================================================================="

    fi
  done

#-------------------------------------------------------------------------------
# Online stage out

  if ((ONLINE_STGOUT == 1)); then
    if ((MACHINE_TYPE == 11)); then
      touch $TMP/loop.${loop}.done
    fi
    if ((builtin_staging && $(datetime $time $((lcycles * CYCLE)) s) <= ETIME)); then
      if ((MACHINE_TYPE == 12)); then
        echo "[$(datetime_now)] ${stimes[1]}: Online stage out"
        bash $SCRP_DIR/src/stage_out.sh s $loop
        pdbash node all $SCRP_DIR/src/stage_out.sh $loop
      else
        echo "[$(datetime_now)] ${stimes[1]}: Online stage out (background job)"
        ( bash $SCRP_DIR/src/stage_out.sh s $loop ;
          pdbash node all $SCRP_DIR/src/stage_out.sh $loop ) &
      fi
    fi
  fi

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo
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

if ((builtin_staging)); then
  echo "[$(datetime_now)] Finalization (stage out)" >&2

  if ((TMPOUT_MODE >= 2)); then
    if ((ONLINE_STGOUT == 1)); then
      wait
      bash $SCRP_DIR/src/stage_out.sh s $loop
      pdbash node all $SCRP_DIR/src/stage_out.sh $loop
    else
      bash $SCRP_DIR/src/stage_out.sh s
      pdbash node all $SCRP_DIR/src/stage_out.sh
    fi
  fi

  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
    safe_rm_tmpdir $TMP
  fi
  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
    safe_rm_tmpdir $TMPL
  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
