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
#    config.main
#    config.fcst
#    config.nml.scale_pp_topo
#    config.nml.scale_pp_landuse
#    config.nml.scale_init
#    config.nml.scale
#
#===============================================================================

cd "$(dirname "$0")"
myname='fcst.sh'
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.main
res=$? && ((res != 0)) && exit $res
. config.$myname1
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

if ((USE_RANKDIR == 1)); then
  SCRP_DIR="."
  if ((TMPDAT_MODE <= 2)); then
    TMPDAT="../dat"
  else
    TMPDAT="./dat"
  fi
  if ((TMPRUN_MODE <= 2)); then
    TMPRUN="../run"
  else
    TMPRUN="./run"
  fi
fi

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9"

#-------------------------------------------------------------------------------

mkdir -p $LOGDIR
exec 3>&1 4>&2
#exec 2>> $LOGDIR/${myname1}.err
exec 2> >(tee -a $LOGDIR/${myname1}.err >&2)

echo "[$(datetime_now)] Start $myname $@" >&2

for vname in DIR OUTDIR DATA_TOPO DATA_LANDUSE DATA_BDY DATA_BDY_WRF OBS OBSNCEP MEMBER NNODES PPN THREADS \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT OUT_OPT LOG_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '                      %-10s = %s\n' $vname "${!vname}" >&2
done

#-------------------------------------------------------------------------------

if ((BUILTIN_STAGING && ISTEP == 1)); then
  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
    safe_init_tmpdir $TMP
  fi
  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
    safe_init_tmpdir $TMPL
  fi
fi

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

#if ((BUILTIN_STAGING && ISTEP == 1)); then
if ((BUILTIN_STAGING)); then
  safe_init_tmpdir $NODEFILE_DIR
  distribute_fcst "$MEMBERS" $CYCLE machinefile $NODEFILE_DIR
else
  distribute_fcst "$MEMBERS" $CYCLE - - $NODEFILE_DIR/distr
fi

#===============================================================================
# Determine the staging list and then stage in

if ((BUILTIN_STAGING)); then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR
  staging_list
  if ((TMPDAT_MODE >= 2 || TMPOUT_MODE >= 2)); then
    pdbash node all $SCRP_DIR/src/stage_in.sh
  fi
fi

#===============================================================================
# Run initialization scripts on all nodes

if ((TMPRUN_MODE <= 2)); then
  pdbash node one $SCRP_DIR/src/init_all_node.sh $myname1
else
  pdbash node all $SCRP_DIR/src/init_all_node.sh $myname1
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
      echo

      ######
      if ((s == 1)); then
        if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
          echo "  ... skip this step (use prepared topo and landuse files)"
          echo
          echo "===================================================================="
          continue
        elif ((BDY_FORMAT == 0)); then
          echo "  ... skip this step (use prepared boundaries)"
          echo
          echo "===================================================================="
          continue
        fi
      fi
      ######
      if ((s == 2)); then
        if ((BDY_FORMAT == 0)); then
          echo "  ... skip this step (use prepared boundaries)"
          continue
        fi
      fi
      ######

      enable_iter=0
      if ((s == 2 && BDY_ENS == 1)); then
        enable_iter=1
      elif ((s == 3)); then
        enable_iter=1
      fi

      if ((enable_iter == 1)); then
        for it in $(seq $nitmax); do
          if ((USE_RANKDIR == 1)); then
            mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf ${stepexecdir[$s]} \
                    "$(rev_path ${stepexecdir[$s]})/fcst_step.sh" $loop $it # > /dev/null
          else
            mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf . \
                    "$SCRP_DIR/fcst_step.sh" $loop $it # > /dev/null
          fi
        done
      else
        if ((USE_RANKDIR == 1)); then
          mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf ${stepexecdir[$s]} \
                  "$(rev_path ${stepexecdir[$s]})/fcst_step.sh" $loop # > /dev/null
        else
          mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf . \
                  "$SCRP_DIR/fcst_step.sh" $loop # > /dev/null
        fi
      fi


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
    if ((BUILTIN_STAGING && $(datetime $time $((lcycles * CYCLE)) s) <= ETIME)); then
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

  exec 1>&3

#-------------------------------------------------------------------------------

  time=$(datetime $time $((lcycles * CYCLE)) s)
  s_flag=0

#-------------------------------------------------------------------------------
done
#-------------------------------------------------------------------------------

#===============================================================================
# Stage out

if ((BUILTIN_STAGING)); then
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

#  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
#    safe_rm_tmpdir $TMP
#  fi
#  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
#    safe_rm_tmpdir $TMPL
#  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
