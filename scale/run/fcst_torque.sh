#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#  August  2014, modified from GFS-LETKF, Guo-Yuan Lien
#  October 2014, modified,                Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP PREP]
#    (see variable explanation in 'config.fcst')
#
#  Use settings:
#    config.all
#    config.fcst (optional)
#
#===============================================================================

cd "$(dirname "$0")"
myname=$(basename "$0")
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.all
(($? != 0)) && exit $?

if [ -f config.$myname1 ]; then
  . config.$myname1
fi

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh

. src/func_$myname1.sh

#-------------------------------------------------------------------------------

if ((MACHINE_TYPE == 10)); then
  PREP=${PREP:-1}
else
  PREP=${PREP:-0}
fi

#-------------------------------------------------------------------------------

mkdir -p $LOGDIR

if ((MACHINE_TYPE == 10 && PREP == 0)); then
  mkdir -p $TMP/runlog
  sleep 0.01s
#  exec 2>> $TMP/runlog/${myname1}.err
else
  sleep 0.01s
####exec 2>> $LOGDIR/${myname1}.err
fi

echo "[$(datetime_now)] Start $myname" >&2

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT FOUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '                      %-10s = %s\n' $vname ${!vname} >&2
done

#===============================================================================
# More configuration

if ((MACHINE_TYPE == 10)); then
  if ((PREP == 1)); then
    safe_init_tmpdir $TMPS
  fi
else
  safe_init_tmpdir $TMP
fi

#-------------------------------------------------------------------------------

# more variable difinition???



#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m


if ((MACHINE_TYPE == 10)); then

  if ((PREP == 1)); then
    NODEFILE_DIR="$TMPS/node"
    safe_init_tmpdir $NODEFILE_DIR
    distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR
  else
    distribute_fcst "$MEMBERS" $CYCLE - -
  fi

else

  safe_init_tmpdir $NODEFILE_DIR
  distribute_fcst "$MEMBERS" $CYCLE machinefile $NODEFILE_DIR

fi

#===============================================================================

if ((MACHINE_TYPE == 10)); then

  if ((PREP == 1)); then

    STAGING_DIR="$TMPS/staging"

    init


    cp $SCRP_DIR/config.all $TMPS
    echo "SCRP_DIR='./runscp'" >> $TMPS/config.all

    cp $SCRP_DIR/config.$myname1 $TMPS
    echo "PREP=0" >> $TMPS/config.$myname1

#### walltime limit as a variable!!!
#### rscgrp automatically determined!!!
#### OMP_NUM_THREADS, PARALLEL as a variable!!!
#### ./runscp ./runlog as a variable

    cat > ${myname1}_pj.sh << EOF
#!/bin/sh
##PJM -g ra000015
#PJM --rsc-list "node=$NNODES"
#PJM --rsc-list "elapse=00:01:00"
#PJM --rsc-list "rscgrp=small"
##PJM --rsc-list "node-quota=29GB"
#PJM --mpi "shape=$NNODES"
#PJM --mpi "proc=$((NNODES*PPN))"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all
EOF

if [ ! -z "$TMPL" ]; then
  echo "#PJM --mpi \"use-rankdir\"" >> ${myname1}_pj.sh
fi

    bash $SCRP_DIR/src/stage_K.sh $STAGING_DIR >> ${myname1}_pj.sh

#########################
    cat >> ${myname1}_pj.sh << EOF
#PJM --stgout "./* /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/"
#PJM --stgout-dir "./node /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/node"
#PJM --stgout-dir "./dat /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/dat"
#PJM --stgout-dir "./run /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/run"
#PJM --stgout-dir "./out /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/out"
#PJM --stgout-dir "./runscp /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/runscp"
#PJM --stgout-dir "./runlog /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/runlog"
EOF
#########################

    cat >> ${myname1}_pj.sh << EOF
#PJM -s
. /work/system/Env_base
export OMP_NUM_THREADS=1
export PARALLEL=1

ls -l .

cd runscp
./${myname}

ls -l .

EOF


    pjstgchk ${myname1}_pj.sh
    (($? != 0)) && exit $?

    ## submit job

    ## wait for job to finish
    
  fi  

else
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  init

####
#safe_init_tmpdir $TMPDAT
#safe_init_tmpdir $TMPOUT
#safe_init_tmpdir $TMPRUN
####

  pdbash node all $SCRP_DIR/src/stage_in.sh
fi

#===============================================================================

if ((PREP == 0)); then

#===============================================================================
# Run cycles of forecasts


declare -a stimes
declare -a stimesfmt
lcycles=$((LCYCLE * CYCLE_SKIP))
s_flag=1
e_flag=0
time=$STIME

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

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

if ((MACHINE_TYPE == 10)); then
  sleep 0.01s
#  exec > $TMP/runlog/${myname1}_${stimes[1]}.log
else
  sleep 0.01s
####  exec > $LOGDIR/${myname1}_${stimes[1]}.log
fi

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
    for n in `seq $NNODES`; do
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

      echo "[$(datetime_now)] ${stepname[$s]}" >&2
      echo
      printf " %2d. %-55s\n" $s "${stepname[$s]}"

      ${stepfunc[$s]}

      echo
      echo "===================================================================="

    fi
  done

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

fi # ((PREP == 0))

#===============================================================================
# Finalization

#safe_rm_tmpdir $TMP
#safe_rm_tmpdir $TMPS

echo "[$(datetime_now)] Finish fcst.sh $@" >&2

#===============================================================================

exit 0
