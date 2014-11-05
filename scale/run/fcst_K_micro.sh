#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a K-computer job script (micro) and run it.
#
#  November 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_K_micro.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP]
#
#===============================================================================

cd "$(dirname "$0")"

#--------------

TIME_LIMIT='00:30:00'

#===============================================================================
# Configuration

. config.all
(($? != 0)) && exit $?
. config.fcst
(($? != 0)) && exit $?

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_fcst.sh

#-------------------------------------------------------------------------------

if ((TMPDAT_MODE != 2 || TMPRUN_MODE != 2 || TMPOUT_MODE != 2)); then
  echo "[Error] $0: When using 'micro' resource group," >&2
  echo "        \$TMPDAT_MODE, \$TMPRUN_MODE, \$TMPOUT_MODE all need to be 2." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

setting

jobscrp="$TMP/fcst_job.sh"

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"
echo

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT FOUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '  %-10s = %s\n' $vname "${!vname}"
done

echo

#-------------------------------------------------------------------------------

safe_init_tmpdir $TMP

#===============================================================================
# Determine the distibution schemes

# K computer
NNODES_real=$NNODES
PPN_real=$PPN
NNODES=$((NNODES*PPN))
PPN=1

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m

safe_init_tmpdir $NODEFILE_DIR
distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR

#===============================================================================
# Determine the staging list and then stage in

echo "[$(datetime_now)] Initialization (stage in)"

safe_init_tmpdir $STAGING_DIR
staging_list
bash $SCRP_DIR/src/stage_in.sh a

#-------------------------------------------------------------------------------
# stage-in: scripts

cp -L -r $SCRP_DIR/config.all $TMP/config.all
cp -L -r $SCRP_DIR/config.fcst $TMP/config.fcst
cp -L -r $SCRP_DIR/fcst.sh $TMP/fcst.sh
mkdir -p $TMP/src
cp -L -r $SCRP_DIR/src/* $TMP/src

echo "SCRP_DIR=\"$TMP\"" >> $TMP/config.all
echo "LOGDIR=\"$TMP/log\"" >> $TMP/config.all

echo "NNODES=$NNODES" >> $TMP/config.all
echo "PPN=$PPN" >> $TMP/config.all
echo "NNODES_real=$NNODES_real" >> $TMP/config.all
echo "PPN_real=$PPN_real" >> $TMP/config.all

#===============================================================================
# Creat a job script

echo "[$(datetime_now)] Create a job script '$jobscrp'"

rscgrp="micro"

cat > $jobscrp << EOF
#!/bin/sh
##PJM -g ra000015
#PJM -N fcst_${SYSNAME}
#PJM -s
#PJM --rsc-list "node=${NNODES_real}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --mpi "shape=${NNODES_real}"
#PJM --mpi "proc=$NNODES"
#PJM --mpi assign-online-node

. /work/system/Env_base
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./fcst.sh
EOF

echo "[$(datetime_now)] Run fcst job on PJM"
echo

job_submit_PJM $jobscrp
echo

if ((ONLINE_STGOUT != 1)); then

  job_end_check_PJM $jobid

else # when using online stage-out, check the joub status in a special way.

  loop=1
  while (($(pjstat $jobid | sed -n '2p' | awk '{print $10}') >= 1)); do
    if [ -e "$TMP/loop.${loop}.done" ]; then
      echo "[$(datetime_now)] Online stage out: Loop # $loop"
      bash $SCRP_DIR/src/stage_out.sh a $loop &
      loop=$((loop+1))
    fi
    sleep 5s
  done
  wait

  while [ -e "$TMP/loop.${loop}.done" ]; do
    echo "[$(datetime_now)] Online stage out: Loop # $loop"
    bash $SCRP_DIR/src/stage_out.sh a $loop
    loop=$((loop+1))
  done

fi

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

if ((ONLINE_STGOUT != 1)); then
  bash $SCRP_DIR/src/stage_out.sh a
fi

mkdir -p $LOGDIR
cp -f $TMP/log/fcst_*.log $LOGDIR
if [ -f "$TMP/log/fcst.err" ]; then
  cat $TMP/log/fcst.err >> $LOGDIR/fcst.err
fi

safe_rm_tmpdir $TMP

echo "[$(datetime_now)] Finish $(basename $0) $@"

#===============================================================================

exit 0
