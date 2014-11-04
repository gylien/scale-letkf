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

echo

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT FOUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '  %-10s = %s\n' $vname "${!vname}"
done

echo
echo "Create a job script '$jobscrp'..."
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

echo "[$(datetime_now)] Initialization (stage in)" >&2

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

#===============================================================================




cp $SCRP_DIR/config.all $TMPS


echo "SCRP_DIR=\"\$TMP\"" >> $TMP/config.all
#echo "NODEFILE_DIR=\"\$(pwd)/node\"" >> $TMPS/config.all
echo "LOGDIR=\"\$TMP/log\"" >> $TMP/config.all

#echo "NNODES=$NNODES" >> $TMPS/config.all
#echo "PPN=$PPN" >> $TMPS/config.all
#echo "NNODES_real=$NNODES_real" >> $TMPS/config.all
#echo "PPN_real=$PPN_real" >> $TMPS/config.all

rscgrp="micro"

#------

cat > $jobscrp << EOF
#!/bin/sh
##PJM -g ra000015
#PJM -N fcst_${SYSNAME}
#PJM -s
##PJM -j
#PJM --rsc-list "node=${NNODES_real}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --rsc-list "rscgrp=${rscgrp}"
##PJM --rsc-list "node-quota=29GB"
#PJM --mpi "shape=${NNODES_real}"
#PJM --mpi "proc=$NNODES"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all

. /work/system/Env_base
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./fcst.sh

EOF

echo

## submit job

## wait for job to finish

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)" >&2

#bash $SCRP_DIR/src/stage_out.sh a

#mkdir -p $LOGDIR
#cp -f $TMP/log/fcst_*.log $LOGDIR
#if [ -f "$TMP/log/fcst.err" ]; then
#  cat $TMP/log/fcst.err >> $LOGDIR/fcst.err
#fi

#safe_rm_tmpdir $TMP
#safe_rm_tmpdir $TMPS

echo "[$(datetime_now)] Finish fcst.sh $@" >&2

#===============================================================================

exit 0
