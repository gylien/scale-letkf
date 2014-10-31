#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a K-computer job script and run it.
#
#  October 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_K.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP]
#
#===============================================================================

cd "$(dirname "$0")"

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

setting

#===============================================================================

#mkdir -p $LOGDIR

#if ((MACHINE_TYPE == 10 && PREP == 0)); then
#  mkdir -p $TMP/runlog
#  sleep 0.01s
##  exec 2>> $TMP/runlog/${myname1}.err
#else
#  sleep 0.01s
#####exec 2>> $LOGDIR/${myname1}.err
#fi

#===============================================================================
# Determine the distibution schemes

safe_init_tmpdir $TMPS

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m

NODEFILE_DIR="$TMPS/node"
safe_init_tmpdir $NODEFILE_DIR
distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR

#===============================================================================

STAGING_DIR="$TMPS/staging"

staging_list



cp $SCRP_DIR/config.all $TMPS
echo "SCRP_DIR='./runscp'" >> $TMPS/config.all

USE_RANKDIR=0
if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
  USE_RANKDIR=1
  echo "USE_RANKDIR=1" >> $TMPS/config.all
  echo "TMP='..'" >> $TMPS/config.all
  echo "TMPL='.'" >> $TMPS/config.all
  echo "SCRP_DIR='.'" >> $TMPS/config.all
else
  USE_RANKDIR=0
  echo "USE_RANKDIR=0" >> $TMPS/config.all
  echo "TMP=/work/\$PJM_JOBDIR" >> $TMPS/config.all
  echo "SCRP_DIR='.'" >> $TMPS/config.all
#  echo "SCRP_DIR=\$TMP" >> $TMPS/config.all
fi

echo "NODEFILE_DIR=\"\$TMP/node\"" >> $TMPS/config.all

if ((TMPDAT_MODE <= 2)); then
  echo "TMPDAT=\"\$TMP/dat\"" >> $TMPS/config.all
else
  echo "TMPDAT=\"\$TMPL/dat\"" >> $TMPS/config.all
fi
if ((TMPRUN_MODE <= 2)); then
  echo "TMPRUN=\"\$TMP/run\"" >> $TMPS/config.all
else
  echo "TMPRUN=\"\$TMPL/run\"" >> $TMPS/config.all
fi
if ((TMPOUT_MODE <= 2)); then
  echo "TMPOUT=\"\$TMP/out\"" >> $TMPS/config.all
else
  echo "TMPOUT=\"\$TMPL/out\"" >> $TMPS/config.all
fi



#### walltime limit as a variable!!!
#### rscgrp automatically determined!!!
#### OMP_NUM_THREADS, PARALLEL as a variable!!!
#### ./runscp ./runlog as a variable

cat > fcst_job.sh << EOF
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

if ((USE_RANKDIR == 1)); then
  echo "#PJM --mpi \"use-rankdir\"" >> fcst_job.sh
fi

bash $SCRP_DIR/src/stage_K.sh $STAGING_DIR >> fcst_job.sh

#########################
cat >> fcst_job.sh << EOF
#PJM --stgout "./* /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/"
#PJM --stgout-dir "./node /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/node"
#PJM --stgout-dir "./dat /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/dat"
#PJM --stgout-dir "./run /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/run"
#PJM --stgout-dir "./out /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/out"
#PJM --stgout-dir "./runscp /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/runscp"
#PJM --stgout-dir "./runlog /volume63/data/ra000015/gylien/scale-letkf/scale/run/tmp/runlog"
EOF
#########################

cat >> fcst_job.sh << EOF
#PJM -s
. /work/system/Env_base
export OMP_NUM_THREADS=1
export PARALLEL=1

ls -l .

cd runscp
./fcst.sh

ls -l .

EOF


pjstgchk fcst_job.sh
(($? != 0)) && exit $?

## submit job

## wait for job to finish

#===============================================================================

#safe_rm_tmpdir $TMPS

#===============================================================================

exit 0
