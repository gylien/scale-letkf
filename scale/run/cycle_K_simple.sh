#!/bin/bash
#===============================================================================
#
#  Wrap cycle.sh in a K-computer job script (micro) and run it.
#
#  February 2015, created,                Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_K_micro.sh [STIME ETIME ISTEP FSTEP TIME_LIMIT]
#
#===============================================================================

cd "$(dirname "$0")"
myname1='cycle'

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
. src/func_${myname1}_simple.sh

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"

setting "$1" "$2" "$3" "$4" "$5"

echo
print_setting
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"

if [ ${TMP:0:8} != '/scratch' ]; then
  echo "[Error] $0: When using 'micro' resource group, \$TMP will be completely removed." >&2
  echo "        Wrong setting detected:" >&2
  echo "        \$TMP = '$TMP'" >&2
  exit 1
fi
safe_init_tmpdir $TMP

#===============================================================================
# Determine the distibution schemes

echo "[$(datetime_now)] Determine the distibution schemes"

declare -a procs
declare -a mem2node
declare -a node
declare -a name_m
declare -a node_m

safe_init_tmpdir $NODEFILE_DIR
if ((IO_ARB == 1)); then                   ##
  distribute_da_cycle_set - $NODEFILE_DIR  ##
else                                       ##
  distribute_da_cycle - $NODEFILE_DIR
fi                                         ##

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

safe_init_tmpdir $STAGING_DIR
staging_list_simple

#===============================================================================
# Generate configuration files

######safe_init_tmpdir $CONFIG_DIR
config_file_list

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

bash $SCRP_DIR/src/stage_in.sh a

#-------------------------------------------------------------------------------

cp -L -r $SCRP_DIR/config.main $TMP/config.main
cp -L -r $SCRP_DIR/config.rc $TMP/config.rc
cp -L -r $SCRP_DIR/config.${myname1} $TMP/config.${myname1}
cp -L -r $SCRP_DIR/${myname1}_simple.sh $TMP/${myname1}_simple.sh
mkdir -p $TMP/src
cp -L -r $SCRP_DIR/src/* $TMP/src

echo "SCRP_DIR=\"$TMP\"" >> $TMP/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMP/config.main

echo "RUN_LEVEL='K_micro'" >> $TMP/config.main

if ((IO_ARB == 1)); then                                        ##
  NNODES=$((NNODES*2))                                          ##
  NNODES_APPAR=$((NNODES_APPAR*2))                              ##
fi                                                              ##

#===============================================================================
# Creat a job script

jobscrp="$TMP/${myname1}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

rscgrp="micro"

cat > $jobscrp << EOF
#!/bin/sh
#PJM -N ${myname1}_${SYSNAME}
#PJM -s
#PJM --rsc-list "node=${NNODES}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --rsc-list "rscgrp=${rscgrp}"
##PJM --mpi "shape=${NNODES}"
#PJM --mpi "proc=${totalnp}"
#PJM --mpi assign-online-node

. /work/system/Env_base_1.2.0-20-1
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./${myname1}_simple.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" || exit \$?
EOF




exit



#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${myname1} job on PJM"
echo

job_submit_PJM $jobscrp
echo

job_end_check_PJM_K $jobid
res=$?

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

bash $SCRP_DIR/src/stage_out.sh a

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

mkdir -p $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.main $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.${myname1} $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.nml.* $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $TMP/${myname1}_job.sh $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $TMP/${myname1}_${SYSNAME}.o${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.o
cp -f $TMP/${myname1}_${SYSNAME}.e${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.e
cp -f $TMP/${myname1}_${SYSNAME}.i${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.i
( cd $SCRP_DIR ; git log -1 --format="SCALE-LETKF version %h (%ai)" > $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )
( cd $MODELDIR ; git log -1 --format="SCALE       version %h (%ai)" >> $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )

finalization

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMP
fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res