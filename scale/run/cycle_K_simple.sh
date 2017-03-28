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
#myname="$(basename "$0")"
job='cycle'

#===============================================================================
# Configuration

. config.main || exit $?
. config.$job || exit $?

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_${job}.sh
. src/func_${job}_simple.sh

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"

setting "$1" "$2" "$3" "$4" "$5"

echo
print_setting || exit $?
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
safe_init_tmpdir $TMP || exit $?

#===============================================================================
# Determine the distibution schemes

echo "[$(datetime_now)] Determine the distibution schemes"

declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

safe_init_tmpdir $NODEFILE_DIR
if ((IO_ARB == 1)); then                              ##
  distribute_da_cycle_set - $NODEFILE_DIR || exit $?  ##
else                                                  ##
  distribute_da_cycle - $NODEFILE_DIR
fi                                                    ##

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

safe_init_tmpdir $STAGING_DIR || exit $?
staging_list_simple || exit $?

#===============================================================================
# Generate configuration files

config_file_list || exit $?

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

stage_in server || exit $?

#-------------------------------------------------------------------------------

cp -L -r $SCRP_DIR/config.main $TMP/config.main
cp -L -r $SCRP_DIR/config.rc $TMP/config.rc
cp -L -r $SCRP_DIR/config.${job} $TMP/config.${job}
cp -L -r $SCRP_DIR/${job}_simple.sh $TMP/${job}_simple.sh
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

jobscrp="$TMP/${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

rscgrp="micro"

cat > $jobscrp << EOF
#!/bin/sh
#PJM -N ${job}_${SYSNAME}
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

./${job}_simple.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" || exit \$?
EOF

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

job_submit_PJM $jobscrp
echo

job_end_check_PJM_K $jobid
res=$?

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_${SYSNAME} 'o e i s' i

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMP
fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res
