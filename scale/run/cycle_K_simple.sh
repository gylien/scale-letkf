#!/bin/bash
#===============================================================================
#
#  Wrap cycle.sh in a K-computer job script and run it.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_K_simple.sh [STIME ETIME ISTEP FSTEP TIME_LIMIT]
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

STAGING_DIR="$TMPS/staging"
NODEFILE_DIR="$TMPS/node"

#-------------------------------------------------------------------------------

if ((USE_TMP_LINK == 1)); then
  echo "[Error] $0: Wrong disk mode for K computer regular jobs." >&2
  exit 1
fi

if [ "$PRESET" = 'K_rankdir' ] && ((PNETCDF == 1)); then
  echo "[Error] When PNETCDF is enabled, 'K_rankdir' preset cannot be used." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"

setting "$1" "$2" "$3" "$4" "$5"

echo
print_setting || exit $?
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"

safe_init_tmpdir $TMPS || exit $?

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

safe_init_tmpdir $TMPS/node
if ((IO_ARB == 1)); then                           ##
  distribute_da_cycle_set - $TMPS/node || exit $?  ##
else                                               ##
  distribute_da_cycle - $TMPS/node || exit $?
fi                                                 ##

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cp -L $SCRP_DIR/config.main $TMPS/config.main

if [ "$PRESET" = 'K_rankdir' ]; then
  echo "TMP='..'" >> $TMPS/config.main
  echo "TMPL='.'" >> $TMPS/config.main
else
  echo "TMP='.'" >> $TMPS/config.main
fi
if ((DISK_MODE == 2)); then
  echo "SCRP_DIR=\"\$TMP\"" >> $TMPS/config.main
elif ((DISK_MODE == 3)); then
  echo "SCRP_DIR=\"\$TMPL\"" >> $TMPS/config.main
fi

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMPS/config.main

echo "RUN_LEVEL=2" >> $TMPS/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
staging_list_simple || exit $?

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${TMPS}/config.main|config.main
${SCRP_DIR}/config.rc|config.rc
${SCRP_DIR}/config.${job}|config.${job}
${SCRP_DIR}/${job}_simple.sh|${job}_simple.sh
${SCRP_DIR}/src/|src/
${TMPS}/node/|node/
EOF

#===============================================================================
# Generate configuration files

config_file_list $TMPS/config || exit $?

#===============================================================================

if ((IO_ARB == 1)); then                                        ##
  NNODES=$((NNODES*2))                                          ##
  NNODES_APPAR=$((NNODES_APPAR*2))                              ##
fi                                                              ##

#===============================================================================
# Creat a job script

jobscrp="${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

if ((NNODES > 36864)); then
  rscgrp="huge"
elif ((NNODES > 384)); then
  rscgrp="large"
else
  rscgrp="small"
fi

cat > $jobscrp << EOF
#!/bin/sh
#PJM -N ${job}_${SYSNAME}
#PJM -s
#PJM --rsc-list "node=${NNODES}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --rsc-list "rscgrp=${rscgrp}"
##PJM --rsc-list "node-quota=29G"
##PJM --mpi "shape=${NNODES}"
#PJM --mpi "proc=${totalnp}"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all
EOF

if [ "$PRESET" = 'K_rankdir' ]; then
  echo "#PJM --mpi \"use-rankdir\"" >> $jobscrp
  stage_K_inout 1
else
  stage_K_inout 0
fi

cat >> $jobscrp << EOF

. /work/system/Env_base_1.2.0-20-1
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" || exit \$?
EOF

#===============================================================================
# Check the staging list

echo "[$(datetime_now)] Run pjstgchk"
echo

pjstgchk $jobscrp || exit $?
echo

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

job_submit_PJM $jobscrp
echo

job_end_check_PJM_K $jobid
res=$?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $SCRP_DIR $jobid ${job}_${SYSNAME} 'o e i s' i

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMPS
fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res
