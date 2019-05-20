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
#    fcst_K.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='fcst'

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

STAGING_DIR="$TMPSL/staging"
NODEFILE_DIR="$TMPS/node"

#-------------------------------------------------------------------------------

if ((USE_TMP_LINK == 1)); then
  echo "[Error] $0: Wrong disk mode for K computer staged jobs." >&2
  exit 1
fi

if [ "$PRESET" = 'K_rankdir' ] && ((PNETCDF == 1)); then
  echo "[Error] When PNETCDF is enabled, 'K_rankdir' preset cannot be used." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@"

setting "$@" || exit $?

if [ "$CONF_MODE" = 'static' ]; then
  . src/func_common_static.sh || exit $?
  . src/func_${job}_static.sh || exit $?
fi

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

safe_init_tmpdir $NODEFILE_DIR || exit $?
distribute_fcst "$MEMBERS" $CYCLE "$NODELIST_TYPE" $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cp -L $SCRP_DIR/config.main $TMPS/config.main

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMPS/config.main
echo "NODEFILE_DIR=\"\$TMPROOT/node\"" >> $TMPS/config.main
echo "RUN_LEVEL=4" >> $TMPS/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMPS/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
if [ "$CONF_MODE" = 'static' ]; then
  staging_list_static || exit $?
  config_file_list $TMPS/config || exit $?
else
  staging_list || exit $?
fi

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${TMPS}/config.main|config.main
${SCRP_DIR}/config.rc|config.rc
${SCRP_DIR}/config.${job}|config.${job}
${SCRP_DIR}/${job}.sh|${job}.sh
${SCRP_DIR}/src/|src/
${NODEFILE_DIR}/|node/
EOF

if [ "$CONF_MODE" != 'static' ]; then
  echo "${SCRP_DIR}/${job}_step.sh|${job}_step.sh" >> ${STAGING_DIR}/${STGINLIST}
fi

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

. /work/system/Env_base_1.2.0-25
export LD_LIBRARY_PATH=/opt/klocal/zlib-1.2.11-gnu/lib:\$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
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

###if [ "$CONF_MODE" = 'static' ]; then
###  config_file_save $TMPS/config || exit $?
###fi

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMPS
  safe_rm_tmpdir $TMPSL
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
