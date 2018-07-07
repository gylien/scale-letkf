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
#    cycle_K_micro.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='cycle'

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

#-------------------------------------------------------------------------------

if ((USE_TMP_LINK == 1 || USE_TMPL == 1)); then
  echo "[Error] $0: Wrong disk mode for K computer micro jobs." >&2
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

safe_init_tmpdir $NODEFILE_DIR || exit $?
if ((DTF_MODE >= 1)); then                                           ##
  distribute_da_cycle_set "$NODELIST_TYPE" $NODEFILE_DIR || exit $?  ##
else                                                                 ##
  distribute_da_cycle "$NODELIST_TYPE" $NODEFILE_DIR || exit $?
fi                                                                   ##

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cat $SCRP_DIR/config.main | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

if ((DTF_MODE >= 1)); then
  cp -L $LIBDTF_PATH/libdtf.so $TMP/libdtf.so
#  cp -L $SCRP_DIR/dtf.ini $TMP/dtf.ini
fi

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMP/config.main
echo "RUN_LEVEL=4" >> $TMP/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMP/config.main

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
${SCRP_DIR}/config.rc|config.rc
${SCRP_DIR}/config.${job}|config.${job}
${SCRP_DIR}/${job}.sh|${job}.sh
${SCRP_DIR}/src/|src/
EOF

if [ "$CONF_MODE" != 'static' ]; then
  echo "${SCRP_DIR}/${job}_step.sh|${job}_step.sh" >> ${STAGING_DIR}/${STGINLIST}
fi

#===============================================================================

if ((DTF_MODE >= 1)); then                                            ##
  NNODES=$((NNODES*2))                                                ##
  NNODES_APPAR=$((NNODES_APPAR*2))                                    ##
fi                                                                    ##

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

stage_in server || exit $?

#===============================================================================
# Creat a job script

jobscrp="$TMP/${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

rscgrp="micro"

if ((SCALE_NP >= 8)); then
  MAX_WORKGROUP_SIZE=$((SCALE_NP / 8))
else
  MAX_WORKGROUP_SIZE=1
fi


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

# clean up LD_LIBRARY_PATH
export LD_LIBRARY_PATH=

. /work/system/Env_base_1.2.0-23
#export LD_LIBRARY_PATH=/opt/klocal/zlib-1.2.11-gnu/lib:\$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE
EOF

if ((DTF_MODE >= 1)); then
  cat >> $jobscrp << EOF

#export DTF_VERBOSE_LEVEL=2 # DEBUG
export DTF_VERBOSE_LEVEL=0
export DTF_SCALE=1
export DTF_INI_FILE=$TMP/dtf.ini
export SCALE_ENSEMBLE_SZ=$((SCALE_NP))
#export DTF_IGNORE_ITER=    #set if necessary
export MAX_WORKGROUP_SIZE=$MAX_WORKGROUP_SIZE

#export LD_LIBRARY_PATH=$TMP:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$TMP:/opt/klocal/zlib-1.2.11-gnu/lib:/opt/aics/netcdf/k-serial-noszip/lib:/opt/FJSVtclang/GM-1.2.0-23/lib64
export LD_LIBRARY_PATH=$TMP:/opt/FJSVtclang/GM-1.2.0-23/lib64
export DTF_GLOBAL_PATH=$TMP
EOF
fi

cat >> $jobscrp << EOF

./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

job_submit_PJM $jobscrp
echo

if ((OBS_USE_JITDT == 1)); then
  rm -f ${SCRP_DIR}/jitdt_send_testdata_on_start.log
  ${SCRP_DIR}/misc/jitdt_send_testdata_on_start.sh > ${SCRP_DIR}/jitdt_send_testdata_on_start.log 2>&1 &
fi

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

backup_exp_setting $job $TMP $jobid ${job}_${SYSNAME} 'o e i' i

if [ "$CONF_MODE" = 'static' ]; then
  config_file_save $TMPS/config || exit $?
fi

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMP
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
