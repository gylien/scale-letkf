#!/bin/bash
#===============================================================================
#
#  Wrap cycle.sh in a K-computer job script and run it.
#
#  October 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_K.sh [..]
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
if ((DTF_MODE >= 1)); then                                           ##
  distribute_da_cycle_set "$NODELIST_TYPE" $NODEFILE_DIR || exit $?  ##
else                                                                 ##
  distribute_da_cycle "$NODELIST_TYPE" $NODEFILE_DIR || exit $?
fi                                                                   ##

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

if [ "$CONF_MODE" = 'static' ] && ((DTF_MODE >= 1)); then
  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${NODEFILE_DIR}/|node/
EOF
else
  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${TMPS}/config.main|config.main
${SCRP_DIR}/config.rc|config.rc
${SCRP_DIR}/config.${job}|config.${job}
${SCRP_DIR}/${job}.sh|${job}.sh
${SCRP_DIR}/src/|src/
${NODEFILE_DIR}/|node/
EOF
fi

#if ((DTF_MODE >= 1)); then
#  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#${LIBDTF_PATH}/libdtf.so|libdtf.so
#EOF
#  if ((DTF_MPMD >= 1)); then
#    cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#${SPLIT_WRAP}|libsplitworld.so 
#EOF
#  fi
#fi

if [ "$CONF_MODE" != 'static' ]; then
  echo "${SCRP_DIR}/${job}_step.sh|${job}_step.sh" >> ${STAGING_DIR}/${STGINLIST}
fi

#===============================================================================

if ((DTF_MODE >= 1)); then                                            ##
  NNODES_ORIG=$NNODES                                                 ##
  NNODES=$((NNODES*2))                                                ##
  NNODES_APPAR=$((NNODES_APPAR*2))                                    ##
fi                                                                    ##

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


if ((DTF_MPMD >= 1)); then
  TOTAL_PROC=$((NNODES / PPN))
else
  TOTAL_PROC=${totalnp}
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
#PJM --mpi "proc=${TOTAL_PROC}"
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

# clean up LD_LIBRARY_PATH
export LD_LIBRARY_PATH=.

#. /work/system/Env_base_1.2.0-23
. /work/system/Env_base
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE
EOF

if ((DTF_MODE >= 1)); then
  if ((SCALE_NP >= 8)); then
    MAX_WORKGROUP_SIZE=$((SCALE_NP / 8))
  else
    MAX_WORKGROUP_SIZE=1
  fi

  cat >> $jobscrp << EOF

#export DTF_VERBOSE_LEVEL=2 # DEBUG
export DTF_VERBOSE_LEVEL=0
export DTF_SCALE=1
export DTF_INI_FILE=./dtf.ini
export SCALE_ENSEMBLE_SZ=$((SCALE_NP))
#export DTF_IGNORE_ITER=    #set if necessary
export MAX_WORKGROUP_SIZE=$MAX_WORKGROUP_SIZE

export DTF_GLOBAL_PATH=.
EOF
fi

if  [ "$CONF_MODE" = 'static' ] && ((DTF_MODE >= 1)) && ((DTF_MPMD >= 1)); then
  cat >> $jobscrp << EOF
#export LD_PRELOAD=./libsplitworld.so

# CREATE OF-PROC FILE
mkdir log_1 log_2
lfs setstripe -c 1 log_1
lfs setstripe -c 1 log_2
lfs setstripe -c 1 log

#for i in \`seq 0 $((NNODES_ORIG*PPN-1))\`; do echo log_1/scale-rm_ens.NOUT_${STIME}.\${i} ; done | xargs -n 200 -P 8 touch
#for i in \`seq 0 $((NNODES_ORIG*PPN-1))\`; do echo log_2/letkf.NOUT_${STIME}.\${i} ; done | xargs -n 200 -P 8 touch
for i in \`seq 0 $((NNODES*PPN-1))\`; do echo log/NOUT.\${i} ; done | xargs -n 200 -P 8 touch

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] ${STIME}: Submitted to background: ensemble forecasts & LETKF by using MPMD" >&2
mpiexec -of-proc log/NOUT \\
    -n $((NNODES_ORIG*PPN)) -x MPMD_COMP=0 ./scale-rm_ens scale-rm_ens_${STIME}.conf :\\
    -n $((NNODES_ORIG*PPN)) -x MPMD_COMP=1 ./letkf letkf_${STIME}.conf
#    -n $((NNODES_ORIG*PPN)) -of-proc log_1/scale-rm_ens.NOUT_${STIME} -x DTF_COMP=0 ./scale-rm_ens scale-rm_ens_${STIME}.conf :\\
#    -n $((NNODES_ORIG*PPN)) -of-proc log_2/letkf.NOUT_${STIME} -x DTF_COMP=1 ./letkf letkf_${STIME}.conf
#    -n $((NNODES_ORIG*PPN)) -vcoordfile ${TMPROOT}/node/set1.proc -of-proc log_1/scale-rm_ens.NOUT_${STIME} -x DTF_COMP=0 ./scale-rm_ens scale-rm_ens_${STIME}.conf :\\
#    -n $((NNODES_ORIG*PPN)) -vcoordfile ${TMPROOT}/node/set2.proc -of-proc log_2/letkf.NOUT_${STIME} -x DTF_COMP=1 ./letkf letkf_${STIME}.conf

jobids=\$(jobs -p)
while [ -n "\$jobids" ]; do
  for job in \$jobids; do
    if ! (kill -0 \$job 2> /dev/null); then
      rcode=0
      wait \$job || rcode=\$?
      if ((rcode != 0)); then
        echo "[\$(date +'%Y-%m-%d %H:%M:%S')] ${STIME}: One of the background programs crashed..." >&2
        exit \$rcode
      fi
    fi
  done
  jobids=\$(jobs -p)
  sleep 1s
done

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] Finish $STIME $ETIME" >&2
EOF
fi

if [ "$CONF_MODE" = 'static' ] && ((DTF_MODE >= 1)) && ((DTF_MPMD < 1)); then
  cat >> $jobscrp << EOF

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] Start $STIME $ETIME" >&2

# CREATE OF-PROC FILE
mkdir log_1 log_2
lfs setstripe -c 1 log_1
lfs setstripe -c 1 log_2

for i in \`seq 0 $((NNODES_ORIG*PPN-1))\`; do echo log_1/scale-rm_ens.NOUT_${STIME}.\${i} ; done | xargs -n 200 -P 8 touch
for i in \`seq 0 $((NNODES_ORIG*PPN-1))\`; do echo log_2/letkf.NOUT_${STIME}.\${i} ; done | xargs -n 200 -P 8 touch


mpiexec -n $((NNODES_ORIG*PPN)) -vcoordfile ${TMPROOT}/node/set1.proc -of-proc log_1/scale-rm_ens.NOUT_${STIME} ./scale-rm_ens scale-rm_ens_${STIME}.conf &

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] ${STIME}: Submitted to background: ensemble forecasts" >&2

mpiexec -n $((NNODES_ORIG*PPN)) -vcoordfile ${TMPROOT}/node/set2.proc -of-proc log_2/letkf.NOUT_${STIME} ./letkf letkf_${STIME}.conf &

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] ${STIME}: Submitted to background: LETKF" >&2

jobids=\$(jobs -p)
while [ -n "\$jobids" ]; do
  for job in \$jobids; do
    if ! (kill -0 \$job 2> /dev/null); then
      rcode=0
      wait \$job || rcode=\$?
      if ((rcode != 0)); then
        echo "[\$(date +'%Y-%m-%d %H:%M:%S')] ${STIME}: One of the background programs crashed..." >&2
        exit \$rcode
      fi
    fi
  done
  jobids=\$(jobs -p)
  sleep 1s
done

echo "[\$(date +'%Y-%m-%d %H:%M:%S')] Finish $STIME $ETIME" >&2

exit 0
EOF
#else
#  cat >> $jobscrp << EOF
#
#./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
#EOF
fi

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

if ((OBS_USE_JITDT == 1)); then
  rm -f ${SCRP_DIR}/jitdt_send_testdata_on_start.log
  ${SCRP_DIR}/misc/jitdt_send_testdata_on_start.sh > ${SCRP_DIR}/jitdt_send_testdata_on_start.log 2>&1 &
fi

job_end_check_PJM_K $jobid
res=$?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $SCRP_DIR $jobid ${job}_${SYSNAME} 'o e i s' i

if [ "$CONF_MODE" = 'static' ]; then
  config_file_save $TMPS/config || exit $?
fi

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMPS
  safe_rm_tmpdir $TMPSL
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
