#!/bin/bash
#===============================================================================
#
#  Wrap cycle.sh in an OFP job script and run it.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_ofp.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='cycle'

RSCGRP=${RSCGRP:-"regular-flat"}
GNAME=${GNAME:-`id -ng`}

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

#-------------------------------------------------------------------------------

statfile=${myname%.*}.stat.$STIME

echo "prep" > $statfile

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

#if [ -e "${TMP}" ]; then
#  echo "[Error] $0: \$TMP will be completely removed." >&2
#  echo "        \$TMP = '$TMP'" >&2
#  exit 1
#fi
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
if ((IO_ARB == 1)); then                              ##
  distribute_da_cycle_set - $NODEFILE_DIR || exit $?  ##
else                                                  ##
  distribute_da_cycle - $NODEFILE_DIR || exit $?
fi                                                    ##

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cat $SCRP_DIR/config.main | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

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

if ((IO_ARB == 1)); then                                              ##
  echo "${SCRP_DIR}/sleep.sh|sleep.sh" >> ${STAGING_DIR}/${STGINLIST} ##
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

NPIN=`expr 255 / \( $PPN \) + 1`

cat > $jobscrp << EOF
#!/bin/sh -l
#PJM -L rscgrp=${RSCGRP}
#PJM -L node=${NNODES}
#PJM -L elapse=${TIME_LIMIT}
#PJM --mpi proc=$((NNODES*PPN))
##PJM --mpi proc=${totalnp}
#PJM --omp thread=${THREADS}
#PJM -g ${GNAME}
##PJM -j
#PJM -S

rm -f machinefile
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
  for ippn in \$(seq $PPN); do
    echo "\$inode" >> machinefile
  done
done

ulimit -s unlimited

#export OMP_STACKSIZE=128m
export OMP_NUM_THREADS=1

export I_MPI_PIN_PROCESSER_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_PJOLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi

export I_MPI_PERHOST=${PPN}
export I_MPI_PIN_DOMAIN=${NPIN}

export KMP_HW_SUBSET=1t

export HFI_NO_CPUAFFINITY=1
unset KMP_AFFINITY

export I_MPI_DEBUG=5

source /work/opt/local/cores/intel/performance_snapshots_2019.6.0.602217/apsvars.sh
source /work/opt/local/cores/intel/parallel_studio_xe_2019/bin/psxevars.sh
source /work/opt/local/cores/intel/itac_2019/bin/itacvars.sh
 
export MPS_STAT_LEVEL=4


module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5


./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

echo 'submit' > $statfile

job_submit_PJM $jobscrp
echo

job_end_check_PJM $jobid
res=$?

echo 'plot' > $statfile

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_job.sh 'o e'

if [ "$CONF_MODE" = 'static' ]; then
  config_file_save $TMPS/config || exit $?
fi

archive_log

#if ((CLEAR_TMP == 1)); then
#  safe_rm_tmpdir $TMP
#fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
