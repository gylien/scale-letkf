#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a OFP job script and run it.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_ofp.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='fcst'

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_common_static.sh || exit $?
. src/func_${job}_static.sh || exit $?

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@"

setting "$@" || exit $?


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

safe_init_tmpdir $NODEFILE_DIR || exit $?
#distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cat $SCRP_DIR/config.main | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMP/config.main
echo "NODEFILE_DIR=\"\$TMPROOT/node\"" >> $TMPS/config.main
echo "RUN_LEVEL=4" >> $TMP/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMP/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
staging_list_static || exit $?
config_file_list $TMPS/config || exit $?

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

cp ${SCRP_DIR}/config.rc $TMP/config.rc
cp ${SCRP_DIR}/config.${job} $TMP/config.${job}
cp ${SCRP_DIR}/${job}.sh $TMP/${job}.sh
cp -r ${SCRP_DIR}/src $TMP/src

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

stage_in server || exit $?

#===============================================================================
# Creat a job script

NPIN=`expr 255 / \( $PPN \) + 1`
jobscrp="$TMP/${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"


# OFP
if [ "$PRESET" = 'OFP' ]; then

  if [ "$RSCGRP" == "" ] ; then
    RSCGRP="regular-cache"
  fi

cat > $jobscrp << EOF
#!/bin/sh
#PJM -L rscgrp=${RSCGRP}
#PJM -L node=${NNODES}
#PJM -L elapse=${TIME_LIMIT}
#PJM --mpi proc=$((NNODES*PPN))
##PJM --mpi proc=${totalnp}
#PJM --omp thread=${THREADS}

#PJM -g $(echo $(id -ng))
# HPC
##PJM -g gx14  

#PJM -s

module unload impi
module unload intel
module load intel/2019.5.281

source /work/opt/local/cores/intel/performance_snapshots_2019.6.0.602217/apsvars.sh
export MPS_STAT_LEVEL=4
 
module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5

export FORT_FMT_RECL=400

export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY
#export KMP_AFFINITY=verbose
#export I_MPI_DEBUG=5

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t

export PSM2_CONNECT_WARN_INTERVAL=2400
export TMI_PSM2_CONNECT_TIMEOUT=2000


#export OMP_STACKSIZE=128m
ulimit -s unlimited

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo
  
  job_submit_PJM $jobscrp
  echo
  
  job_end_check_PJM $jobid
  res=$?

else

# qsub
cat > $jobscrp << EOF
#!/bin/sh
#PBS -N $job
#PBS -q s
#PBS -l nodes=${NNODES}:ppn=${PPN}
#PBS -l walltime=${TIME_LIMIT}
#
#

cd \${PBS_O_WORKDIR}

export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

source /etc/profile.d/modules.sh 
module unload mpt/2.12
module load intelmpi/5.1.2.150


export OMP_NUM_THREADS=${THREADS}
export KMP_AFFINITY=compact


ulimit -s unlimited

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo
  
  job_submit_torque $jobscrp
  echo
  
  job_end_check_torque $jobid
  res=$?

fi

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_job.sh 'o e'

archive_log


#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
