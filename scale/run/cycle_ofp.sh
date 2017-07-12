#!/bin/bash
#===============================================================================
#
#  Wrap cycle.sh in a torque job script and run it.
#
#  November 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_torque.sh [STIME ETIME CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]
#
#===============================================================================

cd "$(dirname "$0")"
job='cycle'

#===============================================================================
# Configuration

. config.main
res=$? && ((res != 0)) && exit $res
. config.$job
res=$? && ((res != 0)) && exit $res

#. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_$job.sh

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"
echo

setting "$1" "$2" "$3" "$4" "$5"

echo
print_setting
echo

#===============================================================================
# Creat a job script

jobscrp="${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

cat > $jobscrp << EOF
#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L node=${NNODES}
#PJM -L elapse=${TIME_LIMIT}
#PJM --mpi proc=${totalnp}
#PJM --omp thread=${THREADS}
#PJM -g gh51
##PJM -j

rm -f machinefile
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
  for ippn in \$(seq $PPN); do
    echo "\$inode" >> machinefile
  done
done

module load hdf5/1.8.17
module load netcdf/4.4.1
module load netcdf-fortran/4.4.3

ulimit -s unlimited
export OMP_STACKSIZE=128m

export RUN_LEVEL=1

./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

job_submit_PJM $jobscrp
echo

job_end_check_PJM $jobid
res=$?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $SCRP_DIR $jobid $jobscrp 'o e'

archive_log

#if ((CLEAR_TMP == 1)); then
#  safe_rm_tmpdir $TMP
#  safe_rm_tmpdir $TMPSL
#fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res
