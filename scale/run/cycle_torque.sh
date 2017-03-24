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
##PBS -N ${job}_${SYSNAME}
#PBS -l nodes=${NNODES}:ppn=${PPN}
##PBS -l walltime=${TIME_LIMIT}
#PBS -W umask=027
##PBS -k oe

ulimit -s unlimited

HOSTLIST=\$(cat \$PBS_NODEFILE | sort | uniq)
HOSTLIST=\$(echo \$HOSTLIST | sed 's/  */,/g')
export MPI_UNIVERSE="\$HOSTLIST $((PPN*THREADS))"

export OMP_NUM_THREADS=${THREADS}
#export PARALLEL=${THREADS}

export FORT_FMT_RECL=400

cd \$PBS_O_WORKDIR

rm -f machinefile
cp -f \$PBS_NODEFILE machinefile

./${job}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" || exit \$?
EOF

echo "[$(datetime_now)] Run ${job} job on PBS"
echo

job_submit_torque $jobscrp
echo

job_end_check_torque $jobid
res=$?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $SCRP_DIR $jobid $jobscrp 'o e'

archive_log

#if ((CLEAR_TMP == 1)); then
#  safe_rm_tmpdir $TMP
#fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res
