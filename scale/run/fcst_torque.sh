#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a torque job script and run it.
#
#  November 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_torque.sh [..]
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
. src/func_${job}.sh || exit $?

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@"
echo

setting "$@" || exit $?

if [ "$CONF_MODE" = 'static' ]; then
  . src/func_common_static.sh || exit $?
  . src/func_${job}_static.sh || exit $?
fi

echo
print_setting || exit $?
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
export MPI_XPMEM_ENABLED=disabled

export OMP_NUM_THREADS=${THREADS}
#export PARALLEL=${THREADS}

export FORT_FMT_RECL=400

cd \$PBS_O_WORKDIR

rm -f machinefile
cp -f \$PBS_NODEFILE machinefile

export RUN_LEVEL=1

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
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

echo "[$(datetime_now)] Finish $myname $@"

exit $res
