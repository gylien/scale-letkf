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
myname1='cycle'

#===============================================================================
# Configuration

. config.main
res=$? && ((res != 0)) && exit $res
. config.$myname1
res=$? && ((res != 0)) && exit $res

#. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"
echo

setting "$1" "$2" "$3" "$4" "$5"

echo
print_setting
echo

#===============================================================================
# Creat a job script

jobscrp="${myname1}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

cat > $jobscrp << EOF
#!/bin/sh
##PBS -N ${myname1}_${SYSNAME}
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

./${myname1}.sh "$STIME" "$ETIME" "$ISTEP" "$FSTEP" || exit \$?
EOF

echo "[$(datetime_now)] Run ${myname1} job on PBS"
echo

job_submit_torque $jobscrp
echo

job_end_check_torque $jobid
res=$?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

n=0
nmax=12
while [ ! -s "${jobscrp}.o${jobid}" ] && ((n < nmax)); do
  n=$((n+1))
  sleep 5s
done

mkdir -p $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.main $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.${myname1} $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/config.nml.* $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f $SCRP_DIR/${myname1}_job.sh $OUTDIR/exp/${jobid}_${myname1}_${STIME}
cp -f ${jobscrp}.o${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.o
cp -f ${jobscrp}.e${jobid} $OUTDIR/exp/${jobid}_${myname1}_${STIME}/job.e
( cd $SCRP_DIR ; git log -1 --format="SCALE-LETKF version %h (%ai)" > $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )
( cd $MODELDIR ; git log -1 --format="SCALE       version %h (%ai)" >> $OUTDIR/exp/${jobid}_${myname1}_${STIME}/version )

finalization

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMPS
  safe_rm_tmpdir $TMPSL
fi

#===============================================================================

echo "[$(datetime_now)] Finish $(basename $0) $@"

exit $res
