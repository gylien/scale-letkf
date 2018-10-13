#!/bin/bash

. config.main



myname1=makeEns

TMPDIR=../tmp/$TMPSUBDIR
mkdir -p $TMPDIR

TIME_LIMIT="00:10:00"
RUNCONF_COMMON=$TMPDIR/makeEns.conf

rm -f $RUNCONF_COMMON
cat << EOF >> $RUNCONF_COMMON
&PARAM_ENSEMBLE
 MEMBER = $MEMBER,
!--MEMBER_RUN--
!--MEMBER_ITER--
/

&PARAM_LETKF_PRC
 NNODES = ${NNODES},
 PPN = ${PPN},
 MEM_NODES = ${SCALE_NP},
 MEM_NP = ${SCALE_NP},
/

EOF


# Creat a job script

jobscrp="${myname1}_job.sh"

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

mpirun -np $((PPN*NNODES)) ../letkf/${myname1} $RUNCONF_COMMON
EOF



