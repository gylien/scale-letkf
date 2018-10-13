#!/bin/sh
##PBS -N makeEns_8km_sc
#PBS -l nodes=1:ppn=4
##PBS -l walltime=00:10:00
#PBS -W umask=027
##PBS -k oe

ulimit -s unlimited

HOSTLIST=$(cat $PBS_NODEFILE | sort | uniq)
HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')
export MPI_UNIVERSE="$HOSTLIST 4"

export OMP_NUM_THREADS=1
#export PARALLEL=1

export FORT_FMT_RECL=400

cd $PBS_O_WORKDIR

rm -f machinefile
cp -f $PBS_NODEFILE machinefile

mpirun -np 4 ../letkf/makeEns ../tmp/scale-letkf_8km_sc/makeEns.conf
