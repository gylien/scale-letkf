#!/bin/sh

#PBS -q s
#PBS -l nodes=1:ppn=16
#PBS -N OBSSIM_2km_CZ2003
#PBS -W umask=027
#PBS -k oe

ulimit -s unlimited

HOSTLIST=$(cat $PBS_NODEFILE | sort | uniq)
HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')
export MPI_XPMEM_ENABLED=disabled
export MPI_UNIVERSE="$HOSTLIST 16"
export MPI_XPMEM_ENABLED=disabled

export OMP_NUM_THREADS=1
#export PARALLEL=1

export FORT_FMT_RECL=400

cd $PBS_O_WORKDIR

rm -f machinefile
cp -f $PBS_NODEFILE machinefile

export RUN_LEVEL=1

mpirun -n 16 /home/honda/work/scale_lt_devel_20181002/scale-letkf/scale/run/../obs/obssim /home/honda/work/scale_lt_devel_20181002/scale-letkf/scale/run/OBSSIM_001.conf &
wait
