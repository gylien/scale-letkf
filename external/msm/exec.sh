#!/bin/sh -l

#PJM -g hp150019
#PJM -L rscgrp=debug-flat
#PJM -L node=1
#PJM -L elapse=00:20:00
#PJM --mpi proc=1
#PJM --omp thread=8
#PJM -X


PPN=1

#ulimit -s unlimited

#export FORT_FMT_RECL=400

#export HFI_NO_CPUAFFINITY=1
#export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
#export I_MPI_HBW_POLICY=hbw_preferred,,
#export I_MPI_FABRICS_LIST=tmi
#unset KMP_AFFINITY

export OMP_NUM_THREADS=8
#export I_MPI_PERHOST=${PPN}
#export KMP_HW_SUBSET=1t

./get_msm.sh 

