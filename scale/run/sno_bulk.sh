#!/bin/bash

#
# SNO is executed by a bulk job
#


EXP=D4_500M_mercator_50mem_verify
STIME="20190610080000"

# Load config.main or specify some shell parameters
#. config/${EXP}/config.main.ofp || exit $?
#RUNDIR="${TMP}_sno"

OUTDIR=/work/hp150019/share/honda/SCALE-LETKF/TEST_INPUT/D4_500m20191206_np256
RUNDIR="$(pwd)/../tmp/sno"
PPN=64


SCALEDIR=/work/hp150019/share/honda/SCALE-LETKF/scale-5.3.3/develop_scale/scale-develop
#PPN=$PPN # Process per node

TYPE=fcst
#TYPE=hist

## Which domain do you want to convert?
#DOM=2 

# Output file (X & Y process number) for each member
NP_OFILE_X=4
NP_OFILE_Y=4
NP_OFILE_X=1
NP_OFILE_Y=1

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

#SNO_MEMBERS=${MEMBER}
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean mdet" # All members + mean + mdet
SNO_MEMBERS=1
SNO_MEM_L=" mean"

# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE}))

# Convert variables (Other variables will NOT be included in converted files)
VARS='"Umet", "Vmet", "W", "T", "QV", "QHYD", "PRES", "RAIN", "MSLP"'
#VARS=' "MSLP"'
VARS=' "W", "QHYD", "T", "Umet", "Vmet", "QV"'

TOPO=1 # Process topography file? # 1: Yes, 0: No
if (( TOPO > 0 )) ; then
  VARS='"TOPO"'
  SNO_MEM_L="mean"
fi

# Get total SNO NODEs
if (( NP_TOTAL < PPN )) ; then
  SNO_NODE=1
else
  SNO_NODE=$((${NP_TOTAL} / ${PPN}))
fi

###############################

# Path for SNO binary
SNOBIN_ORG=${SCALEDIR}/bin/sno
SNOBIN=${RUNDIR}/sno
if [ ! -e ${SNOBIN_ORG} ] ; then
  echo "No SNO binary!"
  exit
fi

###############################


rm -rf ${RUNDIR}

mkdir -p ${RUNDIR}/conf
mkdir -p ${RUNDIR}/log

# copy binary 
cp ${SNOBIN_ORG} ${SNOBIN}

conf_bulk="${RUNDIR}/conf/balk_sno.conf"

cnt=0
for mem in  ${SNO_MEM_L} # member loop
do 
  cnt=$((${cnt} + 1))
  echo $mem

  SNO_OUTPUT=${OUTDIR}/${STIME}/${TYPE}_sno_np$(printf %05d ${NP_OFILE})/${mem}
  SNO_BASENAME_IN="${OUTDIR}/${STIME}/${TYPE}/${mem}/history"
  SNO_BASENAME_OUT="${SNO_OUTPUT}/history"

  if (( TOPO > 0 )) ; then
    SNO_OUTPUT=${OUTDIR}/const/topo_sno_np$(printf %05d ${NP_OFILE})
    SNO_BASENAME_IN="${OUTDIR}/const/topo/topo"
    SNO_BASENAME_OUT="${SNO_OUTPUT}/topo"
  fi
  SNO_DIRPATH_OUT=""

  if [ ! -e ${SNO_OUTPUT} ] ; then
    mkdir -p ${SNO_OUTPUT}
  fi

  conf="${RUNDIR}/conf/sno_${mem}.conf"

cat << EOF >> $conf
&PARAM_IO
 IO_LOG_BASENAME = "log/LOG_${mem}",
 IO_LOG_ALLNODE = .false.,
 IO_LOG_SUPPRESS = .false.,
 IO_LOG_NML_SUPPRESS = .false.,
/
&PARAM_SNO
 dirpath_out = "${SNO_DIRPATH_OUT}",
 basename_in  = "${SNO_BASENAME_IN}",
 basename_out  = "${SNO_BASENAME_OUT}",
 vars         = ${VARS},
 output_single = .true.,
 nprocs_x_out = ${NP_OFILE_X},
 nprocs_y_out = ${NP_OFILE_Y},
 output_gradsctl = .false.,
! output_gradsctl = .true.,
 output_grads = .false.,
! output_grads = .true.,
/
EOF

  ln -s ${conf} ${conf_bulk}.${cnt}

done # member loop



jobsh="${RUNDIR}/job_sno.sh"

cat << EOF >> $jobsh
#!/bin/sh
#PJM -L rscgrp=debug-cache
#PJM -L node=${SNO_NODE}
#PJM -L elapse="00:30:00"
#PJM --mpi proc=${NP_TOTAL}
#PJM --omp thread=1
#PJM -g $(echo $(id -ng))
#PJM -s

module load hdf5/1.8.17
module load netcdf/4.4.1
module load netcdf-fortran/4.4.3

export FORT_FMT_RECL=400
export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY
export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t
ulimit -s unlimited
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] Start SNO"
mpirun -np ${NP_OFILE} ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"
EOF

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd -
