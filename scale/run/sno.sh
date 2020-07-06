#!/bin/bash

#
# SNO is executed by a bulk job
#

# Based on 
#  https://github.com/gylien/scale-letkf/blob/realtime_fcst_D1-3_ope/scale/run_d1-2/config/online_NRT_5.3.X/sno_bulk.sh
#

#STIME=<STIME>
STIME='20180701060000'

. config.main || exit $?
RUNDIR="${TMP}/../run_sno_<STIME>"
OUTDIR=${OUTDIR[1]}

PPN=64 # Process per node

# Which domain do you want to convert?
DOM=1

# Output file (X & Y process number) for each member
NP_OFILE_X=1
NP_OFILE_Y=1

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

SNO_MEMBERS=$((${MEMBER} + 2))
SNO_MEM_L=$(seq -f %04g ${MEMBER})" mean mdet" # All members + mean + mdet
#SNO_MEMBERS=1
#SNO_MEM_L="mdet"
echo $SNO_MEMBERS
echo $SNO_MEM_L

# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE}))

# Convert variables (Other variables will NOT be included in converted files)
#VARS='"U", "V", "T", "QV", "QHYD", "DENS", "PREC", "MSLP", "Gprs", "PW", "ENGT", "ENGP", "ENGI", "ENGK", "MSE" , "Uprs", "Vprs", "Tprs", "Gprs", "QVprs","QHYDprs"'
VARS='"U", "V", "T", "QV", "DENS", "PREC", "MSLP"'

TOPO=0 # Process topography file? # 1: Yes, 0: No
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

conf_bulk="${RUNDIR}/conf/bulk_sno.conf"

cnt=0
for mem in  ${SNO_MEM_L} # member loop
do 
  cnt=$((${cnt} + 1))
  echo $mem

  SNO_OUTPUT=${OUTDIR}/${STIME}/fcst_sno_np$(printf %05d ${NP_OFILE})/${mem}
  SNO_BASENAME_IN="${OUTDIR}/${STIME}/fcst/${mem}/history"
  SNO_BASENAME_OUT="${SNO_OUTPUT}/history"

  if (( TOPO > 0 )) ; then
    SNO_OUTPUT=${OUTDIR}/const/topo_sno_np$(printf %05d ${NP_OFILE})
    SNO_BASENAME_IN="${OUTDIR}/const/topo/topo"
    SNO_BASENAME_OUT="${SNO_OUTPUT}/topo"
  fi

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
 basename_in  = "${SNO_BASENAME_IN}",
 basename_out  = "${SNO_BASENAME_OUT}",
 vars         = ${VARS},
 nprocs_x_out = ${NP_OFILE_X},
 nprocs_y_out = ${NP_OFILE_Y},
 output_grads = .false.,
/
&PARAM_SNOPLGIN_HGRIDOPE
 SNOPLGIN_hgridope_type="LATLON",
 SNOPLGIN_hgridope_lat_start = 12.0,
 SNOPLGIN_hgridope_lat_end = 54.0,
 SNOPLGIN_hgridope_dlat = 0.5,
 SNOPLGIN_hgridope_lon_start = 96.0,
 SNOPLGIN_hgridope_lon_end = 174.0,
 SNOPLGIN_hgridope_dlon = 0.5,
/
EOF

  ln -s ${conf} ${conf_bulk}.${cnt}

done # member loop


NPIN=`expr 255 / \( $PPN \) + 1`
jobsh="${RUNDIR}/job_sno.sh"

cat << EOF >> $jobsh
#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L node=${SNO_NODE}
#PJM -L elapse="00:30:00"
#PJM --mpi proc=${NP_TOTAL}
#PJM --omp thread=1
#PJM -g $(echo $(id -ng))

module unload impi
module unload intel
module load intel/2019.5.281

module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5

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
mpiexec.hydra -n ${NP_OFILE} ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"
EOF

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd -

