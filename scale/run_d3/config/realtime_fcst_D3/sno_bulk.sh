#!/bin/bash

#
# SNO is executed by a bulk job
#

STIME=<STIME>

. config.main || exit $?



PARENT_REF_TIME=<PARENT_REF_TIME>
FCSTLEN=<FCSTLEN>
RUNDIR="${TMP}/../run_sno_ref_${PARENT_REF_TIME}_${STIME}"
OUTDIR=${OUTPUT}/${EXP3}/ref_${PARENT_REF_TIME}
PLOTDIR_SRC=${DIR}/run_d3/plot
PLOTDIR=${RUNDIR}/plot

PPN=64 # Process per node

# Which domain do you want to convert?
DOM=3 

# Output file (X & Y process number) for each member
NP_OFILE_X=1
NP_OFILE_Y=1

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

SNO_MEMBERS=${MEMBER}
#SNO_MEMBERS=1
SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean mdet" # All members + mean + mdet
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean" # All members + mean
#SNO_MEM_L="mdet"

# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE}))

# Convert variables (Other variables will NOT be included in converted files)
#VARS='"U", "V", "W", "T", "QV", "QHYD", "PRES","RAIN", "CAPE"'
#VARS='"U", "V", "W", "T", "QV", "QHYD", "PRES","RAIN"'
#VARS='"U", "V", "W", "T", "RH", "QV", "QC", "QR", "QI", "QS", "QG", "QHYD", "DENS", "RHOT", "PRES", "RAIN", "MSLP", "T2", "U10", "V10"'
VARS='"U", "V", "W", "QR", "QS", "QG", "DENS", "RAIN"'

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
SNOBIN_ORG=${SCALEDIR}/../scale_develop/bin/sno
#SNOBIN_ORG=${SCALEDIR}/bin/sno
SNOBIN=${RUNDIR}/sno
if [ ! -e ${SNOBIN_ORG} ] ; then
  echo "No SNO binary!"
  exit
fi

###############################


rm -rf ${RUNDIR}

mkdir -p ${RUNDIR}/conf
mkdir -p ${RUNDIR}/log

cp -r ${PLOTDIR_SRC} ${PLOTDIR}
echo "character*120,parameter::cdir_base_fcst= \"${OUTPUT}/${EXP3}/\"" >>${PLOTDIR}/common_d3.h  
sed -i -e "s/<MEMBER>/${MEMBER}/g" ${PLOTDIR}/auto.sh

# copy binary 
cp ${SNOBIN_ORG} ${SNOBIN}

conf_bulk="${RUNDIR}/conf/balk_sno.conf"

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


#echo $SNO_OUTPUT
#SNO_DIR=`dirname $SNO_OUTPUT`
#echo `ls $SNO_DIR`
#error_out=`mkdir -p ${SNO_OUTPUT}`

#echo `ls $SNO_DIR`
#echo $error_out
#if [ ! $error_out ] ;then
# echo 'error'
#else
# echo 'OK'
#fi
#exit 0


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
 SNOPLGIN_hgridope_lat_start = 33.8,
 SNOPLGIN_hgridope_lat_end = 35.6,
 SNOPLGIN_hgridope_dlat = 0.01,
 SNOPLGIN_hgridope_lon_start = 133.85,
 SNOPLGIN_hgridope_lon_end = 136.05,
 SNOPLGIN_hgridope_dlon = 0.01,
/

EOF

  ln -s ${conf} ${conf_bulk}.${cnt}

done # member loop



jobsh="${RUNDIR}/job_sno.sh"

cat << EOF >> $jobsh
#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L node=${SNO_NODE}
#PJM -L elapse="01:30:00"
#PJM --mpi proc=${NP_TOTAL}
#PJM --omp thread=1
#PJM -g $(echo $(id -ng))

rm -f machinefile
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
  for ippn in \$(seq ${PPN}); do
    echo "\$inode" >> machinefile
  done
done
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

cd $PLOTDIR
./auto.sh "${PARENT_REF_TIME}" "${STIME}" "${FCSTLEN}" \${PJM_BULKNUM}

EOF

cd ${PLOTDIR}
./compile_rain.sh
./compile_dbz.sh
./compile_uvw.sh
cd -

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd -


