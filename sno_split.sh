#!/bin/bash

#
# SNO is executed by a bulk job
#

# Based on 
#  https://github.com/gylien/scale-letkf/blob/realtime_fcst_D1-3_ope/scale/run_d1-2/config/online_NRT_5.3.X/sno_bulk.sh
#

#STIME=<STIME>
STIME='20180624000000'

. config.main || exit $?
RUNDIR="${TMP}/../run_sno_$STIME"
OUTDIR=${OUTDIR[1]}

PPN=48 # Process per node

# Which domain do you want to convert?
DOM=1

# Output file (X & Y process number) for each member
NP_OFILE_X=16
NP_OFILE_Y=12

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

MEMBER=0
SNO_MEMBERS=$((${MEMBER} + 1))
SNO_MEM_L=$(seq -f %04g ${MEMBER})" mean " # All members + mean + mdet

MEMBER=50
SNO_MEMBERS=20 #$((${MEMBER} + 0))
SNO_MEM_L=$(seq -f %04g 31 ${MEMBER}) 

#SNO_MEMBERS=1
#SNO_MEM_L="mdet"
echo $SNO_MEMBERS
echo $SNO_MEM_L


# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE}))

# Convert variables (Other variables will NOT be included in converted files)
#VARS='"U", "V", "T", "QV", "QHYD", "DENS", "PREC", "MSLP", "Gprs", "PW", "ENGT", "ENGP", "ENGI", "ENGK", "MSE" , "Uprs", "Vprs", "Tprs", "Gprs", "QVprs","QHYDprs"'
#VARS='"U", "V", "T", "QV", "DENS", "PREC", "MSLP"'
VARS='"MOMX", '

VARS=" 'OCEAN_TEMP', 'OCEAN_OCN_Z0M', 'OCEAN_ICE_TEMP', 'OCEAN_ICE_MASS', 'OCEAN_SFC_TEMP', 'OCEAN_SFC_ALB_IR_dir', 'OCEAN_SFC_ALB_IR_dif', 'OCEAN_SFC_ALB_NIR_dir', 'OCEAN_SFC_ALB_NIR_dif', 'OCEAN_SFC_ALB_VIS_dir', 'OCEAN_SFC_ALB_VIS_dif', 'OCEAN_SFC_Z0M', 'OCEAN_SFC_Z0H', 'OCEAN_SFC_Z0E', 'LAND_TEMP', 'LAND_WATER', 'LAND_SFC_TEMP', 'LAND_SFC_ALB_IR_dir', 'LAND_SFC_ALB_IR_dif', 'LAND_SFC_ALB_NIR_dir', 'LAND_SFC_ALB_NIR_dif', 'LAND_SFC_ALB_VIS_dir', 'LAND_SFC_ALB_VIS_dif', 'URBAN_TR', 'URBAN_TB', 'URBAN_TG', 'URBAN_TC', 'URBAN_QC', 'URBAN_UC', 'URBAN_TRL', 'URBAN_TBL', 'URBAN_TGL', 'URBAN_RAINR', 'URBAN_RAINB', 'URBAN_RAING', 'URBAN_ROFF', 'URBAN_SFC_TEMP', 'URBAN_SFC_ALB_IR_dir', 'URBAN_SFC_ALB_IR_dif', 'URBAN_SFC_ALB_NIR_dir', 'URBAN_SFC_ALB_NIR_dif', 'URBAN_SFC_ALB_VIS_dir', 'URBAN_SFC_ALB_VIS_dif', 'DENS', 'MOMZ', 'MOMX', 'MOMY', 'RHOT', 'QV', 'QC', 'QR', 'QI', 'QS', 'QG', 'TKE_MYNN', 'SFLX_rain', 'SFLX_snow', 'SFLX_LW_up', 'SFLX_LW_dn', 'SFLX_SW_up', 'SFLX_SW_dn', 'TOAFLX_LW_up', 'TOAFLX_LW_dn', 'TOAFLX_SW_up', 'TOAFLX_SW_dn', 'SFLX_IR_dn_dir', 'SFLX_IR_dn_dif', 'SFLX_NIR_dn_dir', 'SFLX_NIR_dn_dif', 'SFLX_VIS_dn_dir', 'SFLX_VIS_dn_dif', 'MFLX_cloudbase', 'SFLX_convrain', 'cloudtop', 'cloudbase', 'cldfrac_dp', 'cldfrac_sh', 'w0mean', 'kf_nca', 'DENS_t_CP', 'RHOT_t_CP', 'QV_t_CP', 'QC_t_CP', 'QR_t_CP', 'QI_t_CP', 'QS_t_CP', 'QG_t_CP', 'QH_t_CP' "


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

#  SNO_OUTPUT=${OUTDIR}/${STIME}/fcst_sno_np$(printf %05d ${NP_OFILE})/${mem}
#  SNO_BASENAME_IN="${OUTDIR}/${STIME}/fcst/${mem}/history"
#  SNO_BASENAME_OUT="${SNO_OUTPUT}/history"

  mkdir -p ${OUTDIR}/${STIME}/anal_192p_overwriten/${mem}
  SNO_BASENAME_IN="${OUTDIR}/${STIME}/anal_48p_overwriten/${mem}/init"
  SNO_BASENAME_OUT="${OUTDIR}/${STIME}/anal_192p_overwriten/${mem}/init"

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
!&PARAM_SNOPLGIN_HGRIDOPE
! SNOPLGIN_hgridope_type="LATLON",
! SNOPLGIN_hgridope_lat_start = 12.0,
! SNOPLGIN_hgridope_lat_end = 54.0,
! SNOPLGIN_hgridope_dlat = 0.5,
! SNOPLGIN_hgridope_lon_start = 96.0,
! SNOPLGIN_hgridope_lon_end = 174.0,
! SNOPLGIN_hgridope_dlon = 0.5,
!/
EOF

  ln -s ${conf} ${conf_bulk}.${cnt}

done # member loop



jobsh="${RUNDIR}/job_sno.sh"

cat << EOF >> $jobsh
#!/bin/sh
#PJM -L rscgrp=debug-flat
#PJM -L node=${SNO_NODE}
#PJM -L elapse="00:30:00"
#PJM --mpi proc=${NP_TOTAL}
#PJM --omp thread=1
#PJM -g $(echo $(id -ng))

module unload impi
module unload intel
module load intel/2019.5.281

source /work/opt/local/cores/intel/performance_snapshots_2019.6.0.602217/apsvars.sh
export MPS_STAT_LEVEL=4
 
module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5

export FORT_FMT_RECL=400


export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY
#export KMP_AFFINITY=verbose
#export I_MPI_DEBUG=5

export OMP_NUM_THREADS=1
#export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t

export PSM2_CONNECT_WARN_INTERVAL=2400
export TMI_PSM2_CONNECT_TIMEOUT=2000


#export OMP_STACKSIZE=128m
ulimit -s unlimited

echo "[\$(date "+%Y/%m/%d %H:%M:%S")] Start SNO"
mpiexec.hydra -n ${NP_OFILE} ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"
EOF

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd -

