#!/bin/bash

#
# SNO is executed by a bulk job
#

GRADS=F
PLEV=T


INPUT_FROM_SNOW=F
INPUT_SNOW_NP=16

SEXP=TAIWAN201808_D2_NOHIM8
SEXP=TAIWAN201808_D2_HIM8_1209
SEXP=TAIWAN201808_D2_HIM8_1209_HLOC120
SEXP=TAIWAN201808_D2_HIM8_1209_UVW0.0
SEXP=TAIWAN201808_D2_HIM8_1209_NOPREP

SEXP=TAIWAN201808_D2_HIM8_1218
#SEXP=TAIWAN201808_D2_NOHIM8_1203
SEXP=TAIWAN201808_D2_HIM8_1218_HLOC120km

#SEXP=TAIWAN201808_D2_HIM8_1209_NOVLOC
#SEXP=TAIWAN201808_D2_NOHIM8_1203

SEXP=TAIWAN201808_D2_HIM8_0116_VLOC0.5_NOGC
SEXP=TAIWAN201808_D2_NOHIM8_0116


SEXP=TAIWAN201808_D2_HIM8_0116_VLOC0.5_GC30K

SEXP=TAIWAN201808_D1_1118


tint=3600 # [second]

tstart='2018-08-21 7:00:00'
tend='2018-08-21 12:00:00'
tstart='2018-08-21 12:00:00'
tend=$tstart

. config/${SEXP}/config.main.ofp || exit $?
EXP=${SEXP}
RUNDIR="${TMP}_sno"


SCALEDIR=/work/hp150019/share/honda/SCALE-LETKF/scale-5.3.3/develop_scale/scale-develop
#PPN=$PPN # Process per node

TYPE=fcst
#TYPE=anal
#TYPE=gues
TYPE=hist

## Which domain do you want to convert?
#DOM=2 

# Output file (X & Y process number) for each member
NP_OFILE_X=32
NP_OFILE_Y=32

NP_OFILE_X=4
NP_OFILE_Y=4

if [ "$INPUT_FROM_SNOW" == "T" ] ; then
  NP_OFILE_X=1
  NP_OFILE_Y=1
fi

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

#SNO_MEMBERS=${MEMBER}
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean mdet" # All members + mean + mdet
SNO_MEM_L=" mean"
SNO_MEMBERS=1

#SNO_MEMBERS=50
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})


# Convert variables (Other variables will NOT be included in converted files)
#VARS='"Umet", "Vmet", "W", "T", "QV", "QHYD", "PRES", "RAIN", "MSLP"'
VARS='"MSLP",'
#NVAR='MSLP'
VARS="'PRES', 'SFC_PRES'"

if [ "$PLEV" == "T" ] ; then
  VARS="'QV', 'Umet', 'Vmet', 'W', 'T', 'QHYD', 'MSLP', 'PRES', 'SFC_PRES','PREC', 'PW'"
  #VARS="'QV', 'W', 'T', 'QHYD', 'MSLP', 'PRES', 'SFC_PRES'" # tentative
  if [ "$INPUT_FROM_SNOW" == "T" ] ; then
    VARS="'QV', 'Umet', 'Vmet', 'W', 'T', 'QHYD', 'MSLP', 'SFC_PRES', 'GPH', 'PREC', 'PW'"
  fi
fi

if [ "$TYPE" != "fcst" ] && [ "$TYPE" != "hist" ] ; then
  VARS="'RHOT', 'DENS', 'MOMX', 'MOMY', 'MOMZ'"
  PLEV=F
fi

TOPO=0 # Process topography file? # 1: Yes, 0: No
if (( TOPO > 0 )) ; then
  VARS='"TOPO"'
  SNO_MEM_L="mean"
fi


if [ "$GRADS" == "T" ] ; then
  NP_OFILE_X=1
  NP_OFILE_Y=1
  OUTPUT_GRADS=".true."
  OUTPUT_GRADSCTL=".true."
else
  OUTPUT_GRADS=".false."
  OUTPUT_GRADSCTL=".false."
fi
OUTPUT_SINGLE=".true."
OUTPUT_SINGLE=".false."

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

#exec_l="${RUNDIR}/tmp_exec"
#rm -f $exec_l
#touch $exec_l

time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do # time loop

  timef=$(date -ud "$time" '+%Y-%m-%d %H:%M:%S')

  YYYYs=$(date -ud "$timef" '+%Y')
  MMs=$(date -ud "$timef" '+%m')
  DDs=$(date -ud "$timef" '+%d')
  HHs=$(date -ud "$timef" '+%H')
  MNs=$(date -ud "$timef" '+%M')
  SSs=$(date -ud "$timef" '+%S')

  DTIME=${YYYYs}${MMs}${DDs}${HHs}${MNs}${SSs}
  echo $DTIME


  for mem in  ${SNO_MEM_L} # member loop
  do 
    cnt=$((${cnt} + 1))
    echo $mem
  
    SNO_BASENAME_IN="${OUTDIR}/${DTIME}/${TYPE}/${mem}/history"
  
    if [ "$PLEV" == "T" ] ; then
      SNO_BASENAME_OUT="p_history"
  #    SNO_BASENAME_OUT="$NVAR"
    else
      SNO_BASENAME_OUT="history"
    fi
  
    if [ "$TYPE" != "fcst" ] && [ "$TYPE" != "hist" ]  ; then
      SNO_BASENAME_OUT="$TYPE"
      SNO_BASENAME_IN="${OUTDIR}/${DTIME}/${TYPE}/${mem}/init"
    fi
  
  
    if [ "$INPUT_FROM_SNOW" == "T" ] ; then
      SNO_BASENAME_IN=${OUTDIR}/${DTIME}/${TYPE}_sno_np$(printf %05d ${INPUT_SNOW_NP})/${mem}/${SNO_BASENAME_OUT}
    fi
  
    if [ "$GRADS" == "T" ] ; then
      SNO_OUTPUT_PATH=${OUTDIR}/${DTIME}/${TYPE}_sno_grads/${mem}
    else
      SNO_OUTPUT_PATH=${OUTDIR}/${DTIME}/${TYPE}_sno_np$(printf %05d ${NP_OFILE})/${mem}
    fi
  
    if (( TOPO > 0 )) ; then
      SNO_OUTPUT_PATH=${OUTDIR}/const/topo_sno_np$(printf %05d ${NP_OFILE})
      SNO_BASENAME_IN="${OUTDIR}/const/topo/topo"
      SNO_BASENAME_OUT="${SNO_OUTPUT_PATH}/topo"
    fi
  
    if [ ! -e ${SNO_OUTPUT_PATH} ] ; then
      mkdir -p ${SNO_OUTPUT_PATH}
    fi
  
  
    conf="${RUNDIR}/conf/sno_${mem}_${DTIME}.conf"
  
cat << EOF >> $conf
&PARAM_IO
 IO_LOG_BASENAME = "log/LOG_${mem}_${time}",
 IO_LOG_ALLNODE = .false.,
 IO_LOG_SUPPRESS = .true.,
 IO_LOG_NML_SUPPRESS = .false.,
/

&PARAM_SNO
 basename_in  = "${SNO_BASENAME_IN}",
 basename_out  = "${SNO_BASENAME_OUT}",
 dirpath_out = "${SNO_OUTPUT_PATH}",
 vars         = ${VARS},
 output_single = ${OUTPUT_SINGLE},
 nprocs_x_out = ${NP_OFILE_X},
 nprocs_y_out = ${NP_OFILE_Y},
 output_gradsctl = ${OUTPUT_GRADSCTL},
 output_grads = ${OUTPUT_GRADS},
/

EOF

  if [ "$PLEV" == "T" ] &&  [ "$INPUT_FROM_SNOW" != "T" ] ; then
  
cat << EOF >> $conf
&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type        = 'PLEV', 
 SNOPLGIN_vgridope_lev_num     = 15,
 SNOPLGIN_vgridope_lev_data    = 1000.e+2, 950.e+2, 925.e+2, 900.e+2, 850.e+2, 800.e+2, 700.e+2, 600.e+2, 500.e+2, 400.e+2, 300.e+2, 200.e+2, 100.e+2, 70.e+2, 50.e+2, 
! SNOPLGIN_vgridope_lev_num     = 13,
! SNOPLGIN_vgridope_lev_data    = 1000.e+2, 950.e+2, 925.e+2, 900.e+2, 850.e+2, 800.e+2, 700.e+2, 600.e+2, 500.e+2, 400.e+2, 300.e+2, 200.e+2, 100.e+2, 
/
EOF

  fi
  
    ln -s ${conf} ${conf_bulk}.${cnt}

#cat << EOF >> $exec_l
#mpirun -np ${NP_OFILE} ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
#EOF

  done # member loop

  time=$(date -ud "${tint} second $time" '+%Y-%m-%d %H:%M:%S')
done # time loop


# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE} * ${cnt}))

# Get total SNO NODEs
if (( NP_TOTAL < PPN )) ; then
  SNO_NODE=1
else
  SNO_NODE=$((${NP_TOTAL} / ${PPN}))
  if (( SNO_NODE*PPN < NP_TOTAL )) ; then
    SNO_NODE=$((SNO_NODE + 1))
  fi
fi

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
mpirun -np $((NP_OFILE)) ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"
EOF

cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd - > /dev/null

