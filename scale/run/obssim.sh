#!/bin/bash

MICRO=0 # 1: micro, /=1: large

USER=honda

#EXP=test15km
EXP=test03km
. config/${EXP}/config.main

RTTOV_THREADS=4
RTTOV_ITMAX=3

#
LETKF_RUN="$(pwd)"
if ((MICRO == 1)) ; then
  SWDIR="/scratch/$(id -ng)/${USER}/obssim"
  MODE="micro"
else
  SWDIR="$(pwd)/../tmp_obssim" # Not scratch
  MODE="large"
  rm -rf $SWDIR
fi

OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/obssim.conf_common
SCALE_CONF=${LETKF_RUN}/config/${EXP}/config.nml.scale
TOPO=${OUTDIR}/const/topo


# 15 km
#tstart='2016-06-08 0:00:00'
#tend='2016-06-08 0:00:00'

# 3 km
tstart='2016-06-01 12:00:00'
#tstart='2016-06-03 12:00:00'
#tend='2016-06-03 12:00:00'

tend=$tstart

OBSSIM_TIME_START=2
#OBSSIM_TIME_END=5
#OBSSIM_TIME_START=3
OBSSIM_TIME_START=5
#OBSSIM_TIME_START=1
OBSSIM_TIME_END=$OBSSIM_TIME_START


ELAPSE="0:30:00"
#OBSSIM_TIME_START=4
#OBSSIM_TIME_END=6

#OBSSIM_TIME_START=7
#OBSSIM_TIME_END=9

#OBSSIM_TIME_START=10
#OBSSIM_TIME_END=12

#ctint=21600 # obssim interval 
ctint=600 # obssim interval 
tint=600 # analysis interval (Do not modify!)
ctint=150 # obssim interval 
tint=150 # analysis interval (Do not modify!)

# -- SCALE setting --
MEM_NP=${SCALE_NP}
MEM=mean
TYPE=anal
TYPE=fcst
#TYPE=gues

#OBSSIM_IN_TYPE=restart
#FHEAD=init
OBSSIM_IN_TYPE=history
FHEAD=history

SMEM=1 # 
EMEM=${SMEM} # mean # DEBUG
#EMEM=1 
SMEM=1
EMEM=1
MEM_L=`seq ${SMEM} ${EMEM}`



NNODES=$(( SCALE_NP / PPN))
if (( NNODES < 1 )) ; then
  NNODES=1
fi

MEM_TOTAL=$((EMEM - SMEM + 1))
echo $MEM_TOTAL

NP_TOTAL=$((MEM_TOTAL * SCALE_NP))
NNODE_TOTAL=$((NP_TOTAL / PPN))
if [ $NNODE_TOTAL -lt 1 ] ; then
  NNODE_TOTAL=1
fi
echo $NP_TOTAL" "$NNODE_TOTAL


# -- Him8 DA (RTTOV) setting --

if [ ! -e ${SWDIR} ] ; then
   mkdir -p $SWDIR
fi


ctime="$tstart"
YYYYh=$(date -ud "$ctime" '+%Y')
MMh=$(date -ud "$ctime" '+%m')
DDh=$(date -ud "$ctime" '+%d')
HHh=$(date -ud "$ctime" '+%H')
MNh=$(date -ud "$ctime" '+%M')
SEh=$(date -ud "$ctime" '+%S')

rm -f $RUNCONF_COMMON
cat << EOF >> $RUNCONF_COMMON
&PARAM_ENSEMBLE
! MEMBER = 100,
/
&PARAM_LETKF_PRC
 NNODES = ${NNODES},
 PPN = ${PPN},
 MEM_NODES = 1,
 MEM_NP = ${SCALE_NP},
/

&PARAM_LETKF_H08
 H08_RTTOV_CLD = .true.,
 H08_RTTOV_MINQ = 0.0d0,
 H08_RTTOV_CFRAC_CNST = -0.1d0,
!--H08_RTTOV_COEF_PATH--!
!--H08_NOWDATE--!
!--H08_RTTOV_ITMAX--
!--H08_RTTOV_NTHREAD--
/

EOF


if [ ! -e ${SWDIR}/dat/topo ] ; then
  mkdir -p ${SWDIR}/dat/topo
fi

#-- copy topo

#cp ${TOPO}/topo.pe*[0,1].nc ${SWDIR}/dat/topo/ &
#cp ${TOPO}/topo.pe*[2,3].nc ${SWDIR}/dat/topo/ &
#cp ${TOPO}/topo.pe*[4,5].nc ${SWDIR}/dat/topo/ &
#cp ${TOPO}/topo.pe*[6,7].nc ${SWDIR}/dat/topo/ &

#RAD_DAT=${SWDIR}/dat/rad
#rm -rf $RAD_DAT
#cp -r ${SCALEDIR}/scale-rm/test/data/rad ${RAD_DAT}

# -- Add header for run sh --
TS=$(expr `date -ud "$tstart" '+%s'` - `date -ud "$tint" '+%s'` )
TE=$(expr `date -ud "$tend" '+%s'` - `date -ud "$tint" '+%s'` )
TS=$(expr $TS / $tint + 1 )
TE=$(expr $TE / $tint + 1 )
TLEV=$(($TE - TS + 1))
echo $TLEV
#exit
TNODE=`expr ${MEM_NP} \* $TLEV`

TNODE_CNT=0
VCODE_CNT=0

NPIN=`expr 255 / \( $PPN \) + 1`
rm -f $RUNSH
cat << EOF >> $RUNSH
#!/bin/sh
#PJM -L rscgrp=${RSCGRP}
#PJM -N VIS_OBSSIM
#PJM -L node=$NNODE_TOTAL
#PJM -L elapse=00:30:00
#PJM --mpi proc=$NP_TOTAL
#PJM --omp thread=${RTTOV_THREADS}
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
#export KMP_AFFINITY=verbose
#export I_MPI_DEBUG=5

#export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t


#export OMP_STACKSIZE=128m
ulimit -s unlimited


EOF


ctime="$tstart"
STIME=$(date -ud "$ctime" '+%Y%m%d%H%M%S')

if (( MICRO != 1 )) ; then
  echo "#PJM --stgin-dir \"${SWDIR} ./ recursive=10\"" >> ${RUNSH}
  echo "#PJM --stgout-dir \"./out/${EXP}/${TYPE} $OUTDIR/${STIME}/obssim recursive=10\"" >> ${RUNSH}
#  echo "#PJM --stgout-dir \"./out ${SWDIR}/out recursive=10\"" >> ${RUNSH}
fi


#-- copy init file


ctime="$tstart"
while (($(date -ud "$ctime" '+%s') <= $(date -ud "$tend" '+%s'))); do # -- time

  YYYYh=$(date -ud "$ctime" '+%Y')
  MMh=$(date -ud "$ctime" '+%m')
  DDh=$(date -ud "$ctime" '+%d')
  HHh=$(date -ud "$ctime" '+%H')
  MNh=$(date -ud "$ctime" '+%M')
  SEh=$(date -ud "$ctime" '+%S')

  HTIME=${YYYYh}${MMh}${DDh}${HHh}${MNh}${SEh}

  for MEM in $MEM_L # MEM
  do

    TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
    VCODE_CNT=$(expr ${VCODE_CNT} + 1)   


  if [ $MEM == 0 ] ; then
     MEM=mean
  else
     MEM=$(printf '%04d' $MEM)
  fi

  if (( MICRO == 1 )) ; then
#    ONAME=${SWDIR}/out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}.dat
    ONAME=${SWDIR}/out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}
  else
#    ONAME=./out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}.dat
    ONAME=./out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}
  fi

  #-- copy bin & RTTOV coef files

  ORG_DIR=${INDIR}/${HTIME}/${TYPE}/${MEM}
  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
  mkdir -p $DAT_DIR
  echo $HTIME" "$MEM

  if [ ! -e ${DAT_DIR}/vis ] ; then 
    mkdir -p ${DAT_DIR}/vis
  fi

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
    cp ${RTTOV_COEF} ${DAT_DIR}/
    cp ${RTTOV_SCCOEF} ${DAT_DIR}/

    cp ${RTTOV_COEF_VIS} ${DAT_DIR}/vis
    cp ${RTTOV_SCCOEF_VIS} ${DAT_DIR}/vis
    cp ${RTTOV_MFCOEF} ${DAT_DIR}/vis
  fi

  if [ ! -e ${SWDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
    mkdir -p ${SWDIR}/out/${EXP}/${TYPE}
#    touch ${SWDIR}/out/${EXP}/${TYPE}/tmp
   
  fi
  if [ ! -e ${DAT_DIR} ] ; then
    mkdir -p ${DAT_DIR}
  fi

  cp ${ORG_DIR}/${FHEAD}.pe*[0,1].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/${FHEAD}.pe*[2,3].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/${FHEAD}.pe*[4,5].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/${FHEAD}.pe*[6,7].nc ${DAT_DIR}/ 
  wait

  # copy common parts of obssim.conf 
#  RUNCONF=${SWDIR}/obssim_$(printf %03d $VCODE_CNT).conf
  RUNCONF=${SWDIR}/obssim.conf.$VCODE_CNT
  if (( MICRO == 1 )) ; then
    LRUNCONF=${SWDIR}/obssim_$(printf %03d $VCODE_CNT).conf
    LDAT_DIR=$DAT_DIR
    LDAT_TOPO=$DAT_DIR
  else
    LRUNCONF=./obssim_$(printf %03d $VCODE_CNT).conf
    LDAT_DIR=./dat/${EXP}/${HTIME}/${TYPE}/${MEM}
    LDAT_TOPO=./dat
  fi

  rm -f $RUNCONF
  cat ${RUNCONF_COMMON} | \
     sed -e "/!--H08_NOWDATE--/a H08_NOWDATE = $YYYYh, $MMh, $DDh, $HHh, $MNh, $SEh," \
         -e "/!--H08_RTTOV_COEF_PATH--/a H08_RTTOV_COEF_PATH = \"${LDAT_DIR}\"," \
         -e "/!--H08_RTTOV_ITMAX--/a H08_RTTOV_ITMAX = ${RTTOV_ITMAX}," \
         -e "/!--H08_RTTOV_NTHREAD--/a H08_RTTOV_NTHREAD = ${RTTOV_THREADS}," \
    >> ${RUNCONF}

cat << EOF >> $RUNCONF
&PARAM_OBSSIM
 OBSSIM_IN_TYPE = "${OBSSIM_IN_TYPE}",
 OBSSIM_RESTART_IN_BASENAME = "${LDAT_DIR}/init",
 OBSSIM_HISTORY_IN_BASENAME = "${LDAT_DIR}/history",
 OBSSIM_TOPO_IN_BASENAME = "${LDAT_TOPO}/topo/topo",
 OBSSIM_TIME_START = ${OBSSIM_TIME_START},
 OBSSIM_TIME_END = ${OBSSIM_TIME_END},
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
! OBSSIM_NUM_3D_VARS = 0,
! OBSSIM_3D_VARS_LIST = 4001, 4002, 2819, 2820
! OBSSIM_NUM_2D_VARS = 3,
! OBSSIM_2D_VARS_LIST = 8800,8800,8800,
! OBSSIM_NUM_2D_VARS = 6,
! OBSSIM_2D_VARS_LIST = 8800,8800,8800,8800,8800,8800,
 OBSSIM_NUM_2D_VARS = 16,
 OBSSIM_2D_VARS_LIST = 8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,
/
EOF


  # Add SCALE config
  cat $SCALE_CONF >> $RUNCONF
#  cat $SCALE_CONF |\
#     sed -e "/!--H08_NOWDATE--/a H08_NOWDATE = $YYYYh, $MMh, $DDh, $HHh, $MNh, $SSh," \
#         -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${RAD_DAT}/PARAG.29\"," \
#         -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${RAD_DAT}/PARAPC.29\"," \
#         -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${RAD_DAT}/VARDATA.RM29\"," \
#         -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${RAD_DAT}/cira.nc\"," \
#         -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${RAD_DAT}/MIPAS\"," \
#     >> $RUNCONF



#${SWDIR}/dat


  VCODE_NUM=$(printf %03d $VCODE_CNT)
  if (( MICRO == 1 )) ; then
    BIN=${SWDIR}/obssim
    LVCODE=${SWDIR}/vcode${VCODE_NUM}
  else
    BIN=./obssim
    LVCODE=./vcode${VCODE_NUM}
  fi
  VCODE=${SWDIR}/vcode${VCODE_NUM}

#  rm -f ${VCODE}
#  touch $VCODE

  VNODE_MAX=$(expr ${TNODE_CNT} + ${MEM_NP} - 1)

#  for VN in `seq ${TNODE_CNT} ${VNODE_MAX}`
#  do
#    echo "("${VN}")" >> $VCODE
#  done

#  echo "mpiexec.hydra -np ${MEM_NP} ${BIN} ${LRUNCONF} &" >> $RUNSH


  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "mpiexec.hydra -np ${MEM_NP} ${BIN} obssim.conf.\$PJM_BULKNUM " >> $RUNSH

echo "mv ./out/${EXP}/${TYPE}/Him8_*.nc $OUTDIR/ " >> $RUNSH

sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH

cd $SWDIR > /dev/null
echo pjsub --bulk --sparam 1-${VCODE_CNT} OBSSIM.sh 
pjsub --bulk --sparam 1-${VCODE_CNT} OBSSIM.sh 
cd - > /dev/null
echo ${SWDIR}

exit

