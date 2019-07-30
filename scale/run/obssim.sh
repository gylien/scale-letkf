#!/bin/bash

MICRO=0 # 1: micro, /=1: large

USER=honda

#EXP=test15km
EXP=test03km
. config/${EXP}/config.main

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
SCALE_CONF=${LETKF_RUN}/config.nml.scale
TOPO=${OUTDIR}/const/topo


# 15 km
#tstart='2016-06-08 0:00:00'
#tend='2016-06-08 0:00:00'

# 3 km
tstart='2016-06-03 12:00:00'
tend='2016-06-03 12:00:00'

OBSSIM_TIME_START=2
OBSSIM_TIME_END=5

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

SMEM=0 # 
EMEM=${SMEM} # mean # DEBUG
EMEM=1 
SMEM=2
EMEM=50
MEM_L=`seq ${SMEM} ${EMEM}`



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
 NNODES = ${MEM_NP},
 PPN = 1,
 MEM_NODES = ${MEM_NP},
 MEM_NP = ${MEM_NP},
/

&PARAM_LETKF_H08
 H08_RTTOV_CLD = .true.,
 H08_RTTOV_MINQ = 0.10d0,
!--H08_RTTOV_COEF_PATH--!
!--H08_NOWDATE--!
/

!&PARAM_LETKF_H08
! H08_FORMAT_NC = .true.,
! H08_SIM_ALLG = .true.,
! H08_OBS_THIN_LEV = 2,
!!--H08_NOWDATE--!
! H08_RTTOV_CLD = .true.,
! H08_RTTOV_MINQ_CTOP = 0.10d0,
! H08_RTTOV_CFRAC_CNST = 0.1d0,
! H08_RTTOV_CFRAC = 1, ! scale method for cldfrac
! H08_LIMIT_LEV = 200.0d2,
! H08_VLOCAL_CTOP = .true.,
! H08_RTTOV_KADD = 5,
! H08_RTTOV_PROF_SHIFT = .true.,
!/
EOF


if [ ! -e ${SWDIR}/dat/topo ] ; then
  mkdir -p ${SWDIR}/dat/topo
fi

#-- copy topo

cp ${TOPO}/topo.pe*[0,1].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[2,3].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[4,5].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[6,7].nc ${SWDIR}/dat/topo/ &

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
VCODE_CNT=1

rm -f $RUNSH
cat << EOF >> $RUNSH
#!/bin/sh
#PJM -N Him8_OBSSIM
#PJM -s
#PJM --rsc-list "node=<TNODE_CNT>"
#PJM --rsc-list "elapse=${ELAPSE}"
#PJM --rsc-list "rscgrp=${MODE}"
#PJM --stg-transfiles all
EOF


ctime="$tstart"
STIME=$(date -ud "$ctime" '+%Y%m%d%H%M%S')

if (( MICRO != 1 )) ; then
  echo "#PJM --stgin-dir \"${SWDIR} ./ recursive=10\"" >> ${RUNSH}
  echo "#PJM --stgout-dir \"./out/${EXP}/${TYPE} $OUTDIR/${STIME}/obssim recursive=10\"" >> ${RUNSH}
#  echo "#PJM --stgout-dir \"./out ${SWDIR}/out recursive=10\"" >> ${RUNSH}
fi

cat << EOF >> $RUNSH
. /work/system/Env_base
export F_UFMTENDIAN=big
export OMP_NUM_THREADS=8
export PARALLEL=8

EOF

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

  ORG_DIR=${OUTDIR}/${HTIME}/${TYPE}/${MEM}
  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
  mkdir -p $DAT_DIR
  echo $HTIME" "$MEM

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
    cp ${RTTOV_COEF} ${DAT_DIR}/
    cp ${RTTOV_SCCOEF} ${DAT_DIR}/

    cp ${RTTOV_COEF_VIS} ${DAT_DIR}/
    cp ${RTTOV_SCCOEF_VIS} ${DAT_DIR}/
    cp ${RTTOV_MFCOEF} ${DAT_DIR}/
  fi

  if [ ! -e ${SWDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
    mkdir -p ${SWDIR}/out/${EXP}/${TYPE}
    touch ${SWDIR}/out/${EXP}/${TYPE}/tmp
   
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
  RUNCONF=${SWDIR}/obssim_$(printf %03d $VCODE_CNT).conf
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

  rm -f ${VCODE}
  touch $VCODE

  VNODE_MAX=$(expr ${TNODE_CNT} + ${MEM_NP} - 1)

  for VN in `seq ${TNODE_CNT} ${VNODE_MAX}`
  do
    echo "("${VN}")" >> $VCODE
  done

  echo "mpirun -n ${MEM_NP} --vcoordfile ${LVCODE}  ${BIN} ${LRUNCONF} &" >> $RUNSH

  TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
  VCODE_CNT=$(expr ${VCODE_CNT} + 1)   


  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "wait" >> $RUNSH


sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH

echo ${SWDIR}

exit

