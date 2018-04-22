#!/bin/bash


USER=honda

EXP=KYUSHU2017_D2_NoHim8_2.5min
. config/${EXP}/config.main

#
SWDIR="/scratch/$(id -ng)/${USER}/obssim"
LETKF_RUN="$(pwd)"
OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/obssim.conf_common
SCALE_CONF=${LETKF_RUN}/config.nml.scale
TOPO=${OUTDIR}/const/topo

OBSDIR=$OBS

tstart='2017-07-04 12:00:00'
tend='2017-07-04 12:00:00'


#ctint=21600 # obssim interval 
#ctint=10800 # obssim interval 
ctint=600 # obssim interval 
tint=600 # analysis interval (Do not modify!)
ctint=150 # obssim interval 
tint=150 # analysis interval (Do not modify!)

# -- SCALE setting --
MEM_NP=${SCALE_NP}
MEM=mean
TYPE=anal
#TYPE=gues

SMEM=0 # 
EMEM=${SMEM} # mean
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


EOF


if [ ! -e ${SWDIR}/dat/topo ] ; then
  mkdir -p ${SWDIR}/dat/topo
fi

#-- copy topo
cp ${TOPO}/topo.pe*[0,1].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[2,3].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[4,5].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[6,7].nc ${SWDIR}/dat/topo/ &
cp ${TOPO}/topo.pe*[8,9].nc ${SWDIR}/dat/topo/ 
wait

RAD_DAT=${SWDIR}/dat/rad
rm -rf $RAD_DAT
cp -r ${SCALEDIR}/scale-rm/test/data/rad ${RAD_DAT}

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
#PJM -N FSS_OBSSIM
#PJM -s
#PJM --rsc-list "node=<TNODE_CNT>"
#PJM --rsc-list "elapse=0:10:00"
#PJM --rsc-list "rscgrp=micro"
#PJM --stg-transfiles all

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
  ONAME=${SWDIR}/out/${EXP}/${TYPE}/JMAFSS_${HTIME}_${MEM}.dat

  #-- copy bin & RTTOV coef files

  ORG_DIR=${OUTDIR}/${HTIME}/${TYPE}/${MEM}
  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
  mkdir -p $DAT_DIR
  echo $HTIME

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
  fi

  if [ ! -e ${SWDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
    mkdir -p ${SWDIR}/out/${EXP}/${TYPE}
  fi
  if [ ! -e ${DAT_DIR} ] ; then
    mkdir -p ${DAT_DIR}
  fi

  cp ${ORG_DIR}/history.pe*[0,1].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/history.pe*[2,3].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/history.pe*[4,5].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/history.pe*[6,7].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/history.pe*[8,9].nc ${DAT_DIR}/ &
  wait

  cp ${OBSDIR}/jmaradar_${HTIME}.dat ${DAT_DIR}/ &


  # copy common parts of obssim.conf 
  RUNCONF=${SWDIR}/obssim_$(printf %03d $VCODE_CNT).conf
  rm -f $RUNCONF
#  cat ${RUNCONF_COMMON} | \
#     sed -e "/!--H08_NOWDATE--/a H08_NOWDATE = $YYYYh, $MMh, $DDh, $HHh, $MNh, $SSh," \
#         -e "/!--H08_RTTOV_COEF_PATH--/a H08_RTTOV_COEF_PATH = \"${DAT_DIR}\"," \
#    >> ${RUNCONF}

cat << EOF >> $RUNCONF
&PARAM_OBSSIM
 OBSSIM_IN_TYPE = "history",
! OBSSIM_RESTART_IN_BASENAME = "${DAT_DIR}/init",
 OBSSIM_HISTORY_IN_BASENAME = "${DAT_DIR}/history",
 OBSSIM_TOPO_IN_BASENAME = "${SWDIR}/dat/topo/topo",
 OBSSIM_TIME_START = 1,
 OBSSIM_TIME_END = 1,
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
 OBSSIM_NUM_2D_VARS = 1,
 OBSSIM_2D_VARS_LIST = 5001,
/

&PARAM_OBS_RADAR_JMA
 JMA_RADAR_FILE = "jmaradar_${HTIME}.dat"
/

EOF


  # Add SCALE config
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
  VCODE=${SWDIR}/vcode${VCODE_NUM}
  rm -f ${VCODE}
  touch $VCODE

  VNODE_MAX=$(expr ${TNODE_CNT} + ${MEM_NP} - 1)

  for VN in `seq ${TNODE_CNT} ${VNODE_MAX}`
  do
    echo "("${VN}")" >> $VCODE
  done

  echo "mpirun -n ${MEM_NP} --vcoordfile ${VCODE}  ${SWDIR}/obssim ${RUNCONF} &" >> $RUNSH

  TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
  VCODE_CNT=$(expr ${VCODE_CNT} + 1)   


  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "wait" >> $RUNSH


sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH

echo ${SWDIR}

exit



