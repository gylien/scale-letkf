#!/bin/sh
set -e

SCALE_V=scale-5.1.1
EXP="JAXA2_D2_BDYENS_GPR" #Finish

# Reference
EXP_HIM8=TC201513_D2_PREP_0519
# !! [!--IO_LOG_BASENAME--] should no be specified !!
#

#
SWDIR=/scratch/hp120282/k03117/obssim
OTOP=/volume61/data/hp120282/k03117/scale
DTOP=/volume92/data/hp120282_A/k03117/scale
LETKF_RUN=${OTOP}/${SCALE_V}/letkf/scale/run
OBSSIM_BIN=${OTOP}/${SCALE_V}/letkf/scale_z_181109/obs/obssim
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/obssim.conf_common
SCALE_CONF=${DTOP}/out/${EXP}/20150801130000/log/scale/mean_run.conf
TOPO=${DTOP}/out/JAXA2_D2_BDYENS_TC/const/topo

#
tstart='2015-08-01 18:00:00'
tend='2015-08-01 18:00:00'

ctint=3600 # obssim interval 
tint=600 # analysis interval (Do not modify!)

# -- SCALE setting --
MEM_NP=96
MEM=mean
#TYPE=anal
TYPE=gues
#TYPE=hist

#SMEM=0 # 
#EMEM=${SMEM} # mean
SMEM=1
EMEM=50
#SMEM=13
#EMEM=24
#SMEM=25
#EMEM=36
#SMEM=37
#EMEM=48
#SMEM=49
#EMEM=50
MEM_L=`seq ${SMEM} ${EMEM}`


# -- Him8 DA (RTTOV) setting --

rm -rf $SWDIR
if [ ! -e ${SWDIR} ] ; then
   mkdir -p $SWDIR
fi

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

&PARAM_LETKF_GPR
 MIN_GPR_REF_MEMBER = 1,
 MIN_GPR_REF_MEMBER_OBSREF = 1,
 MIN_GPR_REF_DBZ = 20.0D0,
 LOW_GPR_REF_SHIFT = -5.0D0,
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

# -- Add header for run sh --
###tmp TS=$(expr `date -ud "$tstart" '+%s'` - `date -ud "$tint" '+%s'` )
###tmp TE=$(expr `date -ud "$tend" '+%s'` - `date -ud "$tint" '+%s'` )
###tmp TS=$(expr $TS / $tint + 1 )
###tmp TE=$(expr $TE / $tint + 1 )
###tmp TLEV=$(($TE - TS + 1))
###tmp echo $TS" "$TE" "$TLEV
###tmp TNODE=`expr ${MEM_NP} \* $TLEV`

TNODE_CNT=0
VCODE_CNT=1

rm -f $RUNSH
cat << EOF >> $RUNSH
#!/bin/sh
#PJM -N Him8_OBSSIM
#PJM -s
#PJM --rsc-list "node=<TNODE_CNT>"
#PJM --rsc-list "elapse=1:00:00"
#PJM --rsc-list "rscgrp=small"
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
  SEh=00

  HTIME=${YYYYh}${MMh}${DDh}${HHh}${MNh}00

  for MEM in $MEM_L # MEM
  do


  ONAME=${SWDIR}/out/${EXP}/${TYPE}/gpr_${HTIME}.dat
  #ONAME=${SWDIR}/out/${EXP}/${TYPE}/${MEM}/gpr_${HTIME}.dat
  
  if [ $MEM == 0 ] ; then
     MEM=mean
  else
     MEM=$(printf '%04d' $MEM)
     ONAME=${SWDIR}/out/${EXP}/${TYPE}/gpr_${HTIME}_${MEM}.dat
  fi

  #-- copy bin & RTTOV coef files

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
  fi


  ORG_DIR=${DTOP}/out/${EXP}/${HTIME}/${TYPE}/${MEM}
  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
  echo $HTIME

  if [ ! -e ${SWDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
    mkdir -p ${SWDIR}/out/${EXP}/${TYPE}
  fi
  if [ ! -e ${DAT_DIR} ] ; then
    mkdir -p ${DAT_DIR}
  fi

  cp ${ORG_DIR}/init.pe*[0,1].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/init.pe*[2,3].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/init.pe*[4,5].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/init.pe*[6,7].nc ${DAT_DIR}/ &
  cp ${ORG_DIR}/init.pe*[8,9].nc ${DAT_DIR}/ 
  wait

  # copy common parts of obssim.conf 
  RUNCONF=${SWDIR}/obssim_$(printf %03d $VCODE_CNT).conf
  cp $RUNCONF_COMMON $RUNCONF

cat << EOF >> $RUNCONF
&PARAM_OBSSIM
 OBSSIM_IN_TYPE = "restart",
 OBSSIM_RESTART_IN_BASENAME = "${DAT_DIR}/init",
 OBSSIM_HISTORY_IN_BASENAME = "${DAT_DIR}/history",
 OBSSIM_TOPO_IN_BASENAME = "${SWDIR}/dat/topo/topo",
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
 OBSSIM_TIME_START=1,
 OBSSIM_TIME_END=1,
 OBSSIM_NUM_3D_VARS = 1,
 OBSSIM_3D_VARS_LIST = 5001
/

EOF

  # Add SCALE config
  cat $SCALE_CONF >> $RUNCONF

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

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H%M')
done # -- time
echo "wait" >> $RUNSH


sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH

echo ${SWDIR}

exit



