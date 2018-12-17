#!/bin/sh

SCALE_V=scale-5.1.2

#EXP=TC201513_D2_H08_0416_3

#EXP=TC201513_D2_H08_0519_36
#EXP=TC201513_D2_H08_0519_36_30km
#EXP=TC201513_D2_H08_0519_36_30km_BC
#EXP=TC201513_D2_H08_0519_36_50km
#EXP=TC201513_D2_H08_0519_36_50km_BC
#EXP=TC201513_D2_H08_0519_36_60km_BC
#EXP=TC201513_D2_PREP_0519

#EXP=TC201513_D2_H08_0519_36_60km_30Him8
EXP=TC201513_D2_H08_0519_36_60km #Finish
#EXP=TC201513_D2_PREP_0519 # Finish

#EXP=TC201513_D2_H08_0519_DERO_36_60km_B891013

# Reference
EXP_HIM8=TC201513_D2_PREP_0519
# !! [!--IO_LOG_BASENAME--] should no be specified !!
#

#
SWDIR=/scratch/hp150019/honda/obssim
OTOP=/volume63/data/hp150019/honda/SCALE-LETKF/${SCALE_V}
LETKF_RUN=${OTOP}/${SCALE_V}/letkf/scale/run
OBSSIM_BIN=${OTOP}/${SCALE_V}/letkf/scale/obs/obssim
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/obssim.conf_common
SCALE_CONF=${OTOP}/OUTPUT/${EXP_HIM8}/run.conf
TOPO=${OTOP}/OUTPUT/${EXP_HIM8}/const/topo

# -- RTTOV_DIR --
DIR_RTTOV=/data/share005/honda/RTTOV
RTTOV_COEF=${DIR_RTTOV}/rtcoef_rttov11/rttov7pred54L/rtcoef_himawari_8_ahi.dat
RTTOV_SCCOEF=${DIR_RTTOV}/rtcoef_rttov11/cldaer/sccldcoef_himawari_8_ahi.dat



#tstart='2015-08-02 3:00:00'
#tend='2015-08-02 18:00:00'
tstart='2015-08-02 12:10:00'
tend='2015-08-02 12:10:00'
#tstart='2015-08-02 22:40:00'
#tend='2015-08-02 23:50:00'


#ctint=21600 # obssim interval 
#ctint=10800 # obssim interval 
ctint=600 # obssim interval 
tint=600 # analysis interval (Do not modify!)

# -- SCALE setting --
MEM_NP=72
MEM=mean
TYPE=anal
#TYPE=gues

SMEM=0 # 
EMEM=${SMEM} # mean
#SMEM=1
#EMEM=13
#SMEM=14
#EMEM=26
#SMEM=27
#EMEM=39
#SMEM=40
#EMEM=50
MEM_L=`seq ${SMEM} ${EMEM}`



# -- Him8 DA (RTTOV) setting --

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

&PARAM_LETKF_H08
 H08_RTTOV_MINQ = 0.10D0,
 H08_RTTOV_CLD = .true.,
 H08_RTTOV_MINQ = 0.10d0,
 H08_RTTOV_MINQ_CTOP = 0.10d0,
 H08_RTTOV_CFRAC_CNST = 0.1d0,
 H08_LIMIT_LEV = 200.0d2,
 H08_RTTOV_EXTRA_US76 = .false.,
 H08_VLOCAL_CTOP = .true.,
 H08_CH_USE = 0,0,1,0,0,0,0,0,0,0,
! H08_BT_MIN = 200.0d0, ! turn off? Apply only first 1 hour (6 cycles)?
 H08_BIAS_SIMPLE = .true.,
 H08_CLDERR_SIMPLE = .true.,
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
TS=$(expr `date -ud "$tstart" '+%s'` - `date -ud "$tint" '+%s'` )
TE=$(expr `date -ud "$tend" '+%s'` - `date -ud "$tint" '+%s'` )
TS=$(expr $TS / $tint + 1 )
TE=$(expr $TE / $tint + 1 )
TLEV=$(($TE - TS + 1))
echo $TS" "$TE" "$TLEV
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
  SEh=00

  HTIME=${YYYYh}${MMh}${DDh}${HHh}${MNh}00

  for MEM in $MEM_L # MEM
  do


  ONAME=${SWDIR}/out/${EXP}/${TYPE}/${MEM}/Him8_${HTIME}.dat
  

  if [ $MEM == 0 ] ; then
     MEM=mean
  else
     MEM=$(printf '%04d' $MEM)
     ONAME=${SWDIR}/out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}.dat
  fi

  #-- copy bin & RTTOV coef files

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
    cp ${RTTOV_SCCOEF} ${SWDIR}/
    cp ${RTTOV_COEF} ${SWDIR}/
  fi


  ORG_DIR=${OTOP}/OUTPUT/${EXP}/${HTIME}/${TYPE}/${MEM}
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
! OBSSIM_HISTORY_IN_BASENAME = "history",
 OBSSIM_TOPO_IN_BASENAME = "${SWDIR}/dat/topo/topo",
! OBSSIM_TIME_START = 1,
! OBSSIM_TIME_END = 1,
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
! OBSSIM_NUM_3D_VARS = 0,
! OBSSIM_3D_VARS_LIST = 4001, 4002, 2819, 2820
 OBSSIM_NUM_2D_VARS = 10,
 OBSSIM_2D_VARS_LIST = 8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,
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



