#!/bin/bash

TIME_LIMIT="00:30:00"
USER=honda

EXP=test_Him8_hist_ref
EXP=TAIWAN201808_D2_NOHIM8_0116_RTTOV_TEST
. config/${EXP}/config.main.ofp

#
LETKF_RUN="$(pwd)"
SWDIR="$LETKF_RUN/../tmp/obssim"
OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/obssim.conf_common
SCALE_CONF=${LETKF_RUN}/config.nml.scale
TOPO=${OUTDIR}/const/topo



#tstart='2018-06-24 0:10:00'
#tend='2018-06-24 0:10:00'

tstart='2018-08-21 12:00:00'
tend=$tstart

#ctint=21600 # obssim interval 
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

!&PARAM_LETKF_PRC
! NNODES = ${MEM_NP},
! PPN = ${PPN},
! MEM_NODES = ${MEM_NP},
! MEM_NP = ${MEM_NP},
!/

&PARAM_LOG
 LOG_LEVEL = 2,
/

&PARAM_PROCESS
 PPN = ${PPN},
 MEM_NODES = ${MEM_NP},
/

&PARAM_LETKF_H08
 H08_FORMAT_NC = .true.,
 H08_SIM_ALLG = .true.,
 H08_OBS_THIN_LEV = 2,
!--H08_NOWDATE--!
!--H08_RTTOV_COEF_PATH--!
 H08_RTTOV_CLD = .true.,
 H08_RTTOV_MINQ_CTOP = 0.10d0,
 H08_RTTOV_CFRAC_CNST = 0.1d0,
 H08_RTTOV_CFRAC = 1, ! scale method for cldfrac
 H08_LIMIT_LEV = 200.0d2,
 H08_VLOCAL_CTOP = .true.,
 H08_RTTOV_KADD = 5,
 H08_RTTOV_PROF_SHIFT = .true.,
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
#PJM -L rscgrp=${RSCGRP}
#PJM -L node=${NNODES}
#PJM -L elapse=${TIME_LIMIT}
#PJM --mpi proc=$((PPN))
#PJM --omp thread=${THREADS}
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

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t


#export OMP_STACKSIZE=128m
ulimit -s unlimited


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
  ONAME=${SWDIR}/out/${EXP}/${TYPE}/Him8_${HTIME}_${MEM}.dat

  #-- copy bin & RTTOV coef files

  ORG_DIR=${OUTDIR}/${HTIME}/${TYPE}/${MEM}
  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
  mkdir -p $DAT_DIR
  echo $HTIME

  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
    cp ${OBSSIM_BIN} ${SWDIR}/
    cp ${RTTOV_SCCOEF} ${DAT_DIR}/
    cp ${RTTOV_COEF} ${DAT_DIR}/
  fi

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
  rm -f $RUNCONF
  cat ${RUNCONF_COMMON} | \
     sed -e "/!--H08_NOWDATE--/a H08_NOWDATE = $YYYYh, $MMh, $DDh, $HHh, $MNh, $SSh," \
         -e "/!--H08_RTTOV_COEF_PATH--/a H08_RTTOV_COEF_PATH = \"${DAT_DIR}\"," \
    >> ${RUNCONF}

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
  cat $SCALE_CONF |\
     sed -e "/!--H08_NOWDATE--/a H08_NOWDATE = $YYYYh, $MMh, $DDh, $HHh, $MNh, $SSh," \
         -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${RAD_DAT}/PARAG.29\"," \
         -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${RAD_DAT}/PARAPC.29\"," \
         -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${RAD_DAT}/VARDATA.RM29\"," \
         -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${RAD_DAT}/cira.nc\"," \
         -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${RAD_DAT}/MIPAS\"," \
     >> $RUNCONF

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

  echo "mpiexec.hydra -n ${MEM_NP}  ${SWDIR}/obssim ${RUNCONF} &" >> $RUNSH

  TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
  VCODE_CNT=$(expr ${VCODE_CNT} + 1)   


  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "wait" >> $RUNSH

TNODE_CNT=$((SCALE_NP/PPN))
echo $TNODE_CNT
sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH

echo ${SWDIR}

cd $SWDIR > /dev/null
pjsub OBSSIM.sh
cd - > /dev/null

exit



