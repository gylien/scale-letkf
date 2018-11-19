#!/bin/bash


USER=honda

EXP=8km_sc
. config/${EXP}/config.main.hakushu

#
LETKF_RUN="$(pwd)"

#SWDIR="/scratch/$(id -ng)/${USER}/obssim"
#SWDIR=${TMPL}
SWDIR=${LETKF_RUN}
OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$SWDIR/OBSSIM.sh
RUNCONF_COMMON=$SWDIR/OBSSIM.conf_common
SCALE_CONF=${LETKF_RUN}/config.nml.scale
TOPO=${OUTDIR}/const/topo



tstart='2000-01-01 0:00:00'
tend='2000-01-01 0:00:00'


ctint=600 # obssim interval 
tint=600 # analysis interval (Do not modify!)

# -- SCALE setting --
MEM_NP=${SCALE_NP}
MEM=mean
TYPE=fcst

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


&PARAM_LETKF_RADAR
 USE_OBSERR_RADAR_REF = .true.
 USE_OBSERR_RADAR_VR = .true.
 RADAR_REF_THRES_DBZ = 10.0D0,
 MIN_RADAR_REF_MEMBER = 20,
 MIN_RADAR_REF_MEMBER_OBSREF = 1,
 MIN_RADAR_REF_DBZ = 10.0D0,
 LOW_REF_SHIFT = -5.0D0,
 RADAR_ZMAX = 11.0D3,
/

&PARAM_OBS_ERROR
 OBSERR_RADAR_REF = 5.0D0,
 OBSERR_RADAR_VR = 3.0D0,
/

EOF

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

echo $RUNSH
rm -f $RUNSH
cat > $RUNSH << EOF
#!/bin/sh
#PBS -l nodes=1:ppn=${SCALE_NP}
#PBS -W umask=027
##PBS -k oe

ulimit -s unlimited

HOSTLIST=\$(cat \$PBS_NODEFILE | sort | uniq)
HOSTLIST=\$(echo \$HOSTLIST | sed 's/  */,/g')
export MPI_XPMEM_ENABLED=disabled
export MPI_UNIVERSE="\$HOSTLIST $((PPN*THREADS))"
export MPI_XPMEM_ENABLED=disabled

export OMP_NUM_THREADS=${THREADS}
#export PARALLEL=${THREADS}

export FORT_FMT_RECL=400

cd \$PBS_O_WORKDIR

rm -f machinefile
cp -f \$PBS_NODEFILE machinefile

export RUN_LEVEL=1

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
  ORG_DIR=${OUTDIR}/${HTIME}/${TYPE}/${MEM}

  ONAME=${ORG_DIR}/radar_${HTIME}_${MEM}.dat

  #-- copy bin & RTTOV coef files

#  DAT_DIR=${SWDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
#  mkdir -p $DAT_DIR
#  echo $HTIME
#
#  if [ ! -e ${SWDIR}/${OBSSIM_BIN} ] ; then 
#    cp ${OBSSIM_BIN} ${SWDIR}/
#  fi
#
#  if [ ! -e ${SWDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
#    mkdir -p ${SWDIR}/out/${EXP}/${TYPE}
#  fi
#  if [ ! -e ${DAT_DIR} ] ; then
#    mkdir -p ${DAT_DIR}
#  fi
#
#  cp ${ORG_DIR}/history.pe*[0,1].nc ${DAT_DIR}/ &
#  cp ${ORG_DIR}/history.pe*[2,3].nc ${DAT_DIR}/ &
#  cp ${ORG_DIR}/history.pe*[4,5].nc ${DAT_DIR}/ &
#  cp ${ORG_DIR}/history.pe*[6,7].nc ${DAT_DIR}/ &
#  cp ${ORG_DIR}/history.pe*[8,9].nc ${DAT_DIR}/ 

#  echo $DAT_DIR
#  wait

  # copy common parts of obssim.conf 
  RUNCONF=${SWDIR}/OBSSIM_$(printf %03d $VCODE_CNT).conf
  rm -f $RUNCONF

cat ${RUNCONF_COMMON} > ${RUNCONF}

cat << EOF >> $RUNCONF

&PARAM_OBSSIM
 OBSSIM_IN_TYPE = "history",
! OBSSIM_RESTART_IN_BASENAME = "${DAT_DIR}/init",
 OBSSIM_HISTORY_IN_BASENAME = "${ORG_DIR}/history",
 OBSSIM_TOPO_IN_BASENAME = "${SWDIR}/dat/topo/topo",
 OBSSIM_TIME_START = 1,
 OBSSIM_TIME_END = 7,
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
 OBSSIM_NUM_3D_VARS = 2,
 OBSSIM_3D_VARS_LIST = 4001, 4002, !2819, 2820,
! OBSSIM_NUM_2D_VARS = 10,
! OBSSIM_2D_VARS_LIST = 8800,8800,8800,8800,8800,8800,8800,8800,8800,8800,
 OBSSIM_RADAR_LON = 0.0d0,
 OBSSIM_RADAR_LAT = 0.0d0,
 OBSSIM_RADAR_Z = 0.0d0,

/
EOF

  # Add SCALE config
  cat $SCALE_CONF >> $RUNCONF



#  echo "mpirun -n ${MEM_NP} ${SWDIR}/obssim ${RUNCONF} &" >> $RUNSH
  echo "mpirun -n ${MEM_NP} ${LETKF_RUN}/../obs/obssim ${RUNCONF} &" >> $RUNSH

  TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
  VCODE_CNT=$(expr ${VCODE_CNT} + 1)   


  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "wait" >> $RUNSH


#sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH
#
#echo ${SWDIR}

qsub $RUNSH

exit



