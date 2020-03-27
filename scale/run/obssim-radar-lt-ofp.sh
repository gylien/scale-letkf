#!/bin/bash


USER=honda
SYS=ofp

OBSTYPE="RADAR"
OBSTYPE="H08"
#OBSTYPE="LT"
#OBSTYPE="FP" # Flash point
#OBSTYPE="CONV"
#OBSTYPE="ALL"
#OBSTYPE="FP" # Flash point

# Generate new obs format file
OBSSIM_OBSOUT=".false." # anal/gues
#OBSSIM_OBSOUT=".true." # fcst

H08_RTTOV_CFRAC=1

TYPE=fcst
#TYPE=hist
TYPE=anal
TYPE=gues



EXP=2000m_DA_0306_TEST_Him8
#EXP=2000m_DA_0306

. config/${EXP}/config.main.$SYS
. config/${EXP}/config.fcst

OBSSIM_RADAR_LON=180
OBSSIM_RADAR_LAT=180




tstart='2001-01-01 1:05:00'
tstart='2001-01-01 1:10:00'
tstart='2001-01-01 1:15:00'
tstart='2001-01-01 1:20:00'
tstart='2001-01-01 1:25:00'
tstart='2001-01-01 1:30:00'
#tstart='2001-01-01 1:00:00'
tend='2001-01-01 2:00:00'
tstart='2001-01-01 1:00:00'
tstart='2001-01-01 1:40:00'
tend=$tstart
#tend='2001-01-01 1:25:00'

if [ "$TYPE" == "fcst" ] || [ "$TYPE" == "hist" ]; then
  tstart='2001-01-01 1:00:00'

  FCSTLEN=1800 
  FCSTLEN=3600 
  TS=1
  TE=121


  tstart='2001-01-01 1:10:00'
  #tstart='2001-01-01 1:20:00'
  tstart='2001-01-01 1:30:00'
  #tstart='2001-01-01 1:05:00'
  tstart='2001-01-01 1:00:00'



  tend=$tstart
  TS=1
  FCSTLEN=3600 
  TE=13

  #FCSTLEN=1800 
  #TE=7

  if [ "$TYPE" == "hist" ] ; then
    TE=2
  fi

elif [ "$TYPE" == "anal" ] || [ "$TYPE" == "gues" ] ; then
  FCSTOUT=30
  FCSTLEN=30 
  FCSTLEN=600

  TS=1
  TE=1
fi

# -- SCALE setting --
MEM_NP=${SCALE_NP}
MEM=mean

SMEM=0 # 
#SMEM=320 # 
#SMEM=252 # 
EMEM=${SMEM} # mean





#--SMEM--
#--EMEM--

#
LETKF_RUN="$(pwd)"

#WDIR="/scratch/$(id -ng)/${USER}/obssim"
#WDIR=${TMPL}
WDIR=${LETKF_RUN}/../tmp_obssim_${EXP}_${SMEM}_${EMEM}
OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$WDIR/OBSSIM.sh
RUNCONF_COMMON=$WDIR/OBSSIM.conf_common
SCALE_CONF=${LETKF_RUN}/config/${EXP}/config.nml.scale
TOPO=${OUTDIR}/const/topo



#tend=$(date -ud "${FCSTLEN} second $tstart" '+%Y-%m-%d %H:%M:%S')

if [ "$OBSTYPE" = 'H08' ]; then
  PPN=32
  PPN=16
fi

#
MEM_L=`seq ${SMEM} ${EMEM}`

#
LETKF_RUN="$(pwd)"

#WDIR="/scratch/$(id -ng)/${USER}/obssim"
#WDIR=${TMPL}
WDIR=${LETKF_RUN}/../tmp_obssim_${EXP}_${SMEM}_${EMEM}
OBSSIM_BIN="${LETKF_RUN}/../obs/obssim"
RUNSH=$WDIR/OBSSIM.sh
RUNCONF_COMMON=$WDIR/OBSSIM.conf_common
SCALE_CONF=${LETKF_RUN}/config/${EXP}/config.nml.scale
TOPO=${OUTDIR}/const/topo

# -- RTTOV_DIR --
DIR_RTTOV=/work/hp150019/share/honda/RTTOV12.3
RTTOV_COEF=${DIR_RTTOV}/rtcoef_rttov12/rttov7pred54L/rtcoef_himawari_8_ahi.dat
RTTOV_SCCOEF=${DIR_RTTOV}/rtcoef_rttov12/cldaer_ir/sccldcoef_himawari_8_ahi.dat


#tend=$(date -ud "${FCSTLEN} second $tstart" '+%Y-%m-%d %H:%M:%S')

if [ "$TYPE" == "fcst" ] || [ "$TYPE" == "hist" ] ; then
  OBSSIM_IN_TYPE="history"
  ctint=$(( FCSTLEN * 2 )) # obssim interval  # initial time loop
elif [ "$TYPE" == "anal" ] || [ "$TYPE" == "gues" ] ; then
  OBSSIM_IN_TYPE="restart"
  ctint=$LCYCLE # analysis interval (Do not modify!)
#  FCSTLEN=30 #
fi


#if [ ! -e ${WDIR} ] ; then
#fi

# clean up
rm -rf $WDIR
mkdir -p $WDIR
mkdir -p $WDIR/log

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
 NNODES = $((SCALE_NP / PPN )),
 PPN = ${PPN},
 MEM_NODES = $((SCALE_NP / PPN)),
 MEM_NP = ${SCALE_NP},
/


&PARAM_LETKF_RADAR
 USE_OBSERR_RADAR_REF = .true.
 USE_OBSERR_RADAR_VR = .true.
 RADAR_REF_THRES_DBZ = 15.0D0,
 MIN_RADAR_REF_MEMBER = 20,
 MIN_RADAR_REF_MEMBER_OBSREF = 1,
 MIN_RADAR_REF_DBZ = 5.0D0,
 LOW_REF_SHIFT = 0.0D0, ! Do not use
 RADAR_ZMAX = 31.0D3, ! Entire the domain
/

&PARAM_OBS_ERROR
 OBSERR_RADAR_REF = 5.0D0,
 OBSERR_RADAR_VR = 3.0D0,
/

&PARAM_LETKF_H08
 H08_RTTOV_COEF_PATH = "./"
 H08_RTTOV_CFRAC = ${H08_RTTOV_CFRAC},
/

EOF

#RAD_DAT=${WDIR}/dat/rad
#rm -rf $RAD_DAT
#cp -r ${SCALEDIR}/scale-rm/test/data/rad ${RAD_DAT}




# -- Add header for run sh --
#TS=1
#TE=$((FCSTLEN / tint + 1))
#if [ "$TYPE" == "anal" ] || [ "$TYPE" == "gues" ] ; then
#  TE=1
#fi


echo "TIME LEVELS: "$TS" "$TE


TNODE_CNT=0
VCODE_CNT=1

echo $RUNSH
rm -f $RUNSH
cat > $RUNSH << EOF
#!/bin/sh

#PJM -N OBSSIM
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L node=<TNODE_CNT>
##PJM -L node=$((SCALE_NP/PPN))
#PJM -L elapse=00:30:00
#PJM --mpi proc=<TPRC>
#PJM --omp thread=${THREADS}
#PJM -g hp150019

ulimit -s unlimited

export FORT_FMT_RECL=400

#rm -f machinefile
#for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
#  for ippn in \$(seq $PPN); do
#    echo "\$inode" >> machinefile
#  done
#done

module load hdf5/1.8.17
module load netcdf/4.4.1
module load netcdf-fortran/4.4.3
ulimit -s unlimited
export OMP_STACKSIZE=128m

EOF

echo ""
echo $MEM_L
echo $tstart
echo $tend

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
  OBSSIM_TOPO_IN_BASENAME=${OUTDIR}/const/topo/topo

  if [ "$OBSTYPE" = 'RADAR' ]; then
    ONAME=${ORG_DIR}/radar_${HTIME}_${MEM}.dat
    OBSSIM_NUM_3D_VARS="2"
    OBSSIM_3D_VARS_LIST="4001, 4002" 
    OBSSIM_NUM_2D_VARS="0"
    OBSSIM_2D_VARS_LIST="4001" 
    OHEAD="radar3d" 

  elif [ "$OBSTYPE" = 'FP' ]; then
    ONAME=${ORG_DIR}/fp_${HTIME}_${MEM}.dat
    OBSSIM_NUM_3D_VARS="1"
    OBSSIM_3D_VARS_LIST="5003" 
    OBSSIM_NUM_2D_VARS="1"
    OBSSIM_2D_VARS_LIST="5004" 
    OHEAD="fp" 

  elif [ "$OBSTYPE" = 'LT' ]; then
    ONAME=${ORG_DIR}/lt_${HTIME}_${MEM}.dat
    OBSSIM_NUM_3D_VARS="3"
    OBSSIM_3D_VARS_LIST="5001, 5003, 5005" 
    OBSSIM_NUM_2D_VARS="0"
    OBSSIM_2D_VARS_LIST="5002, 5004"
    OHEAD="lt" 
#  elif [ "$OBSTYPE" = 'CONV' ]; then
  elif [ "$OBSTYPE" = 'ALL' ]; then
#    ONAME=${ORG_DIR}/conv_${HTIME}_${MEM}.dat
    ONAME=${ORG_DIR}/all3d_${HTIME}_${MEM}.dat
#    OBSSIM_NUM_3D_VARS="8"
#    # u, v, w, t, p, qv, qh, qcrg
#    OBSSIM_3D_VARS_LIST="-999, -999, -999, -999, -999, -999, -999, -999"
    OBSSIM_NUM_3D_VARS="23"
    OBSSIM_3D_VARS_LIST="-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -14, -15, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25"
#    OBSSIM_3D_VARS_LIST="-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999"
#, -999, -999, -999, -999"
#, -999, -999, -999, -999" 
    OBSSIM_NUM_2D_VARS="0"
    OBSSIM_2D_VARS_LIST="4001" 
#    OHEAD="conv3d" 
    OHEAD="all3d" 
  elif [ "$OBSTYPE" = 'H08' ]; then
    ONAME=${ORG_DIR}/Him8_${HTIME}_${MEM}.dat
    OBSSIM_NUM_3D_VARS="0"
    OBSSIM_3D_VARS_LIST="1"
    OBSSIM_NUM_2D_VARS="10"
    OBSSIM_2D_VARS_LIST="8800, 8800, 8800, 8800, 8800, 8800, 8800, 8800, 8800, 8800" 
    OHEAD="Him8" 
  fi
  ONAME_OBS=${ORG_DIR}/${OHEAD}_i${HTIME}_${MEM}_${OBSSIM_RADAR_LON}_${OBSSIM_RADAR_LAT}

  #-- copy bin & RTTOV coef files

#  DAT_DIR=${WDIR}/dat/${EXP}/${HTIME}/${TYPE}/${MEM}
#  mkdir -p $DAT_DIR
#  echo $HTIME
#
#  if [ ! -e ${WDIR}/${OBSSIM_BIN} ] ; then 
#    cp ${OBSSIM_BIN} ${WDIR}/
#  fi
#
#  if [ ! -e ${WDIR}/out/${EXP}/${TYPE}/${MEM} ] ; then
#    mkdir -p ${WDIR}/out/${EXP}/${TYPE}
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

  cp $RTTOV_COEF ${WDIR}
  cp $RTTOV_SCCOEF ${WDIR}

  # copy common parts of obssim.conf 
  RUNCONF=${WDIR}/OBSSIM_$(printf %03d $VCODE_CNT).conf
  rm -f $RUNCONF

cat ${RUNCONF_COMMON} > ${RUNCONF}

cat << EOF >> $RUNCONF


&PARAM_OBSERR
 OBSERR_RADAR_REF = 5.0d0,
 OBSERR_RADAR_VR = 3.0d0,
/

&PARAM_OBSSIM
 OBSSIM_OBSOUT = ${OBSSIM_OBSOUT},
 OBSSIM_IN_TYPE = "${OBSSIM_IN_TYPE}",
 OBSSIM_RESTART_IN_BASENAME = "${ORG_DIR}/init",
 OBSSIM_HISTORY_IN_BASENAME = "${ORG_DIR}/history",
 OBSSIM_TOPO_IN_BASENAME = "${OBSSIM_TOPO_IN_BASENAME}",
 OBSSIM_TIME_START = ${TS},
 OBSSIM_TIME_END = ${TE},
 OBSSIM_TIME_INT = ${FCSTOUT},
 OBSSIM_GRADS_OUT_NAME = "${ONAME}",
 OBSSIM_NUM_3D_VARS = ${OBSSIM_NUM_3D_VARS},
 OBSSIM_3D_VARS_LIST = ${OBSSIM_3D_VARS_LIST}, 
 OBSSIM_NUM_2D_VARS = ${OBSSIM_NUM_2D_VARS},
 OBSSIM_2D_VARS_LIST = ${OBSSIM_2D_VARS_LIST},
!! ! About sqrt(20**2+20**2) km away from the storm (Similar to Zhang et al. 2004MWR)
! About sqrt(10**2+10**2) km from the storm (similar to Aksoy et al. 2009's Fig. 1a)
 OBSSIM_RADAR_LON = ${OBSSIM_RADAR_LON}.0d3,
 OBSSIM_RADAR_LAT = ${OBSSIM_RADAR_LAT}.0d3,
 OBSSIM_RADAR_Z = 0.0d0,
 OBSSIM_RADAR_ERR_10 = .true., ! Similar to Maejima et al. and Xue et al. GRL
 OBSSIM_RADAR_CLR_THIN = 2, ! Similar to Aksoy et al. 
 OBSSIM_RADAR_RANGE = 200.0d3, ! everywhere in the domain
 OBSSIM_OBSOUT_FNAME = "${ONAME_OBS}",
/
EOF

  # Add SCALE config
  cat $SCALE_CONF >> $RUNCONF



  echo "mpiexec.hydra -n ${SCALE_NP} ${LETKF_RUN}/../obs/obssim ${RUNCONF} ${WDIR}/log/NOUT-${TNODE_CNT} &" >> $RUNSH
#  echo "mpiexec -n ${SCALE_NP} ${LETKF_RUN}/../obs/obssim ${RUNCONF} ${WDIR}/NOUT" >> $RUNSH

#  TNODE_CNT=$(expr ${TNODE_CNT} + ${MEM_NP})   
  TNODE_CNT=$(( TNODE_CNT + MEM_NP / PPN  ))   
  VCODE_CNT=$(expr ${VCODE_CNT} + 1)   

  done # -- MEM

  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H:%M:%S')
done # -- time
echo "wait" >> $RUNSH


sed -i -e  's/<TNODE_CNT>/'${TNODE_CNT}'/g' $RUNSH
sed -i -e  's/<TPRC>/'$((TNODE_CNT*SCALE_NP))'/g' $RUNSH
#
#echo ${WDIR}


cd $WDIR
pjsub $RUNSH
cd -

exit



