#/bin/bash

CP_THRDMAX=8

# Load functions
. src/func_datetime.sh

# Load basic setting
#. config.main
#. config.fcst

EXP=2000m_InSnd_LT_SN14_Mac_0605
OUTDIR="/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/$EXP"
MEMBER=80
SCALE_NP=64

SCALE_NP_TMP=$((SCALE_NP - 1))

# CHECK
STIME='20000101000000'
ETIME='20000101000000'
FCSTLEN=2400 # 40min

ATIME=$(datetime $STIME ${FCSTLEN} s)

# prepare init (restart) directory
AMEAN=$OUTDIR/$ATIME/anal/mean
mkdir -p $AMEAN

FMEAN=$OUTDIR/$STIME/fcst/mean


for m in `seq 0 $((MEMBER))`
do
  mem=$(printf %04d $m)
  if ((m > MEMBER)) ; then
    mem='mdet'
  elif ((m == 0)) ; then
    mem='mean'
  fi
  AMEM=$OUTDIR/$ATIME/anal/$mem
 
  mkdir -p ${AMEM}
  echo $mem
  
  for p in `seq 0 $((SCALE_NP-1))`
  do
    PRC=$(printf %06d $p)

    ORG=$OUTDIR/$STIME/fcst/${mem}/init_${ATIME}.pe${PRC}.nc
    NEW=${AMEM}/init.pe${PRC}.nc
    ln -sf ${ORG} ${NEW}

  done

done

