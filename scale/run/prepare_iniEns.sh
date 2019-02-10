#/bin/bash

CP_THRDMAX=8

# Load functions
. src/func_datetime.sh

# Load basic setting
. config.main
. config.fcst

SCALE_NP_TMP=$((SCALE_NP - 1))

# CHECK
STIME='20000101000000'
ETIME='20000101000000'
FCSTLEN=2100 # 35min
FCSTLEN=1800 # 30min

ATIME=$(datetime $STIME ${FCSTLEN} s)

# prepare init (restart) directory
AMEAN=$OUTDIR/$ATIME/anal/mean
mkdir -p $AMEAN

FMEAN=$OUTDIR/$STIME/fcst/mean


# copy original restart files for "mean" from a nature run
if [ ! -e $AMEAN/init.pe$(printf %06d $SCALE_NP_TMP).nc ] ; then

  echo "copy mean"
  CP_THRD=0

  for p in `seq 0 $((SCALE_NP-1))`
  do
    if (( CP_THRD >= CP_THRDMAX )) ; then
      wait
      CP_THRD=0
    fi
    ORG_INITFILE=$FMEAN/init_${ATIME}.pe$(printf %06d $p).nc
    NEW_INITFILE=$AMEAN/init.pe$(printf %06d $p).nc
    cp  $ORG_INITFILE $NEW_INITFILE &

    CP_THRD=$((CP_THRD + 1))
  done
  wait
else
  echo "mean is ready (not copy)"
fi

CP_THRD=0
# copy the mean directory including restart files for ensemble members
for m in `seq 1 $((MEMBER))`
do
  echo $m
  if (( CP_THRD >= CP_THRDMAX )) ; then
    wait
    CP_THRD=0
  fi

  mem=$(printf %04d $m)
  if ((m > MEMBER)) ; then
    mem='mdet'
  fi

  AMEM=$OUTDIR/$ATIME/anal/$mem
  if [ ! -e $AMEM/init.pe$(printf %06d $SCALE_NP_TMP).nc ] ; then
    rm -rf $AMEM
    echo "copy mean "$mem
    cp -r $AMEAN $AMEM &

    CP_THRD=$((CP_THRD + 1))
  fi
done



