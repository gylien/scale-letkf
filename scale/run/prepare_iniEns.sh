#/bin/bash

# Load functions
. src/func_datetime.sh

# Load basic setting
. config.main
. config.fcst

# CHECK
FCSTLEN=0

ATIME=$(datetime $STIME ${FCSTLEN} s)

# prepare init (restart) directory
AMEAN=$OUTDIR/$ATIME/anal/mean
mkdir -p $AMEAN

FMEAN=$OUTDIR/$STIME/fcst/mean

# copy original restart files for "mean" from a nature run
if [ ! -e $AMEAN/init.pe$(printf %06d 0).nc ] ; then
  for p in `seq 0 $((SCALE_NP-1))`
  do
    ORG_INITFILE=$FMEAN/init_${ATIME}.pe$(printf %06d $p).nc
    NEW_INITFILE=$AMEAN/init.pe$(printf %06d $p).nc
    cp  $ORG_INITFILE $NEW_INITFILE
  done
else
  echo "mean is ready (not copy)"
fi

# copy the mean directory including restart files for ensemble members
for m in `seq 1 $((MEMBER))`
do
  mem=$(printf %04d $m)
  if ((m > MEMBER)) ; then
    mem='mdet'
  fi

  AMEM=$OUTDIR/$ATIME/anal/$mem
  if [ ! -e $AMEM/init.pe$(printf %06d 0).nc ] ; then
    rm -rf $AMEM
    echo "copy mean "$mem
    cp -r $AMEAN $AMEM
  fi
done



