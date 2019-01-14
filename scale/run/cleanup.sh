#/bin/bash

. config.main

STIME='20000101003500'
ETIME='20000101013500'

Y4S=${STIME:0:4}
M2S=${STIME:4:2}
D2S=${STIME:6:2}
H2S=${STIME:8:2}
N2S=${STIME:10:2}
S2S=${STIME:12:2}

Y4E=${ETIME:0:4}
M2E=${ETIME:4:2}
D2E=${ETIME:6:2}
H2E=${ETIME:8:2}
N2E=${ETIME:10:2}
S2E=${ETIME:12:2}


GUES=0 # Only gues/mean & gues/sprd
HIST=0 # No hist
ANAL=1 # All member stored only every 10 min
OBSG=0 # No obsgues

LOG=0 # head node only

tstart=''$Y4S'-'$M2S'-'$D2S' '$H2S':'$N2S':'$S2S''
tend=''$Y4E'-'$M2E'-'$D2E' '$H2E':'$N2E':'$S2E''
tint=300 # sec

MV=1 # mv data from tmp to OUTDIR
THRD=1
THRD_MAX=8

time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do

#  timef=$(date -ud "$time" '+%Y-%m-%d %H:%M:%S')
  timef=$(date -ud "${tint} second $time" '+%Y-%m-%d %H:%M:%S')

  YYYY=$(date -ud "$timef" '+%Y')
  MM=$(date -ud "$timef" '+%m')
  DD=$(date -ud "$timef" '+%d')
  HH=$(date -ud "$timef" '+%H')
  MN=$(date -ud "$timef" '+%M')
  SS=$(date -ud "$timef" '+%S')

  DTIME=${YYYY}${MM}${DD}${HH}${MN}${SS}
#  echo $DTIME

  TMPOUT=${TMP}/out/$DTIME
  #echo $TMPOUT

  if (( GUES < 1 )); then
#    echo "Remove GUES"
    GUES_DIR=$TMPOUT/gues  
    rm -rf $GUES_DIR/0*
    rm -rf $GUES_DIR/meanf
  fi

  if (( HIST == 0 )); then
    HIST_DIR=$TMPOUT/hist
    rm -rf $HIST_DIR
  fi

  if (( ANAL == 1 )) && [ "${MN:1:1}" != '0' ]; then
    ANAL_DIR=$TMPOUT/anal
    rm -rf $ANAL_DIR/0*
#echo "ANAL_DIR "$ANAL_DIR" "${MN:1:1}
  fi

  if (( OBSG == 0 )); then
    OBSG_DIR=$TMPOUT/obsgues
    rm -rf $OBSG_DIR
  fi

  if (( LOG == 0 )); then
    LOG_DIR=$TMPOUT/log
    LOG_OUT=$OUTDIR/$DTIME/log

    if [ -e "$LOG_DIR/letkf/NOUT-000000" ] ; then 
      mkdir -p $LOG_OUT/letkf
      cp $LOG_DIR/letkf/NOUT-000000 $LOG_OUT/letkf/
      rm -rf $LOG_DIR
    fi
  fi

  if (( MV == 1 )) && (( THRD < THRD_MAX )) ; then
    cp -r $TMPOUT $OUTDIR/ &
    THRD=$(( THRD + 1))
  elif (( THRD == THRD_MAX )) ; then
    wait
    THRD=1
  fi
echo "CHK "$THRD

  time=$(date -ud "${tint} second $time" '+%Y-%m-%d %H:%M:%S')
done




