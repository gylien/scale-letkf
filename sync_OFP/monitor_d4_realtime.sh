#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc


FCSTBASE="$DATADIR/d4/realtime"
zref="z01127m"
end_timef="2020-12-31 00:00:00"


lastmin=202001010000
while [ `date -u +%s` -le `date -ud "$end_timef" +%s` ] ;do
  nowmin=`date -u +%Y%m%d%H%M`
  if [ $nowmin != $lastmin ] ; then
    if [ ! -z "`ls $FCSTBASE/$zref/anal*.png`" ] ; then
      analmin=`ls $FCSTBASE/$zref/anal*.png | tail -n 1 | grep -o '[0-9]\{14\}' | cut -c 1-12`
    else
      analmin='-'
    fi
    if [ ! -z "`ls $FCSTBASE/$zref/fcst*.png`" ] ; then
      fcstmin=`ls $FCSTBASE/$zref/fcst*.png | tail -n 1 | grep -o '[0-9]\{14\}' | cut -c 1-12`
      fcsts=`date -ud "${fcstmin:0:4}-${fcstmin:4:2}-${fcstmin:6:2} ${fcstmin:8:2}:${fcstmin:10:2}" +%s`
      nows=`date -u +%s`
      fcstdif=`expr $fcsts - $nows`
    else
      fcstmin='-'
      fcstdif='-'
    fi
    echo "$nowmin $analmin $fcstmin $fcstdif" >> monitor.log
    lastmin=$nowmin
  fi
  sleep 10s
done
