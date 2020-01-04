#!/bin/bash -l

sec_sleep=600
LCYCLE=21600

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd $wkdir

running='running_fcst_d1_ext'

if [ -f $running ] ;then
  echo 'already running.'
  exit 1
else
  echo $HOSTNAME $$ > $running
fi

function unlock () {
 [ `cat $running | awk '{print $2}'` == $$ ] && rm -f $running
}
trap unlock EXIT

. ${wkdir}/admin.rc || exit $?

INITSRC=$realtimebase/result/$expname/d1


wkdir="$(cd "$( dirname "$0" )" && pwd)"

source ~/.bashrc

START_TIME=$1
RUN_TIME=$START_TIME
RUN_TIME_END="2019-01-01 00:00:00"
[ -s "$timestopfile" ] && RUN_TIME_END="$(cat $timestopfile)"


while [ "$RUN_TIME" != "$RUN_TIME_END" ] ;do
  CYCLE_TIME="$(cat admin_cycle.time)"
  CTIMEf="$(date -ud "${CYCLE_TIME}" +'%Y%m%d%H%M%S')"
  TIMEf="$(date -ud "${RUN_TIME}" +'%Y%m%d%H%M%S')"
   if [ $CTIMEf -ge $TIMEf ] ;then
    if [ -s $INITSRC/$TIMEf/anal/mean/init.pe000000.nc ] ;then
     now="$(date -u +'%Y-%m-%d %H:%M:%S')"
     echo "$now $RUN_TIME start " 
     nohup ./admin_fcst_d1_ext.sh "$RUN_TIME" &> admin_fcst.log.${TIMEf} &
     RUN_TIME="$(date -ud "$LCYCLE second $RUN_TIME" +'%Y-%m-%d %H:%M:%S')"
    sleep 30
    fi
   else
    now="$(date -u +'%Y-%m-%d %H:%M:%S')"
    echo "$now wait "
   sleep $sec_sleep
   fi
done

