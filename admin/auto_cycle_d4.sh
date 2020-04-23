#!/bin/sh

source ~/.bashrc

wkdir="$(cd "$( dirname "$0" )" && pwd)"
. $wkdir/admin.rc || exit $1

#-----------------------------
running=$wkdir'/running_dacycle_d4'

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
#-----------------------------

time_offset=0
[ -f "time_offset.txt" ] && time_offset=`cat time_offset.txt` 

isec=0
cyclesec=60
limitsec=86400

END_TIME="2020-04-23 18:00:00"

rundir=${realtimebase}/scale_${scale_ver}/scale-letkf_${letkf_ver}_d4/scale/run

INIT_TIME="$(date -ud "$time_offset second now" +'%Y-%m-%d %H:00:00')"
INIT_TIMEf="$(date -ud "$INIT_TIME" +'%Y%m%d%H0000')"

cd $rundir

#$rundir/prep.sh init 

  echo " $INIT_TIME start "

  nohup ./start.sh "$INIT_TIMEf" &> admin_cycle_d4.log.${INIT_TIMEf}
 res=$?
 [ "$res" != "0" ] && exit 99
  echo " $INIT_TIME complete "

while [ `date -ud "1 hour $INIT_TIME" +%s` -le `date -ud "$END_TIME" +%s` ] ; do
 INIT_TIME="$(date -ud "1 hour $INIT_TIME" +'%Y-%m-%d %H:00:00')"
 INIT_TIMEf="$(date -ud "$INIT_TIME" +'%Y%m%d%H0000')"
  echo " $INIT_TIME start "
 nohup ./restart.sh "$INIT_TIMEf" &> admin_cycle_d4.log.${INIT_TIMEf}
 [ "$res" != "0" ] && exit 99
  echo " $INIT_TIME complete "
done 

