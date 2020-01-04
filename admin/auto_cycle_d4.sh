#!/bin/sh

isec=0
cyclesec=60
limitsec=21600


. `pwd`/admin.rc || exit $1

FCSTLEN=600

EXPBASE=${realtimebase}/r0051_nest/exp_d4_1km

wkdir="$(cd "$( dirname "$0" )" && pwd)"

source ~/.bashrc


while [ $isec -le $limitsec ] ;do


START_TIME="$(date -u +'%Y-%m-%d %H:%M:%S')"
START_TIMEf="$(date -ud "${START_TIME}" +'%Y%m%d%H%M%S')"
LIMIT_TIMEf="$(date -ud "-1 hour ${START_TIME}" +'%Y%m%d%H%M%S')"


INIT_TIME="$(date -ud "${START_TIME}" +'%Y-%m-%d %H:00:00')"
INIT_TIMEf="$(date -ud "${INIT_TIME}" +'%Y%m%d%H%M%S')"

while [ ! -d $EXPBASE/$INIT_TIMEf/bdy/0050 ];do
 INIT_TIME="$(date -ud "-1 hour ${INIT_TIME}" +'%Y-%m-%d %H:00:00')"
 INIT_TIMEf="$(date -ud "${INIT_TIME}" +'%Y%m%d%H%M%S')"
done
 if [ $INIT_TIMEf -ge $LIMIT_TIMEf ] ; then
 nohup ./admin_cycle_d4.sh "$INIT_TIME" $FCSTLEN &> admin_cycle_d4.log&
 fi

 sleep ${cyclesec}s 
 isec=`expr $isec + $cyclesec`
 echo "$isec / $limitsec wait..."
done


