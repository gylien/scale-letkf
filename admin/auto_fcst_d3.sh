#!/bin/bash -l

running='running_fcst_d3'

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

wkdir="$(cd "$( dirname "$0" )" && pwd)"

. ${wkdir}/admin.rc || exit $1

isec=0
cyclesec=300
limitsec=43200

PARENT_TIME_B=`date -d "-1 year" +"%F %T"`

spinup_hour=2
fcst_hour=7
start_hour_fromnow=0

FCSTHOUR_DEF=`expr $fcst_hour + $spinup_hour` ### 6+2


####PARENT_FCSTLEN=${FCSTLEN_d2}
PARENT_FCSTLEN=64800

ofp_parentdir=$realtimebase/result/ope/d2

source ~/.bashrc

time_offset=0
[ -f "time_offset.txt" ] && time_offset=`cat time_offset.txt` 



while [ $isec -le $limitsec ] ;do

echo "$isec / $limitsec"

START_TIME="$(date -ud "$time_offset second now" +'%Y-%m-%d %H:%M:%S')"
START_TIMEf="$(date -ud "${START_TIME}" +'%Y%m%d%H%M%S')"

[ -s "$timestopfile" ] && RUN_TIME_END="$(cat $timestopfile)" || RUN_TIME_END="$(date -ud "1 day $START_TIME" +'%Y-%m-%d %H:%M:%S')"
RUN_TIME_ENDf="$(date -ud "${RUN_TIME_END}" +'%Y%m%d%H%M%S')"

PARENT_TIME="$(date -ud "${START_TIME}" +'%Y-%m-%d %H:00:00')"
PARENT_TIMEf="$(date -ud "${PARENT_TIME}" +'%Y%m%d%H%M%S')"


cmem=`printf %04d $nmem_d3`

if [ ! -z "`ls $ofp_parentdir/*/fcst/mean/history.pe000000.nc`" ] ;then

while [ ! -f $ofp_parentdir/$PARENT_TIMEf/fcst/$cmem/history.pe000000.nc ] || [ -f ./admin_fcst_d1-2.lock.${PARENT_TIMEf} ] ;do
 PARENT_TIME="$(date -ud "-1 hour ${PARENT_TIME}" +'%Y-%m-%d %H:00:00')"
 PARENT_TIMEf="$(date -ud "${PARENT_TIME}" +'%Y%m%d%H%M%S')"
done

 if [ `date -d "$PARENT_TIME" +%s` -gt `date -d "$PARENT_TIME_B" +%s` ] ;then

# echo "PARENT_TIME" $PARENT_TIME $PARENT_TIME_B
 PARENT_TIME_B=$PARENT_TIME

 INIT_LIMITf="$(date -ud "$PARENT_FCSTLEN second -2 hour -${spinup_hour} hour  ${PARENT_TIME}"  +'%Y%m%d%H%M%S')"
 PARENT_LIMITf="$(date -ud "$PARENT_FCSTLEN second ${PARENT_TIME}"  +'%Y%m%d%H%M%S')"


 starth="$(date -ud "${START_TIME}" +%H)"
 hdif=`expr \( $starth + ${spinup_hour} \) \% 6` ### set aviable forecast hour to 0-6, 6-12, 12-18, or 18-24

 INIT_START="$(date -ud " - $hdif hour ${START_TIME}" +'%Y-%m-%d %H:00:00')"
 INIT_STARTf="$(date -ud "${INIT_START}" +'%Y%m%d%H%M%S')"

 while [ $INIT_STARTf -le $INIT_LIMITf ]; do
 echo "INIT_STARTf" $INIT_STARTf $INIT_LIMITf
    now="$(date -u +'%Y-%m-%d %H:%M:%S')"
    echo "$now ${PARENT_TIMEf}.${INIT_STARTf} start "
    FCSTHOUR=$FCSTHOUR_DEF
    DESTFCSTf="$(date -ud "$FCSTHOUR hour ${INIT_START}" +'%Y%m%d%H%M%S')"    
    while [ $DESTFCSTf -gt $PARENT_LIMITf ] ;do
     FCSTHOUR=`expr $FCSTHOUR - 1`
     DESTFCSTf="$(date -ud "$FCSTHOUR hour ${INIT_START}" +'%Y%m%d%H%M%S')"    
    done
    FCSTLEN=`expr $FCSTHOUR \* 3600`
    nohup ./admin_fcst_d3.sh "$PARENT_TIME" "$INIT_START" "$FCSTLEN" &> /dev/null &
#    echo "start" "$PARENT_TIME" "$INIT_START" "$FCSTLEN" 
    INIT_START="$(date -ud " $FCSTHOUR hour -${spinup_hour} hour -1 hour ${INIT_START}" +'%Y-%m-%d %H:00:00')"
    INIT_STARTf="$(date -ud "${INIT_START}" +'%Y%m%d%H%M%S')"
 sleep 31
 isec=0
 done

 fi

 fi

 sleep ${cyclesec}s 
 isec=`expr $isec + $cyclesec`

 echo "$isec / $limitsec wait..."
done


