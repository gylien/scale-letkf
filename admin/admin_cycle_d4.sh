#!/bin/bash -lx
#-------------------------------------------------------------------------------

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd ${wkdir}
myname=$(basename "$0")
myname1=${myname%.*}

. ${wkdir}/admin.rc || exit $?

TIME=$1
FCSTLEN_ext_d4=$2

if  [ -z $FCSTLEN_ext_d4 ] ; then
 echo "specify forecast length (second)"
 exit 0
fi

TIMEf="$(date -ud "${TIME}" +'%Y%m%d%H%M%S')"

lockfile_time="${lockfile}.${TIMEf}"

rundir=${realtimebase}/scale_${scale_ver}/scale-letkf-rt_dacycle/scale/run

#-------------------------------------------------------------------------------

function unlock () {
  if (($(cat $lockfile_time) == $$)); then
    rm -f $lockfile_time
  fi
}
trap unlock EXIT

#-------------------------------------------------------------------------------

now="$(date -u +'%Y-%m-%d %H:%M:%S')"

if [ -e "$lockfile_time" ]; then
  echo "$now [PREV]" >> $logfile
  exit
else
  echo $$ >> $lockfile_time
fi


### 
### START
###


#-------------------------------------------------------------------------------

now="$(date -u +'%Y-%m-%d %H:%M:%S')"
echo "$now [TRY ] $TIMEf" >> $logfile


success=0
  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [RUN ] $TIMEf" >> $logfile
  cd ${rundir} 
  ./admin.sh ${TIMEf} 10 "00:30:00" ${n_mem} ${FCSTLEN_ext_d4} &> admin.log
  res=$?
  echo $res
  if ((res == 0)); then
    success=1
  else
    now="$(date -u +'%Y-%m-%d %H:%M:%S')"
    echo "$now [ERR ] $TIMEf - Exit code: $res" >> $logfile
    sleep 10s
  fi

if ((success == 1)); then

  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [DONE] $TIMEf" >> $logfile
else
  echo "$now [FAIL] $TIMEf - All trials failed" >> $logfile
fi
