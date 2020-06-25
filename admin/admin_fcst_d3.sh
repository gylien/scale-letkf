#!/bin/bash -l
#-------------------------------------------------------------------------------

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd ${wkdir}
myname=$(basename "$0")
myname1=${myname%.*}

. ./admin.rc || exit $?

#if (($(ssh-add -L | cut -d " " -f 3) != ${HOME}/.ssh/${rsa_key})) ;then
# echo "prepare ssh-agent first."
# exit $?
#fi

REF_TIME=$1
TIME=$2
FCSTLEN=$3

if  [ -z $FCSTLEN ] ; then
 echo "specify forecast length (second)"
 exit 0
fi

REFTIMEf="$(date -ud "${REF_TIME}" +'%Y%m%d%H%M%S')"
TIMEf="$(date -ud "${TIME}" +'%Y%m%d%H%M%S')"

lockfile_time="${lockfile}.${REFTIMEf}.${TIMEf}"

rundir=${realtimebase}/scale_${scale_ver}/scale-letkf_${letkf_ver}/scale/run_d3

#web_url=daweb.r-ccs27.riken.jp
#web_outdir=/home/amemiya/public_html/scale/data/ens/fcst_d3/ref_${REFTIMEf}
#rsakey_web='/home/amemiya/.ssh/id_rsa_mac'

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

istime="$TIME"
istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"

success=0
  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [RUN ] $TIMEf" >> $logfile
  cd ${rundir} 
  ./admin.sh ${REFTIMEf} ${TIMEf} ${FCSTLEN} "02:00:00" ${nmem_d3} &> admin.log.${REFTIMEf}.${TIMEf}
  res=$?

  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  if ((res == 0)); then
  echo "$now [DONE] $TIMEf" >> $logfile
  else
    now="$(date -u +'%Y-%m-%d %H:%M:%S')"
    echo "$now [ERR ] $TIMEf - Exit code: $res" >> $logfile
  echo "$now [FAIL] $TIMEf - All trials failed" >> $logfile
  fi
