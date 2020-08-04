#/bin/bash -l
#-------------------------------------------------------------------------------

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd ${wkdir}
myname=$(basename "$0")
myname1=${myname%.*}

TIME=$1
TIMEf="$(date -ud "${TIME}" +'%Y%m%d%H%M%S')"


. ./admin.rc || exit $?

rundir="${realtimebase}/scale_${scale_ver}/scale-letkf_${letkf_ver}/scale/run_d1-2"

initdir=${realtimebase}/result/$expname/d1
bdytmpdir=${realtimebase}/result/$expname/ncepgfs_grads


lockfile_time="${lockfile}.${TIMEf}"
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


TIMEstop="$(date -ud "$TIME" +'%Y-%m-%d %H:%M:%S')"

if [ -s "$timestopfile" ]; then
  TIMEstoptmp="$(cat $timestopfile)"
  if (($(date -ud "$TIMEstoptmp" +'%s') < $(date -ud "$TIMEstop" +'%s'))); then
    TIMEstop="$TIMEstoptmp"
  fi
fi

if [ -s "$wkdir/../admin_r0051_AA/admin_cycle.time" ]; then
  CYCLE_TIME="$(cat "$wkdir/../admin_r0051_AA/admin_cycle.time")"
  if (($(date -ud "$CYCLE_TIME" +'%s') < $(date -ud "$TIMEstop" +'%s'))); then
    TIMEstop="$CYCLE_TIME"
  fi
fi


#-------------------------------------------------------------------------------


#if (($(date -ud "$TIME" +'%s') > $(date -ud "$TIMEstop" +'%s'))); then
#  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
#  echo "$now [STOP] $TIMEf" >> $logfile
#  exit
#fi

now="$(date -u +'%Y-%m-%d %H:%M:%S')"
echo "$now [TRY ] $TIMEf" >> $logfile

n=0
istime="$TIME"
istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
while ((istimef <= $(date -ud "$TIMEstop" +'%Y%m%d%H%M%S') || n == 0)); do
  ready=1

  tfcst=0
  while ((tfcst <= FCSTLEN_d2)); do
    itimef="$(date -ud "$tfcst second $istime" +'%Y%m%d%H%M%S')"
#    if [ ! -s "$wrfdir/$(date -ud "$istime" +'%Y%m%d%H')/mean/wrfout_${itimef}" ]; then
    if [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/atm_${itimef}.grd" ] ||
       [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/sfc_${itimef}.grd" ] ||
       [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/land_${itimef}.grd" ]; then
      ready=0
      break
    fi
    tfcst=$((tfcst+LCYCLE))
  done
  if ((ready == 0)); then
    if ((n == 0)); then
      echo "$now [WAIT] $istimef - Model files are not ready." >> $logfile
      exit
    else
      break
    fi
  fi

  if [ ! -s "$initdir/${istimef}/anal/mean/init.pe000000.nc" ]; then
    if ((n == 0)); then
      echo "$now [WAIT] $istimef - LETKF analyses are not ready." >> $logfile
      exit
    else
      break
    fi
  fi

  n=$((n+1))
  istime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
done

ETIME="$(date -ud "- $LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
ETIMEf="$(date -ud "${ETIME}" +'%Y%m%d%H%M%S')"
NCYCLE=$n

if ((NCYCLE < F_MIN_CYCLE)); then
  echo "$now [WAIT] $TIMEf - Available cycles = $NCYCLE; minimum cycles to run = $F_MIN_CYCLE" >> $logfile
  exit
fi

if [ "$TIMEf" = "$ETIMEf" ]; then
  TIMEPRINT="$TIMEf"
else
  TIMEPRINT="${TIMEf}-${ETIMEf}"
fi

#-------------------------------------------------------------------------------

rm -f $outfile

cd $RUNDIR


now="$(date -u +'%Y-%m-%d %H:%M:%S')"
echo "$now [TRAN] $TIMEPRINT - Upload files" >> $logfile

istime="$TIME"
istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
while ((istimef <= ETIMEf)); do
  mkdir -p $bdytmpdir/${istimef}/mean
  rm -fr $bdytmpdir/${istimef}/mean/*
  tfcst=0
  while ((tfcst <= FCSTLEN_d2+LCYCLE)); do
    itimef="$(date -ud "$tfcst second $istime" +'%Y%m%d%H%M%S')"
    ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/atm_${itimef}.grd $bdytmpdir/${istimef}/mean
    ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/sfc_${itimef}.grd $bdytmpdir/${istimef}/mean
    ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/land_${itimef}.grd $bdytmpdir/${istimef}/mean
    tfcst=$((tfcst+LCYCLE))
  done


  istime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
done

  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [RUN ] $TIMEPRINT " >> $logfile
  cd ${rundir} 
  ./admin.sh ${TIMEf} ${FCSTLEN_d2} '03:00:00' $nmem_d2 
  res=$?

if ((res == 0)); then
  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [DONE] $TIMEPRINT" >> $logfile
else
  echo "$now [FAIL] $TIMEPRINT - All trials failed" >> $logfile
fi

cd $wkdir
