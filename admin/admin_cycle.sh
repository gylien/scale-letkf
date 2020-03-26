#!/bin/bash -l
#-------------------------------------------------------------------------------

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd ${wkdir}

myname=$(basename "$0")
myname1=${myname%.*}

. admin.rc || exit $?

rundir=${realtimebase}/scale_${scale_ver}/scale-letkf_${letkf_ver}/scale/run

outdir="${realtimebase}/result/${expname}" 

bdytmpdir="${outdir}/ncepgfs_grads_da"
obstmpdir="${outdir}/ncepobs"

#-------------------------------------------------------------------------------

function unlock () {
  if (($(cat $lockfile) == $$)); then
    rm -f $lockfile
  fi
}
trap unlock EXIT

#-------------------------------------------------------------------------------

now="$(date -u +'%Y-%m-%d %H:%M:%S')"
if [ ! -s "$timefile" ]; then
  echo "$now [ERR ] Cannot find previous model time." >> $logfile
  exit
fi
if [ -e "$lockfile" ]; then
  echo "$now [PREV]" >> $logfile
  exit
else
  echo $$ >> $lockfile
fi

PREVIOUS_TIME="$(cat $timefile)"
TIME="$(date -ud "${PREVIOUS_TIME}" +'%Y-%m-%d %H:%M:%S')"
TIMEf="$(date -ud "$TIME" +'%Y%m%d%H%M%S')"
ATIME="$(date -ud "$LCYCLE second $TIME" +'%Y-%m-%d %H:%M:%S')"
ATIMEf="$(date -ud "$ATIME" +'%Y%m%d%H%M%S')"

TIMEstop="$(date -ud "$TIME" +'%Y-%m-%d %H:%M:%S')"

if [ -s "$timestopfile" ]; then
  TIMEstoptmp="$(cat $timestopfile)"
  if (($(date -ud "$TIMEstoptmp" +'%s') < $(date -ud "$TIMEstop" +'%s'))); then
    TIMEstop="$TIMEstoptmp"
  fi
fi

if (($(date -ud "$TIME" +'%s') > $(date -ud "$TIMEstop" +'%s'))); then
  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [STOP] $ATIMEf" >> $logfile
  exit
fi

#-------------------------------------------------------------------------------

now="$(date -u +'%Y-%m-%d %H:%M:%S')"
echo "$now [TRY ] $ATIMEf" >> $logfile

n=0
istime="$TIME"
istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
iatime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
iatimef="$(date -ud "$iatime" +'%Y%m%d%H%M%S')"
while ((istimef <= $(date -ud "$TIMEstop" +'%Y%m%d%H%M%S') || n == 0)); do

  if [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/atm_${istimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/atm_${iatimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/atm_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/sfc_${istimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/sfc_${iatimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/sfc_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/land_${istimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/land_${iatimef}.grd" ] ||
     [ ! -s "$gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/land_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd" ]; then
    if ((n == 0)); then
      echo "$now [WAIT] $iatimef - Model files are not ready." >> $logfile
      exit
    else
      break
    fi
  fi


  if [ ! -s "$obsdir/$(date -ud "$iatime" +'%Y%m%d%H')/obs_${iatimef}.dat" ]; then
echo "$obsdir/$(date -ud "$iatime" +'%Y%m%d%H')/obs_${iatimef}.dat" # DEBUG
    if ((n == 0)); then
      echo "$now [WAIT] $iatimef - Observation files are not ready." >> $logfile
      exit
    else
      break
    fi
  fi

  n=$((n+1))
  istime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
  iatime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  iatimef="$(date -ud "$iatime" +'%Y%m%d%H%M%S')"
done

ETIME="$(date -ud "- $LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
ETIMEf="$(date -ud "$ETIME" +'%Y%m%d%H%M%S')"
EATIME="$(date -ud "$LCYCLE second $ETIME" +'%Y-%m-%d %H:%M:%S')"
EATIMEf="$(date -ud "$EATIME" +'%Y%m%d%H%M%S')"
NCYCLE=$n

if [ "$TIMEf" = "$ETIMEf" ]; then
  TIMEPRINT="$ATIMEf"
else
  TIMEPRINT="${ATIMEf}-${EATIMEf}"
fi

#-------------------------------------------------------------------------------

rm -f $outfile

now="$(date -u +'%Y-%m-%d %H:%M:%S')"
echo "$now [TRAN] $TIMEPRINT - Prepare external files" >> $logfile

istime="$TIME"
istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
iatime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
iatimef="$(date -ud "$iatime" +'%Y%m%d%H%M%S')"
while ((istimef <= ETIMEf)); do
  mkdir -p $bdytmpdir/${istimef}/mean
  rm -fr $bdytmpdir/${istimef}/mean/*

  ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/atm_${istimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/atm_${iatimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/atm_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/sfc_${istimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/sfc_${iatimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/sfc_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$istime" +'%Y%m%d%H')/land_${istimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/land_${iatimef}.grd $bdytmpdir/${istimef}/mean
  ln -s $gradsdir/$(date -ud "$iatime" +'%Y%m%d%H')/land_$(date -ud "$LCYCLE second $iatime" +'%Y%m%d%H%M%S').grd $bdytmpdir/${istimef}/mean

  ln -s $obsdir/$(date -ud "$iatime" +'%Y%m%d%H')/obs_${iatimef}.dat ${obstmpdir}

  istime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  istimef="$(date -ud "$istime" +'%Y%m%d%H%M%S')"
  iatime="$(date -ud "$LCYCLE second $istime" +'%Y-%m-%d %H:%M:%S')"
  iatimef="$(date -ud "$iatime" +'%Y%m%d%H%M%S')"
done



success=0
  total_sec=$(date -ud "1970-01-01 ${WTIME_L[$ntry]}" +'%s')
  WTIME_L_use=$(date -ud "$((total_sec * NCYCLE)) second 1970-01-01 00:00:00" +'%H:%M:%S')
  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [RUN ] $TIMEPRINT/$ntry" >> $logfile

  cd ${rundir}
  ./admin.sh ${TIMEf} "00:30:00" $nmem

  res=$?

  if ((res == 0)); then
    success=1
  else
    now="$(date -u +'%Y-%m-%d %H:%M:%S')"
    echo "$now [ERR ] $TIMEPRINT/$ntry - Exit code: $res" >> $logfile
    sleep 10s
  fi

if ((success == 1)); then
  cd ${plotdir}

  now="$(date -u +'%Y-%m-%d %H:%M:%S')"
  echo "$now [DONE] $TIMEPRINT" >> $logfile
  echo "$EATIME" > $timefile
else
  echo "$now [FAIL] $TIMEPRINT - All trials failed" >> $logfile
fi
