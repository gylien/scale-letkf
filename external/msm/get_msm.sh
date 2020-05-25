#!/bin/bash -l

wkdir="$( cd "$( dirname "$0" )" && pwd )"

MSMSRCDIR=http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/latest
datadir=$wkdir/..  ### EDIT ME

cd ${wkdir}
now=`date -u +'%Y-%m-%d %H:%M:%S'`

if [ ! -s "${wkdir}/mtime" ]; then
  echo "$now [ERR ] Cannot find previous model time."
  exit
fi

if [ -e "${wkdir}/running" ]; then
  echo "$now [PREV]" >> ${wkdir}/get_msm.log
  exit
else
  touch ${wkdir}/running
fi

PREVIOUS_TIME=`cat ${wkdir}/mtime`
GET_TIME=`date -u -d "+ 6 hour ${PREVIOUS_TIME}" +'%Y-%m-%d %H'`
YYYY=`date -u -d "$GET_TIME" +'%Y'`
MM=`date -u -d "$GET_TIME" +'%m'`
DD=`date -u -d "$GET_TIME" +'%d'`
HH=`date -u -d "$GET_TIME" +'%H'`
YYYYMMDDHH="$YYYY$MM$DD$HH"
YYYYMMDD=${YYYYMMDDHH:0:8}

echo "$now [TRY ] $YYYYMMDDHH" >> ${wkdir}/get_msm.log

mkdir -p ${datadir}/msm/$YYYYMMDDHH
cd ${datadir}/msm/$YYYYMMDDHH

cfiles=MSM${YYYYMMDDHH}S.nc
cfilep=MSM${YYYYMMDDHH}P.nc


if [ -s $cfiles ] && [ -s $cfilep ]; then
  allget=1
else
wget --cache=off $MSMSRCDIR/$YYYYMMDD/$cfiles
wget --cache=off $MSMSRCDIR/$YYYYMMDD/$cfilep
if [ -s $cfiles ] && [ -s $cfilep ]; then
  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [GET ] $cfiles and $cfilep " >> ${wkdir}/get_msm.log
  allget=1
else
  allget=0
fi
fi



  if ((allget == 1)); then
      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [CONV] $YYYYMMDDHH -> msm_{s/p}.$YYYYMMDDHHMMSS - GrADS" >> ${wkdir}/get_msm.log
      bash $wkdir/run/grads/convert.sh "$datadir/msm/${YYYYMMDDHH}" \
       > ${wkdir}/convert_grads.log 2>&1




### hourly plot
[ -f convert_plot.log ] && rm convert_plot.log
t=0
while ((t <= 33)); do ### extended ?
tf=`printf '%03d' $t`
TIME_fcst=`date -ud "${t} hour $YYYY-$MM-$DD $HH" +'%Y-%m-%d %H:%M:%S'`
YYYYMMDDHHMMSS_fcst=`date -ud "$TIME_fcst" +'%Y%m%d%H%M%S'`
      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [PLOT] $YYYYMMDDHH - Plot (background job)" >> ${wkdir}/get_msm.log
      bash $wkdir/run/plot/plot.sh "$TIME_fcst" "$GET_TIME" $((t+1)) "$datadir/msm/${YYYYMMDDHH}" \
       >> ${wkdir}/convert_plot.log 2>&1
      bash $wkdir/run/plot/plot_d3.sh "$TIME_fcst" "$GET_TIME" $((t+1)) "$datadir/msm/${YYYYMMDDHH}" \
       >> ${wkdir}/convert_plot_d3.log 2>&1
t=$((t+1))
done

  fi



if ((allget == 1)); then
  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [DONE] $YYYYMMDDHH" >> ${wkdir}/get_msm.log
  echo "$GET_TIME" > ${wkdir}/mtime
fi

rm -f ${wkdir}/running
