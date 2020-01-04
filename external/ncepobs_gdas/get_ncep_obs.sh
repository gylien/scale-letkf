#!/bin/bash -l

source ~/.bashrc

wkdir="$( cd "$( dirname "$0" )" && pwd )"
datadir=$wkdir/..


cd ${wkdir}
now=`date -u +'%Y-%m-%d %H:%M:%S'`

if [ ! -s "${wkdir}/mtime" ]; then
  echo "$now [ERR ] Cannot find previous model time."
  exit
fi

if [ -e "${wkdir}/running" ]; then
  echo "$now [PREV]" >> ${wkdir}/get_ncep_obs.log
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

echo "$now [TRY ] $YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log

mkdir -p ${wkdir}/$YYYYMMDDHH
cd ${wkdir}/$YYYYMMDDHH

allget=1
if [ ! -s "prepbufr.$YYYYMMDDHH" ]; then
  rm -f gdas.t${HH}z.prepbufr.nr
#  wget --limit-rate=500k ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.${YYYY}${MM}${DD}/gdas.t${HH}z.adpupa.tm00.bufr_d.unblok
    echo "$now [GET ] start" >> ${wkdir}/get_ncep_obs.log
  wget --cache=off https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.${YYYY}${MM}${DD}/${HH}/gdas.t${HH}z.prepbufr.nr
  if [ -s "gdas.t${HH}z.prepbufr.nr" ]; then
    mv -f gdas.t${HH}z.prepbufr.nr prepbufr.$YYYYMMDDHH
    now=`date -u +'%Y-%m-%d %H:%M:%S'`
    echo "$now [GET ] prepbufr.$YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log
  else
    allget=0
  fi
fi

#if [ ! -s "prepbufr.$YYYYMMDDHH.unblok" ]; then
#  rm -f gdas.t${HH}z.prepbufr.unblok.nr
##  wget --limit-rate=500k ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.${YYYY}${MM}${DD}/gdas.t${HH}z.adpupa.tm00.bufr_d.unblok
#  wget --cache=off http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.${YYYY}${MM}${DD}/gdas.t${HH}z.prepbufr.unblok.nr
#  if [ -s "gdas.t${HH}z.prepbufr.unblok.nr" ]; then
#    mv -f gdas.t${HH}z.prepbufr.unblok.nr prepbufr.$YYYYMMDDHH.unblok
#    now=`date -u +'%Y-%m-%d %H:%M:%S'`
#    echo "$now [GET ] prepbufr.$YYYYMMDDHH.unblok" >> ${wkdir}/get_ncep_obs.log
#  else
#    allget=0
#  fi
#fi

if [ "$allget" -eq 1 ]; then
  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [CONV] $YYYYMMDDHH: dec_prepbufr" >> ${wkdir}/get_ncep_obs.log
  bash $wkdir/run/dec_prepbufr/convert.sh "${YYYY}-${MM}-${DD} ${HH}" "$wkdir/${YYYYMMDDHH}" "$datadir/ncepobs_gdas_letkf/${YYYYMMDDHH}" \
   > ${wkdir}/convert_dec_prepbufr.log 2>&1

#  now=`date -u +'%Y-%m-%d %H:%M:%S'`
#  echo "$now [PLOT] $YYYYMMDDHH: map" >> ${wkdir}/get_ncep_obs.log
#  $wkdir/run/plot/map $YYYYMMDDHH 
#  mkdir $wkdir/../ncepobs_gdas_letkf/${YYYYMMDDHH}/map
#  mv $wkdir/run/plot/*.png $wkdir/../ncepobs_gdas_letkf/${YYYYMMDDHH}/map

  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [DONE] $YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log
  echo "$GET_TIME" > ${wkdir}/mtime
fi

rm -f ${wkdir}/running
