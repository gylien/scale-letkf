#!/bin/bash -l

tmax=36
nmem=20

wkdir="$( cd "$( dirname "$0" )" && pwd )"
datadir=$wkdir/../../external

cd ${wkdir}
now=`date -u +'%Y-%m-%d %H:%M:%S'`

if [ ! -s "${wkdir}/mtime" ]; then
  echo "$now [ERR ] Cannot find previous model time."
  exit
fi

if [ -e "${wkdir}/running" ]; then
  echo "$now [PREV]" >> ${wkdir}/get_ncep_gefs.log
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

for imem in `seq $nmem` mdet mean ;do ### mean must come later than mdet 
#for imem in mdet ;do
  if [ "$imem" == "mean" ];then
    cmem=$imem
    cmem_src='avg'
  elif [ "$imem" == "mdet" ]; then
    cmem=$imem
    cmem_src='c00' 
  else
    cmem=`printf %04d $imem`
    cmem_src='p'`printf %02d $imem`
  fi


echo "$now [TRY ] $YYYYMMDDHH $cmem" >> ${wkdir}/get_ncep_gefs.log

mkdir -p ${wkdir}/$YYYYMMDDHH/$cmem
cd ${wkdir}/$YYYYMMDDHH/$cmem

allget=1
### download using grib filter
t=0
while ((t <= $tmax)); do

  tf=`printf '%03d' $t`
  TIME_fcst=`date -ud "${t} hour $YYYY-$MM-$DD $HH" +'%Y-%m-%d %H:%M:%S'`
  YYYYMMDDHHMMSS_fcst=`date -ud "$TIME_fcst" +'%Y%m%d%H%M%S'`

  if [ ! -s "gfs_a.$YYYYMMDDHHMMSS_fcst" ]; then
    rm -f gefs_${cmem}.${YYYYMMDDHH}_f${tf}_a
    rm -f gefs_${cmem}.${YYYYMMDDHH}_f${tf}_b
#
curl "https://nomads.ncep.noaa.gov/cgi-bin/filter_gens_0p50.pl?file=ge${cmem_src}.t${HH}z.pgrb2a.0p50.f${tf}&all_lev=on&var_HGT=on&var_PRES=on&var_PRMSL=on&var_RH=on&var_SOILW=on&var_TMP=on&var_TSOIL=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=90&rightlon=180&toplat=60&bottomlat=10&dir=%2Fgefs.${YYYY}${MM}${DD}%2F${HH}%2Fpgrb2ap5" -o gefs_${cmem}_${YYYYMMDDHH}_f${tf}_a

#echo "https://nomads.ncep.noaa.gov/cgi-bin/filter_gens_0p50.pl?file=ge${cmem_src}.t${HH}z.pgrb2a.0p50.f${tf}&all_lev=on&var_HGT=on&var_PRES=on&var_PRMSL=on&var_RH=on&var_SOILW=on&var_TMP=on&var_TSOIL=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=90&rightlon=180&toplat=60&bottomlat=10&dir=%2Fgefs.${YYYY}${MM}${DD}%2F${HH}%2Fpgrb2ap5" -o gefs_${cmem}_${YYYYMMDDHH}_f${tf}_a
 
#exit 0

[ "$imem" == "mean" ] || curl "https://nomads.ncep.noaa.gov/cgi-bin/filter_gens_0p50.pl?file=ge${cmem_src}.t${HH}z.pgrb2b.0p50.f${tf}&all_lev=on&var_TSOIL=on&var_SOILW=on&var_UGRD=on&var_VGRD=on&var_TMP=on&var_HGT=on&var_RH=on&subregion=&leftlon=90&rightlon=180&toplat=60&bottomlat=10&dir=%2Fgefs.${YYYY}${MM}${DD}%2F${HH}%2Fpgrb2bp5" -o gefs_${cmem}_${YYYYMMDDHH}_f${tf}_b
    if [ -s "gefs_${cmem}_${YYYYMMDDHH}_f${tf}_a" ] ; then
      mv gefs_${cmem}_${YYYYMMDDHH}_f${tf}_a gfs_a.${YYYYMMDDHHMMSS_fcst} 
      [ "$imem" == "mean" ] || mv gefs_${cmem}_${YYYYMMDDHH}_f${tf}_b gfs_b.${YYYYMMDDHHMMSS_fcst} 
      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [GET ] $YYYYMMDDHH $cmem-> gfs.$YYYYMMDDHHMMSS_fcst" >> ${wkdir}/get_ncep_gefs.log
    else
      allget=0
    fi
  fi


### generate avg for b with wgrib2
#cat gefs_0*_${inittime}_f${fhour}_b >temp_merge
#$wgrib2_exe temp_merge -set_ensm_derived_fcst 1 $nmem -grib gefs_mean_${inittime}_f${fhour}_b 

  if ((allget == 1)); then

#      now=`date -u +'%Y-%m-%d %H:%M:%S'`
#      echo "$now [CONV] $YYYYMMDDHH $cmem -> gfs.$YYYYMMDDHHMMSS_fcst - GrADS" >> ${wkdir}/get_ncep_gefs.log
#      bash $wkdir/run/grads/convert.sh "$TIME_fcst" "$TIME_fcst" "$datadir/ncepgefs/${YYYYMMDDHH}/$cmem/gfs" \
#       > ${wkdir}/convert_grads.log 2>&1

  ###    now=`date -u +'%Y-%m-%d %H:%M:%S'`
  ###    echo "$now [CONV] $YYYYMMDDHH $cmem -> gfs.$YYYYMMDDHHMMSS_fcst - Plot (background job)" >> ${wkdir}/get_ncep_gefs.log
  ###    bash $wkdir/run/plot/plot.sh "$TIME_fcst" "$GET_TIME" $((t/6+1)) "$datadir/ncepgefs/${YYYYMMDDHH}/$cmem/gfs" \
  ###     > ${wkdir}/convert_plot.log 2>&1 &

      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [CONV] $YYYYMMDDHH $cmem -> gfs.$YYYYMMDDHHMMSS_fcst - GrADS for SCALE" >> ${wkdir}/get_ncep_gefs.log
      bash $wkdir/run/grads_scale/convert.sh "$TIME_fcst" "$TIME_fcst" "$datadir/ncepgefs/${YYYYMMDDHH}/$cmem/gfs" "$datadir/ncepgefs_grads/${YYYYMMDDHH}/$cmem" \
       > ${wkdir}/convert_grads_scale.log 2>&1
  fi

t=$((t+6))

done ### tf

done ### imem

if ((allget == 1)); then
  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [DONE] $YYYYMMDDHH $cmem" >> ${wkdir}/get_ncep_gefs.log
  echo "$GET_TIME" > ${wkdir}/mtime
fi

rm -f ${wkdir}/running
