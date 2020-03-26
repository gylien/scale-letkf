#!/bin/bash

tstart="$1"
tend="$2"
gribfile_prefix="$3"
outdir="$4"

#tstart='2016-08-10 06'
#tend='2016-08-10 06'
#gribfile_prefix='/data7/gylien/realtime/ncepgfs/2016081000/gfs'
#outdir='/data7/gylien/realtime/ncepgfs_grads/2016081000'

gribfile_dir="$(dirname ${gribfile_prefix})"
gribfile_base="$(basename ${gribfile_prefix})"
gribfile_dir_gfs=`echo $gribfile_dir | sed -e "s/ncepgefs/ncepgfs/g" | rev | cut -c 6- | rev`

tint=21600

region='90:180 10:60'
nlev_gfs=31
nlev_gefs=26
levels_mb=(   1000  975  950  925  900  850  800  750  700  650  600  550  500  450  400  350  300  250  200  150  100   70   50   30   20   10    7    5    3    2    1 ) 
levels_uv_a=(    1    0    0    1    0    1    0    0    1    0    0    0    1    0    1    0    1    1    1    0    1    0    1    0    0    1    0    0    0    0    0) 
levels_g_a=(     1    0    0    1    0    1    0    0    1    0    0    0    1    0    0    0    1    1    1    0    1    0    1    0    0    1    0    0    0    0    0) 
levels_tr_a=(    1    0    0    1    0    1    0    0    1    0    0    0    1    0    0    0    0    1    1    0    1    0    1    0    0    1    0    0    0    0    0) 


#region_levs='1000 925 850 700 500 300 200 100 50 20 10'
#region_latlon='95.5:174.5 12.5:54.5'


#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p ${outdir}

export LANG=en
stimefgrads=$(date -ud "$tstart" '+%H:%MZ%d%b%Y')

cat ${wkdir}/ctl/land.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/land.ctl
cat ${wkdir}/ctl/sfc.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/sfc.ctl
cat ${wkdir}/ctl/atm.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/atm.ctl
#cat ${wkdir}/ctl/reg.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/reg.ctl

time="$tstart"
###timef2s=$(date -ud "$tstart" '+%Y%m%d%H%M%S')
timef2s=`echo $gribfile_dir | grep -o '[0-9]\{10\}'`0000

while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do

  timef=$(date -ud "$time" '+%Y-%m-%d %H')
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  echo "[$timef]"

  gfile_topo=${gribfile_prefix}_a.${timef2s}
  gfile_a=${gribfile_prefix}_a.${timef2}
  gfile_b=${gribfile_prefix}_b.${timef2}
 

###
  gfile_b=`echo $gfile_b | sed "s/mean/mdet/g"` ### TORI AEZU 
###

  ofile_land=${outdir}/land_${timef2}.grd
  ofile_sfc=${outdir}/sfc_${timef2}.grd
  ofile_atm=${outdir}/atm_${timef2}.grd
  ofile_reg=${outdir}/reg_${timef2}.grd
  
  gfile_gfs=${gribfile_dir_gfs}/${gribfile_base}.${timef2}

#----
  tfile="tmp.grb2"
  rm -f ${tfile} ${ofile_land}
  #--Land data
#  wgrib2 ${gfile_b} -match ":LAND:"  -match ":surface:"   -grib ${tfile}
  wgrib2 ${gfile_gfs} -match ":LAND:"  -match ":surface:"   -small_grib ${region} ${tfile}
  wgrib2 ${gfile_b} -match ":TMP:"   -match ":surface:"     -append -grib ${tfile}

  wgrib2 ${gfile_a} -match ":TSOIL:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile_b} -match ":TSOIL:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile_a} -match ":SOILW:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile_b} -match ":SOILW:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_land} -g2clib 0

  rm -f ${tfile} ${ofile_sfc}
  #--Surface data
#  wgrib2 ${gfile} -match ":PRES:"  -match ":mean sea level:"            -grib ${tfile}
#  wgrib2 ${gfile} -match ":MSLET:"  -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile_a} -match ":PRMSL:" -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile_a} -match ":PRES:"  -match ":surface:"           -append -grib ${tfile}
  wgrib2 ${gfile_a} -match ":UGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile_a} -match ":VGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile_a} -match ":TMP:"   -match ":2 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile_a} -match ":RH:"    -match ":2 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile_topo} -match ":HGT:"   -match ":surface:"        -append -grib ${tfile} ### available only in analysis
  wgrib2 ${tfile} -no_header -ieee ${ofile_sfc} -g2clib 0

  rm -f ${tfile} ${ofile_atm}
  #--Upper data
  lev_indx=0
  while [ $lev_indx -lt $nlev_gfs ] ;do
    level="${levels_mb[$lev_indx]}"
    if [ $lev_indx -lt $nlev_gefs ] ;then
      [ ${levels_g_a[$lev_indx]} == 1 ] && gfile="$gfile_a" || gfile="$gfile_b"
      wgrib2 ${gfile} -match ":HGT:" -match ":${level} mb" -append -grib ${tfile}
    else
#      wgrib2 ${gfile_gfs} -match ":HGT:" -match ":${level} mb" -append -small_grib ${region} ${tfile}
      wgrib2 ${gfile_gfs} -match ":HGT:" -match ":${level} mb" -small_grib ${region} tmp_small.grb2 ### small_grib cannot into append (why?)
      wgrib2 tmp_small.grb2 -append -grib $tfile
    fi
    lev_indx=`expr $lev_indx + 1`
  done

  lev_indx=0
  while [ $lev_indx -lt $nlev_gfs ] ;do
    level="${levels_mb[$lev_indx]}"
    if [ $lev_indx -lt $nlev_gefs ] ;then
      [ ${levels_uv_a[$lev_indx]} == 1 ] && gfile="$gfile_a" || gfile="$gfile_b"
      wgrib2 ${gfile} -match ":UGRD:"  -match "mb" -not "mb above ground" -match ":${level} mb" -append -grib ${tfile}
    else
#      wgrib2 ${gfile_gfs} -match ":UGRD:"   -match "mb" -not "mb above ground" -match ":$level mb" -append -small_grib ${region} ${tfile}
      wgrib2 ${gfile_gfs} -match ":UGRD:" -match ":${level} mb" -small_grib ${region} tmp_small.grb2
      wgrib2 tmp_small.grb2 -append -grib $tfile
    fi
    lev_indx=`expr $lev_indx + 1`
  done

  lev_indx=0
  while [ $lev_indx -lt $nlev_gfs ] ;do
    level="${levels_mb[$lev_indx]}"
    if [ $lev_indx -lt $nlev_gefs ] ;then
      [ ${levels_uv_a[$lev_indx]} == 1 ] && gfile="$gfile_a" || gfile="$gfile_b"
      wgrib2 ${gfile} -match ":VGRD:" -match ":${level} mb" -append -grib ${tfile}
    else
      wgrib2 ${gfile_gfs} -match ":VGRD:" -match ":${level} mb" -small_grib ${region} tmp_small.grb2
      wgrib2 tmp_small.grb2 -append -grib $tfile
     fi
    lev_indx=`expr $lev_indx + 1`
  done
 
  lev_indx=0
  while [ $lev_indx -lt $nlev_gfs ] ;do
    level="${levels_mb[$lev_indx]}"
    if [ $lev_indx -lt $nlev_gefs ] ;then
      [ ${levels_tr_a[$lev_indx]} == 1 ] && gfile="$gfile_a" || gfile="$gfile_b"
      wgrib2 ${gfile} -match ":TMP:" -match ":${level} mb" -append -grib ${tfile}
    else
      wgrib2 ${gfile_gfs} -match ":TMP:" -match ":${level} mb" -small_grib ${region} tmp_small.grb2
      wgrib2 tmp_small.grb2 -append -grib $tfile
    fi
    lev_indx=`expr $lev_indx + 1`
  done
 
  lev_indx=0
  while [ $lev_indx -lt $nlev_gfs ] ;do
    level="${levels_mb[$lev_indx]}"
    if [ $lev_indx -lt $nlev_gefs ] ;then
      [ ${levels_tr_a[$lev_indx]} == 1 ] && gfile="$gfile_a" || gfile="$gfile_b"
      wgrib2 ${gfile} -match ":RH:"  -match ":${level} mb" -append -grib ${tfile}
    else
      wgrib2 ${gfile_gfs} -match ":RH:" -match ":${level} mb" -small_grib ${region} tmp_small.grb2
      wgrib2 tmp_small.grb2 -append -grib $tfile
    fi
    lev_indx=`expr $lev_indx + 1`
  done

  wgrib2 ${tfile} -no_header -ieee ${ofile_atm} -g2clib 0

time=$(date -ud "$tint second $time" '+%Y-%m-%d %H')

done
