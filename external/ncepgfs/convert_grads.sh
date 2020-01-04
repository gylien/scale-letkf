#!/bin/bash

tstart="2016-04-02 18:00:00"
tend="2016-04-02 18:00:00"
tint=21600

time=`date -ud "$tstart" +'%Y-%m-%d %H:%M:%S'`
while [ `date -ud "${time}" +'%s'` -le `date -ud "$tend" +'%s'` ]; do

  yyyymmddhhiiss=`date -ud "${time}" +'%Y%m%d%H%M%S'`
  yyyymmddhh=`date -ud "${time}" +'%Y%m%d%H'`

  echo "[${yyyymmddhh}]"

  cd $yyyymmddhh
  grads -blc "../convert_grads.gs gfs.${yyyymmddhhiiss}.ctl"
  mv grads.fwrite gfs.${yyyymmddhhiiss}.grd
  cd ..

time=`date -ud "$tint second ${time}" +'%Y-%m-%d %H:%M:%S'`
done
