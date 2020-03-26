#!/bin/bash

tstart="$1"
tend="$2"
gribfile_prefix="$3"

#tstart='2015-06-01 00'
#tend='2015-06-06 00'
#gribfile_prefix='/data7/gylien/realtime/ncepgfs/2015060100/gfs'

gribfile_dir="$(dirname ${gribfile_prefix})"
gribfile_base="$(basename ${gribfile_prefix})"
tint=21600

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do

  timef=$(date -ud "$time" '+%Y-%m-%d %H')
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  echo "[$timef]"

  cd $gribfile_dir
  ${wkdir}/g2ctl ${gribfile_base}.${timef2} > ${gribfile_base}.${timef2}.ctl
  gribmap -i ${gribfile_base}.${timef2}.ctl

  if [ "$(basename ${gribfile_dir})" = "${timef2:0:10}" ]; then
    grads -blc "${wkdir}/convert_grads.gs ${gribfile_base}.${timef2}.ctl"
    mv -f grads.fwrite ${gribfile_base}.${timef2}.grd
  fi

time=$(date -ud "$tint second $time" '+%Y-%m-%d %H')
done
