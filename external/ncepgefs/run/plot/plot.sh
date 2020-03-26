#!/bin/bash -l

tplot="$1"
tfcstbase="$2"
tfcst="$3"
gribfile_prefix="$4"

#tplot='2015-06-01 00'
#tfcstbase='2015-06-01 00'
#tfcst=1
#gribfile_prefix='/data7/gylien/realtime/ncepgfs/2015060100/gfs'

#figlist='sfc_prcp sfc_wind sfc_2mtemp sfc_temp 850_temp 700_vvel 500_vort 300_wspd olr_ir max_ref'
figlist='sfc_prcp sfc_wind sfc_2mtemp sfc_temp 850_temp 850_theq 700_vvel 500_vort 300_wspd'

gribfile_dir="$(dirname ${gribfile_prefix})"
gribfile_base="$(basename ${gribfile_prefix})"
tint=21600

#web_host="daweb"
#web_remote_dir="/home/amemiya/public_html/scale/data/gfs"

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p $wkdir/out

timef=$(date -ud "$tplot" '+%Y-%m-%d %H')
timef2=$(date -ud "$tplot" '+%Y%m%d%H%M%S')
echo "[$timef]"

###
###timebasef=$(date -ud "$tfcstbase" '+%HZ%d%b%Y' | tr '[a-z]' '[A-Z]')
###
timebasef=$(date -ud "$tfcstbase" '+%HZ%d%m%Y' | tr '[a-z]' '[A-Z]')

timebasef2=$(date -ud "$tfcstbase" '+%Y%m%d%H%M%S')


#echo   ${gribfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}
#exit




###grads -blc "plot_driver_6h.gs ${gribfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}" #\

### Use GrADS 2.0.2 -- not 2.1.1
echo "plot_driver_6h.gs ${gribfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}"
grads -blc "plot_driver_6h.gs ${gribfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}" 
#/apps/SLES11/opt/grads/2.0.2/bin/grads -blc "plot_driver_6h.gs ${gribfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}" 
#      > /dev/null 2>&1
#      > plot_driver_6h.log 2>&1

#----

outdir="${gribfile_dir}/plot"
mkdir -p $outdir

itfcstf=$(printf '%06d' $(((tfcst-1) * tint)))



for prefix in $figlist; do
  if [ -s "$wkdir/out/${prefix}_f${itfcstf}.png" ]; then
#    ssh -i /home/amemiya/.ssh/id_rsa_mac $web_host "mkdir -p ${web_remote_dir}/${timebasef2}/${prefix}"

### TORI AEZU 
    mv $wkdir/out/${prefix}_f${itfcstf}.png $wkdir/out/${prefix}_f${itfcstf}_tmp.png
    convert  $wkdir/out/${prefix}_f${itfcstf}_tmp.png -rotate 90 -resize 757x520! $wkdir/out/${prefix}_f${itfcstf}.png

#    rsync  -e "ssh -i /home/amemiya/.ssh/id_rsa_mac" -avz $wkdir/out/${prefix}_f${itfcstf}.png ${web_host}:${web_remote_dir}/${timebasef2}/${prefix}
    mv -f $wkdir/out/${prefix}_f${itfcstf}.png $outdir
  fi


  if [ -s "$wkdir/out/${prefix}_f${itfcstf}.eps" ]; then
#    ssh $web_host "mkdir -p ${web_remote_dir}/${timebasef2}/${prefix}"
#    rsync -avz $wkdir/out/${prefix}_f${itfcstf}.eps ${web_host}:${web_remote_dir}/${timebasef2}/${prefix}
    mv -f $wkdir/out/${prefix}_f${itfcstf}.eps $outdir
  fi
done
