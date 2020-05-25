#!/bin/bash -l

tplot="$1"
tfcstbase="$2"
tfcst="$3"
gradsfile_dir="$4"

#tplot='2019-03-14 00'
#tfcstbase='2019-03-14 00'
#tfcst=1
#gradsfile_dir='/data9/amemiya/realtime_OFP/data/msm/2019031400'

#figlist='sfc_prcp sfc_wind sfc_2mtemp sfc_temp 850_temp 700_vvel 500_vort 300_wspd olr_ir max_ref'

figlist='sfc_prcp sfc_wind sfc_2mtemp 925_temp 925_theq 850_temp 850_theq 700_temp 500_temp'


gradsfile_base_s=msm_s
gradsfile_base_p=msm_p
tint=3600

web_host="daweb"
web_remote_dir="/home/amemiya/public_html/scale/data/msm"

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p $wkdir/out

timef=$(date -ud "$tplot" '+%Y-%m-%d %H')
timef2=$(date -ud "$tplot" '+%Y%m%d%H%M%S')
echo "[$timef]"


export LANG=C
###


#timebasef=$(date -ud "$tfcstbase" '+%HZ%d%m%Y' | tr '[a-z]' '[A-Z]')
#timeplotf=$(date -ud "$tplot" '+%HZ%d%m%Y' | tr '[a-z]' '[A-Z]')

timebasef=$(date -ud "$tfcstbase" '+%HZ%d%b%Y' | tr '[a-z]' '[A-Z]')
timeplotf=$(date -ud "$tplot" '+%HZ%d%b%Y' | tr '[a-z]' '[A-Z]')

timebasef2=$(date -ud "$tfcstbase" '+%Y%m%d%H%M%S')


#echo   ${gradsfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}
#exit




###grads -blc "plot_driver_6h.gs ${gradsfile_prefix}.${timef2}.ctl ${timebasef} ${tfcst}" #\
#echo ${gradsfile_prefix}.${timef2}.ctl

#exit 0

### Use GrADS 2.0.2 -- not 2.1.1

#/apps/SLES11/opt/grads/2.0.2/bin/grads -blc "plot_driver_6h_s.gs ${gradsfile_dir}/${gradsfile_base_s}.${timebasef2}.ctl ${timebasef} ${timeplotf} ${tfcst}" 
#/apps/SLES11/opt/grads/2.0.2/bin/grads -blc "plot_driver_6h_p.gs ${gradsfile_dir}/${gradsfile_base_p}.${timebasef2}.ctl ${timebasef} ${timeplotf} ${tfcst}" 

grads -blc "plot_driver_1h_s.gs ${gradsfile_dir}/${gradsfile_base_s}.${timebasef2}.ctl ${timebasef} ${timeplotf} ${tfcst}"   > plot_driver_1h_s.log 2>&1

if [ `expr $tfcst % 3` == 1 ];then
 tfcstp=`expr $tfcst \/ 3 + 1`
grads -blc "plot_driver_3h_p.gs ${gradsfile_dir}/${gradsfile_base_p}.${timebasef2}.ctl ${timebasef} ${timeplotf} ${tfcstp}"   > plot_driver_3h_p.log 2>&1
fi


#      > /dev/null 2>&1

#----

outdir="${gradsfile_dir}/plot"
mkdir -p $outdir

itfcstf=$(printf '%06d' $(((tfcst-1) * tint)))

for prefix in $figlist; do
  if [ -s "$wkdir/out/${prefix}_f${itfcstf}.png" ]; then
#    ssh -i /home/amemiya/.ssh/id_rsa_mac $web_host "mkdir -p ${web_remote_dir}/${timebasef2}/${prefix}"

### TORI AEZU 
###    mv $wkdir/out/${prefix}_f${itfcstf}.png $wkdir/out/${prefix}_f${itfcstf}_tmp.png
###    convert  -resize 757x520! $wkdir/out/${prefix}_f${itfcstf}_tmp.png $wkdir/out/${prefix}_f${itfcstf}.png

#    rsync  -e "ssh -i /home/amemiya/.ssh/id_rsa_mac" -avz $wkdir/out/${prefix}_f${itfcstf}.png ${web_host}:${web_remote_dir}/${timebasef2}/${prefix}
    mv -f $wkdir/out/${prefix}_f${itfcstf}.png $outdir
  fi
  if [ -s "$wkdir/out/${prefix}_f${itfcstf}.eps" ]; then
#    ssh $web_host "mkdir -p ${web_remote_dir}/${timebasef2}/${prefix}"
#    rsync -avz $wkdir/out/${prefix}_f${itfcstf}.eps ${web_host}:${web_remote_dir}/${timebasef2}/${prefix}
    mv -f $wkdir/out/${prefix}_f${itfcstf}.eps $outdir
  fi
done
