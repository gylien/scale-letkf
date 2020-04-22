#!/bin/sh


myname=$0
mydir=`dirname $myname`
cd $mydir

time_offset=0
[ -f "../../time_offset.txt" ] && time_offset=`cat ../../time_offset.txt` 
nowtime="$(date -ud "$time_offset second now" +'%Y-%m-%d %H:%M:%S')"


monitor_dir="$mydir/.."
OUTBASE="${mydir}/../../../result/ope"

nowsec=`date -ud "$nowtime" +%s` 
nowsec_t=`expr $nowsec / 600 \* 600 - 32400`

timeref=`date -d @$nowsec_t +%Y%m%d%H%M%S`

timeobs=`tail -n 1 ${monitor_dir}/monitor_obs.txt | awk '{print $1}'`0000
timegfs=`tail -n 1 ${monitor_dir}/monitor_gfs.txt | awk '{print $1}'`0000
timed1=`tail -n 1 ${monitor_dir}/monitor_cycle_temp.txt | awk '{print $1}'`0000

[ "$timed1" -eq "0000" ] && timed1=`date -d " -3 days $nowtime" +%Y%m%d%H%M%S` ### TORI AEZU


latest=`ls -1 ${OUTBASE}/d2/*/fcstgpi/mdet/sfc_prcp_* | tail -n 1`
based2=`echo $latest | egrep --only-matching [0-9]{14}`
lend2=`basename $latest | sed -e 's/[^0-9]//g' `


yyyy=`echo $based2 | cut -c 1-4`
mm=`echo $based2 | cut -c 5-6`
dd=`echo $based2 | cut -c 7-8`
hh=`echo $based2 | cut -c 9-10`
mon=`echo $based2 | cut -c 11-12`
sec=`echo $based2 | cut -c 13-14`


timed2=`date -d "$lend2 sec ${yyyy}-${mm}-${dd} ${hh}:${mon}:${sec}" +%Y%m%d%H%M%S`


latest=`ls -1t ${OUTBASE}/d3/ref_*/*/plot/mean/rain_*.png | head -n 1`
based3=`echo $latest | egrep --only-matching /[0-9]{14} | cut -c 2-15`
lend3=`basename $latest | egrep --only-matching f[0-9]{6} | cut -c 2-7`


yyyy=`echo $based3 | cut -c 1-4`
mm=`echo $based3 | cut -c 5-6`
dd=`echo $based3 | cut -c 7-8`
hh=`echo $based3 | cut -c 9-10`
mon=`echo $based3 | cut -c 11-12`
sec=`echo $based3 | cut -c 13-14`
timed3=`date -d "$lend3 sec ${yyyy}-${mm}-${dd} ${hh}:${mon}:${sec}" +%Y%m%d%H%M%S`


echo $timeref $timeobs $timegfs $timed1 $timed2 $timed3 >> data/d1-3.txt
tac data/d1-3.txt > data_inv/d1-3.txt 

./draw_d1-3
mv ./figs/d1-3_0001.png ./figs/realtime_d1-3.png 
mogrify -trim ./figs/realtime_d1-3.png 

if [ ${timeref:8:6} -eq "000000" ]; then ### archive
 mv figs/realtime_d1-3.png save/realtime_d1-3_${timeref:0:8}.png 
fi
