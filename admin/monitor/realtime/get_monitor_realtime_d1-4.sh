#!/bin/bash -l


myname=$0
mydir=`dirname $myname`
cd $mydir

. $mydir/../../admin.rc || exit $?

time_offset=0
[ -f "$mydir/../../time_offset.txt" ] && time_offset=`cat $mydir/../../time_offset.txt` 
nowtime="$(date -ud "$time_offset second now" +'%Y-%m-%d %H:%M:%S')"

monitor_dir="$mydir/.."
OUTBASE="${mydir}/../../../result/ope"

prev=`cat data/d1-4.txt | tail -n 1`
if [ ! -z "$prev" ] ;then 
timeref_p=`echo $prev | awk '{print $1}'`
timeobs_p=`echo $prev | awk '{print $2}'`
timegfs_p=`echo $prev | awk '{print $3}'`
timed1_p=`echo $prev | awk '{print $4}'`
timed2_p=`echo $prev | awk '{print $5}'`
timed3_p=`echo $prev | awk '{print $6}'`
timed4_p=`echo $prev | awk '{print $7}'`
else
dummy=`date -ud "-5 days $nowtime" +%Y%m%d%H%M%S`
timeref_p=$dummy
timeobs_p=$dummy
timegfs_p=$dummy
timed1_p=$dummy
timed2_p=$dummy
timed3_p=$dummy
timed4_p=$dummy
fi


nowsec=`date -ud "$nowtime" +%s` 
nowsec_t=`expr $nowsec / 600 \* 600 - 32400`

timeref=`date -d @$nowsec_t +%Y%m%d%H%M%S`

if [ "$timeref_p" == "$timeref" ];then
 mv data/d1-4.txt temp
 sed '$d' temp > data/d1-4.txt
 rm temp
fi


timeobs=`tail -n 1 ${monitor_dir}/monitor_obs.txt | awk '{print $1}'`0000
timegfs=`tail -n 1 ${monitor_dir}/monitor_gfs.txt | awk '{print $1}'`0000
timed1=`tail -n 1 ${monitor_dir}/monitor_cycle_temp.txt | awk '{print $1}'`0000

[ "$timeobs" -eq "0000" ] && timeobs=$timeobs_p ### TORI AEZU
[ "$timegfs" -eq "0000" ] && timegfs=$timegfs_p ### TORI AEZU
[ "$timed1" -eq "0000" ] && timed1=$timed1_p ### TORI AEZU

timed2=$timed2_p

latest=`ls -1 ${OUTBASE}/d2/20*/fcstgpi/mdet/sfc_prcp_f000000.png | tail -n 1`

if [ ! -z "$latest" ];then
based2=`echo $latest | egrep --only-matching [0-9]{14}`
testdir=`dirname $latest`
latest=`ls -1 $testdir/sfc_prcp_f*.png | tail -n 1`
lend2=`basename $latest | sed -e 's/[^0-9]//g' `

yyyy=`echo $based2 | cut -c 1-4`
mm=`echo $based2 | cut -c 5-6`
dd=`echo $based2 | cut -c 7-8`
hh=`echo $based2 | cut -c 9-10`
mon=`echo $based2 | cut -c 11-12`
sec=`echo $based2 | cut -c 13-14`

timed2=`date -d "$lend2 sec ${yyyy}-${mm}-${dd} ${hh}:${mon}:${sec}" +%Y%m%d%H%M%S`
fi

timed3=$timed3_p

list=`ls -1td ${OUTBASE}/d3/ref_20*/20*/fcstgpi/mdet | head -n 10`

if [ ! -z "$list" ];then
for path in $list;do
latest=`ls -1 $path/sfc_prcp_*.png | tail -n 1` 

based3=`echo $latest | egrep --only-matching /[0-9]{14} | cut -c 2-15`
lend3=`basename $latest | egrep --only-matching f[0-9]{6} | cut -c 2-7`

yyyy=`echo $based3 | cut -c 1-4`
mm=`echo $based3 | cut -c 5-6`
dd=`echo $based3 | cut -c 7-8`
hh=`echo $based3 | cut -c 9-10`
mon=`echo $based3 | cut -c 11-12`
sec=`echo $based3 | cut -c 13-14`
timed3_new=`date -d "$lend3 sec ${yyyy}-${mm}-${dd} ${hh}:${mon}:${sec}" +%Y%m%d%H%M%S`
[ $timed3_new -gt $timed3 ] && timed3=$timed3_new
done
fi

timed4=$timed4_p

list=`ls -1t ${OUTBASE}/d3/ref_20*/20*/d4_${dx_d4}/anal/mean/init_*.pe000000.nc | head -n 10`
if [ ! -z "$list" ];then
for file in $list;do
based3=`echo $file | egrep --only-matching "init_[0-9]{8}-[0-9]{6}" | cut -c 6-20`
yyyy=`echo $based3 | cut -c 1-4`
mm=`echo $based3 | cut -c 5-6`
dd=`echo $based3 | cut -c 7-8`
hh=`echo $based3 | cut -c 10-11`
min=`echo $based3 | cut -c 12-13`
sec=`echo $based3 | cut -c 14-15`
timed4_new=`date -ud "1 hour ${yyyy}-${mm}-${dd} ${hh}:${min}:${sec}" +%Y%m%d%H%M%S` ### hourly data

[ $timed4_new -gt $timed4 ] && timed4=$timed4_new
done
fi


#echo $timeref $timeobs $timegfs $timed1 $timed2 $timed3 $timed4
#exit 0
echo $timeref $timeobs $timegfs $timed1 $timed2 $timed3 $timed4 >> data/d1-4.txt
tac data/d1-4.txt > data_inv/d1-4.txt 

[ -f draw_d1-4 ] || dclfrt draw_d1-4.f90 -o draw_d1-4
./draw_d1-4 `date -ud "9 hour $nowtime" +%Y%m%d%H%M%S`
mv ./figs/d1-4_0001.png ./figs/realtime_d1-4.png 
mogrify -trim ./figs/realtime_d1-4.png 

if [ ${timeref:8:6} -eq "000000" ]; then ### archive
 cp figs/realtime_d1-4.png save/realtime_d1-4_${timeref:0:8}.png 
fi
