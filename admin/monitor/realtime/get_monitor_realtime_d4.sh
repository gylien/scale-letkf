#!/bin/sh

myname=$0
mydir=`dirname $myname`
cd $mydir

time_offset=0
[ -f "../../time_offset.txt" ] && time_offset=`cat ../../time_offset.txt` 
nowtime="$(date -ud "$time_offset second now" +'%Y-%m-%d %H:%M:%S')"


nowsec=`date -ud "$nowtime" +%s` 
nowsec_t=`expr $nowsec / 5 \* 5 - 32400`

timeref=`date -d @$nowsec_t +%Y%m%d%H%M%S`

timeobs=`ssh -i ~/.ssh/id_rsa_mac daweb "cd public_html/HPCC_scale/data/d4/realtime ; ls -1t obs_* | tail -n 1 | sed -e 's/[^0-9]//g'" `
timeanal=`ssh -i ~/.ssh/id_rsa_mac daweb "cd public_html/HPCC_scale/data/d4/realtime ; ls -1t anal_* | tail -n 1 | sed -e 's/[^0-9]//g'" `
timefcst=`ssh -i ~/.ssh/id_rsa_mac daweb "cd public_html/HPCC_scale/data/d4/realtime ; ls -1t fcst_* | tail -n 1 | sed -e 's/[^0-9]//g'" `


[ "$timeobs" == "" ] && timeobs=`date -d " -1 days $nowtime" +%Y%m%d%H%M%S` ### TORI AEZU
[ "$timeanal" == "" ] && timeanal=`date -d " -1 days $nowtime" +%Y%m%d%H%M%S` ### TORI AEZU
[ "$timefcst" == "" ] && timefcst=`date -d " -1 days $nowtime" +%Y%m%d%H%M%S` ### TORI AEZU


echo $timeref $timeobs $timeanal $timefcst >> data/d4.txt
tac data/d4.txt > data_inv/d4.txt 

./draw_d4
mv ./figs/d4_0001.png ./figs/realtime_d4.png 
mogrify -trim ./figs/realtime_d4.png 

scp -i ~/.ssh/id_rsa_mac figs/realtime_d4.png daweb:public_html/scale/

#if [ ${timeref:8:6} -eq "000000" ]; then ### archive
# mv fig/realtime_d1-3.png save/realtime_d4_${timeref:0:8}.png 
#fi
