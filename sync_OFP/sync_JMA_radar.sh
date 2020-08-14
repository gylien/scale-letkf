#!/bin/sh

if [ -f /home/amemiya/ssh-agent-log.sh ] ;then
source /home/amemiya/ssh-agent-log.sh

sshcommand="ssh -i $HOME/.ssh/id_rsa_daweb"
hostname=macallan
JMADIR_SRC=/data9/amemiya/JMA_precip/images/
JMADIR=$HOME/public_html/scale/data/JMA_precip

oldest=`date -ud "-2 day" +%s`
YMD=`date -ud "0 day" +%Y%m%d`
YMDb=`date -ud "-1 day" +%Y%m%d`
YMDbb=`date -ud "-2 day" +%Y%m%d`

for d in d2 d3 d3_grads d4_Tokyo; do
 rsync -e "$sshcommand" $hostname:$JMADIR_SRC/$d/$YMDbb/radar_*.png $JMADIR/nowcast_$d/realtime/  
 rsync -e "$sshcommand" $hostname:$JMADIR_SRC/$d/$YMDb/radar_*.png $JMADIR/nowcast_$d/realtime/  
 rsync -e "$sshcommand" $hostname:$JMADIR_SRC/$d/$YMD/radar_*.png $JMADIR/nowcast_$d/realtime/  
# rsync -e "$sshcommand" $hostname:$JMADIR_SRC/$d/nowcast_*.png $JMADIR/nowcast_$d/realtime/  
# rm $JMADIR/nowcast_$d/realtime/nowcast_??_0001.png
# rm $JMADIR/nowcast_$d/realtime/*_tmp.png
 chmod 644 $JMADIR/nowcast_$d/realtime/*.png
# list=`ls -x $JMADIR/nowcast_$d/realtime/nowcast_??????????????.png`
# for fname in $list ;do
#  file=`basename $fname`
#  ftime=${file:8:14}
#  ftimes=`date -d "${ftime:0:4}-${ftime:4:2}-${ftime:6:2} ${ftime:8:2}:${ftime:10:2}:${ftime:12:2}" +%s`
#  [ $ftimes -le $oldest ] && rm $fname 
#  [ -f $JMADIR/nowcast_$d/realtime/radar_${ftime}.png ] && rm $fname
# done
 list=`ls -x $JMADIR/nowcast_$d/realtime/radar_??????????????.png`
 for fname in $list ;do
  file=`basename $fname`
  ftime=${file:6:14}
  ftimes=`date -d "${ftime:0:4}-${ftime:4:2}-${ftime:6:2} ${ftime:8:2}:${ftime:10:2}:${ftime:12:2}" +%s`
  if [ $ftimes -le $oldest ] then
    rm $fname 
  fi
 done
done

cd $HOME/public_html/scale/
php ./monitor_plot_d2.php
php ./monitor_plot_d3.php
#php ./monitor_plot_d4.php
php ./monitor_plot_JMAprecip.php


else
  echo 'ssh-agent is not running.'
fi
