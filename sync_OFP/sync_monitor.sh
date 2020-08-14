#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc


echo "check ssh-add ..."
echo `ssh-add -l`


TXTDIR_OFP=$SCALEDIR_OFP/admin/monitor
TXTDIR=$mydir/../monitor

cycle_old=`tail -n 1 $TXTDIR/monitor_cycle_temp.txt | cut -c 9-10`

items="monitor_cycle.txt monitor_cycle_temp.txt monitor_fcst.txt monitor_fcst_d2.txt monitor_fcst_d3.txt monitor_gfs.txt monitor_obs.txt"

for item in $items ;do
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/$item $TXTDIR/  
done

### LETKF monitoring
for var in ps q t u v ;do
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/letkf/figs/${var}_letkf_0001.png $TXTDIR/letkf_monitor/all
done
 chmod 644 $TXTDIR/letkf_monitor/all/*.png 

cycle_new=`tail -n 1 $TXTDIR/monitor_cycle_temp.txt | cut -c 9-10`

### LETKF detailed monitoring


if [ "$cycle_old" != "$cycle_new" ] ;then
  if [ "$cycle_new" == "00" ] || [ "$cycle_new" == "12" ];then
for imem in `seq 0 191` ;do
 cmem=`printf %03d $imem`
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/letkf/subdomain/figs/$cmem/ $TXTDIR/letkf_monitor/sub/$cmem/
done
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/letkf/subdomain/figs/map/ $TXTDIR/letkf_monitor/sub/map/
 chmod 644 $TXTDIR/letkf_monitor/sub/*/*.png 
  fi
fi

### realtime progress monitoring
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/realtime/figs/realtime_d1-3.png $TXTDIR/
 rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/realtime/figs/realtime_d1-4.png $TXTDIR/

chmod 644 $TXTDIR/realtime_d1-4.png 

cd $mydir/..
php ./monitor_plot_d1.php
php ./monitor_plot_d2.php
php ./monitor_plot_d3.php

res=`$sshcommand $hostname "ls $TXTDIR_OFP/../time_offset.txt"`
if [ ! -z "$res" ] ;then
   rsync -e "$sshcommand" $hostname:$TXTDIR_OFP/../time_offset.txt $TXTDIR/../
fi
