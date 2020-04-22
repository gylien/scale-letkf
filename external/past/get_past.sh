#!/bin/sh

myname=$0
wkdir=`dirname $myname`
txtfile=$wkdir/../../admin/time_offset.txt


#nows=`date -ud "2020-04-23 10:00:00" +%s`
#pasts=`date -ud "2019-09-10 18:00:00" +%s`
#time_offset=`expr $pasts - $nows`
#echo $time_offset > $txtfile
#exit 0



lag_obs="6 hour"
lag_gfs="4 hour"
lag_msm="6 hour"

if [ -f $txtfile ]; then
  time_offset=`cat $time_offset`
  pasttime=`date -ud "$time_offset now" +"%Y-%m-%d %H:%M:%S"`
  for obsdir in `/bin/ls -xd $wkdir/ncepobs_gdas_letkf/*`;do
    obstime=`date -ud "${obsdir:0:4}-${obsdir:4:2}-${obsdir:6:2} ${obsdir:8:2}"`
    if [ `date -ud "$obstime" +%s` -lt `date -ud "- $lag_obs $pasttime" +%s` ];then
      mv $wkdir/ncepobs_gdas_letkf/$obstime $wkdir/../ncepobs_gdas_letkf/
    fi    
  done
  for gfsdir in `/bin/ls -xd $wkdir/ncepgfs_grads/*`;do
    gfstime=`date -ud "${gfsdir:0:4}-${gfsdir:4:2}-${gfsdir:6:2} ${gfsdir:8:2}"`
    if [ `date -ud "$gfstime" +%s` -lt `date -ud "- $lag_gfs $pasttime" +%s` ];then
      mv $wkdir/ncepgfs_grads/$gfstime $wkdir/../ncepgfs_grads/
    fi    
  done
  for msmdir in `/bin/ls -xd $wkdir/msm/*`;do
    msmtime=`date -ud "${msmdir:0:4}-${msmdir:4:2}-${msmdir:6:2} ${msmdir:8:2}"`
    if [ `date -ud "$msmtime" +%s` -lt `date -ud "- $lag_msm $pasttime" +%s` ];then
      mv $wkdir/msm/$msmtime $wkdir/../msm/
    fi    
  done
fi
