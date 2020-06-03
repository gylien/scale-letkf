#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc

<<<<<<< HEAD
GFSDIR_OFP=$SCALEDIR_OFP/external/ncepgfs
=======
GFSDIR_OFP=$SCALEDIR/external/ncepgfs
>>>>>>> 504b3ad... initial commit
GFSDIR=$DATADIR/gfs

testfile="sfc_wind_f432000.png"

timeget_OFP=`$sshcommand $hostname "cat $GFSDIR_OFP/mtime"`
timeget=`date -d "$timeget_OFP" +%Y%m%d%H%M%S`

timeget_in=$1
[ ! -z $timeget_in ] && timeget=$timeget_in

ymdh=${timeget:0:10}

[ -f $GFSDIR/$timeget/sfc_wind/$testfile ] && exit


### TEST
res=''
res=`$sshcommand $hostname "ls $GFSDIR_OFP/$ymdh/plot/$testfile"`
###
if [ ! -z $res ] ;then

mkdir $GFSDIR/$timeget
items="300_wspd  500_vort  700_vvel  850_temp  850_theq  sfc_2mtemp  sfc_prcp  sfc_temp  sfc_wind"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $GFSDIR_OFP/$ymdh/plot/${item}_f432000.png"`
 if [ ! -z $res ] ;then 
 mkdir $GFSDIR/$timeget/$item
 rsync -e "$sshcommand" $hostname:$GFSDIR_OFP/$ymdh/plot/$item*.png $GFSDIR/$timeget/$item/  
 chmod 644 $GFSDIR/$timeget/$item/*.png
 fi
done
fi
cd $mydir/..
php ./monitor_plot_d1.php

