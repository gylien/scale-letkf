#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc

nftime_msm=33
tint=3600

MSMDIR_OFP=$SCALEDIR_OFP/external/msm
MSMDIR=$DATADIR/msm

testfile="sfc_wind_f086400.png"

timeget_OFP=`$sshcommand $hostname "cat $MSMDIR_OFP/mtime"`
timeget=`date -d "$timeget_OFP" +%Y%m%d%H%M%S`

timeget_in=$1
[ ! -z $timeget_in ] && timeget=$timeget_in

ymdh=${timeget:0:10}

[ -f $MSMDIR/$timeget/sfc_wind/$testfile ] && exit

### TEST
res=''
res=`$sshcommand $hostname "ls $MSMDIR_OFP/$ymdh/plot/$testfile"`
###

if [ ! -z $res ] ;then
 mkdir $MSMDIR/$timeget

items="sfc_2mtemp  sfc_prcp  sfc_wind"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $MSMDIR_OFP/$ymdh/plot/${item}_f086400.png"`
 if [ ! -z $res ] ;then 
 mkdir $MSMDIR/$timeget/$item
 rsync -e "$sshcommand" $hostname:$MSMDIR_OFP/$ymdh/plot/$item*.png $MSMDIR/$timeget/$item/  
 chmod 644 $MSMDIR/$timeget/$item/*.png
 fi
done

items="d3_sfc_uvt  d3_sfc_prcp"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $MSMDIR_OFP/$ymdh/plot/${item}_f086400.png"`
 if [ ! -z $res ] ;then 
 mkdir -p ${MSMDIR}_d3/$timeget/$item
 rsync -e "$sshcommand" $hostname:$MSMDIR_OFP/$ymdh/plot/$item*.png ${MSMDIR}_d3/$timeget/$item/  
 chmod 644 ${MSMDIR}_d3/$timeget/$item/*.png
 files=`ls -x ${MSMDIR}_d3/$timeget/$item/d3_*.png`
 for file in $files ;do
   fname=`basename $file`
   file_new=`echo $fname | cut -c 4-`
   mv $file ${MSMDIR}_d3/$timeget/$item/$file_new 
 done
 item_new=`echo $item | cut -c 4-`
 mv ${MSMDIR}_d3/$timeget/$item ${MSMDIR}_d3/$timeget/$item_new
 fi
done


items="925_theq 925_temp 850_theq 850_temp 700_temp 500_temp"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $MSMDIR_OFP/$ymdh/plot/${item}_f086400.png"`
 if [ ! -z $res ] ;then 
 mkdir $MSMDIR/$timeget/$item
 rsync -e "$sshcommand" $hostname:$MSMDIR_OFP/$ymdh/plot/$item*.png $MSMDIR/$timeget/$item/  
 chmod 644 $MSMDIR/$timeget/$item/*.png
 for iftime in `seq $nftime_msm`;do
  ftime=`expr $iftime \* $tint`
  ftimep=`expr \( $iftime / 3 \) \* 3 \* $tint`
  ftimef=`printf %06d $ftime`
  ftimepf=`printf %06d $ftimep`
  [ $ftimef != $ftimepf ] && ln -s $MSMDIR/$timeget/$item/${item}_f${ftimepf}.png  $MSMDIR/$timeget/$item/${item}_f${ftimef}.png 
 done
 fi
done

fi
cd $mydir/..
php ./monitor_plot_d2.php
<<<<<<< HEAD
php ./monitor_plot_d3.php
=======
>>>>>>> 504b3ad... initial commit

