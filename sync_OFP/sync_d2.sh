#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc

timeget_in=$1

FCSTDIR_OFP=$SCALEDIR_OFP/result/ope/d2
FCSTDIR=$DATADIR/d2

testfile="sfc_wind_f086400.png"
testdir='fcstgpi'
res=`$sshcommand $hostname "find $FCSTDIR_OFP -type d -maxdepth 2 -name $testdir -mmin -10"`


[ -z "$res" ] && [ -z "$timeget_in" ] && exit

timegets=`echo $res | grep -o '[0-9]\{14\}'`

[ ! -z $timeget_in ] && timegets=$timeget_in

for timeget in $timegets; do
echo 'get ' $timeget ' ...'

ymdh=${timeget:0:10}

[ -f $FCSTDIR/$timeget/sfc_wind/$testfile ] && exit


### TEST
res=''
res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/fcstgpi/mdet/$testfile"`
###
if [ ! -z $res ] ;then

mkdir $FCSTDIR/$timeget
items=" sfc_2mtemp  sfc_prcp  sfc_wind  925_theq 925_temp 850_theq 850_temp 700_temp 500_temp"
#items=" sfc_2mtemp  sfc_prcp  sfc_wind"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/fcstgpi/mdet/${item}_f086400.png"`
 if [ ! -z $res ] ;then 
 mkdir $FCSTDIR/$timeget/$item
 rsync -e "$sshcommand" $hostname:$FCSTDIR_OFP/$timeget/fcstgpi/mdet/$item*.png $FCSTDIR/$timeget/$item/  
 chmod 644 $FCSTDIR/$timeget/$item/*.png
 fi
done


items="d3_sfc_uvt  d3_sfc_prcp"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/fcstgpi/mdet/${item}_f086400.png"`
 if [ ! -z $res ] ;then 
 mkdir -p ${FCSTDIR}_d3/$timeget/$item
 rsync -e "$sshcommand" $hostname:$FCSTDIR_OFP/$timeget/fcstgpi/mdet/$item*.png ${FCSTDIR}_d3/$timeget/$item/  
 chmod 644 ${FCSTDIR}_d3/$timeget/$item/*.png
 files=`ls -x ${FCSTDIR}_d3/$timeget/$item/d3_*.png`
 for file in $files ;do
   fname=`basename $file`
   file_new=`echo $fname | cut -c 4-`
   mv $file ${FCSTDIR}_d3/$timeget/$item/$file_new 
 done
 item_new=`echo $item | cut -c 4-`
 mv ${FCSTDIR}_d3/$timeget/$item ${FCSTDIR}_d3/$timeget/$item_new
 fi
done


fi
done

cd $mydir/..
php ./monitor_plot_d2.php
php ./monitor_plot_d3.php


