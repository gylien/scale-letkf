#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc

FCSTDIR_OFP=$SCALEDIR_OFP/result/ope/d1_ext
FCSTDIR=$DATADIR/d1

testfile="sfc_wind_f432000.png"
testdir='fcstgpi'
res=`$sshcommand $hostname "find $FCSTDIR_OFP -type d -maxdepth 2 -name $testdir -mmin -10"`

[ -z "$res" ] && [ ! -z "$timeget_in" ] && exit

timegets=`echo $res | grep -o '[0-9]\{14\}'`

timeget_in=$1
[ ! -z "$timeget_in" ] && timegets=$timeget_in

for timeget in $timegets; do
echo 'get ' $timeget ' ...'
ymdh=${timeget:0:10}

[ -f $FCSTDIR/$timeget/sfc_wind/$testfile ] && continue


### TEST
res=''
res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/fcstgpi/mdet/$testfile"`
###
if [ ! -z $res ] ;then

mkdir $FCSTDIR/$timeget
items="300_wspd  500_vort  700_vvel  850_temp  850_theq  sfc_2mtemp  sfc_prcp  sfc_temp  sfc_wind"

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/fcstgpi/mdet/${item}_f432000.png"`
 if [ ! -z $res ] ;then 
 mkdir $FCSTDIR/$timeget/$item
 rsync -e "$sshcommand" $hostname:$FCSTDIR_OFP/$timeget/fcstgpi/mdet/$item*.png $FCSTDIR/$timeget/$item/  
 chmod 644 $FCSTDIR/$timeget/$item/*.png
 fi
done
fi

done

cd $mydir/..
php ./monitor_plot_d1.php


