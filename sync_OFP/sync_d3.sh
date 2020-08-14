#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc


timeget_in_ref=$1
timeget_in=$2

#FCSTDIR_OFP=/work/hp150019/c24140/HPCC_SCALE-LETKF-rt2/result/ope/d3
FCSTDIR_OFP=$SCALEDIR_OFP/result/ope/d3
FCSTDIR=$DATADIR/d3

testfile="sfc_prcp_f000000.png"
testdir='fcstgpi'
res=`$sshcommand $hostname "find $FCSTDIR_OFP -type d -maxdepth 3 -name $testdir -mmin -10"`

if [ -z "$res" ] && [ -z "$timeget_in" ] ;then
   exit
fi

timegets=`echo $res | grep -o 'ref_[0-9]\{14\}\/[0-9]\{14\}'`

if [ ! -z $timeget_in ] ;then
   timegets=ref_$timeget_in_ref/$timeget_in
fi

for timeget in $timegets; do
echo 'get ' $timeget ' ...'

if [ -f $FCSTDIR/$timeget/$testdir/mdet/$testfile ] ;then
   exit
fi

for cmem in mdet mean; do

### TEST
res=''
res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/$testdir/$cmem/$testfile"`
###
if [ ! -z $res ] ;then

mkdir -p $FCSTDIR/$timeget
items=" sfc_prcp  sfc_uvt 1000m_uvw 1000m_dbz 3000m_uvw 3000m_dbz 5000m_uvw 5000m_dbz "

for item in $items ;do
 res=''
 res=`$sshcommand $hostname "ls $FCSTDIR_OFP/$timeget/$testdir/$cmem/${item}_f000000.png"`
 if [ ! -z $res ] ;then 
 mkdir -p $FCSTDIR/$timeget/$cmem/$item
 rsync -e "$sshcommand" $hostname:$FCSTDIR_OFP/$timeget/$testdir/$cmem/$item*.png $FCSTDIR/$timeget/$cmem/$item  
 chmod 644 $FCSTDIR/$timeget/$cmem/$item/*.png
 fi
done
fi

done #cmem
done #timeget

cd $mydir/..
php ./monitor_plot_d3.php

