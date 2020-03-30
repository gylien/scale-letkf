#!/bin/bash -l

myname=$0
cd `dirname $myname`
wkdir=`pwd`

. $wkdir/config.main
. $wkdir/src/func_datetime.sh


TMP=`pwd`/bgdata 


PARENT=$OUTPUT/$EXP3

INITRUNDIR=$SCALEDIR/scale-letkf_verify/scale/run_d4_init
statfile=fcst_ofp.stat

mode=$1

YMDH=`date -u +%Y%m%d%H0000`

if [ "$mode" == "init" ];then

testmem=`printf %04d $MEMBER`
testfile=init_$(datetime_scale $YMDH).pe000000.nc
#echo $testmem/$testfile
res=`ls $PARENT/ref_*/*/${EXP4}/anal/$testmem/$testfile 2>/dev/null`
#res=`ls -d $PARENT/ref_* 2>/dev/null`

tslatest=`date -ud "-1 day" +%s`
if [ ! -z "$res" ]; then
 for path in $res;do 
  YMDHinitD2=`echo $path | grep -o 'ref_[0-9]\{14\}' | cut -c 5-18`
  YMDHinitD2D3=`echo $path | grep -o 'ref_[0-9]\{14\}\/[0-9]\{14\}' | cut -c 5-33`
#  YMDHinitD3=`echo $path | grep -o '[0-9]\{14\}' | cut -c 15-28`
#  YMDHinitD3=`echo $path | grep -o '\/[0-9]\{14\}' `
# echo $path
#  echo $YMDHinitD2
 # echo $YMDHinitD2D3
#exit
  if [ ! -f ${INITRUNDIR}/${statfile}.${YMDHinitD2}.${YMDH} ] ;then ### not running now
   tsfile=`ls -l ${PARENT}/ref_${YMDHinitD2D3}/${EXP4}/anal/$testmem/$testfile --time-style="+%s" | awk '{print $6}'`
   if [ $tsfile -gt $tslatest ] ;then
     tslatest=$tsfile
     YMDHinitD2D3latest=$YMDHinitD2D3
   fi 
  fi
 done
fi

mkdir -p $TMP
[ -h $TMP/$YMDH ] && rm $TMP/$YMDH
path=$PARENT/ref_${YMDHinitD2D3latest}/$EXP4
ln -s $path $TMP/$YMDH

fi ### init

if [ "$mode" == "init" ] || [ "$mode" == "next" ];then
 YMDHplus=`date -ud "1 hour" +%Y%m%d%H0000`
 tslatest=`date -ud "-1 day" +%s`
elif [ ${#mode} == 14 ] ;then
 YMDHplus=$mode
 tslatest=`date -ud "-10 year" +%s`
else
   echo " usage: ./prep.sh init/next/YYYYMMDDHHMMSS "
   exit
fi

testmem=`printf %04d $MEMBER`
testfile=init_$(datetime_scale $YMDHplus).pe000000.nc
#echo $testmem/$testfile
res=`ls $PARENT/ref_*/*/${EXP4}/anal/$testmem/$testfile 2>/dev/null`
#res=`ls -d $PARENT/ref_* 2>/dev/null`
#echo $testmem/$testfile
#echo $res


if [ ! -z "$res" ]; then
 for path in $res;do 
  YMDHinitD2=`echo $path | grep -o 'ref_[0-9]\{14\}' | cut -c 5-18`
  YMDHinitD2D3=`echo $path | grep -o 'ref_[0-9]\{14\}\/[0-9]\{14\}' | cut -c 5-33`
# echo $path
# echo $YMDHinitD2
# echo $YMDHinitD2D3
# echo ${INITRUNDIR}/${statfile}.${YMDHinitD2}.${YMDHplus}

  if [ ! -f ${INITRUNDIR}/${statfile}.${YMDHinitD2}.${YMDHplus} ] ;then ### not running now
   tsfile=`ls -l ${PARENT}/ref_${YMDHinitD2D3}/${EXP4}/anal/$testmem/$testfile --time-style="+%s" | awk '{print $6}'`
   if [ $tsfile -gt $tslatest ] ;then
     tslatest=$tsfile
     YMDHinitD2D3latest=$YMDHinitD2D3
   fi 
  fi
 done
fi

mkdir -p $TMP
[ -h $TMP/$YMDHplus ] && rm $TMP/$YMDHplus
path=$PARENT/ref_${YMDHinitD2D3latest}/$EXP4
ln -s $path $TMP/$YMDHplus

if [ "$mode" == "init" ] || [ "$mode" == "next" ];then
links=`ls -x $TMP`
for link in $links;do
 [ $link != $YMDH ] && [ $link != $YMDHplus ] && rm $TMP/$link
done
fi

