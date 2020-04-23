#!/bin/bash -l

myname=$0
cd `dirname $myname`
wkdir=`pwd`


. $wkdir/config.main
. $wkdir/src/func_datetime.sh

TIMEIN=$1
TIMEOUT=$2

SRCDIR=$OUTDIR/$TIMEIN


for imem in mean mdet `seq 1 $MEMBER`;do
if [ "$imem" == 'mdet' ] || [ "$imem" == "mean" ];then
 cmem=$imem
else
 cmem=`printf %04d $imem`
fi
echo $cmem"..."
mkdir -p $OUTDIR/$TIMEOUT/anal/$cmem
mv $OUTDIR/$TIMEIN/anal/$cmem/init_${TIMEOUT:0:8}-${TIMEOUT:8:6}.000.pe*.nc $OUTDIR/$TIMEOUT/anal/$cmem
#echo $OUTDIR/$TIMEIN/anal/$cmem/init_${TIMEOUT:0:8}-${TIMEOUT:8:6}.000.pe000000.nc $OUTDIR/$TIMEOUT/anal/$cmem
#exit 0
done
