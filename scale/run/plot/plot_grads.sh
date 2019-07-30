#!/bin/sh

export GASCRP=/work/hp150019/share/honda/local/grads/lib

myname=$0
mydir=`dirname $myname`
cd $mydir

BASETIME=$1
STIME=$2
istart=$3
iend=$4
iskip=$5

OUTDIR=/work/hp150019/c24140/scale-letkf-rt/r0051_nest/exp_d4/plot/$STIME

YMD=`echo $STIME | cut -c 1-8`
HMS=`echo $STIME | cut -c 9-14`

rm out/*
cp fcst_template.ctl fcst.ctl
sed -i -e "s/<BASETIME>/${BASETIME}/g" fcst.ctl
sed -i -e "s/<YYYYMMDD>/${YMD}/g" fcst.ctl
sed -i -e "s/<HHMMSS>/${HMS}/g" fcst.ctl


grads -bcl "1500m_ref.gs fcst.ctl $istart $iend $istip"
grads -bcl "5000m_ref.gs fcst.ctl $istart $iend $istip"

[ -d $OUTDIR ] || mkdir -p $OUTDIR
rm $OUTDIR/*
cp -r out/* $OUTDIR/
