#!/bin/bash -l

source ~/.bashrc

MEMBER=<MEMBER>
OUTDIR=<OUTDIR>

myname=$0

PARENT_REF_TIME=$1
TIME=$2
FCSTLEN=$3
inum=$4

ntime=`expr $FCSTLEN \/ 600 + 1`


if [ ${#inum} -eq 0 ] ;then
 echo "specify ref_time and time and fcstlen and cmem"
 exit 
fi

mydir=`dirname $myname`

inum_mean=`expr $MEMBER + 1`
inum_mdet=`expr $MEMBER + 2`
if [ $inum -eq $inum_mean ];then
cmem=mean
elif [ $inum -eq $inum_mdet ];then
cmem=mdet
else
cmem=`printf %04d $inum`
fi

srcfile=${OUTDIR}/${TIME}/fcst_sno_np00001/mean/history.pe000000.nc
if [ -f $srcfile ] ;then

ulimit -s unlimited 

module load hdf5
module load netcdf
module load netcdf-fortran

echo "draw "$cmem

mkdir -p $mydir/plot_temp
rm $mydir/plot_temp/* ###links

echo 'plot rain...'
./draw_rain ${PARENT_REF_TIME} ${TIME} $cmem 1 $ntime 1 
echo 'plot uvw1500...'
./draw_uvw ${PARENT_REF_TIME} ${TIME} $cmem 1 $ntime 1 1500 
echo 'plot uvw5000...'
./draw_uvw ${PARENT_REF_TIME} ${TIME} $cmem 1 $ntime 1 5000 
echo 'plot dbz1500...'
./draw_dbz ${PARENT_REF_TIME} ${TIME} $cmem 1 $ntime 1 1500 
echo 'plot dbz5000...'
./draw_dbz ${PARENT_REF_TIME} ${TIME} $cmem 1 $ntime 1 5000 
echo 'trim...'
mogrify -trim $mydir/plot_temp/$cmem/*.png
echo '=== complete ==='
