#!/bin/sh

source ~/.bashrc

r_url="c24140@ofp.jcahpc.jp"
web_url="amemiya@daweb.r-ccs27.riken.jp"

WEBDIRBASE="/home/amemiya/public_html/scale/data/ens/fcst_d3"

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
if [ $inum -eq 51 ];then
cmem=mean
elif [ $inum -eq 52 ];then
cmem=mdet
else
cmem=`printf %04d $inum`
fi

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
#ssh $web_url "mkdir -p ${WEBDIRBASE}/ref_${PARENT_REF_TIME}/${TIME}"

#for prod in rain uvw1500m uvw5000m dbz1500m ;do
#for imem in `seq $nmem`;do
#mem=`printf %04d $imem`
#rsync -av ./figure/d3/${prod}_${mem}_f*_0001.png ${web_url}:${WEBDIRBASE}/ref_${PARENT_REF_TIME}/${TIME}
#done
#rsync -av ./figure/d3/${prod}_mean_f*_0001.png ${web_url}:${WEBDIRBASE}/ref_${PARENT_REF_TIME}/${TIME}
#done
