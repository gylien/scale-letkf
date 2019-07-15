#!/bin/sh

r_url="c24140@ofp.jcahpc.jp"
web_url="amemiya@daweb.r-ccs27.riken.jp"

nmem=36

R_SRCDIRBASE="/work/hp150019/c24140/scale-letkf-rt/r0051_nest/exp_d3"
SRCDIRBASE="/data_ballantine01/miyoshi-t/amemiya/scale-letkf-rt/r0051_nest_d3"

WEBDIRBASE="/home/amemiya/public_html/scale/data/ens/fcst_d3"

PARENT_REF_TIME=$1
TIME=$2

if [ ${#TIME} -eq 0 ] ;then
 echo "specify ref_time and time"
 exit 
fi

ulimit -s unlimited 

#mv ./figure/*.png ./figure/d3/

#mogrify -trim ./figure/d3/*.png


ssh $web_url "mkdir -p ${WEBDIRBASE}/fcst_d3/${PARENT_REF_TIME}/${TIME}"

for prod in rain uvw1500m uvw5000m dbz1500m ;do
for imem in `seq $nmem`;do
mem=`printf %04d $imem`
rsync -av ./figure/d3/${prod}_${mem}_f*_0001.png ${web_url}:${WEBDIRBASE}/fcst_d3/${PARENT_REF_TIME}/${TIME}
done
rsync -av ./figure/d3/${prod}_mean_f*_0001.png ${web_url}:${WEBDIRBASE}/fcst_d3/${PARENT_REF_TIME}/${TIME}
done
