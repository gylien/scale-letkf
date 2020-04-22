#!/bin/sh


time_end_fmt=`cat ../../../admin_cycle.time` 
time_end=`date -d "$time_end_fmt" +%Y%m%d%H%M%S`
time_start=`date -d "-6 day $time_end_fmt" +%Y%m%d%H%M%S`

wkdir="$( cd "$( dirname "$0" )" && pwd )"

cd ${wkdir}

BASEDIR=$wkdir/../../../../result/ope/d1

vars_u=( X X U V T Q PS ) 
vars_l=( x x u v t q ps ) 

for varnum in `seq 2 6`;do
 cd ./draw
 dclfrt draw_${vars_l[$varnum]}_num.f90 -o draw_${vars_l[$varnum]}
 cd -
done

nsubdom=192
nmem=50

for idom in `seq $nsubdom` ;do
 idomm=`expr $idom - 1`
 cdom=`printf %03d $idomm`

 echo $cdom ' ...'


 mkdir -p data/$cdom
 mkdir -p  data_inv
 mkdir -p figs/$cdom

 RECDIR=data/$cdom

rm $RECDIR/*

time=$time_start
while [ $time -le $time_end ];do
echo $time $time_end

if [ -f $BASEDIR/$time/log/letkf/NOUT-000${nsubdom} ] ;then 
imems=$nsubdom
imeme=`expr \( $nmem + 2 \) \* $nsubdom - 1`
  for imem in `seq $imems $imeme` ; do
    cmem=`printf %06d $imem`
    [ -f $BASEDIR/$time/log/letkf/NOUT-${cmem} ] && rm $BASEDIR/$time/log/letkf/NOUT-${cmem}
  done
fi

[ -f $BASEDIR/$time/log/letkf/NOUT-000001 ] && cdomt=$cdom || cdomt='000'

key_begin='OBSERVATIONAL\ DEPARTURE\ STATISTICS\ \[ANALYSIS\] \(IN\ THIS\ SUBDOMAIN\)'
key_end='NUMBER'
cat $BASEDIR/$time/log/letkf/NOUT-000${cdomt} | awk "/${key_begin}/,/${key_end}/" > anal.txt 

key_begin='OBSERVATIONAL\ DEPARTURE\ STATISTICS\ \[GUESS\] \(IN\ THIS\ SUBDOMAIN\)'
key_end='NUMBER'
cat $BASEDIR/$time/log/letkf/NOUT-000${cdomt} | awk "/${key_begin}/,/${key_end}/" > gues.txt


for step in anal gues ;do
for varnum in `seq 2 6` ;do
 bias=`cat ${step}.txt | grep BIAS   | awk -v "vnum=${varnum}" '{print $vnum}'`
 rmse=`cat ${step}.txt | grep RMSE   | awk -v "vnum=${varnum}" '{print $vnum}'`
 num=`cat ${step}.txt | grep NUMBER | awk -v "vnum=${varnum}" '{print $vnum}'`
 recfile=${vars_l[$varnum]}_${step}.txt
 echo "$time $bias $rmse $num" | awk '{printf "%14d" ,$1} {printf "%13.6f" ,$2} {printf "%13.6f" ,$3} {printf "%6d\n" ,$4}'>> $RECDIR/$recfile
done
done

yyyymmdd=`echo $time | cut -c 1-8`
hh=`echo $time | cut -c 9-10`
time=`date -d "6 hours $yyyymmdd $hh" +%Y%m%d%H0000`
done

for varnum in `seq 2 6` ;do
for step in anal gues ;do
 recfile=${vars_l[$varnum]}_${step}.txt 
 tac ./data/$cdom/$recfile  > ./data_inv/$recfile
done
 cd ./draw
 echo "draw_${vars_l[$varnum]}..."
 ./draw_${vars_l[$varnum]} 
 cd $wkdir
done
  
  mogrify -trim ./temp_figs/*.png
 
  cp ./temp_figs/*.png ./figs/$cdom/


if [ `date -d "$time_end_fmt" +%H` == "00" ];then
  mkdir -p ./save/$cdom/$time_end
  cp ./temp_figs/*.png ./save/$cdom/$time_end/
fi

done ### idom


### draw map
cd ./draw
for varnum in `seq 2 6` ;do
 ./draw_map_subdom ${vars_l[$varnum]} rmse anal
 ./draw_map_subdom ${vars_l[$varnum]} bias anal
 ./draw_map_subdom ${vars_l[$varnum]} nobs anal
done
cd -

  mogrify -trim ./figs/map/*.png

if [ `date -d "$time_end_fmt" +%H` == "00" ];then
  mkdir -p ./save/map/$time_end
  cp ./figs/map/*.png ./save/map/$time_end/
fi

rm anal.txt
rm gues.txt

