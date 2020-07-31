#!/bin/sh


time_end_fmt=`cat ../../admin_cycle.time` 
time_end=`date -d "$time_end_fmt" +%Y%m%d%H%M%S`
time_start=`date -d "-6 day $time_end_fmt" +%Y%m%d%H%M%S`

wkdir="$( cd "$( dirname "$0" )" && pwd )"

cd ${wkdir}

BASEDIR=$wkdir/../../../result/ope/d1
RECDIR=./data

vars_u=( X X U V T Q PS ) 
vars_l=( x x u v t q ps ) 

[ -d data ] || mkdir data
[ -d data_inv ] || mkdir data_inv
[ -d figs ] || mkdir figs

rm $RECDIR/*

time=$time_start
while [ $time -le $time_end ];do
#echo $time $time_end

key_begin='OBSERVATIONAL\ DEPARTURE\ STATISTICS\ \[ANALYSIS\] \(GLOBAL\)'
key_end='NUMBER'
cat $BASEDIR/$time/log/letkf/NOUT-000000 | awk "/${key_begin}/,/${key_end}/" > anal.txt 

key_begin='OBSERVATIONAL\ DEPARTURE\ STATISTICS\ \[GUESS\] \(GLOBAL\)'
key_end='NUMBER'
cat $BASEDIR/$time/log/letkf/NOUT-000000 | awk "/${key_begin}/,/${key_end}/" > gues.txt

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
 tac ./data/$recfile  > ./data_inv/$recfile
done
 cd ./draw
 dclfrt draw_${vars_l[$varnum]}_num.f90
 ./a.out
 cd -
done

mogrify -trim ./figs/*.png 

if [ `date -d "$time_end_fmt" +%H` == "00" ];then
  mkdir -p ./save/$time_end
  cp ./figs/*.png ./save/$time_end/
fi

### TORI AEZU
exit 0 


  cd $wkdir/subdomain
  ./get_monitor_letkf.sh
  cd $wkdir/obsdep/sounding 
  ./drawall.sh $time_end



