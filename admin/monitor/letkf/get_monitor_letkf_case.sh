#!/bin/sh


time_start=20190901000000 
time_end=20190930000000

BASEDIR=/work/hp150019/c24140/scale-letkf-rt/r0051/exp_d1
RECDIR=./data_case

vars_u=( X X U V T Q PS ) 
vars_l=( x x u v t q ps ) 

[ -d data ] || mkdir data
[ -d data_inv ] || mkdir data_inv
[ -d figs ] || mkdir figs

rm $RECDIR/*

time=$time_start
while [ $time -le $time_end ];do
echo $time

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

mkdir $RECDIR/$yyyymmdd
mv $RECDIR/*.txt $RECDIR/$yyyymmdd/

hh=`echo $time | cut -c 9-10`
time=`date -d "24 hours $yyyymmdd $hh" +%Y%m%d%H0000`
done


