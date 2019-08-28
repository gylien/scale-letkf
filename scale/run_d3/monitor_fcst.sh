#!/bin/sh

TMPDIR="../tmp/scale-letkf_d3"
PARENT_REF_TIME=$1
STIME=$2
#logfile=${TMPDIR}_ref_${PARENT_REF_TIME}_$STIME/out/$STIME/log/fcst_scale/mdet_LOG.pe000000
logfile=/work/hp150019/c24140/scale-letkf-rt/r0051_nest/exp_d3/ref_${PARENT_REF_TIME}/$STIME/log/fcst_scale/mdet_LOG.pe000000


stat=`cat fcst_ofp.stat.$PARENT_REF_TIME.$STIME`

if [ `echo $stat| awk '{print $1}'` == "submit" ] ;then
 jobid=`echo $stat | awk '{print $2}'`
 wait=`pjstat | grep $jobid | cut -c 64-76`
 if [ `echo $wait | cut -c 1` == '(' ] ;then
  echo $wait
 else
  if [ -s $logfile ] ;then
  step=`tail -n 100 $logfile | grep STEP: | tail -n 1 | awk '{print $8}'`  
  step=${step:0:-1}
  steptot=`tail -n 100 $logfile | grep STEP: | tail -n 1 | awk '{print $9}'`  
  echo `expr $step \* 100 \/ $steptot`'%'
  else
   echo 'init'
  fi
 fi 
else
 echo $stat
fi
