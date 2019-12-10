#!/bin/bash

. config.main || exit $?
PARENT_REF_TIME=$1
STIME=$2
logfile=$OUTPUT/d3/ref_${PARENT_REF_TIME}/$STIME/log/fcst_scale/mdet_LOG.pe000000


stat=`cat fcst_ofp.stat.$PARENT_REF_TIME.$STIME`
#stat=`cat fcst_obcx.stat.$PARENT_REF_TIME.$STIME`

if [ `echo $stat| awk '{print $1}'` == "submit" ] ;then
 jobid=`echo $stat | awk '{print $2}'`
 wait=`/usr/local/bin/pjstat | grep $jobid | cut -c 64-76`
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
