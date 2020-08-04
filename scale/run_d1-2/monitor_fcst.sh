#!/bin/bash

. config.main || exit $?

STIME=$1

iter=2 ##### run with 182 nodes 

logfile=$OUTPUT/d2/$STIME/log/0001_LOG_${STIME}.pe000000
logfile2=$OUTPUT/d2/$STIME/log/0050_LOG_${STIME}.pe000000

stat=`cat fcst_ofp.stat.$STIME`
#stat=`cat fcst_obcx.stat.$STIME`
if [ `echo $stat| awk '{print $1}'` == "submit" ] ;then
 jobid=`echo $stat | awk '{print $2}'`
 wait=`/usr/local/bin/pjstat | grep $jobid | cut -c 64-76`
 if [ `echo $wait | cut -c 1` == '(' ] ;then
  echo $wait
 else
  if [ -s $logfile ] ;then
  step=`tail -n 200 $logfile | grep STEP: | tail -n 1 | awk '{print $8}'`  
  step=${step:0:-1}
  steptot=`tail -n 200 $logfile | grep STEP: | tail -n 1 | awk '{print $9}'`  
  if [ $iter == 2 ] && [ "$step" == "$steptot" ] ;then
    step2=`tail -n 200 $logfile2 | grep STEP: | tail -n 1 | awk '{print $8}'`  
    step2=${step2:0:-1}
    step=`expr $step + $step2`
    steptot=`expr $steptot \* 2`
  fi  
  echo `expr $step \* 100 \/ $steptot`'%'
  else
   echo 'init'
  fi
 fi 
else
 echo $stat
fi
