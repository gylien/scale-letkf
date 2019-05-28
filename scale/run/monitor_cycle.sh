#!/bin/sh

TMPDIR="../tmp/scale-letkf_exp_d1"
STIME=$1
logfile=$TMPDIR/out/$STIME/log/scale/0001_LOG.pe000000

stat=`cat cycle_ofp.stat.$STIME`

if [ "$stat" == "submit" ] ;then
 pjstat > pjstat.txt
 wait=`cat pjstat.txt | grep cycle_job | cut -c 64-76` 
 rm pjstat.txt
 if [ `echo $wait | cut -c 1` == '(' ] ;then
  echo $wait
 else
  if [ -s $logfile ] ;then
  step=`tail -n 100 $logfile | grep STEP: | tail -n 1 | awk '{print $8}'`  
  step=${step:-0}
  steptot=`tail -n 100 $logfile | grep STEP: | tail -n 1 | awk '{print $9}'`  
  steptot=${steptot:-100}

  step=`echo $step | cut -d "/" -f 1`

  echo `expr $step \* 100 \/ $steptot`'%'
  else
   echo 'init'
  fi
 fi 
else
 echo $stat
fi
