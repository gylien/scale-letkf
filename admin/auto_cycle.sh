#!/bin/bash -l

running='running_cycle'

if [ -f $running ] ;then
  echo 'already running.'
  exit 1
else
  echo $HOSTNAME $$ > $running
fi

function unlock () {
 [ `cat $running | awk '{print $2}'` == $$ ] && rm -f $running
}
trap unlock EXIT


inittime_end='2021-01-01 00:00:00'
sec_sleep=600

source ~/.bashrc

while [ "$inittime_bef" != "$inittime_end" ] ;do
  inittime_now=`cat admin_cycle.time`  
  echo " $inittime_now start "
  nohup ./admin_cycle.sh &> log
  if [ `tail -1 admin_cycle.log | cut -d " " -f 3` == "[DONE]" ] ; then
   echo " $inittime_now complete "
   inittime_bef=$inittime_now
  elif [ `tail -1 admin_cycle.log | cut -d " " -f 3` == "[WAIT]" ] ; then ### wait
   sleep $sec_sleep
  else
   echo " cycle finishes abnormally. something is wrong. "
   echo " auto.sh stopped."
   exit 
  fi
done

