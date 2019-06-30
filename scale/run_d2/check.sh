#!/bin/sh

for n in `seq 10`;do
# echo `pjstat -H | grep cycle_job | awk '{print $1}' | tail -n $n | head -n 1`
 jobid=`pjstat -H | grep 208 | awk '{print $1}' | tail -n $n | head -n 1`
 echo `pjstat -H -s $jobid | grep MAX\ MEMORY` 
done
