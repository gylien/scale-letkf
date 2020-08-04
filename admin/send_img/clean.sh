#!/bin/bash  -l

me=`whoami`
list=`ps aux | grep $me | grep auto_send_img.sh | grep -v "grep" | awk '{print $2}'`

for pid in $list ;do
  kill $pid
done

