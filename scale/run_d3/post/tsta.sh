#!/bin/sh


STIME=20200412180000 

STIMEf="${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}"
export LANG=en_us
tstamp=`date -d "$STIMEf" +%H:%MZ%d%b%Y | tr '[a-z]' '[A-Z]'`

echo $tstamp
