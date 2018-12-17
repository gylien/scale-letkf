#!/bin/sh
set -e


#
tstart='2015-07-29 8:00:00'
tend='2015-07-29 9:00:00'

ctint=3600 # obssim interval 
tint=600 # analysis interval (Do not modify!)

ctime="$tstart"
while (($(date -ud "$ctime" '+%s') <= $(date -ud "$tend" '+%s'))); do # -- time
  ctime=$(date -ud "${ctint} second $ctime" '+%Y-%m-%d %H%M')
	echo $ctime
done # -- time

exit



