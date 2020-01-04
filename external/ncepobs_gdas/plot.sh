#!/bin/bash -l

source ~/.bashrc

wkdir="$( cd "$( dirname "$0" )" && pwd )"

START_TIME=$1
GET_TIME=$START_TIME

END_TIME=`cat mtime`
while  [ "$GET_TIME" != "$END_TIME" ] ;do
YYYY=`date -u -d "$GET_TIME" +'%Y'`
MM=`date -u -d "$GET_TIME" +'%m'`
DD=`date -u -d "$GET_TIME" +'%d'`
HH=`date -u -d "$GET_TIME" +'%H'`
YYYYMMDDHH="$YYYY$MM$DD$HH"

cd $wkdir/run/plot
rm *.png
./map $YYYYMMDDHH 
mkdir $wkdir/../ncepobs_gdas_letkf/${YYYYMMDDHH}/map
[ ! -z "`ls *.png`" ] && mogrify -trim *.png
mv *.png $wkdir/../ncepobs_gdas_letkf/${YYYYMMDDHH}/map

GET_TIME=`date -u -d "+ 6 hour ${GET_TIME}" +'%Y-%m-%d %H'`


done
