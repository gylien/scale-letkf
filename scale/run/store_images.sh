#!/bin/bash

TMP=$1
OUTDIR=$2

  echo
  echo "Start: store images"
    
  cd $TMP

  line=`cat ./cycle_job.sh.e* | grep "Start cycle.sh"`
  TIMEstart=`echo $line | awk '{print $5}'`
  TIMEstop=`echo $line | awk '{print $6}'`

  TIMEf="${TIMEstart:0:4}-${TIMEstart:4:2}-${TIMEstart:6:2} ${TIMEstart:8:2}"
  TIME=`date -ud "$TIMEf" +%Y%m%d%H%M%S`

  while [ $TIME -le $TIMEstop ] ;do
  echo $TIMEf" ..."
  DEST=$OUTDIR/$TIME/dafcst

  [ -d $DEST ] || mkdir -p $DEST
 
#  mkdir -p save

  for z in `seq 0 10`; do
    ZLEV="z"$(printf %02d $z)
    echo ${ZLEV}000...


    timestamp=`date -ud "$TIMEf" +%Y%m%d-%H`
    timestamp2=`date -ud "1 hour $TIMEf" +%Y%m%d-%H%M%S`

    find . -name "anal_dbz_${timestamp}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_anal.txt
    find . -name "obs_dbz_${timestamp}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_obs.txt
    find . -name "fcst_dbz_${timestamp}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_fcst.txt
    find . -name "anal_dbz_${timestamp2}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print >> list_anal.txt
    find . -name "obs_dbz_${timestamp2}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print >> list_obs.txt
    find . -name "fcst_dbz_${timestamp2}*${ZLEV}[0-9][0-9][0-9]m_*.png" -print >> list_fcst.txt
    tar -zcf $DEST/anal_${ZLEV}000m.tar.gz --files-from list_anal.txt    
    tar -zcf $DEST/obs_${ZLEV}000m.tar.gz --files-from list_obs.txt    
    tar -zcf $DEST/fcst_${ZLEV}000m.tar.gz --files-from list_fcst.txt    

    xargs rm < list_anal.txt
    xargs rm < list_obs.txt
    xargs rm < list_fcst.txt

#if [ $z -eq 1 ]; then
#    mv list_[a,f,o]*.txt ./save/
#fi
#####    find [a,o]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs rm -f &> /dev/null

  done ### z
  rm list_*.txt

 
  TIMEf=`date -ud "1 hour $TIMEf" +"%Y-%m-%d %H"`
  TIME=`date -ud "$TIMEf" +%Y%m%d%H%M%S`

  done ### time

  cd - > /dev/null
  

  echo "End: store images"


