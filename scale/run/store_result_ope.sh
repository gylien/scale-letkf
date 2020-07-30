#!/bin/bash -l

myname=$0
mydir=`dirname $myname`

cd $mydir
. config.main

STIME=$1
TMP=$OUTDIR/dafcst_img
TMPNC=$OUTDIR/dafcst_nc
DEST=$OUTDIR/$STIME/dafcst

STIMEf="${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}:${STIME:10:2}:${STIME:12:2}"

  echo
  echo "Start: store images"

  [ -d $DEST ] || mkdir -p $DEST

# for z in `seq 0 10`; do
#    ZLEV="z"$(printf %02d $z)
#    echo ${ZLEV}000...
####    find ${TMP}/[a,o,f]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs mv -t ${OUTDIR}/${STIME}/dafcst/ &> /dev/null

####    cd ${OUTDIR}/${STIME}/dafcst

    TIMEENDf=`date -ud "1 hour $STIMEf" +"%Y-%m-%d %H:%M:%S"`
    TIMEf=`date -ud "30 second $STIMEf" +"%Y-%m-%d %H:%M:%S"`
   
    rm -f $DEST/list_*.txt
    cd $TMP

while [ `date -ud "$TIMEf" +%s` -le `date -ud "$TIMEENDf" +%s` ] ;do
    echo ${TIMEf}...
    tstamp=`date -ud "$TIMEf" +"%Y%m%d-%H%M%S"`
    find . -name "anal_dbz_${tstamp}*.png" -print >> $DEST/list_anal.txt
    find . -name "obs_dbz_${tstamp}*.png"  -print >> $DEST/list_obs.txt
    find . -name "fcst_dbz_${tstamp}*.png" -print >> $DEST/list_fcst.txt
    [ -f $TMPNC/fcst_ref3d_${tstamp}.nc ] && echo "fcst_ref3d_${tstamp}.nc" >> $DEST/list_nc.txt
    TIMEf=`date -ud "30 second $TIMEf" +"%Y-%m-%d %H:%M:%S"`
done

    cd $TMP
    echo "make tar anal..."
    tar -zcf --remove-files $DEST/anal_ope.tar.gz --files-from $DEST/list_anal.txt    
    echo "make tar obs..."
    tar -zcf --remove-files  $DEST/obs_ope.tar.gz  --files-from $DEST/list_obs.txt    
    echo "make tar fcst..."
    tar -zcf --remove-files $DEST/fcst_ope.tar.gz --files-from $DEST/list_fcst.txt    
    cd $TMPNC
    echo "make tar nc..."
    tar -zcf --remove-files $DEST/ncfile_fcst_ref3d.tar.gz --files-from $DEST/list_nc.txt    

    cd $TMP
    xargs rm < $DEST/list_anal.txt
    xargs rm < $DEST/list_obs.txt
    xargs rm < $DEST/list_fcst.txt
    cd $TMPNC
    xargs rm < $DEST/list_nc.txt

#if [ $z -eq 1 ]; then
#    mv list_[a,f,o]*.txt ./save/
#fi
#####    find [a,o]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs rm -f &> /dev/null

#  done

  rm $DEST/list_*.txt

  cd - > /dev/null

  echo "End: store images"


