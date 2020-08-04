#!/bin/bash -l

myname=$0
mydir=`dirname $myname`

cd $mydir
. config.main

STIME=$1
TMP=$OUTDIR/dafcst_img
TMPNC=$OUTDIR/dafcst_nc
TMPGRD=$OUTDIR/dafcst
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
    tstamp_jst=`date -ud "9 hour $TIMEf" +"%Y%m%d-%H%M%S"`
    find . -name "anal_dbz_${tstamp}*.png" -print >> $DEST/list_anal.txt
    find . -name "obs_dbz_${tstamp}*.png"  -print >> $DEST/list_obs.txt
    find . -name "fcst_dbz_${tstamp}*.png" -print >> $DEST/list_fcst.txt
    [ -f $TMPNC/${tstamp_jst}.nc ] && echo "${tstamp_jst}.nc" >> $DEST/list_nc.txt
    [ -f $TMPGRD/${tstamp}.grd ] && echo "${tstamp}.grd" >> $DEST/list_grd.txt
   TIMEf=`date -ud "30 second $TIMEf" +"%Y-%m-%d %H:%M:%S"`
done

    cd $TMP
    echo "make tar anal..."
    tar --remove-files -z -v -c -f $DEST/anal_ope.tar.gz --files-from $DEST/list_anal.txt    
    echo "make tar obs..."
    tar --remove-files -z -v -c -f $DEST/obs_ope.tar.gz  --files-from $DEST/list_obs.txt    
    echo "make tar fcst..."
    tar --remove-files -z -v -c -f $DEST/fcst_ope.tar.gz --files-from $DEST/list_fcst.txt    
    cd $TMPNC
    echo "make tar nc..."
    tar --remove-files -z -v -c -f $DEST/ncfile_fcst_ref3d.tar.gz --files-from $DEST/list_nc.txt    
    cd $TMPGRD
    echo "make tar grads..."
    tar --remove-files -v -c -f $DEST/grads_ref3d.tar --files-from $DEST/list_grd.txt    
###    tar --remove-files -z -c -v -f $DEST/grads_ref3d.tar.gz --files-from $DEST/list_grd.txt   ### Too slow to use

#    cd $TMP
#    xargs rm < $DEST/list_anal.txt
#    xargs rm < $DEST/list_obs.txt
#    xargs rm < $DEST/list_fcst.txt
#    cd $TMPNC
#    xargs rm < $DEST/list_nc.txt
#    cd $TMPGRD
#    xargs rm < $DEST/list_grd.txt


#if [ $z -eq 1 ]; then
#    mv list_[a,f,o]*.txt ./save/
#fi
#####    find [a,o]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs rm -f &> /dev/null

#  done

  rm $DEST/list_*.txt

  cd - > /dev/null

  echo "End: store results"


