#!/bin/bash -l

. config.main

file_list=`ls -xd dacycle_1km_*`

for TMP in $file_list ;do


STIME=${TMP:12:25}

if [ $STIME -lt 20200713150000 ] ;then

DEST=$OUTDIR/$STIME/dafcst

  echo
  echo "Start $TMP"

  [ -d $DEST ] || mkdir -p $DEST
 
  cd $TMP

#  mkdir -p save

  for z in `seq 0 10`; do
    ZLEV="z"$(printf %02d $z)
    echo ${ZLEV}000...
####    find ${TMP}/[a,o,f]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs mv -t ${OUTDIR}/${STIME}/dafcst/ &> /dev/null

####    cd ${OUTDIR}/${STIME}/dafcst
    find . -name "anal_*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_anal.txt
    find . -name "obs_*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_obs.txt
    find . -name "fcst_*${ZLEV}[0-9][0-9][0-9]m_*.png" -print > list_fcst.txt
    tar -zcf anal_${ZLEV}000m.tar.gz --files-from list_anal.txt    
    tar -zcf obs_${ZLEV}000m.tar.gz --files-from list_obs.txt    
    tar -zcf fcst_${ZLEV}000m.tar.gz --files-from list_fcst.txt    

    mv [a,o,f]*.tar.gz $DEST/

    xargs rm < list_anal.txt
    xargs rm < list_obs.txt
    xargs rm < list_fcst.txt

#if [ $z -eq 1 ]; then
#    mv list_[a,f,o]*.txt ./save/
#fi
#####    find [a,o]*_${ZLEV}[0-9][0-9][0-9]m_*.png -type f | xargs rm -f &> /dev/null

  done
  rm list_*.txt

  cd - > /dev/null


  fi

  done

  echo "End: store images"


