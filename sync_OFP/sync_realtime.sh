#!/bin/sh

mydir=`dirname $0`

source $mydir/admin.rc

RUNDIR_OFP=$SCALEDIR_OFP/scale_ope/scale-letkf_ope_d4/scale/run/dacycle_1km ### 1km
FCSTBASE=$DATADIR/d4/realtime
TEMPDIR=$DATADIR/d4/temp

wait_sec=5
limit_sec=600
levels='z01127m z01957m z03048m z04480m z06360m z08829m'
num_fcst_past=20

sec=0

rm -r $TEMPDIR
mkdir -p $TEMPDIR
for level in $levels ;do
FCSTDIR=$FCSTBASE/$level
[ -d $FCSTDIR ] || mkdir -p $FCSTDIR
rm $FCSTDIR/*.png
done

touch previous.time

while [ $sec -lt $limit_sec ]; do

for level in $levels ;do
FCSTDIR=$FCSTBASE/$level
#rsync -e "$sshcommand" $hostname:$RUNDIR_OFP/obs_dbz_[0-9]{8}-[0-9]{6}_${level}_0001.png $TEMPDIR/  
#rsync -e "$sshcommand" $hostname:$RUNDIR_OFP/anal_dbz_[0-9]{8}-[0-9]{6}_${level}_0001.png $TEMPDIR/  
#rsync -e "$sshcommand" $hostname:$RUNDIR_OFP/fcst_dbz_[0-9]{8}-[0-9]{6}_FT[0-9]{4}s_${level}_0001.png $TEMPDIR/  
rsync -e "$sshcommand" --no-t "$hostname:$RUNDIR_OFP/obs_dbz_*_${level}_0001.png" $TEMPDIR/  
rsync -e "$sshcommand" --no-t "$hostname:$RUNDIR_OFP/anal_dbz_*_${level}_0001.png" $TEMPDIR/  
rsync -e "$sshcommand" --no-t "$hostname:$RUNDIR_OFP/fcst_dbz_*_${level}_0001.png" $TEMPDIR/  


##### from new files  
### find latest anal and rename it and copy to FCSTDIR
list=`find $TEMPDIR -name anal_dbz_*_${level}_0001.png -newer previous.time`
if [ ! -z "$list" ];then
  tlabelf_latest=0
  for file in $list ;do
    file_anal=`basename $file`
    file_obs=`echo $file_anal | sed -e "s/anal/obs/g"`
    tlabel=`echo $file_anal | grep -o '[0-9]\{8\}-[0-9]\{6\}'`
    tlabelf=${tlabel:0:8}${tlabel:9:6}
    [ $tlabelf -gt $tlabelf_latest ] && tlabelf_latest=$tlabelf
    cp $TEMPDIR/$file_anal $FCSTDIR/anal_${tlabelf}.png
    cp $TEMPDIR/$file_obs  $FCSTDIR/obs_${tlabelf}.png
    chmod 644 $FCSTDIR/anal_${tlabelf}.png
    chmod 644 $FCSTDIR/obs_${tlabelf}.png
  done
### remove fcst in FCSTDIR behind anal time
  list=`find $FCSTDIR -name fcst_*.png`
  if [ ! -z "$list" ]; then
   for file in $list; do
    file_fcst=`basename $file`
    fcst_base=`echo $file_fcst | grep -o '[0-9]\{14\}'`    
   [ $fcst_base -le $tlabelf_latest ] && rm $FCSTDIR/$file_fcst
   done
  fi
fi

tlabel_latest="${tlabelf_latest:0:8}-${tlabelf_latest:8:6}"
date_latest="${tlabelf_latest:0:4}-${tlabelf_latest:4:2}-${tlabelf_latest:6:2} ${tlabelf_latest:8:2}:${tlabelf_latest:10:2}:${tlabelf_latest:12:2}"

### find fcsts from latest anal and rename /copy to FCSTDIR
list_base=`ls -1t $FCSTDIR/anal_*.png | head -n $num_fcst_past | tac `
if [ ! -z "$list_base" ];then
for file_base in $list_base;do
tlabelf_base=`basename $file_base | grep -o '[0-9]\{14\}'`
date_base="${tlabelf_base:0:4}-${tlabelf_base:4:2}-${tlabelf_base:6:2} ${tlabelf_base:8:2}:${tlabelf_base:10:2}:${tlabelf_base:12:2}"
list=`find $TEMPDIR -name fcst_dbz_${tlabelf_base:0:8}-${tlabelf_base:8:6}_*_${level}_0001.png -newer previous.time`
#echo $file_base $tlabelf_base
#echo fcst_dbz_${tlabelf_base:0:8}-${tlabelf_base:8:6}
#echo $list
#exit 
if [ ! -z "$list" ] ; then
#  rm $FCSTDIR/fcst_*
  for file in $list;do
    file_fcst=`basename $file`
    fcst_sec=`echo $file_fcst | grep -o 'FT[0-9]\{4\}' | cut -c 3-6`
    tlabelf=`date -ud "$fcst_sec sec $date_base" +%Y%m%d%H%M%S`
     [ $tlabelf -gt $tlabelf_latest ] && cp $TEMPDIR/$file_fcst $FCSTDIR/fcst_${tlabelf}.png
    chmod 644 $FCSTDIR/fcst_${tlabelf}.png
  done
fi

done
fi

done ### levels

touch previous.time

cd $mydir/..
php ./monitor_plot_d4.php
cd - 
sec=`expr $sec + $wait_sec`

done ### while

rm previous.time
