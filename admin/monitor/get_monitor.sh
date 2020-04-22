#!/bin/bash -l


wkdir="$( cd "$( dirname "$0" )" && pwd )"

cd ${wkdir}

TOPDIR="${wkdir}/../.."
OBSDIR=$TOPDIR/external/ncepobs_gdas_letkf
GFSDIR=$TOPDIR/external/ncepgfs_grads
OUTDIR=$TOPDIR


### OFP runtime directories 
RUNBASE=$TOPDIR/scale_ope/scale-letkf_ope/scale
rundir_cycle=$RUNBASE/run
rundir_fcst=$RUNBASE/run_d1_ext
rundir_fcst_d2=$RUNBASE/run_d1-2
rundir_fcst_d3=$RUNBASE/run_d3


now="$(date -u +'%Y-%m-%d %H:%M:%S')"
nowf="$(date -u  +'%Y%m%d%H')"
ystf="$(date -ud "-30 hour $now" +'%Y%m%d%H')"

oldest=$ystf
HH=`echo $oldest | cut -c 9-10` 
while [ $HH != 00  ] && [ $HH != 06 ] && [ $HH != 12 ] && [ $HH != 18 ] ;do
 oldest=`expr $oldest - 1`
 HH=`echo $oldest | cut -c 9-10` 
done
 HH=`echo $oldest | cut -c 9-10` 
 YMD=`echo $oldest | cut -c 1-8` 


### NCEP prepbufr and GFS

echo 'monitor prepbufr and GFS ...'

[ -e monitor_obs.txt ] && rm monitor_obs.txt
[ -e monitor_gfs.txt ] && rm monitor_gfs.txt

iftime=$oldest
while [ $iftime -le $nowf ] ;do
# echo $iftime $nowf
  if [ -d $OBSDIR/$iftime ]; then 
    tstmp=`ls -ld $OBSDIR/$iftime | awk '{print $8}'`
    echo $iftime $tstmp >> monitor_obs.txt  
  fi

### file size check
atmsize=161150400
sfcsize=7277760
landsize=10396800
regsize=3189540

 HH=`echo $iftime | cut -c 9-10` 
 YMD=`echo $iftime | cut -c 1-8` 
 iftime_next=`date -ud "6 hour $YMD $HH" +'%Y%m%d%H'`

  if [ -d $GFSDIR/$iftime ]; then 
    tstmp=`ls -ld $GFSDIR/$iftime | awk '{print $8}'`
    atmtest=`ls -l $GFSDIR/$iftime/atm_${iftime}0000.grd | awk '{print $5}'`
    sfctest=`ls -l $GFSDIR/$iftime/sfc_${iftime}0000.grd | awk '{print $5}'`
    landtest=`ls -l $GFSDIR/$iftime/land_${iftime}0000.grd | awk '{print $5}'`
    regtest=`ls -l $GFSDIR/$iftime/reg_${iftime}0000.grd | awk '{print $5}'`
 if [ $atmtest != $atmsize ] || [ $sfctest != $sfcsize ] || [ $landtest != $landsize ] || [ $regtest != $regsize ] ;then
    echo $iftime $tstmp C >> monitor_gfs.txt  
  else
    atmtest=`ls -l $GFSDIR/$iftime/atm_${iftime_next}0000.grd | awk '{print $5}'`
    sfctest=`ls -l $GFSDIR/$iftime/sfc_${iftime_next}0000.grd | awk '{print $5}'`
    landtest=`ls -l $GFSDIR/$iftime/land_${iftime_next}0000.grd | awk '{print $5}'`
    regtest=`ls -l $GFSDIR/$iftime/reg_${iftime_next}0000.grd | awk '{print $5}'`
    if [ $atmtest != $atmsize ] || [ $sfctest != $sfcsize ] || [ $landtest != $landsize ] || [ $regtest != $regsize ] ;then
     echo $iftime $tstmp C >> monitor_gfs.txt  
    else
     echo $iftime $tstmp >> monitor_gfs.txt  
    fi
  fi
 fi 

 HH=`echo $iftime | cut -c 9-10` 
 YMD=`echo $iftime | cut -c 1-8` 
 iftime=`date -ud "6 hour $YMD $HH" +'%Y%m%d%H'`
done

if [ "$noOFP" == "T" ] ;then
 scp -i /home/amemiya/.ssh/id_rsa_mac monitor_*.txt ${web_url}:${web_basedir}
 exit 
fi

### ANALYSIS and FORECAST

echo 'monitor D1 analysis ...'

rm monitor_cycle.txt
if [ -s ../admin_cycle.lock ] ; then
 itime=`cat ../admin_cycle.time`
 itimef=`date -ud "$itime" +'%Y%m%d%H%M%S'`
 cd $rundir_cycle
 stat=`./monitor_cycle.sh ${itimef}`
 cd $wkdir
 echo $itimef $stat
 echo $itimef "$stat" > monitor_cycle.txt
else
 touch monitor_cycle.txt
fi

echo 'monitor D1 fcst ...'

lockfiles=`ls -x ../admin_fcst_d1_ext.lock.* 2>/dev/null`
 rm monitor_fcst.txt
if [ "$lockfiles" != "" ] ; then
 for lockfile in $lockfiles;do
 itime=`echo $lockfile | rev | cut -c 1-14 | rev`
  cd $rundir_fcst 
  stat=`./monitor_fcst.sh ${itime}`
  echo $itime $stat
  cd $wkdir
 echo $itime $stat >> monitor_fcst.txt
 done
else
 touch monitor_fcst.txt
fi

echo 'monitor D2 fcst ...'

lockfiles=`ls -x ../admin_fcst_d1-2.lock.* 2>/dev/null`
 rm monitor_fcst_d2.txt
if [ "$lockfiles" != "" ] ; then
 for lockfile in $lockfiles;do
  itime=`echo $lockfile | rev | cut -c 1-14 | rev`
  stat=`cd $rundir_fcst_d2 && ./monitor_fcst.sh ${itime}`
  cd $wkdir
  echo $itime $stat >> monitor_fcst_d2.txt
 done
else
 touch monitor_fcst_d2.txt
fi

echo 'monitor D3 fcst ...'

lockfiles=`ls -x ../admin_fcst_d3.lock.* 2>/dev/null`
 rm monitor_fcst_d3.txt
if [ "$lockfiles" != "" ] ; then
 for lockfile in $lockfiles;do
  itime=`echo $lockfile | rev | cut -c 1-14 | rev`
  itimeref=`echo $lockfile | rev | cut -c 16-29 | rev`
#  echo 'itimeref itime' $itimeref $itime 
#  exit 0
  stat=`cd $rundir_fcst_d3 && ./monitor_fcst.sh ${itimeref} ${itime}`
#  echo 'stat = ' $stat
  echo $itimeref $itime $stat >> monitor_fcst_d3.txt
#  echo $stat >> stat.txt
 done
else
 touch monitor_fcst_d3.txt
fi
cd $wkdir


### ANALYSIS (TORI AEZU)

time_cycle_old=`tail -n 1 monitor_cycle_temp.txt | cut -c 1-14`

[ -e monitor_cycle_temp.txt ] && rm monitor_cycle_temp.txt
touch monitor_cycle_temp.txt

iftime=$oldest
running=''

while [ $iftime -le $nowf ] ;do
  if [ -s $OUTDIR/result/ope/d1/${iftime}0000/anal/mean/init.pe000000.nc ]; then 
    tstmp=`ls -l $OUTDIR/result/ope/d1/${iftime}0000/anal/mean/init.pe000000.nc | awk '{print $8}'`
    echo $iftime $tstmp >> monitor_cycle_temp.txt
  fi
 HH=`echo $iftime | cut -c 9-10` 
 YMD=`echo $iftime | cut -c 1-8` 
 iftime=`date -ud "6 hour $YMD $HH" +'%Y%m%d%H'`
done

time_cycle=`tail -n 1 monitor_cycle_temp.txt | cut -c 1-14`
if [ "$time_cycle" != "$time_cycle_old" ]; then
  echo 'monitor D1 LETKF performance ...'
  cd ./letkf
  ./get_monitor_letkf.sh
  cd ..
  cd ./realtime
  ./get_monitor_realtime_d1-4.sh
  cd ..
fi
