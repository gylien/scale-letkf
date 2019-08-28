#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

#if (($# < 5)); then
#  echo "$0: Insufficient arguments" >&2
#  exit 1
#fi

#SCPNAME="$1"; shift
PARENT_REF_TIME_D2="$1"; shift
PARENT_REF_TIME="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
#ETIME="$1"; shift
#TIME_DT="$1"; shift
#TIME_DT_DYN="$1"; shift
#NNODES="$1"; shift
WTIME_L="$1"; shift
NMEM="$1"

#PARENT_REF_TIME=20190608180000
#STIME=20190608180000
#FCSTLEN=600
#WTIME_L=00:10:00
#NMEM=2 
#NP_OFILE_X=2
#NP_OFILE_Y=2


SCPNAME=fcst
ETIME="$STIME"

#NODE=`expr \( $NMEM  + 2 \) \* 16` ### D3 
NODE=`expr \( $NMEM  + 2 \) \* 4` ### D3 
###WTIME_L="00:20:00"

CONFIG='realtime_fcst_D4_1km'
PRESET='OFP'

#-------------------------------------------------------------------------------

 config_suffix='ofp'
 script_suffix='_ofp'

if [ "$SCPNAME" = 'cycle' ]; then
  DATA_BDY_WRF="ncepgfs_wrf_da"
  DATA_BDY_GRADS="ncepgfs_grads_da"
else
  DATA_BDY_WRF="ncepgfs_wrf"
  DATA_BDY_GRADS="ncepgfs_grads"
fi

ntry=1
while [ $ntry -le 3 ] ;do


#-------------------------------------------------------------------------------

### wait until the submittion of previous jobs are completed

iwait=1
while [ $iwait == 1 ];do
iwait=0
running_jobs=`ls -x fcst_ofp.stat.*`
if [ "$running_jobs" != "" ];then
for statfile in $running_jobs ;do
if [ `cat $statfile` == "prep" ] ;then
 iwait=1
 sleep 13s
fi
done
fi
[ -s temp.lock ] && iwait=1
done

echo $PARENT_REF_TIME_D2 $PARENT_REF_TIME $STIME > temp.lock

#-------------------------------------------------------------------------------

###rm -f config.*

cp config/${CONFIG}/config.* .

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<FCSTLEN>/${FCSTLEN}/g" \
    > config.${SCPNAME}

cat config.main.ofp | \
   sed -e "s/<MEMBER>/${NMEM}/g" | \
   sed -e "s/<NNODES>/${NODE}/g" | \
   sed -e "s/<STIME>/${STIME}/g" | \
   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" | \
   sed -e "s/<PARENT_REF_TIME_D2>/${PARENT_REF_TIME_D2}/g" \
 > config.main
rm config.main.ofp

### NCDF file merge
#cat config/${CONFIG}/sno_bulk.sh | \
#   sed -e "s/<STIME>/${STIME}/g" | \
#   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" | \
#   sed -e "s/<NP_OFILE_X>/${NP_OFILE_X}/g" | \
#   sed -e "s/<NP_OFILE_Y>/${NP_OFILE_Y}/g" \
# > sno_bulk_ref_${PARENT_REF_TIME}_${STIME}.sh

#chmod 750 sno_bulk_ref_${PARENT_REF_TIME}_${STIME}.sh

. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${PARENT_REF_TIME}.${STIME} 2>&1


#-------------------------------------------------------------------------------

res=$?
if [ $res -eq 77 ] ; then
 NODE=`expr $NODE \/ 2` ### D3 

 sec1=`date -d "$WTIME_L" +%s`
 sec0=`date -d "00:00:00" +%s`
 wtime_sec=`expr \( $sec1 - $sec0 \) \* 2`
 WTIME_L=`date -d "2000/1/1 $wtime_sec second" +%H:%M:%S`

 ntry=`expr $ntry + 1`
elif [ $res -ne 0 ] ; then
 echo 'abort : res=' $res 
 rm ${SCPNAME}${script_suffix}.stat.${PARENT_REF_TIME}.${STIME} 
exit $res
else
 echo 'done.'
 break
fi

done

if [ $ntry -eq 4 ];then
 echo 'abort : OFP is too crowded!'
 exit 1
 rm ${SCPNAME}${script_suffix}.stat.${PARENT_REF_TIME}.${STIME} 
fi
#-------------------------------------------------------------------------------

  jobname="${SCPNAME}_${SYSNAME}"
  jobid=$(grep 'pjsub Job' ${SCPNAME}${script_suffix}.log.${PARENT_REF_TIME}.${STIME} | cut -d ' ' -f6)
  logdir="$OUTDIR/exp/${jobid}_${SCPNAME}_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  jobinfo="$logdir/job.i"


#if [ ! -e "$stdout" ] || [ ! -e "$stderr" ]; then
#  exit 101
#fi
#if [ -z "$(tail -n 1 $stderr | grep "Finish ${SCPNAME}.sh")" ]; then
#  exit 102
#fi

#-------------------------------------------------------------------------------

#rm -f ${SCPNAME}_job.sh
#rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

##### ./sno_bulk_ref_${PATENT_REF_TIME}_${STIME}.sh > ./sno_bulk_${PARENT_REF_TIME}_${STIME}.log 2>&1 || exit $?

#-------------------------------------------------------------------------------

exit 0
