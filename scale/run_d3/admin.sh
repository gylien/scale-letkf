#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

#if (($# < 5)); then
#  echo "$0: Insufficient arguments" >&2
#  exit 1
#fi

#SCPNAME="$1"; shift
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

if [ $NMEM == 'mdet' ];then
 MEMBERS='mdet'
 NMEM=1
else
 MEMBERS='all'
fi

CONFIG='realtime_fcst_D3'
PRESET='OFP'

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'OFP' ]; then
 if [ "$MEMBERS" = "mdet" ];then
   NNODES=`expr \( $NMEM  \) \* 16` ### D2
 else
   NNODES=`expr \( $NMEM + 2 \) \* 16` ### D2
 fi
 config_suffix='ofp'
 script_suffix='_ofp'
elif [ "$PRESET" = 'OBCX' ]; then
 if [ "$MEMBERS" = "mdet" ];then
   NNODES=`expr \( $NMEM  \) \* 64` ### D2
 else
   NNODES=`expr \( $NMEM + 2 \) \* 64` ### D2
    while [ $NNODES -gt 256 ] ;do
      NNODES=`expr $NNODES \/ 2`
    done
  fi
 config_suffix='obcx'
 script_suffix='_obcx'
fi

ntry=1
while [ $ntry -le 3 ] ;do

if [ "$SCPNAME" = 'cycle' ]; then
  DATA_BDY_WRF="ncepgfs_wrf_da"
  DATA_BDY_GRADS="ncepgfs_grads_da"
else
  DATA_BDY_WRF="ncepgfs_wrf"
  DATA_BDY_GRADS="ncepgfs_grads"
fi


#-------------------------------------------------------------------------------

### wait until the submittion of previous jobs are completed

iwait=1
while [ $iwait == 1 ];do
iwait=0
running_jobs=`ls -x fcst${script_suffix}.stat.* 2>/dev/null`
if [ ! -z "$running_jobs" ];then
for statfile in $running_jobs ;do
if [ "`cat $statfile`" == "prep" ] ;then
 iwait=1
 sleep 23s
fi
done
fi
[ -s temp.lock ] && iwait=1 
done

#-------------------------------------------------------------------------------

###rm -f config.*

echo ${PARENT_REF_TIME}.${STIME} > temp.lock

cp config/${CONFIG}/config.* .

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<MEMBERS>/${MEMBERS}/g" | \
    sed -e "s/<FCSTLEN>/${FCSTLEN}/g" \
    > config.${SCPNAME}

cat config.main.${config_suffix} | \
   sed -e "s/<MEMBER>/${NMEM}/g" | \
   sed -e "s/<NNODES>/${NNODES}/g" | \
   sed -e "s/<STIME>/${STIME}/g" | \
   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" \
 > config.main
rm config.main.${config_suffix}

### NCDF file merge
cat config/${CONFIG}/sno_bulk.sh | \
   sed -e "s/<STIME>/${STIME}/g" | \
   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" | \
   sed -e "s/<FCSTLEN>/${FCSTLEN}/g" | \
   sed -e "s/<NP_OFILE_X>/${NP_OFILE_X}/g" | \
   sed -e "s/<NP_OFILE_Y>/${NP_OFILE_Y}/g" \
 > sno_bulk_ref_${PARENT_REF_TIME}_${STIME}.sh

chmod 750 sno_bulk_ref_${PARENT_REF_TIME}_${STIME}.sh

. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

cat config/${CONFIG}/common_d3.h > plot/common_d3.h
cat "character*120,parameter::cdir_base_fcst= \"${OUTPUT}/${EXP3}/\"," >>plot/common_d3.h  


#-------------------------------------------------------------------------------

echo 'exec job'

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${PARENT_REF_TIME}.${STIME} 2>&1

res=$?

if [ $res -eq 77 ] ; then
 NNODES=`expr $NNODES \/ 2` ### D3 

 sec1=`date -d "$WTIME_L" +%s`
 sec0=`date -d "00:00:00" +%s`
 wtime_sec=`expr \( $sec1 - $sec0 \) \* 2`
 WTIME_L=`date -d "2000/1/1 $wtime_sec second" +%H:%M:%S`

 ntry=`expr $ntry + 1`
 echo "retry :: NNODES="$NNODES" WTIME_L="$WTIME_L
elif [ $res -ne 0 ] ; then
 echo 'abort : res='$res
 exit $res
else
 echo 'exec finish successfully.'
 break
fi

done

if [ $ntry -eq 4 ] ;then
 echo "abort : OFP is too crowded."
 exit 1 
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

echo 'exec sno'
 ./sno_bulk_ref_${PARENT_REF_TIME}_${STIME}.sh > ./sno_bulk_${PARENT_REF_TIME}_${STIME}.log 2>&1 &

#-------------------------------------------------------------------------------

echo 'exec d4 init...'
spinup_d3=3600
intv_d4=3600

cd ../run_d4_init
 STIME_f=`date -d "${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}:${STIME:10:2}:${STIME:12:2}" +"%F %T"`
 STIME_D4_f=`date -d "${spinup_d3} second ${STIME_f}" +"%F %T"`
 while [ `date -d "$STIME_D4_f" +%s` -le `date -d "${FCSTLEN} second -${intv_d4} second ${STIME_f}" +%s` ] ;do
   STIME_D4=`date -d "${STIME_D4_f}" +%Y%m%d%H%M%S`
   sleep 17s
   ./admin.sh ${PARENT_REF_TIME} ${STIME} ${STIME_D4} ${intv_d4} "00:25:00" $NMEM &>admin.log.${PARENT_REF_TIME}.${STIME}.${STIME_D4} &
 STIME_D4_f=`date -d "${intv_d4} second ${STIME_D4_f}" +"%F %T"`
 done
cd -

#-------------------------------------------------------------------------------

exit 0
