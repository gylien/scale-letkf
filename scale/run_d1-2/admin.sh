#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 4)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

#SCPNAME="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
#ETIME="$1"; shift
#TIME_DT="$1"; shift
#TIME_DT_DYN="$1"; shift
#NNODES="$1"; shift
WTIME_L="$1"; shift
NMEM="$1"


SCPNAME=fcst
ETIME="$STIME"

if [ $NMEM == 'mdet' ];then
 MEMBERS='mean mdet'
 NMEM=2
else
 MEMBERS='all'
fi

CONFIG='online_NRT_5.3.X'
PRESET='OFP'

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'OFP' ]; then
 if [ "$MEMBERS" = "mean mdet" ];then
   NNODES=`expr \( $NMEM  \) \* 7` ### D2
 else
   NNODES=`expr \( $NMEM + 2 \) \* 7` ### D2
 fi
 config_suffix='ofp'
 script_suffix='_ofp'
elif [ "$PRESET" = 'OBCX' ]; then
 if [ "$MEMBERS" = "mean mdet" ];then
   NNODES=`expr \( $NMEM  \) \* 16` ### D2
 else
   NNODES=`expr \( $NMEM + 2 \) \* 16` ### D2
    while [ $NNODES -gt 256 ] ;do
      NNODES=`expr $NNODES \/ 2`
    done
  fi
 config_suffix='obcx'
 script_suffix='_obcx'
fi

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
running_jobs=`ls -x fcst${config_suffix}.stat.*`
if [ "$running_jobs" != "" ];then
for statfile in $running_jobs ;do
if [ `cat $statfile` == "prep" ] ;then
 iwait=1
 sleep 10s
fi
done
fi
done

#-------------------------------------------------------------------------------

###rm -f config.*


cp config/${CONFIG}/config.* .

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<MEMBERS>/${MEMBERS}/g" | \
    sed -e "s/<FCSTLEN>/${FCSTLEN}/g" \
    > config.${SCPNAME}

cat config.main.ofp | \
   sed -e "s/<MEMBER>/${NMEM}/g" | \
   sed -e "s/<NNODES>/${NNODES}/g" | \
   sed -e "s/<STIME>/${STIME}/g" \
 > config.main
rm config.main.ofp


. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${STIME} 2>&1 || exit $?

#-------------------------------------------------------------------------------

  jobname="${SCPNAME}_${SYSNAME}"
  jobid=$(grep 'pjsub Job' ${SCPNAME}${script_suffix}.log.${STIME} | cut -d ' ' -f6)
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
cd post
./post${script_suffix}.sh ${STIME} &> post.log.${STIME}
cd -
##-------------------------------------------------------------------------------

#rm -f ${SCPNAME}_job.sh
#rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
