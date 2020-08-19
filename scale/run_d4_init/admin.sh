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


SCPNAME=fcst
ETIME="$STIME"

CONFIG='realtime_fcst_D4_500m'
PRESET=`hostname | cut -d '.' -f 2 | tr '[a-z]' '[A-Z]'`

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'OFP' ]; then
  if [ "$CONFIG" ==  "realtime_fcst_D4_1km" ];then
    NNODES=`expr \( $NMEM + 2 \) ` 
  elif [ "$CONFIG" ==  "realtime_fcst_D4_500m" ];then
    NNODES=`expr \( $NMEM + 2 \) \* 16`
  elif [ "$CONFIG" ==  "realtime_fcst_D4_250m" ];then
    NNODES=`expr \( $NMEM + 2 \) \* 64` 
  else
    echo "CONFIG="$CONFIG" not supported."
    exit 1 
  fi
    while [ $NNODES -gt 256 ] ;do
      NNODES=`expr $NNODES \/ 2`
    done
 config_suffix='ofp'
 script_suffix='_ofp'
elif [ "$PRESET" = 'OBCX' ]; then
  if [ "$CONFIG" ==  "realtime_fcst_D4_1km" ];then
    NNODES=`expr \( $NMEM + 2 \) \* 4` 
  elif [ "$CONFIG" ==  "realtime_fcst_D4_500m" ];then
    NNODES=`expr \( $NMEM + 2 \) \* 64`
  else
    echo "CONFIG="$CONFIG" not supported."
    exit 1 
  fi
  NNODES=`expr \( $NMEM + 2 \) \* 4` ### D4 1km
    while [ $NNODES -gt 256 ] ;do
      NNODES=`expr $NNODES \/ 2`
    done
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

ntry=1
while [ $ntry -le 3 ] ;do


#-------------------------------------------------------------------------------

### wait until the submittion of previous jobs are completed

iwait=1
res=`grep prep fcst${script_suffix}.stat.*`
mytime=${PARENT_REF_TIME}.${STIME}
if [ -s waiting_list ] || [ ! -z "$res" ];then
 echo $mytime >> waiting_list
else
iwait=0
fi

while [ $iwait == 1 ];do
if [ -s waiting_list ] ;then
 next=`cat waiting_list | head -n 1`
 res=`grep prep fcst${script_suffix}.stat.*`
 if [ "$next" == "$mytime" ] && [ -z "$res" ] ;then 
   iwait=0 
   wcl=`cat waiting_list| wc -l`
   wcl=`expr $wcl - 1`
   cp waiting_list temp_list
   cat temp_list | tail -n $wcl > waiting_list 
   rm temp_list
 else
  sleep 23s
 fi
else
 iwait=0 ### never occur
 echo 'ERROR'
 exit 1
fi
done

[ -f waiting_list ] && [ ! -s waiting_list ] && rm waiting_list ### clean

echo 'prep' > fcst${script_suffix}.stat.$mytime


#-------------------------------------------------------------------------------

###rm -f config.*

cp config/${CONFIG}/config.* .

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<FCSTLEN>/${FCSTLEN}/g" \
    > config.${SCPNAME}

cat config.main.${config_suffix} | \
   sed -e "s/<MEMBER>/${NMEM}/g" | \
   sed -e "s/<NNODES>/${NNODES}/g" | \
   sed -e "s/<STIME>/${STIME}/g" | \
   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" | \
   sed -e "s/<PARENT_REF_TIME_D2>/${PARENT_REF_TIME_D2}/g" \
 > config.main
rm config.main.${config_suffix}


. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${PARENT_REF_TIME}.${STIME} 2>&1


#-------------------------------------------------------------------------------

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
cd $TOPDIR/scale_ope/scale-letkf_ope_d4/scale/run
./prep.sh $STIME
cd -

#-------------------------------------------------------------------------------

exit 0
