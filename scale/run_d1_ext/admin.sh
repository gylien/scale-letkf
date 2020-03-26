#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 3)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

#SCPNAME="$1"; shift
SCPNAME=fcst
STIME="$1"; shift
WTIME_L="$1"; shift
MEMBER="$1"

ETIME="$STIME"
TIME_DT=40.0D0
TIME_DT_DYN=8.0D0

CONFIG='realtime_ope_d1'
PRESET=`hostname | cut -d '.' -f 2 | tr '[a-z]' '[A-Z]'`

FP_SUFFIX='_single'

if [ $MEMBER == 'mdet' ] ;then
  MEMBERS='mdet'
  MEMBER=1
  NNODES=3
else
  MEMBERS='all'
fi

#-------------------------------------------------------------------------------

#if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ]; then
#  config_suffix='K'
#  script_suffix='_K'
#elif [ "$PRESET" = 'K_micro' ]; then
#  config_suffix='K'
#  script_suffix='_K_micro'
if [ "$PRESET" = "OFP" ]; then
  config_suffix='ofp'
  script_suffix='_ofp'
  if [ $MEMBERS='mdet' ];then
    NNODES=3
  else
    NNODES=`\( $MEMBER + 2 \) * 3`
  fi
elif [ "$PRESET" = "OBCX" ]; then
  config_suffix='obcx'
  script_suffix='_obcx'
  if [ $MEMBERS='mdet' ];then
    NNODES=8
  else
    NNODES=`\( $MEMBER + 2 \) * 8`
    while [ $NNODES -gt 256 ] ;do
      NNODES=`expr $NNODES \/ 2`
    done
   fi
else
  echo "[Error] Unsupported \$PRESET" >&2
  exit 1
fi

if [ "$SCPNAME" = 'cycle' ]; then
  DATA_BDY_WRF="ncepgfs_wrf_da"
  DATA_BDY_GRADS="ncepgfs_grads_da"
else
  DATA_BDY_WRF="ncepgfs_wrf"
#  DATA_BDY_GRADS="ncepgfs_grads"
  DATA_BDY_GRADS="ncepgfs_grads_ext"
fi

#-------------------------------------------------------------------------------
 
### wait until the submittion of previous jobs are completed

iwait=1
res=`grep prep fcst${script_suffix}.stat.*`
mytime=${STIME}
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

rm -f config.main
rm -f config.${SCPNAME}
rm -f config.nml.*

cat config/${CONFIG}/config.main.${config_suffix} | \
    sed -e "s/<PRESET>/${PRESET}/g" | \
    sed -e "s/<DATA_BDY_WRF>/${DATA_BDY_WRF}/g" | \
    sed -e "s/<DATA_BDY_GRADS>/${DATA_BDY_GRADS}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" | \
    sed -e "s/<MEMBER>/${MEMBER}/g" | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<FP_SUFFIX>/${FP_SUFFIX}/g" \
    > config.main

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" \
    > config.${SCPNAME}

cat config/${CONFIG}/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

ln -fs config/${CONFIG}/config.nml.ensmodel .
ln -fs config/${CONFIG}/config.nml.letkf .
ln -fs config/${CONFIG}/config.nml.scale_pp .
ln -fs config/${CONFIG}/config.nml.scale_init .
if [ -e "config/${CONFIG}/config.nml.scale_user" ]; then
  ln -fs config/${CONFIG}/config.nml.scale_user .
fi
if [ -e "config/${CONFIG}/config.nml.obsope" ]; then
  ln -fs config/${CONFIG}/config.nml.obsope .
fi
if [ -e "config/${CONFIG}/config.nml.grads_boundary" ]; then
  ln -fs config/${CONFIG}/config.nml.grads_boundary .
#cat config/${CONFIG}/config.nml.grads_boundary | \
#    sed -e "s/--FNAME_SFC--/--DIR--\/bdysfc/g" | \
#    sed -e "s/--FNAME_ATMOS--/--DIR--\/bdyatm/g" | \
#    sed -e "s/--FNAME_LAND--/--DIR--\/bdyland/g"  \
#    > config.nml.grads_boundary
fi

. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${STIME} 2>&1 || exit $?

#-------------------------------------------------------------------------------

#if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ] || [ "$PRESET" = 'K_micro' ] || [ "$PRESET" = "OFP" ]; then
  jobname="${SCPNAME}_${SYSNAME}"
  jobid=$(grep 'pjsub Job' ${SCPNAME}${script_suffix}.log.${STIME} | cut -d ' ' -f6)
  logdir="$OUTDIR/exp/${jobid}_${SCPNAME}_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  jobinfo="$logdir/job.i"
#  stdout="${jobname}.o${jobid}"
#  stderr="${jobname}.e${jobid}"
#  jobinfo="${jobname}.i${jobid}"
#fi

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
#-------------------------------------------------------------------------------

rm -f ${SCPNAME}_job.sh
rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
