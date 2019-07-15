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


###NODE=`expr \( $NMEM  + 2 \) \* 11` ### D3
NODE=`expr \( $NMEM  + 2 \) \* 7` ### D2

####NODE=`expr \( $NMEM  + 2 \) \* 12` ### MODIFIED D2 DOMAIN DECOMP


CONFIG='online_NRT_5.3.X'
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
    sed -e "s/<FCSTLEN>/${FCSTLEN}/g" \
    > config.${SCPNAME}

cat config.main.ofp | \
   sed -e "s/<MEMBER>/${NMEM}/g" | \
   sed -e "s/<NNODES>/${NODE}/g" | \
   sed -e "s/<STIME>/${STIME}/g" \
 > config.main
rm config.main.ofp

cat config/${CONFIG}/sno_bulk.sh | \
   sed -e "s/<STIME>/${STIME}/g" \
 > sno_bulk_${STIME}.sh
cat config/${CONFIG}/sno_bulk_d2.sh | \
   sed -e "s/<STIME>/${STIME}/g" \
 > sno_bulk_d2_${STIME}.sh

chmod 755 sno_bulk*${STIME}.sh


. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

## shorter DT
#cp config.nml.scale.d3.new config.nml.scale.d3
#cp config.nml.scale.new config.nml.scale

### modified domain decomp
#cp config.nml.scale.d2.new config.nml.scale.d2
#cp config.nml.scale_pp.d2.new config.nml.scale_pp.d2
#cp config.nml.scale_init.d2.new config.nml.scale_init.d2


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
./sno_bulk_${STIME}.sh > sno_bulk_${STIME}.log 2>&1 || exit $?
./sno_bulk_d2_${STIME}.sh > sno_bulk_d2_${STIME}.log 2>&1 || exit $?
#-------------------------------------------------------------------------------

#rm -f ${SCPNAME}_job.sh
#rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
