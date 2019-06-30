#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 5)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

#SCPNAME="$1"; shift
PARENT_REF_TIME="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
#ETIME="$1"; shift
#TIME_DT="$1"; shift
#TIME_DT_DYN="$1"; shift
#NNODES="$1"; shift
WTIME_L="$1"; shift
NMEM="$1" ; shift
NP_OFILE_X="$1" ; shift
NP_OFILE_Y="$1"


SCPNAME=fcst
ETIME="$STIME"

###
NODE=`expr \( $NMEM  + 2 \) \* 16` ### D3 
###NODE=`expr \( $NMEM  + 2 \) \* 4`

###WTIME_L="03:00:00"


CONFIG='realtime_fcst_D3'
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
   sed -e "s/<PARENT_REF_TIME>/${PARENT_REF_TIME}/g" \
 > config.main
rm config.main.ofp

### NCDF file merge
cat config/${CONFIG}/sno_bulk.sh | \
   sed -e "s/<STIME>/${STIME}/g" | \
   sed -e "s/<NP_OFILE_X>/${NP_OFILE_X}/g" | \
   sed -e "s/<NP_OFILE_Y>/${NP_OFILE_Y}/g" \
 > sno_bulk.sh

. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log.${PARENT_REF_TIME}.${STIME} 2>&1 || exit $?

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

./sno_bulk.sh > ./sno_bulk.log 2>&1 || exit $?

#-------------------------------------------------------------------------------

#rm -f ${SCPNAME}_job.sh
#rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
