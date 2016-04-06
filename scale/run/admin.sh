#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

. config.main
res=$? && ((res != 0)) && exit $res

#. src/func_datetime.sh

#-------------------------------------------------------------------------------

if (($# < 5)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

SCPNAME="$1"
STIME="$2"
TIME_DT="$3"
TIME_DT_DYN="$4"
NNODES="$5"
WTIME_L="$6"

#ETIME=$(datetime "$STIME" ${LCYCLE} s)

CONFIG='realtime_v160405_d1'

#-------------------------------------------------------------------------------

if [ "$SCPNAME" = 'cycle' ]; then
  STIME_DIR="${STIME}_da"
else
  STIME_DIR="${STIME}"
fi
cat config/${CONFIG}/config.main.K | \
    sed -e "s/<STIME>/${STIME_DIR}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" \
    > config.main

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" \
    > config.${SCPNAME}

cat config/${CONFIG}/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

#-------------------------------------------------------------------------------

./${SCPNAME}_K.sh > ${SCPNAME}_K.log 2>&1
res=$? && ((res != 0)) && exit $res

jobname="${SCPNAME}_${SYSNAME}"
jobid=$(grep 'pjsub Job' ${SCPNAME}_K.log | cut -d ' ' -f6)

#-------------------------------------------------------------------------------

if [ ! -s "${jobname}.o${jobid}" ] || [ ! -s "${jobname}.e${jobid}" ] || \
   [ ! -s "${jobname}.i${jobid}" ] || [ ! -s "${jobname}.s${jobid}" ]; then
  exit 101
elif [ -n "$(grep 'ERR.' ${jobname}.e${jobid})" ]; then
  exit 102
elif [ -n "$(grep 'terminated' ${jobname}.e${jobid})" ]; then
  exit 103
#elif [ ! -s "${jobname}.s${jobid}" ]; then
#  exit 104
#elif [ "$(tail -n 1 ${jobname}.s${jobid})" != "---(Stage-Out Error Information)---" ]; then
#  exit 105
fi

rm -f ${SCPNAME}_job.sh
rm -f ${jobname}.o${jobid}
rm -f ${jobname}.e${jobid}
rm -f ${jobname}.s${jobid}
rm -f ${jobname}.i${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
