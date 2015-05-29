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

CONFDIR="$OUTDIR/${STIME}/${SCPNAME}_conf"

#-------------------------------------------------------------------------------

if [ "$SCPNAME" = 'cycle' ]; then
  STIME_DIR="${STIME}_da"
else
  STIME_DIR="${STIME}"
fi
cat config/EastAsia_18km_48p/config.main.K | \
    sed -e "s/<STIME>/${STIME_DIR}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" \
    > config.main

cat config/EastAsia_18km_48p/config.${SCPNAME} | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" \
    > config.${SCPNAME}

cat config/EastAsia_18km_48p/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

mkdir -p $CONFDIR
cp config.* $CONFDIR

#-------------------------------------------------------------------------------

./${SCPNAME}_K.sh > ${SCPNAME}_K.log 2>&1
res=$? && ((res != 0)) && exit $res

jobname="${SCPNAME}_${SYSNAME}"
jobid=$(grep 'pjsub Job' ${SCPNAME}_K.log | cut -d ' ' -f6)

n=0
nmax=120
while [ ! -s "${jobname}.o${jobid}" ] || [ ! -s "${jobname}.e${jobid}" ] ||
      [ ! -s "${jobname}.s${jobid}" ] || [ ! -s "${jobname}.i${jobid}" ] && ((n < nmax)); do
  n=$((n+1))
  sleep 5s
done

#-------------------------------------------------------------------------------

res=0
if ((n >= nmax)); then
  res=101
elif [ -n "$(grep 'ERR.' ${jobname}.e${jobid})" ]; then
  res=102
elif [ -n "$(grep 'terminated' ${jobname}.e${jobid})" ]; then
  res=103
elif [ ! -s "${jobname}.s${jobid}" ]; then
  res=104
elif [ "$(tail -n 1 ${jobname}.s${jobid})" != "---(Stage-Out Error Information)---" ]; then
  res=105
fi

mv -f ${SCPNAME}_job.sh $CONFDIR/${jobid}.b
mv -f ${SCPNAME}_K.log $CONFDIR/${jobid}.l

mv -f ${jobname}.o${jobid} $CONFDIR/${jobid}.o
mv -f ${jobname}.e${jobid} $CONFDIR/${jobid}.e
mv -f ${jobname}.s${jobid} $CONFDIR/${jobid}.s
mv -f ${jobname}.i${jobid} $CONFDIR/${jobid}.i

mv -f $LOGDIR/${SCPNAME}_${STIME}.log $CONFDIR/${jobid}.lo
mv -f $LOGDIR/${SCPNAME}.err $CONFDIR/${jobid}.le

#-------------------------------------------------------------------------------

exit $res
