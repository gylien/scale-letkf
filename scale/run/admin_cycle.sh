#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

. config.main
(($? != 0)) && exit $?

#-------------------------------------------------------------------------------

if (($# < 3)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

STIME="$1"
TIME_DT="$2"
TIME_DT_DYN="$3"

CONFDIR="$OUTDIR/${STIME}/anal_conf"

#-------------------------------------------------------------------------------

cat config/EastAsia_18km_48p/config.main.K | \
    sed -e "s/<STIME>/${STIME}/g" \
    > config.main

cat config/EastAsia_18km_48p/config.cycle | \
    sed -e "s/<STIME>/${STIME}/g" \
    > config.cycle

cat config/EastAsia_18km_48p/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

mkdir -p $CONFDIR
cp config.* $CONFDIR

#-------------------------------------------------------------------------------

#./cycle_K.sh > cycle_K.log 2>&1

#jobname=$(head -n 3 cycle_job.sh | grep '#PJM -N' | cut -d ' ' -f3)
jobname="cycle_${SYSNAME}"
jobid=$(grep 'pjsub Job' cycle_K.log | cut -d ' ' -f6)

if [ ! -s "${jobname}.o${jobid}" ]; then
  exit 1
elif [ ! -s "${jobname}.e${jobid}" ]; then
  exit 2
elif [ ! -s "${jobname}.s${jobid}" ]; then
  exit 3
elif [ ! -s "${jobname}.i${jobid}" ]; then
  exit 4
elif [ -n "$(grep 'ERR.' ${jobname}.e${jobid})" ]; then
  exit 5
elif [ -n "$(grep 'terminated' ${jobname}.e${jobid})" ]; then
  exit 6
elif [ ! -s "${jobname}.s${jobid}" ]; then
  exit 7
elif [ "$(tail -n 1 ${jobname}.s${jobid})" != "---(Stage-Out Error Information)---" ]; then
  exit 8
fi

mv -f ${jobname}.o${jobid} $CONFDIR
mv -f ${jobname}.e${jobid} $CONFDIR
mv -f ${jobname}.s${jobid} $CONFDIR
mv -f ${jobname}.i${jobid} $CONFDIR

mv -f $LOGDIR/cycle_${STIME}.log $CONFDIR
mv -f $LOGDIR/cycle.err $CONFDIR

#
# remove unused file......
#


exit 0
