#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 6)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

SCPNAME="cycle"

DX="$1"; shift
STIME="$1"; shift
NCYCLE="$1"; shift
WTIME_L="$1"; shift
NMEM="$1"; shift 
DACYCLE_RUN_FCST_TIME="$1"

intv_sec=`expr \( $NCYCLE - 1 \) \* 30`
STIME_in="${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}:${STIME:10:2}:${STIME:12:2}"
ETIME=`date -d "${intv_sec} second ${STIME_in}" +'%Y%m%d%H%M%S'`

MAX_DACYCLE_RUN_FCST=$NCYCLE

NNODES=`expr \( $NMEM + 2 + $MAX_DACYCLE_RUN_FCST \) / 1` ### 1km / 64domain


#-------------------------------------------------------------------------------

CONFIG="D4_${DX}_mercator"
if [ ! -d config/${CONFIG} ];then
 echo "EXPTYPE "$EXPTYPE" not supported !"
 exit
fi

cp config/${CONFIG}/config.* .

cat config.main.ofp | \
    sed -e "s/<MEMBER>/${NMEM}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" \
    > config.main

cat config/$CONFIG/config.cycle | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<MAX_DACYCLE_RUN_FCST>/${MAX_DACYCLE_RUN_FCST}/g" | \
    sed -e "s/<DACYCLE_RUN_FCST_TIME>/${DACYCLE_RUN_FCST_TIME}/g"  \
    > config.cycle


. config.main || exit $?

#-------------------------------------------------------------------------------

./cycle_ofp.sh > cycle_ofp.log 2>&1 || exit $?

#-------------------------------------------------------------------------------

  jobname="${SCPNAME}_${SYSNAME}"
  jobid=$(grep 'pjsub Job' cycle_ofp.log | cut -d ' ' -f6)
  logdir="$OUTDIR/exp/${jobid}_${SCPNAME}_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  jobinfo="$logdir/job.i"

#-------------------------------------------------------------------------------

rm -f ${SCPNAME}_job.sh
rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
