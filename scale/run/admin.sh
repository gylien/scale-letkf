#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 7)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

SCPNAME="$1"; shift
STIME="$1"; shift
ETIME="$1"; shift
TIME_DT="$1"; shift
TIME_DT_DYN="$1"; shift
NNODES="$1"; shift
WTIME_L="$1"

if [ "$ETIME" = '-' ]; then
  ETIME="$STIME"
fi

CONFIG='realtime_r0051_d1'
PRESET='OFP'

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ]; then
  config_suffix='K'
  script_suffix='_K'
elif [ "$PRESET" = 'K_micro' ]; then
  config_suffix='K'
  script_suffix='_K_micro'
elif [ "$PRESET" = "OFP" ]; then
  config_suffix='ofp'
  script_suffix='_ofp'
else
  echo "[Error] Unsupported \$PRESET" >&2
  exit 1
fi

if [ "$SCPNAME" = 'cycle' ]; then
  DATA_BDY_WRF="ncepgfs_wrf_da"
  DATA_BDY_GRADS="ncepgfs_grads_da"
else
  DATA_BDY_WRF="ncepgfs_wrf"
  DATA_BDY_GRADS="ncepgfs_grads"
fi

#-------------------------------------------------------------------------------

rm -f config.main
rm -f config.${SCPNAME}
rm -f config.nml.*

cat config/${CONFIG}/config.main.${config_suffix} | \
    sed -e "s/<PRESET>/${PRESET}/g" | \
    sed -e "s/<DATA_BDY_WRF>/${DATA_BDY_WRF}/g" | \
    sed -e "s/<DATA_BDY_GRADS>/${DATA_BDY_GRADS}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g" \
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
fi

. config.main || exit $?
#. config.$SCPNAME || exit $?
#. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}${script_suffix}.sh > ${SCPNAME}${script_suffix}.log 2>&1 || exit $?

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ] || [ "$PRESET" = 'K_micro' ] || [ "$PRESET" = "OFP" ]; then
  jobname="${SCPNAME}_${SYSNAME}"
  jobid=$(grep 'pjsub Job' ${SCPNAME}${script_suffix}.log | cut -d ' ' -f6)
  logdir="$OUTDIR/exp/${jobid}_${SCPNAME}_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  jobinfo="$logdir/job.i"
#  stdout="${jobname}.o${jobid}"
#  stderr="${jobname}.e${jobid}"
#  jobinfo="${jobname}.i${jobid}"
fi

if [ ! -e "$stdout" ] || [ ! -e "$stderr" ]; then
  exit 101
fi
if [ -z "$(tail -n 1 $stderr | grep "Finish ${SCPNAME}.sh")" ]; then
  exit 102
fi

#-------------------------------------------------------------------------------

rm -f ${SCPNAME}_job.sh
rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_${SCPNAME}_${STIME} exp

#-------------------------------------------------------------------------------

exit 0
