#!/bin/bash
#===============================================================================

function print_msg () {
#  rm -f config.main
#  rm -f config.${SCPNAME}
#  rm -f config.nml.ensmodel
#  rm -f config.nml.letkf
#  rm -f config.nml.obsope
#  rm -f config.nml.scale
#  rm -f config.nml.scale_pp
#  rm -f config.nml.scale_init

#  rm -f ${SCPNAME}_job.sh
##  rm -f ${jobname}.o${jobid}
##  rm -f ${jobname}.e${jobid}
##  rm -f ${jobname}.s${jobid}
##  rm -f ${jobname}.i${jobid}

#  rm -f test.letkf
#  rm -f test.letkf.ref
#  rm -f test.log

  echo
  echo "========================================"
  echo "JOB STATUS:   $status_job"
  echo "RESULT:       $status_result"
  echo "========================================"
}
trap print_msg EXIT

status_job='Not run'
status_result='No result'

#-------------------------------------------------------------------------------

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if (($# < 5)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

CONFIG="$1"; shift
PRESET="$1"; shift
MPI_TYPE="$1"; shift
SCPNAME="$1"; shift
WTIME_L="$1"

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ]; then
  config_suffix='K'
  script_suffix='K'
elif [ "$PRESET" = 'K_micro' ]; then
  config_suffix='K'
  script_suffix='K_micro'
elif [ "$PRESET" = 'Linux_torque' ]; then
  config_suffix='hakushu'
  script_suffix='torque'
else
  echo "[Error] Unsupported \$PRESET" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

cat config/${CONFIG}/config.main.${config_suffix} | \
    sed -e "s/<PRESET>/${PRESET}/g" | \
    sed -e "s/<MPI_TYPE>/${MPI_TYPE}/g" \
    > config.main

cat config/${CONFIG}/config.${SCPNAME} | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" \
    > config.${SCPNAME}

ln -fs config/${CONFIG}/config.nml.ensmodel .
ln -fs config/${CONFIG}/config.nml.letkf .
ln -fs config/${CONFIG}/config.nml.obsope .
ln -fs config/${CONFIG}/config.nml.scale .
ln -fs config/${CONFIG}/config.nml.scale_pp .
ln -fs config/${CONFIG}/config.nml.scale_init .

. config.main || exit $?
. config.$SCPNAME || exit $?
. src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

./${SCPNAME}_${script_suffix}.sh 2>&1 | tee test.log
rc=${PIPESTATUS[0]}
status_job='Done with errors'
((rc == 0)) || exit $rc

jobname="${SCPNAME}_${CONFIG}"
jobid=$(grep 'pjsub Job' test.log | cut -d ' ' -f6)

#echo "jobname = ${jobname}"
#echo "jobid   = ${jobid}"

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ] || [ "$PRESET" = 'K_micro' ]; then
  logdir="$OUTDIR/exp/${jobid}_${SCPNAME}_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  letkflogname="NOUT.0"
#elif [ "$PRESET" = 'Linux_torque' ]; then
#
fi

#echo $stdout
#echo $stderr

if [ ! -e "$stdout" ] || [ ! -e "$stderr" ]; then
  exit 101
fi

echo
echo "STANDARD ERROR MESSAGE:"
echo "========================================"
cat $stderr
echo "========================================"

if [ -z "$(tail -n 1 $stderr | grep "Finish ${SCPNAME}.sh")" ]; then
  exit 102
fi

status_job='Done'
status_result='Incorrect'

eatime=$(datetime $ETIME $LCYCLE s)
str_search="OBSERVATIONAL DEPARTURE STATISTICS (GLOBAL)"
sed -n "/${str_search}/,+7p" $OUTDIR/${eatime}/log/letkf/${letkflogname} | grep -v "$str_search" > test.letkf
sed -n "/${str_search}/,+7p" $OUTDIR/results/${eatime}/log/letkf/${letkflogname} | grep -v "$str_search" > test.letkf.ref

if [ -n "$(diff test.letkf test.letkf.ref)" ]; then
  exit 201
fi

status_result='Correct'

#-------------------------------------------------------------------------------

exit 0
