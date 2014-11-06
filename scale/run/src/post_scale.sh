#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all
. src/func_datetime.sh

if (($# < 5)); then
  cat >&2 << EOF

[post_scale.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK MEM_NP STIME MEM FCSTLEN TMPDIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member (not used !!!)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM      Name of the ensemble member
  FCSTLEN  Forecast length (second)
  TMPDIR   Temporary directory to run the model

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
FCSTLEN="$1"; shift
TMPDIR="$1"

restartbaselen=23  # 7 + 16

#===============================================================================


#echo
#echo '##### post #####'
#ls -lL $TMPDIR
#echo


if ((OUT_OPT <= 2)); then
  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
#  for ifile in $(ls $TMPDIR/history*.nc); do
#    ifilebase=$(basename $ifile)
#    mv -f $ifile $TMPOUT/${STIME}/fcst/${MEM}/out${ifilebase/history/out}
#  done
  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
fi

if ((OUT_OPT <= 1)); then
  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
  for ifile in $(ls $TMPDIR/restart*.nc); do
    ifilebase=$(basename $ifile)
    mv -f $ifile $TMPOUT/${STIME}/fcst/${MEM}/init_$(datetime ${STIME} $FCSTLEN s)${ifilebase:$restartbaselen}
  done
fi

if ((LOG_OPT <= 3)); then
  mkdir -p $TMPOUT/${STIME}/log/scale
  if [ -f "$TMPDIR/LOG${SCALE_LOG_SFX}" ]; then
    mv -f $TMPDIR/LOG${SCALE_LOG_SFX} $TMPOUT/${STIME}/log/scale/${MEM}_LOG${SCALE_LOG_SFX}
  fi
fi

#===============================================================================

exit 0
