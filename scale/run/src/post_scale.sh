#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main
. src/func_datetime.sh

if (($# < 8)); then
  cat >&2 << EOF

[post_scale.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK MEM_NP STIME MEM FCSTLEN TMPDIR LOG_OPT SCPCALL

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member (not used !!!)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM      Name of the ensemble member
  FCSTLEN  Forecast length (second)
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  SCPCALL  Called from which script? (fcst/cycle)

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
FCSTLEN="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
SCPCALL="$1"

ATIME=$(datetime $STIME $LCYCLE s)

restartbaselen=23  # 7 + 16

#===============================================================================


#echo
#echo '##### post #####'
#ls -lL $TMPDIR
#echo


#if [ "$SCPCALL" == 'fcst' ] && ((OUT_OPT <= 2)); then
#  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
#  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
#elif [ "$SCPCALL" == 'cycle' ]; then
#  mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
#  mv -f $TMPDIR/history*.nc $TMPOUT/${ATIME}/gues/${MEM}
#fi

#if [ "$SCPCALL" == 'fcst' ] && ((OUT_OPT <= 1)); then
#  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
#  for ifile in $(cd $TMPDIR ; ls restart*.nc); do
#    mv -f ${TMPDIR}/${ifile} $TMPOUT/${STIME}/fcst/${MEM}/init_$(datetime ${STIME} $FCSTLEN s)${ifile:$restartbaselen}
#  done
#elif [ "$SCPCALL" == 'cycle' ]; then
#  mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
#  for ifile in $(cd $TMPDIR ; ls restart*.nc); do
#    mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/${MEM}/init${ifile:$restartbaselen}
#  done
#fi

if [ "$SCPCALL" == 'fcst' ]; then
  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
  file_prefix=$(cd $TMPDIR ; ls restart*.nc | head -n 1) # pick up the first restart output. ###### TO DO: explicitly calculate the time string???
  for ifile in $(cd $TMPDIR ; ls ${file_prefix:0:$restartbaselen}*.nc); do
    mv -f ${TMPDIR}/${ifile} $TMPOUT/${STIME}/fcst/${MEM}/init_$(datetime ${STIME} $FCSTLEN s)${ifile:$restartbaselen}
  done
elif [ "$SCPCALL" == 'cycle' ]; then
  MEMtmp=$MEM
  if [ "$MEM" == 'mean' ]; then
    MEMtmp='meanf'
  fi
  mkdir -p $TMPOUT/${ATIME}/gues/${MEMtmp}
  mv -f $TMPDIR/history*.nc $TMPOUT/${ATIME}/gues/${MEMtmp}
  file_prefix=$(cd $TMPDIR ; ls restart*.nc | head -n 1) # pick up the first restart output. ###### TO DO: explicitly calculate the time string???
  for ifile in $(cd $TMPDIR ; ls ${file_prefix:0:$restartbaselen}*.nc); do
    mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/${MEMtmp}/init${ifile:$restartbaselen}
  done
fi

if ((LOG_OPT <= 3)); then
  mkdir -p $TMPOUT/${STIME}/log/scale
  if [ -f "$TMPDIR/LOG${SCALE_LOG_SFX}" ]; then
    mv -f $TMPDIR/LOG${SCALE_LOG_SFX} $TMPOUT/${STIME}/log/scale/${MEM}_LOG${SCALE_LOG_SFX}
  fi
fi

if ((LOG_OPT <= 1)); then
  mkdir -p $TMPOUT/${STIME}/log/scale
  if [ -f "$TMPDIR/run.conf" ]; then
    mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/scale/${MEM}_run.conf
  fi
fi

######
if ((MYRANK == 0)); then
  mkdir -p $TMPOUT/${STIME}/log/scale
  if [ -f "$TMPDIR/../latlon_domain_catalogue.txt" ]; then
    mv -f $TMPDIR/../latlon_domain_catalogue.txt $TMPOUT/${STIME}/log/scale/latlon_domain_catalogue.txt
  fi
fi
######

if ((MYRANK < MEM_NP)); then
  if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK)" ]; then
    mkdir -p $TMPOUT/${STIME}/log/scale
    mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK) $TMPOUT/${STIME}/log/scale
  fi
fi
#if [ "$MEM" == '0001' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
#  mkdir -p $TMPOUT/${STIME}/log/scale
#  for q in $(seq $MEM_NP); do
#    if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1)))" ]; then
#      mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1))) $TMPOUT/${STIME}/log/scale
#    fi
#  done
#fi

######
#mkdir -p $TMPOUT/${STIME}/log/scale/Fprofd
#mv -f $TMPDIR/../Fprofd/* $TMPOUT/${STIME}/log/scale/Fprofd
######

#===============================================================================

exit 0
