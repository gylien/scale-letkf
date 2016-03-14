#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 7)); then
  cat >&2 << EOF

[post_scale_init.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK MEM_NP STIME MKINIT MEM TMPDIR LOG_OPT [SCPCALL]

  MYRANK   My rank number (not used)
  MEM_NP  Number of processes per member
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MKINIT   Make initial condition as well?
            0: No
            1: Yes
  MEM      Name of the ensemble member
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP=$1; shift
STIME="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
SCPCALL="${1:-cycle}"

initbaselen=20

#===============================================================================

#mkdir -p $TMPOUT/${STIME}/bdy/${MEM}
#mv -f $TMPDIR/boundary*.nc $TMPOUT/${STIME}/bdy/${MEM}

if ((MKINIT == 1)); then
  mkdir -p $TMPOUT/${STIME}/anal/${MEM}
  for ifile in $(cd $TMPDIR ; ls init_*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/anal/${MEM}/init${ifile:$initbaselen}
  done
elif ((OCEAN_INPUT == 1 && OCEAN_FORMAT == 99)); then
  mkdir -p $TMPOUT/${STIME}/anal/${MEM}
  for ifile in $(cd $TMPDIR ; ls init_*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/anal/${MEM}/init_ocean${ifile:$initbaselen}
  done
fi

#if ((LOG_OPT <= 2)); then
#  mkdir -p $TMPOUT/${STIME}/log/scale_init
#  if [ -f "$TMPDIR/init_LOG${SCALE_LOG_SFX}" ]; then
#    mv -f $TMPDIR/init_LOG${SCALE_LOG_SFX} $TMPOUT/${STIME}/log/scale_init/${MEM}_init_LOG${SCALE_LOG_SFX}
#  fi
#fi

if [ "$SCPCALL" = 'fcst' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/init.conf" ]; then
      mv -f $TMPDIR/init.conf $TMPOUT/${STIME}/log/scale_init/${MEM}_fcst_init.conf
    fi
  fi
elif [ "$SCPCALL" = 'cycle' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/init.conf" ]; then
      mv -f $TMPDIR/init.conf $TMPOUT/${STIME}/log/scale_init/${MEM}_init.conf
    fi
  fi
fi

#if ((MYRANK < MEM_NP)); then
#  if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK)" ]; then
#    mkdir -p $TMPOUT/${STIME}/log/scale_init
#    mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK) $TMPOUT/${STIME}/log/scale_init
#  fi
#fi
#if [ "$MEM" == '0001' ] || [ "$MEM" == 'mean' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
#  mkdir -p $TMPOUT/${STIME}/log/scale_init
#  for q in $(seq $MEM_NP); do
#    if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1)))" ]; then
#      mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1))) $TMPOUT/${STIME}/log/scale_init
#    fi
#  done
#fi

#===============================================================================

exit 0
