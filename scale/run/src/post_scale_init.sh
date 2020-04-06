#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[post_scale_init.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK STIME MKINIT MEM TMPDIR LOG_OPT [SCPCALL]

  MYRANK   My rank number (not used)
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
STIME="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
SCPCALL="${1:-cycle}"

initbaselen=24

#===============================================================================

if ((MKINIT == 1)); then
  if ((PNETCDF == 1)); then
    mkdir -p $TMPOUT/${STIME}/anal
    ifile="$(cd $TMPDIR ; ls init_*.nc 2> /dev/null)"
    if [ -e "$TMPDIR/${ifile}" ]; then
      mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/anal/${MEM}.init.nc
    fi
  else
    mkdir -p $TMPOUT/${STIME}/anal/${MEM}
    for ifile in $(cd $TMPDIR ; ls init_*.nc 2> /dev/null); do
      mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/anal/${MEM}/init${ifile:$initbaselen}
    done
  fi
elif ((USE_INIT_FROM_BDY == 1)); then
  if ((PNETCDF == 1)); then
    mkdir -p $TMPOUT/${STIME}/anal
    ifile="$(cd $TMPDIR ; ls init_*.nc 2> /dev/null)"
    if [ -e "$TMPDIR/${ifile}" ]; then
      mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/bdy/${MEM}.init_bdy.nc
    fi
  else
    mkdir -p $TMPOUT/${STIME}/anal/${MEM}
    for ifile in $(cd $TMPDIR ; ls init_*.nc 2> /dev/null); do
      #mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/bdy/${MEM}/init_bdy${ifile:$initbaselen}
      mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/bdy/${MEM}/
    done
  fi
fi

if [ "$SCPCALL" = 'cycle' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/init.conf" ]; then
      mv -f $TMPDIR/init.conf $TMPOUT/${STIME}/log/scale_init/${MEM}_init.conf
    fi
    if [ -f "$TMPDIR/gradsbdy.conf" ]; then
      mv -f $TMPDIR/gradsbdy.conf $TMPOUT/${STIME}/log/scale_init/${MEM}_gradsbdy.conf
    fi
  fi
elif [ "$SCPCALL" = 'fcst' ]; then
  if ((LOG_OPT <= 2)); then
    if [ -f "$TMPDIR/init.conf" ]; then
      mv -f $TMPDIR/init.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale_init/${MEM}_init.conf
    fi
    if [ -f "$TMPDIR/gradsbdy.conf" ]; then
      mv -f $TMPDIR/gradsbdy.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale_init/${MEM}_gradsbdy.conf
    fi
  fi
fi

#===============================================================================

exit 0
