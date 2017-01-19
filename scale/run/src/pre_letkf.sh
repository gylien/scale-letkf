#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 4)); then
  cat >&2 << EOF

[pre_letkf.sh]

Usage: $0 MYRANK ATIME MEM OUT_OPT [ADAPTINFL RTPS_INFL_OUT NOBS_OUT]

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  OUT_OPT
  ADAPTINFL
  RTPS_INFL_OUT
  NOBS_OUT

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
OUT_OPT="$1"; shift
ADAPTINFL="${1:-0}"; shift
RTPS_INFL_OUT="${1:-0}"; shift
NOBS_OUT="${1:-0}"

#===============================================================================

if [ "$MEM" == 'mean' ]; then ###### using a variable for 'meanf', 'mean', 'sprd'
  for ifile in $(cd $TMPOUT/${ATIME}/gues/meanf ; ls init*.nc 2> /dev/null); do
    if ((ADAPTINFL == 1)) && [ ! -s "$TMPOUT/${ATIME}/diag/infl" ]; then
      mkdir -p $TMPOUT/${ATIME}/diag/infl
      cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/diag/infl
    fi
  done

  if ((ENABLE_PARAM_USER != 1)); then
    for ifile in $(cd $TMPOUT/${ATIME}/gues/meanf ; ls init*.nc 2> /dev/null); do
      mkdir -p $TMPOUT/${ATIME}/gues/mean
      cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/gues/mean
      mkdir -p $TMPOUT/${ATIME}/anal/mean
      cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/anal/mean
      mkdir -p $TMPOUT/${ATIME}/gues/sprd
      cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/gues/sprd
      mkdir -p $TMPOUT/${ATIME}/anal/sprd
      cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/anal/sprd
      if ((RTPS_INFL_OUT == 1)); then
        mkdir -p $TMPOUT/${ATIME}/diag/rtps
        cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/diag/rtps
      fi
      if ((NOBS_OUT == 1)); then
        mkdir -p $TMPOUT/${ATIME}/diag/nobs
        cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/diag/nobs
      fi
    done
  fi
else
  mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}

  if ((OUT_OPT <= 3 && ENABLE_PARAM_USER != 1)); then
    for ifile in $(cd $TMPOUT/${ATIME}/anal/${MEM} ; ls init*.nc 2> /dev/null); do
      mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
      cp -f $TMPOUT/${ATIME}/anal/${MEM}/${ifile} $TMPOUT/${ATIME}/gues/${MEM}
    done
  fi
fi

#===============================================================================

exit 0
