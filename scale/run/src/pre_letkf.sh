#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[pre_letkf.sh]

Usage: $0 MYRANK ATIME MEM OUT_OPT OBSOUT_OPT [ADAPTINFL SPRD_OUT RTPS_INFL_OUT NOBS_OUT]

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  OUT_OPT
  OBSOUT_OPT
  ADAPTINFL
  SPRD_OUT
  RTPS_INFL_OUT
  NOBS_OUT

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
OUT_OPT="$1"; shift
OBSOUT_OPT="$1"; shift
ADAPTINFL="${1:-0}"; shift
SPRD_OUT="${1:-1}"; shift
RTPS_INFL_OUT="${1:-0}"; shift
NOBS_OUT="${1:-0}"

#===============================================================================

if [ "$MEM" = 'mean' ]; then ###### using a variable for 'mean', 'sprd'
  if ((ADAPTINFL == 1)) && [ ! -s "$TMPOUT/${ATIME}/diag/infl" ]; then
    for ifile in $(cd $TMPOUT/${ATIME}/anal/mean ; ls init*.nc 2> /dev/null); do
      mkdir -p $TMPOUT/${ATIME}/diag/infl
      cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/diag/infl
    done
  fi

  if ((ENABLE_PARAM_USER != 1)); then
    for ifile in $(cd $TMPOUT/${ATIME}/anal/mean ; ls init*.nc 2> /dev/null); do
      mkdir -p $TMPOUT/${ATIME}/gues/mean
      cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/gues/mean
      if ((SPRD_OUT == 1)); then
        mkdir -p $TMPOUT/${ATIME}/anal/sprd
        cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/anal/sprd
        mkdir -p $TMPOUT/${ATIME}/gues/sprd
        cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/gues/sprd
      fi
      if ((RTPS_INFL_OUT == 1)); then
        mkdir -p $TMPOUT/${ATIME}/diag/rtps
        cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/diag/rtps
      fi
      if ((NOBS_OUT == 1)); then
        mkdir -p $TMPOUT/${ATIME}/diag/nobs
        cp -f $TMPOUT/${ATIME}/anal/mean/${ifile} $TMPOUT/${ATIME}/diag/nobs
      fi
    done
  fi
else
  if ((ENABLE_PARAM_USER != 1)); then
    if ((OUT_OPT <= 3)) || [ "$MEM" = 'mdet' ]; then
      for ifile in $(cd $TMPOUT/${ATIME}/anal/${MEM} ; ls init*.nc 2> /dev/null); do
        mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
        cp -f $TMPOUT/${ATIME}/anal/${MEM}/${ifile} $TMPOUT/${ATIME}/gues/${MEM}
      done
    fi
  fi
fi

if ((OBSOUT_OPT <= 2)); then
  mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}
fi

#===============================================================================

exit 0
