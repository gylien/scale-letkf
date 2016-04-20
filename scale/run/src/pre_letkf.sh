#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 3)); then
  cat >&2 << EOF

[pre_letkf.sh]

Usage: $0 MYRANK ATIME MEM

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"

#===============================================================================

if [ "$MEM" == 'mean' ]; then ###### using a variable for 'meanf', 'mean', 'sprd'
#if [ -d "$TMPOUT/${ATIME}/gues/meanf" ]; then  # required....
  for ifile in $(cd $TMPOUT/${ATIME}/gues/meanf ; ls init*.nc 2> /dev/null); do
    mkdir -p $TMPOUT/${ATIME}/gues/mean
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/gues/mean
    mkdir -p $TMPOUT/${ATIME}/anal/mean
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/anal/mean
    mkdir -p $TMPOUT/${ATIME}/gues/sprd
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/gues/sprd
    mkdir -p $TMPOUT/${ATIME}/anal/sprd
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/anal/sprd



    mkdir -p $TMPOUT/${ATIME}/diag/infl
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/diag/infl
    mkdir -p $TMPOUT/${ATIME}/diag/nobs
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPOUT/${ATIME}/diag/nobs



  done
#fi
else
  for ifile in $(cd $TMPOUT/${ATIME}/gues/${MEM} ; ls init*.nc 2> /dev/null); do
    mkdir -p $TMPOUT/${ATIME}/anal/${MEM}
    cp -f $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPOUT/${ATIME}/anal/${MEM}
  done
fi

#===============================================================================

exit 0
