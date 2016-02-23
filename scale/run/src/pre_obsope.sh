#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[pre_obsope.sh]

Usage: $0 MYRANK ATIME MEM TMPDIR

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  MEMSEQ
  TMPDIR  Temporary directory to run the program

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
MEMSEQ="$1"; shift
TMPDIR="$1"

historybaselen=7
#initbaselen=4

#===============================================================================

#if [ -d "$TMPOUT/${ATIME}/gues/${MEM}" ]; then
#  for ifile in $(cd $TMPOUT/${ATIME}/gues/${MEM} ; ls history*.nc 2> /dev/null); do
#    ln -fs $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPDIR/hist.${MEMSEQ}${ifile:$historybaselen}
#  done
#  for ifile in $(cd $TMPOUT/${ATIME}/gues/${MEM} ; ls init*.nc 2> /dev/null); do
#    ln -fs $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPDIR/init.${MEMSEQ}${ifile:$initbaselen}
#  done
#fi

if [ "$MEM" != 'mean' ]; then
  mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}
fi

#===============================================================================

exit 0
