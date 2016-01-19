#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each member.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[pre_letkf.sh]

Usage: $0 MYRANK TOPO ATIME MEM TMPDIR

  MYRANK  My rank number (not used)
  TOPO    Basename of SCALE topography files
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  MEMSEQ
  TMPDIR  Temporary directory to run the program

EOF
  exit 1
fi

MYRANK="$1"; shift
TOPO="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
MEMSEQ="$1"; shift
TMPDIR="$1"

#historybaselen=7
obsdabaselen=10
initbaselen=4

#topodir=$(dirname $TOPO)
#topobase=$(basename $TOPO)
#topobaselen=${#topobase}

#===============================================================================

if [ "$MEM" == 'mean' ]; then ###### using a variable for 'meanf', 'mean', 'sprd'
#if [ -d "$TMPOUT/${ATIME}/gues/meanf" ]; then  # required....
  for ifile in $(cd $TMPOUT/${ATIME}/gues/meanf ; ls init*.nc 2> /dev/null); do
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPDIR/gues.mean${ifile:$initbaselen}
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPDIR/anal.mean${ifile:$initbaselen}
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPDIR/gues.sprd${ifile:$initbaselen}
    cp -f $TMPOUT/${ATIME}/gues/meanf/${ifile} $TMPDIR/anal.sprd${ifile:$initbaselen}
  done
#fi

  ln -fs ${TOPO}*.nc $TMPDIR
#  for ifile in $(cd $topodir ; ls ${topobase}*.nc 2> /dev/null); do
#    ln -fs $topodir/${ifile} $TMPDIR/topo${ifile:$topobaselen}
#  done
else
  if [ -d "$TMPOUT/${ATIME}/obsgues/${MEM}" ]; then
    for ifile in $(cd $TMPOUT/${ATIME}/obsgues/${MEM} ; ls obsda.${MEM}.*.dat 2> /dev/null); do
#      ln -fs $TMPOUT/${ATIME}/obsgues/${MEM}/${ifile} $TMPDIR/${ifile}
      ln -fs $TMPOUT/${ATIME}/obsgues/${MEM}/${ifile} $TMPDIR/obsda.${MEMSEQ}${ifile:$obsdabaselen}
    done
  fi

  if [ -d "$TMPOUT/${ATIME}/gues/${MEM}" ]; then
#    for ifile in $(cd $TMPOUT/${ATIME}/gues/${MEM} ; ls history*.nc 2> /dev/null); do
#      ln -fs $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPDIR/hist.${MEMSEQ}${ifile:$historybaselen}
#    done

    for ifile in $(cd $TMPOUT/${ATIME}/gues/${MEM} ; ls init*.nc 2> /dev/null); do
      ln -fs $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPDIR/gues.${MEMSEQ}${ifile:$initbaselen}
      cp -f $TMPOUT/${ATIME}/gues/${MEM}/${ifile} $TMPDIR/anal.${MEMSEQ}${ifile:$initbaselen}
    done
  fi
fi

#===============================================================================

exit 0
