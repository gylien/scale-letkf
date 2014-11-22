#!/bin/bash
#===============================================================================
#
#  Script to post-process the obsope outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all
. src/func_datetime.sh

if (($# < 3)); then
  cat >&2 << EOF

[post_obsope.sh] Post-process the obsope outputs.

Usage: $0 MYRANK STIME TMPDIR

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  TMPDIR   Temporary directory to run the model

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
TMPDIR="$1"

ATIME=$(datetime $STIME $LCYCLE s)

#===============================================================================

mkdir -p $TMPOUT/${ATIME}/obsgues
if [ -n "$(ls $TMPDIR/obsval*.dat 2> /dev/null)" ]; then  ### to be improved!!
  mv -f $TMPDIR/obsval*.dat $TMPOUT/${ATIME}/obsgues      #
fi                                                        #

if ((LOG_OPT <= 4)); then
  mkdir -p $TMPOUT/${ATIME}/log/obsope
  if [ -f "$TMPDIR/NOUT-000000" ]; then  # using a variable for '000000'
    mv -f $TMPDIR/NOUT-000000 $TMPOUT/${ATIME}/log/obsope/NOUT-000000
  fi
fi

#===============================================================================

exit 0
