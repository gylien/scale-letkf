#!/bin/bash
#===============================================================================
#
#  Script to post-process the obsope outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 4)); then
  cat >&2 << EOF

[post_obsope.sh]

Usage: $0 MYRANK ATIME MEM TMPDIR

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  TMPDIR  Temporary directory to run the program

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"

#===============================================================================

mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}
for ifile in $(cd $TMPDIR ; ls obsda.${MEM}.*.dat 2> /dev/null); do
  mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/obsgues/${MEM}/${ifile}
done

if [ "$MEM" == '0001' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
  mkdir -p $TMPOUT/${ATIME}/log/obsope
  mv -f $TMPDIR/NOUT-000000 $TMPOUT/${ATIME}/log/obsope/NOUT-000000 ###### using a variable for 'NOUT-000000'
fi

#===============================================================================

exit 0
