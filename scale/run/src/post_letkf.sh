#!/bin/bash
#===============================================================================
#
#  Script to post-process the LETKF outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 4)); then
  cat >&2 << EOF

[post_letkf.sh]

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

letkfbaselen=9

#===============================================================================

mkdir -p $TMPOUT/${ATIME}/anal/${MEM}
for ifile in $(cd $TMPDIR ; ls anal.${MEM}.*.nc 2> /dev/null); do
  mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/${MEM}/init${ifile:$letkfbaselen}
done

if [ "$MEM" == 'mean' ]; then ###### using a variable for 'meanf', 'mean', 'sprd'
  mkdir -p $TMPOUT/${ATIME}/gues/mean
  for ifile in $(cd $TMPDIR ; ls gues.mean.*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/mean/init${ifile:$letkfbaselen}
  done
  mkdir -p $TMPOUT/${ATIME}/anal/mean
  for ifile in $(cd $TMPDIR ; ls anal.mean.*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/mean/init${ifile:$letkfbaselen}
  done
  mkdir -p $TMPOUT/${ATIME}/gues/sprd
  for ifile in $(cd $TMPDIR ; ls gues.sprd.*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/sprd/init${ifile:$letkfbaselen}
  done
  mkdir -p $TMPOUT/${ATIME}/anal/sprd
  for ifile in $(cd $TMPDIR ; ls anal.sprd.*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/sprd/init${ifile:$letkfbaselen}
  done
fi

if [ "$MEM" == '0001' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
  mkdir -p $TMPOUT/${ATIME}/log/letkf
  mv -f $TMPDIR/NOUT-000000 $TMPOUT/${ATIME}/log/letkf/NOUT-000000 ###### using a variable for 'NOUT-000000'
fi

#===============================================================================

exit 0
