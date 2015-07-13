#!/bin/bash
#===============================================================================
#
#  Script to post-process the LETKF outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[post_letkf.sh]

Usage: $0 MYRANK MEM_NP ATIME MEM TMPDIR LOG_OPT

  MYRANK  My rank number (not used)
  MEM_NP  Number of processes per member
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  MEM     Name of the ensemble member
  TMPDIR  Temporary directory to run the program
  LOG_OPT

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
ATIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"

letkfbaselen=9

#===============================================================================

mkdir -p $TMPOUT/${ATIME}/anal/${MEM}

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
else
  for ifile in $(cd $TMPDIR ; ls anal.${MEM}.*.nc 2> /dev/null); do
    mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/${MEM}/init${ifile:$letkfbaselen}
  done
fi

if ((MYRANK < MEM_NP)); then
  mkdir -p $TMPOUT/${ATIME}/log/letkf
  if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK)" ]; then
    mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK) $TMPOUT/${ATIME}/log/letkf
  fi
fi
#if [ "$MEM" == '0001' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
#  mkdir -p $TMPOUT/${ATIME}/log/letkf
#  for q in $(seq $MEM_NP); do
#    if [ -e "$TMPDIR/NOUT-$(printf $PROCESS_FMT $((q-1)))" ]; then
#      mv -f $TMPDIR/NOUT-$(printf $PROCESS_FMT $((q-1))) $TMPOUT/${ATIME}/log/letkf
#    fi
#  done
#fi

#if [ "$MEM" == 'mean' ] && ((LOG_OPT <= 4)); then ###### using a variable for 'mean'
#  mkdir -p $TMPOUT/${ATIME}/log/letkf
#  for q in $(seq $MEM_NP); do
#    mv -f $TMPDIR/NOUT-$(printf $PROCESS_FMT $((MEMBER*MEM_NP+q-1))) $TMPOUT/${ATIME}/log/letkf # m=MEMBER+1 (mmean is not declared in this script)
#  done
#fi

#===============================================================================

exit 0
