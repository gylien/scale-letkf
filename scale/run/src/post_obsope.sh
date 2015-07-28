#!/bin/bash
#===============================================================================
#
#  Script to post-process the obsope outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[post_obsope.sh]

Usage: $0 MYRANK MEM_NP ATIME MEM TMPDIR LOG_OPT

  MYRANK  My rank number
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

#===============================================================================

mkdir -p $TMPOUT/${ATIME}/obsgues/${MEM}
for ifile in $(cd $TMPDIR ; ls obsda.${MEM}.*.dat 2> /dev/null); do
  mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/obsgues/${MEM}/${ifile}
done

if ((MYRANK < MEM_NP)); then
  if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK)" ]; then
    mkdir -p $TMPOUT/${ATIME}/log/obsope
    mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK) $TMPOUT/${ATIME}/log/obsope
  fi
fi
#if [ "$MEM" == '0001' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
#  mkdir -p $TMPOUT/${ATIME}/log/obsope
#  for q in $(seq $MEM_NP); do
#    if [ -e "$TMPDIR/NOUT-$(printf $PROCESS_FMT $((q-1)))" ]; then
#      mv -f $TMPDIR/NOUT-$(printf $PROCESS_FMT $((q-1))) $TMPOUT/${ATIME}/log/obsope
#    fi
#  done
#fi

#if [ "$MEM" == 'mean' ] && ((LOG_OPT <= 4)); then ###### using a variable for 'mean'
#  mkdir -p $TMPOUT/${ATIME}/log/obsope
#  for q in $(seq $MEM_NP); do
#    mv -f $TMPDIR/NOUT-$(printf $PROCESS_FMT $((MEMBER*MEM_NP+q-1))) $TMPOUT/${ATIME}/log/obsope # m=MEMBER+1 (mmean is not declared in this script)
#  done
#fi

#===============================================================================

exit 0
