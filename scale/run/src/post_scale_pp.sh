#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[post_scale_pp.sh] Post-process the scale-les_init outputs.

Usage: $0 MYRANK MEM_NP STIME MEM TMPDIR LOG_OPT [SCPCALL]

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
SCPCALL="${1:-cycle}"

#===============================================================================

#if [ "$TOPO_FORMAT" != 'prep' ]; then
#  mkdir -p $TMPOUT/${STIME}/topo
#  mv -f $TMPDIR/topo*.nc $TMPOUT/${STIME}/topo
#fi

#if [ "$LANDUSE_FORMAT" != 'prep' ]; then
#  mkdir -p $TMPOUT/${STIME}/landuse
#  mv -f $TMPDIR/landuse*.nc $TMPOUT/${STIME}/landuse
#fi

#if ((LOG_OPT <= 2)); then
#  mkdir -p $TMPOUT/${STIME}/log/scale_pp
#  if [ -f "$TMPDIR/pp_LOG${SCALE_LOG_SFX}" ]; then
#    mv -f $TMPDIR/pp_LOG${SCALE_LOG_SFX} $TMPOUT/${STIME}/log/scale_pp/${MEM}_pp_LOG${SCALE_LOG_SFX}
#  fi
#fi

if [ "$SCPCALL" = 'fcst' ]; then
  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/pp.conf" ]; then
      mv -f $TMPDIR/pp.conf $TMPOUT/${STIME}/log/scale_pp/${MEM}_fcst_pp.conf
    fi
  fi
elif [ "$SCPCALL" = 'cycle' ]; then
  if ((LOG_OPT <= 4)); then
    if [ -f "$TMPDIR/pp.conf" ]; then
      mv -f $TMPDIR/pp.conf $TMPOUT/${STIME}/log/scale_pp/${MEM}_pp.conf
    fi
  fi
fi

#if ((MYRANK < MEM_NP)); then
#  if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK)" ]; then
#    mkdir -p $TMPOUT/${STIME}/log/scale_pp
#    mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $MYRANK) $TMPOUT/${STIME}/log/scale_pp
#  fi
#fi
#if [ "$MEM" == '0001' ] || [ "$MEM" == 'mean' ] && ((LOG_OPT <= 4)); then ###### using a variable for '0001'
#  mkdir -p $TMPOUT/${STIME}/log/scale_pp
#  for q in $(seq $MEM_NP); do
#    if [ -e "$TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1)))" ]; then
#      mv -f $TMPDIR/../NOUT-$(printf $PROCESS_FMT $((q-1))) $TMPOUT/${STIME}/log/scale_pp
#    fi
#  done
#fi

#===============================================================================

exit 0
