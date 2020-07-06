#!/bin/bash
#===============================================================================
#
#  Script to post-process the LETKF outputs.
#  December 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 4)); then
  cat >&2 << EOF

[post_letkf.sh]

Usage: $0 MYRANK STIME ATIME TMPDIR LOG_OPT OUT_OPT

  MYRANK  My rank number (not used)
  STIME
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR  Temporary directory to run the program
  LOG_OPT
  OUT_OPT

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
OUT_OPT="$1"

#===============================================================================

if ((LOG_OPT <= 4 && MYRANK == 0)); then
  if [ -f "$TMPDIR/letkf.conf" ]; then
    mv -f $TMPDIR/letkf.conf $TMPOUT/${ATIME}/log/letkf/letkf.conf
  fi
fi

if ((OUT_OPT >= 3 && MYRANK == 0)); then
  rm -rf ${OUTDIR}/${STIME}/hist/00??
  rm -rf ${OUTDIR}/${ATIME}/gues/0???
fi

if (( MYRANK == 0)); then
  rm -f $TMPDAT_OBS/*_${ATIME}.*
fi

#if [ -f "${TMPOUT}/vbc/Him8_vbca.dat" ] && [ ${MYRANK} -eq 0 ] ; then
#  mv ${TMPOUT}/vbc/Him8_vbca.dat ${TMPOUT}/vbc/Him8_vbca_${ATIME}.dat
#fi


#===============================================================================

exit 0
