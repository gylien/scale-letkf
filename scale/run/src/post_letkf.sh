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

Usage: $0 MYRANK ATIME TMPDIR LOG_OPT

  MYRANK  My rank number (not used)
  ATIME   Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR  Temporary directory to run the program
  LOG_OPT

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"

#===============================================================================

if ((LOG_OPT <= 4 && MYRANK == 0)); then
  if [ -f "$TMPDIR/letkf.conf" ]; then
    mv -f $TMPDIR/letkf.conf $TMPOUT/${ATIME}/log/letkf/letkf.conf
    mv -f $TMPDIR/LOG.pe000000 $TMPOUT/${ATIME}/log/letkf/LOG.pe000000
  fi
fi

# -- Cloud dependent obs err --

if ((MYRANK == 0)); then

  if ls $TMPDIR/Him8_ERR_CA_B??.txt > /dev/null 2>&1 ; then
    cp ${TMPDIR}/Him8_ERR_CA_B??.txt $TMPOUT/${ATIME}/log/letkf/
  fi

  BB_LIST="07 08 09 10 11 12 13 14 15 16"
  for BB in ${BB_LIST} ; do
    CA_FILE1="${TMPDIR}/Him8_ERR_CA_B${BB}.dat"
    CA_FILE2="${TMPDAT}/Him8/Him8_ERR_CA_B${BB}_${ATIME}.dat"
    CA_FILE3="${TMPOUT}/${ATIME}/log/letkf/Him8_ERR_CA_B${BB}_${ATIME}.dat"
    if [ -e ${CA_FILE1} ] ; then
      cp $CA_FILE1 $CA_FILE2
      cp $CA_FILE1 $CA_FILE3
    fi
  done

fi

#===============================================================================

exit 0
