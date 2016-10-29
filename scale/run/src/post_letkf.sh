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
    CA_FILE_A="${TMPDIR}/Him8_ERR_CA_A_B${BB}.dat"
    CA_FILE_B="${TMPDIR}/Him8_ERR_CA_B_B${BB}.dat"
    CA_FILE2_A="${TMPDAT}/Him8/Him8_ERR_CA_A_B${BB}_${ATIME}.dat"
    CA_FILE3_A="${TMPOUT}/${ATIME}/log/letkf/Him8_ERR_CA_A_B${BB}_${ATIME}.dat"
    CA_FILE2_B="${TMPDAT}/Him8/Him8_ERR_CA_B_B${BB}_${ATIME}.dat"
    CA_FILE3_B="${TMPOUT}/${ATIME}/log/letkf/Him8_ERR_CA_B_B${BB}_${ATIME}.dat"
    if [ -e ${CA_FILE_A} ] && [ -e ${CA_FILE_B} ] ; then
      cp $CA_FILE_A $CA_FILE2_A
      cp $CA_FILE_A $CA_FILE3_A
      cp $CA_FILE_B $CA_FILE2_B
      cp $CA_FILE_B $CA_FILE3_B
    fi
  done

fi

#===============================================================================

exit 0
