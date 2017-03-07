#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main
. src/func_datetime.sh

if (($# < 7)); then
  cat >&2 << EOF

[post_scale.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK STIME MEM FCSTLEN TMPDIR LOG_OPT OUT_OPT [SCPCALL DELETE_MEMBER SPRD_OUT RTPS_INFL_OUT NOBS_OUT]

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM      Name of the ensemble member
  FCSTLEN  Forecast length (second)
  TMPDIR   Temporary directory to run the model
  LOG_OPT
  OUT_OPT
  SCPCALL  Called from which script? (fcst/cycle)
  DELETE_MEMBER
  SPRD_OUT
  RTPS_INFL_OUT
  NOBS_OUT

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
FCSTLEN="$1"; shift
TMPDIR="$1"; shift
LOG_OPT="$1"; shift
OUT_OPT="$1"; shift
SCPCALL="${1:-cycle}"; shift
DELETE_MEMBER="${1:-0}"; shift
SPRD_OUT="${1:-1}"; shift
RTPS_INFL_OUT="${1:-0}"; shift
NOBS_OUT="${1:-0}"

#===============================================================================

if [ "$SCPCALL" = 'cycle' ]; then

  ATIME=$(datetime $STIME $LCYCLE s)

  if [ "$MEM" = 'mean' ]; then ###### using a variable for 'mean', 'mdet', 'sprd'

    mkdir -p $TMPOUT/${STIME}/hist/mean
    mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/hist/mean

    mkdir -p $TMPOUT/${ATIME}/anal/mean
    file_prefix="restart_${ATIME:0:8}-${ATIME:8:6}.000"
    restartbaselen=27
    for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
      mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/anal/mean/init${ifile:$restartbaselen}
    done

    if ((ENABLE_PARAM_USER == 1)); then
      restartbaselen=29

      mkdir -p $TMPOUT/${ATIME}/gues/mean
      file_prefix="restart_2_${ATIME:0:8}-${ATIME:8:6}.000"
      for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
        mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/mean/init${ifile:$restartbaselen}
      done

      icopy=2
      if ((SPRD_OUT == 1)); then
        icopy=$((icopy+1))
        mkdir -p $TMPOUT/${ATIME}/anal/sprd
        file_prefix="restart_${icopy}_${ATIME:0:8}-${ATIME:8:6}.000"
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/anal/sprd/init${ifile:$restartbaselen}
        done

        icopy=$((icopy+1))
        mkdir -p $TMPOUT/${ATIME}/gues/sprd
        file_prefix="restart_${icopy}_${ATIME:0:8}-${ATIME:8:6}.000"
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/sprd/init${ifile:$restartbaselen}
        done
      fi

      if ((RTPS_INFL_OUT == 1)); then
        icopy=$((icopy+1))
        mkdir -p $TMPOUT/${ATIME}/diag/rtps
        file_prefix="restart_${icopy}_${ATIME:0:8}-${ATIME:8:6}.000"
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/diag/rtps/init${ifile:$restartbaselen}
        done
      fi

      if ((NOBS_OUT == 1)); then
        icopy=$((icopy+1))
        mkdir -p $TMPOUT/${ATIME}/diag/nobs
        file_prefix="restart_${icopy}_${ATIME:0:8}-${ATIME:8:6}.000"
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/diag/nobs/init${ifile:$restartbaselen}
        done
      fi
    fi # [ ENABLE_PARAM_USER == 1 ]

  else

    mkdir -p $TMPOUT/${STIME}/hist/${MEM}
    mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/hist/${MEM}

    mkdir -p $TMPOUT/${ATIME}/anal/${MEM}
    file_prefix="restart_${ATIME:0:8}-${ATIME:8:6}.000"
    restartbaselen=27
    for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
      mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/anal/${MEM}/init${ifile:$restartbaselen}
    done

    if ((ENABLE_PARAM_USER == 1)); then
      if ((OUT_OPT <= 3)) || [ "$MEM" = 'mdet' ]; then
        mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
        file_prefix="restart_2_${ATIME:0:8}-${ATIME:8:6}.000"
        restartbaselen=29
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f ${TMPDIR}/${ifile} $TMPOUT/${ATIME}/gues/${MEM}/init${ifile:$restartbaselen}
        done
      fi
    fi

    if ((DELETE_MEMBER == 1)) && [ "$MEM" != 'mdet' ]; then
      if [ -d "$TMPOUT/${STIME}/anal/${MEM}" ]; then
        rm -f $TMPOUT/${STIME}/anal/${MEM}/*
      fi
    fi

  fi

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/scale/${MEM}_run.conf
    fi
  fi

  if ((MYRANK == 0)); then
    if [ -f "$TMPDIR/../latlon_domain_catalogue.txt" ]; then
      mv -f $TMPDIR/../latlon_domain_catalogue.txt $TMPOUT/${STIME}/log/scale/latlon_domain_catalogue.txt
    fi
  fi

elif [ "$SCPCALL" = 'fcst' ]; then

  FTIME=$(datetime $STIME $FCSTLEN s)

  mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
  mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
  if ((OUT_OPT <= 1)); then
    file_prefix="restart_${FTIME:0:8}-${FTIME:8:6}.000"
    restartbaselen=27
    for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
      mv -f ${TMPDIR}/${ifile} $TMPOUT/${STIME}/fcst/${MEM}/init_${FTIME}${ifile:$restartbaselen}
    done
  fi

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale/${MEM}_run.conf
    fi
  fi

  if ((MYRANK == 0)); then
    if [ -f "$TMPDIR/../latlon_domain_catalogue.txt" ]; then
      mv -f $TMPDIR/../latlon_domain_catalogue.txt $TMPOUT/${STIME}/log/${SCPCALL}_scale/latlon_domain_catalogue.txt
    fi
  fi

fi

#===============================================================================

exit 0
