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

    if ((PNETCDF == 1)); then
      mkdir -p $TMPOUT/${STIME}/hist
      mv -f $TMPDIR/history.nc $TMPOUT/${STIME}/hist/mean.history.nc
    else
      mkdir -p $TMPOUT/${STIME}/hist/mean
      mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/hist/mean
    fi

    if ((PNETCDF == 1)); then
      mkdir -p $TMPOUT/${ATIME}/anal
      ifile="restart_$(datetime_scale $ATIME).nc"
      if [ -e "$TMPDIR/${ifile}" ]; then
        mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/mean.init.nc
      fi
    else
      mkdir -p $TMPOUT/${ATIME}/anal/mean
      file_prefix="restart_$(datetime_scale $ATIME)"
      restartbaselen=27
      for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
        mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/mean/init${ifile:$restartbaselen}
      done
    fi

    if ((ENABLE_PARAM_USER == 1)); then
      if ((PNETCDF == 1)); then
        ifile="restart_2_$(datetime_scale $ATIME).nc"
        if [ -e "$TMPDIR/${ifile}" ]; then
          mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/mean.init.nc
        fi

        icopy=2
        if ((SPRD_OUT == 1)); then
          icopy=$((icopy+1))
          mkdir -p $TMPOUT/${ATIME}/anal
          ifile="restart_${icopy}_$(datetime_scale $ATIME).nc"
          if [ -e "$TMPDIR/${ifile}" ]; then
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/sprd.init.nc
          fi

          icopy=$((icopy+1))
          ifile="restart_${icopy}_$(datetime_scale $ATIME).nc"
          if [ -e "$TMPDIR/${ifile}" ]; then
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/sprd.init.nc
          fi
        fi
      else
        restartbaselen=29

        mkdir -p $TMPOUT/${ATIME}/gues/mean
        file_prefix="restart_2_$(datetime_scale $ATIME)"
        for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
          mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/mean/init${ifile:$restartbaselen}
        done

        icopy=2
        if ((SPRD_OUT == 1)); then
          icopy=$((icopy+1))
          mkdir -p $TMPOUT/${ATIME}/anal/sprd
          file_prefix="restart_${icopy}_$(datetime_scale $ATIME)"
          for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/sprd/init${ifile:$restartbaselen}
          done

          icopy=$((icopy+1))
          mkdir -p $TMPOUT/${ATIME}/gues/sprd
          file_prefix="restart_${icopy}_$(datetime_scale $ATIME)"
          for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/sprd/init${ifile:$restartbaselen}
          done
        fi
      fi

      if ((RTPS_INFL_OUT == 1)); then
        icopy=$((icopy+1))
        if ((PNETCDF == 1)); then
          mkdir -p $TMPOUT/${ATIME}/diag
          ifile="restart_${icopy}_$(datetime_scale $ATIME).nc"
          if [ -e "$TMPDIR/${ifile}" ]; then
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/diag/rtps.init.nc
          fi
        else
          mkdir -p $TMPOUT/${ATIME}/diag/rtps
          file_prefix="restart_${icopy}_$(datetime_scale $ATIME)"
          for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/diag/rtps/init${ifile:$restartbaselen}
          done
        fi
      fi

      if ((NOBS_OUT == 1)); then
        icopy=$((icopy+1))
        if ((PNETCDF == 1)); then
          mkdir -p $TMPOUT/${ATIME}/diag
          ifile="restart_${icopy}_$(datetime_scale $ATIME).nc"
          if [ -e "$TMPDIR/${ifile}" ]; then
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/diag/nobs.init.nc
          fi
        else
          mkdir -p $TMPOUT/${ATIME}/diag/nobs
          file_prefix="restart_${icopy}_$(datetime_scale $ATIME)"
          for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/diag/nobs/init${ifile:$restartbaselen}
          done
        fi
      fi
    fi # [ ENABLE_PARAM_USER == 1 ]

  else

    if ((PNETCDF == 1)); then
      mkdir -p $TMPOUT/${STIME}/hist
      mv -f $TMPDIR/history.nc $TMPOUT/${STIME}/hist/${MEM}.history.nc
    else
      mkdir -p $TMPOUT/${STIME}/hist/${MEM}
      mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/hist/${MEM}
    fi

    if ((PNETCDF == 1)); then
      mkdir -p $TMPOUT/${ATIME}/anal
      ifile="restart_$(datetime_scale $ATIME).nc"
      if [ -e "$TMPDIR/${ifile}" ]; then
        mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/${MEM}.init.nc
      fi
    else
      mkdir -p $TMPOUT/${ATIME}/anal/${MEM}
      file_prefix="restart_$(datetime_scale $ATIME)"
      restartbaselen=27
      for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
        mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/anal/${MEM}/init${ifile:$restartbaselen}
      done
    fi

    if ((ENABLE_PARAM_USER == 1)); then
      if ((OUT_OPT <= 3)) || [ "$MEM" = 'mdet' ]; then
        if ((PNETCDF == 1)); then
          mkdir -p $TMPOUT/${ATIME}/gues
          ifile="restart_2_$(datetime_scale $ATIME).nc"
          if [ -e "$TMPDIR/${ifile}" ]; then
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/${MEM}.init.nc
          fi
        else
          mkdir -p $TMPOUT/${ATIME}/gues/${MEM}
          file_prefix="restart_2_$(datetime_scale $ATIME)"
          restartbaselen=29
          for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
            mv -f $TMPDIR/${ifile} $TMPOUT/${ATIME}/gues/${MEM}/init${ifile:$restartbaselen}
          done
        fi
      fi
    fi

    if ((DELETE_MEMBER == 1)) && [ "$MEM" != 'mdet' ]; then
      if ((PNETCDF == 1)); then
        if [ -e "$TMPOUT/${STIME}/anal/${MEM}.init.nc" ]; then
          rm -f $TMPOUT/${STIME}/anal/${MEM}.init.nc
        fi
      else
        if [ -d "$TMPOUT/${STIME}/anal/${MEM}" ]; then
          rm -f $TMPOUT/${STIME}/anal/${MEM}/*
        fi
      fi
    fi

  fi

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/scale/${MEM}_run.conf
    fi
  fi

elif [ "$SCPCALL" = 'fcst' ]; then

  FTIME=$(datetime $STIME $FCSTLEN s)

  if ((PNETCDF == 1)); then
    mkdir -p $TMPOUT/${STIME}/fcst
    mv -f $TMPDIR/history.nc $TMPOUT/${STIME}/fcst/${MEM}.history.nc
  else
    mkdir -p $TMPOUT/${STIME}/fcst/${MEM}
    mv -f $TMPDIR/history*.nc $TMPOUT/${STIME}/fcst/${MEM}
  fi

  if ((OUT_OPT <= 1)); then
    if ((PNETCDF == 1)); then
      ifile="restart_${FTIME:0:8}-${FTIME:8:6}.000.nc"
      if [ -e "$TMPDIR/${ifile}" ]; then
        mv -f $TMPDIR/${ifile} $TMPOUT/${STIME}/fcst/${MEM}.init_${FTIME}.nc
      fi
    else
      file_prefix="restart_${FTIME:0:8}-${FTIME:8:6}.000"
      restartbaselen=27
      for ifile in $(cd $TMPDIR ; ls ${file_prefix}*.nc); do
        mv -f ${TMPDIR}/${ifile} $TMPOUT/${STIME}/fcst/${MEM}/init_${FTIME}${ifile:$restartbaselen}
      done
    fi
  fi

  if ((LOG_OPT <= 3)); then
    if [ -f "$TMPDIR/run.conf" ]; then
      mv -f $TMPDIR/run.conf $TMPOUT/${STIME}/log/${SCPCALL}_scale/${MEM}_run.conf
    fi
  fi

fi

#===============================================================================

exit 0
