#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE model run.
#  August  2014  created,  Guo-Yuan Lien
#  October 2014  modified, Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 14)); then
  cat >&2 << EOF

[pre_scale.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK MEM INIT OCEAN LAND BDY TOPO LANDUSE STIME FCSTLEN FCSTINT HISTINT TMPDIR OUT_OPT [SCPCALL BDY_STIME SPRD_OUT RTPS_INFL_OUT NOBS_OUT]

  MYRANK   My rank number (not used)
  MEM      Name of the ensemble member
  INIT     Basename of SCALE initial files
  OCEAN    Basename of SCALE initial ocean files
  LAND     Basename of SCALE initial land files
  BDY      Basename of SCALE boundary files
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval of restart files (second)
  HISTINT  Output interval of history files (second)
  TMPDIR   Temporary directory to run the model
  OUT_OPT
  SCPCALL
  BDY_STIME  (format: YYYYMMDDHHMMSS)
  SPRD_OUT
  RTPS_INFL_OUT
  NOBS_OUT

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM="$1"; shift
INIT="$1"; shift
OCEAN="$1"; shift
LAND="$1"; shift
BDY="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
HISTINT="$1"; shift
TMPDIR="$1"; shift
OUT_OPT="$1"; shift
SCPCALL="${1:-cycle}"; shift
BDY_STIME="${1:-$STIME}"; shift
SPRD_OUT="${1:-1}"; shift
RTPS_INFL_OUT="${1:-0}"; shift
NOBS_OUT="${1:-0}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

BS_YYYY=${BDY_STIME:0:4}
BS_MM=${BDY_STIME:4:2}
BS_DD=${BDY_STIME:6:2}
BS_HH=${BDY_STIME:8:2}
BS_II=${BDY_STIME:10:2}
BS_SS=${BDY_STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

FILE_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
fi

RESTART_OUT_ADDITIONAL_COPIES=0
RESTART_OUT_ADDITIONAL_BASENAME=
if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale'
  RESTART_OUTPUT='.true.'
  if [ "$MEM" = 'mean' ]; then ###### using a variable for 'mean', 'mdet', 'sprd'
    RESTART_OUT_ADDITIONAL_COPIES=1
    RESTART_OUT_ADDITIONAL_BASENAME="\"${TMPDIR}\/restart_2\", "
    icopy=2
    if ((SPRD_OUT == 1)); then
      RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+2))
      icopy=$((icopy+1))
      RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"${TMPDIR}\/restart_${icopy}\", "
      icopy=$((icopy+1))
      RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"${TMPDIR}\/restart_${icopy}\", "
    fi
    if ((RTPS_INFL_OUT == 1)); then
      RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
      icopy=$((icopy+1))
      RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"${TMPDIR}\/restart_${icopy}\", "
    fi
    if ((NOBS_OUT == 1)); then
      RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
      icopy=$((icopy+1))
      RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"${TMPDIR}\/restart_${icopy}\", "
    fi
  elif [ "$MEM" = 'mdet' ]; then
    RESTART_OUT_ADDITIONAL_COPIES=1
    RESTART_OUT_ADDITIONAL_BASENAME="\"${TMPDIR}\/restart_2\", "
  elif ((OUT_OPT <= 3)); then
    RESTART_OUT_ADDITIONAL_COPIES=1
    RESTART_OUT_ADDITIONAL_BASENAME="\"${TMPDIR}\/restart_2\", "
  fi
else
  IO_LOG_DIR="${SCPCALL}_scale"
  if ((OUT_OPT <= 1)); then
    RESTART_OUTPUT='.true.'
  else
    RESTART_OUTPUT='.false.'
  fi
fi

DOMAIN_CATALOGUE_OUTPUT=".false."
if [ "$(basename $TMPDIR)" == '0001' ]; then ###### using a variable for '0001'
  DOMAIN_CATALOGUE_OUTPUT=".true."
fi

#===============================================================================

conf="$(cat $TMPDAT/conf/config.nml.scale | \
        sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
            -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
            -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
            -e "/!--TIME_DURATION--/a TIME_DURATION = ${FCSTLEN}.D0," \
            -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${FCSTINT}.D0," \
            -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${FCSTINT}.D0," \
            -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${FCSTINT}.D0," \
            -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${FCSTINT}.D0," \
            -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${INIT}\"," \
            -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
            -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${TMPDIR}\/restart\"," \
            -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"${TOPO}\"," \
            -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE}\"," \
            -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${BDY}\"," \
            -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = $BS_YYYY, $BS_MM, $BS_DD, $BS_HH, $BS_II, $BS_SS," \
            -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0," \
            -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${TMPDIR}\/history\"," \
            -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${HISTINT}.D0," \
            -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_monitor\"," \
            -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPDAT_CONSTDB}/land/param.bucket.conf\"," \
            -e "/!--DOMAIN_CATALOGUE_FNAME--/a DOMAIN_CATALOGUE_FNAME = \"$TMPOUT/const/log/latlon_domain_catalogue.txt\"," \
            -e "/!--DOMAIN_CATALOGUE_OUTPUT--/a DOMAIN_CATALOGUE_OUTPUT = ${DOMAIN_CATALOGUE_OUTPUT}," \
            -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/PARAG.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/PARAPC.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/VARDATA.RM29\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/cira.nc\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPDAT_CONSTDB}/rad/MIPAS\",")"
if ((ENABLE_PARAM_USER == 1)); then
  conf="$(echo "$conf" | sed -e "/!--TIME_END_RESTART_OUT--/a TIME_END_RESTART_OUT = .false.,")"
  conf="$(echo "$conf" | sed -e "/!--RESTART_OUT_ADDITIONAL_COPIES--/a RESTART_OUT_ADDITIONAL_COPIES = ${RESTART_OUT_ADDITIONAL_COPIES},")"
  conf="$(echo "$conf" | sed -e "/!--RESTART_OUT_ADDITIONAL_BASENAME--/a RESTART_OUT_ADDITIONAL_BASENAME = ${RESTART_OUT_ADDITIONAL_BASENAME}")"
else
  if [ "$OCEAN" != '-' ]; then
    conf="$(echo "$conf" | sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${OCEAN}\",")"
  fi
  if [ "$LAND" != '-' ]; then
    conf="$(echo "$conf" | sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${LAND}\",")"
  fi
fi
echo "$conf" > $TMPDIR/run.conf

if ((ENABLE_PARAM_USER == 1)); then
  conf="$(cat $TMPDAT/conf/config.nml.scale_user)"
  if [ "$OCEAN" != '-' ]; then
    conf="$(echo "$conf" | sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${OCEAN}\",")"
  fi
  if [ "$LAND" != '-' ]; then
    conf="$(echo "$conf" | sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${LAND}\",")"
  fi
  echo "$conf" >> $TMPDIR/run.conf
fi

#===============================================================================

exit 0
