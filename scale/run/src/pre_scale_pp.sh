#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE topo creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[pre_scale_pp.sh] Prepare a temporary directory for scale-rm_pp run.

Usage: $0 MYRANK STIME MEM TMPDIR DATADIR [SCPCALL COPYTOPO CATALOGUE]

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM
  TMPDIR   Temporary directory to run scale-rm_pp
  DATADIR  Directory of SCALE data files
  SCPCALL  Called from which script? (fcst/cycle)
  COPYTOPO
  CATALOGUE

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
DATADIR="$1"; shift   ###### no use
SCPCALL="${1:-cycle}"; shift
COPYTOPO="$1"; shift
CATALOGUE="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
else
  FILE_AGGREGATE=".false"
fi

if ((PNETCDF == 1)); then
  TOPO_OUT_BASENAME="$TMPOUT/const/topo"
  if ((LANDUSE_UPDATE == 1)); then
    LANDUSE_OUT_BASENAME="$TMPOUT/${STIME}/landuse"
  else
    LANDUSE_OUT_BASENAME="$TMPOUT/const/landuse"
  fi
else
  TOPO_OUT_BASENAME="$TMPOUT/const/topo/topo"
  if ((LANDUSE_UPDATE == 1)); then
    LANDUSE_OUT_BASENAME="$TMPOUT/${STIME}/landuse/landuse"
  else
    LANDUSE_OUT_BASENAME="$TMPOUT/const/landuse/landuse"
  fi
fi

if [ "$TOPO_FORMAT" != 'prep' ]; then
  CONVERT_TOPO='.true.'
else
  CONVERT_TOPO='.false.'
fi

if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  CONVERT_LANDUSE='.true.'
else
  CONVERT_LANDUSE='.false.'
fi

OFFLINE_PARENT_BASENAME=
if ((BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
  OFFLINE_PARENT_BASENAME="$COPYTOPO"
fi

if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale_pp'
else
  IO_LOG_DIR="${SCPCALL}_scale_pp"
fi

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_pp | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
        -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
        -e "/!--TOPO_OUT_BASENAME--/a TOPO_OUT_BASENAME = \"${TOPO_OUT_BASENAME}\"," \
        -e "/!--LANDUSE_OUT_BASENAME--/a LANDUSE_OUT_BASENAME = \"${LANDUSE_OUT_BASENAME}\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--CONVERT_TOPO--/a CONVERT_TOPO = $CONVERT_TOPO," \
        -e "/!--CONVERT_LANDUSE--/a CONVERT_LANDUSE = $CONVERT_LANDUSE," \
        -e "/!--CNVTOPO_name--/a CNVTOPO_name = \"$TOPO_FORMAT\"," \
        -e "/!--GTOPO30_IN_DIR--/a GTOPO30_IN_DIR = \"${TMPDAT_CONSTDB}/topo/GTOPO30/Products\"," \
        -e "/!--DEM50M_IN_DIR--/a DEM50M_IN_DIR = \"${TMPDAT_CONSTDB}/topo/DEM50M/Products\"," \
        -e "/!--CNVLANDUSE_name--/a CNVLANDUSE_name = '$LANDUSE_FORMAT'," \
        -e "/!--GLCCv2_IN_DIR--/a GLCCv2_IN_DIR = \"${TMPDAT_CONSTDB}/landuse/GLCCv2/Products\"," \
        -e "/!--LU100M_IN_DIR--/a LU100M_IN_DIR = \"${TMPDAT_CONSTDB}/landuse/LU100M/Products\"," \
        -e "/!--COPYTOPO_IN_BASENAME--/a COPYTOPO_IN_BASENAME = \"${COPYTOPO}\"," \
        -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${CATALOGUE}\"," \
        -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${OFFLINE_PARENT_BASENAME}\"," \
        -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
        -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
    > $TMPDIR/pp.conf

#===============================================================================

exit 0
