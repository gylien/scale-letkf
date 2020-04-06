#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE boundary creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 12)); then
  cat >&2 << EOF

[pre_scale_init.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK TOPO LANDUSE BDYORG STIME MKINIT MEM MEM_BDY TMPDIR BDY_TIME_LIST NUMBER_OF_TSTEPS NUMBER_OF_SKIP_TSTEPS [SCPCALL]

  MYRANK   My rank number (not used)
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  BDYORG   Path of the source boundary files
           SCALE history: XXX
           WRF: Basename of WRF files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MKINIT   Make initial condition as well?
            0: No
            1: Yes
  MEM      Name of the ensemble member
  MEM_BDY  Name of the ensemble member of the boundary data source
  TMPDIR   Temporary directory to run scale-rm_init
  BDY_TIME_LIST
  NUMBER_OF_TSTEPS
  NUMBER_OF_SKIP_TSTEPS
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
BDYORG="$1"; shift
STIME="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
MEM_BDY="$1"; shift
TMPDIR="$1"; shift
BDY_TIME_LIST="$1"; shift
NUMBER_OF_TSTEPS="$1"; shift
NUMBER_OF_SKIP_TSTEPS="$1"; shift
SCPCALL="${1:-cycle}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

historybaselen=7

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
else
  FILE_AGGREGATE=".false"
fi

if ((MKINIT == 1 || USE_INIT_FROM_BDY == 1)); then
  RESTART_OUTPUT='.true.'
  RESTART_OUT_BASENAME="${TMPDIR}\/init"
  if [ "$SCPCALL" = 'fcst' ] && ((MKINIT == 1)); then
    RESTART_OUT_BASENAME="${TMPDIR}\/init"
  else
    RESTART_OUT_BASENAME="${TMPDIR}\/init_bdy"
  fi
else
  RESTART_OUTPUT='.false.'
fi

if ((PNETCDF == 1)); then
  BASENAME_BOUNDARY="$TMPOUT/${STIME}/bdy/${MEM_BDY}.boundary"
else
  BASENAME_BOUNDARY="$TMPOUT/${STIME}/bdy/${MEM_BDY}/boundary"
fi

if ((BDY_FORMAT == 1)); then
  BASENAME_ORG="${TMPDIR}/bdydata"
  FILETYPE_ORG='SCALE-RM'
  LATLON_CATALOGUE_FNAME="$BDYORG/latlon_domain_catalogue.txt"
elif ((BDY_FORMAT == 2)); then
  BASENAME_ORG="${TMPDIR}/bdydata"
  FILETYPE_ORG='WRF-ARW'
  LATLON_CATALOGUE_FNAME=
elif ((BDY_FORMAT == 4)); then
  BASENAME_ORG="${TMPDIR}/gradsbdy.conf"
  FILETYPE_ORG='GrADS'
  LATLON_CATALOGUE_FNAME=
else
  echo "[Error] $0: Unsupport boundary file types" >&2
  exit 1
fi

NUMBER_OF_FILES=0
for time_bdy in $BDY_TIME_LIST; do
  NUMBER_OF_FILES=$((NUMBER_OF_FILES+1))
done

OFFLINE_PARENT_BASENAME=
if ((BDY_FORMAT == 1)); then
  if ((NUMBER_OF_FILES <= 1)); then
    OFFLINE_PARENT_BASENAME="${TMPDIR}/bdydata"
  else
    OFFLINE_PARENT_BASENAME="${TMPDIR}/bdydata_$(printf %05d 0)"
  fi
fi

i=0
for time_bdy in $BDY_TIME_LIST; do
  if ((NUMBER_OF_FILES <= 1)); then
    file_number=''
  else
    file_number="_$(printf %05d $i)"
  fi
  if ((BDY_ROTATING == 1)); then
    bdyorg_path="$(cd "${BDYORG}/${STIME}" && pwd)"
  else
    bdyorg_path="$(cd "${BDYORG}/const" && pwd)"
  fi
  if ((BDY_FORMAT == 1)); then
    if ((PNETCDF_BDY_SCALE == 1)); then
      if [ -s "${bdyorg_path}/${time_bdy}/${MEM_BDY}.history.nc" ]; then
        ln -fs "${bdyorg_path}/${time_bdy}/${MEM_BDY}.history.nc" $TMPDIR/bdydata${file_number}.nc ############ need to check !!!!
      else
        echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${time_bdy}/${MEM_BDY}.history.nc'."
        exit 1
      fi
    else
      if [ -s "${bdyorg_path}/${time_bdy}/${MEM_BDY}/history${SCALE_SFX_0}" ]; then
        for ifile in $(cd ${bdyorg_path}/${time_bdy}/${MEM_BDY} ; ls history*.nc 2> /dev/null); do
          ln -fs "${bdyorg_path}/${time_bdy}/${MEM_BDY}/${ifile}" $TMPDIR/bdydata${file_number}${ifile:$historybaselen}
        done
      else
        echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${time_bdy}/${MEM_BDY}/history.*.nc'."
        exit 1
      fi
    fi
  elif ((BDY_FORMAT == 2)); then
    if [ -s "${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}" ]; then
      ln -fs "${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}" $TMPDIR/bdydata${file_number}
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}'."
      exit 1
    fi
  elif ((BDY_FORMAT == 4)); then
    if [ -s "${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd" ]; then
      ln -fs "${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd" $TMPDIR/bdyatm${file_number}.grd
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd'."
      exit 1
    fi
    if [ -s "${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd" ]; then
      ln -fs "${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd" $TMPDIR/bdysfc${file_number}.grd
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd'."
      exit 1
    fi
    if [ -s "${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd" ]; then
      ln -fs "${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd" $TMPDIR/bdyland${file_number}.grd
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd'."
      exit 1
    fi
  fi
  i=$((i+1))
done

if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale_init'
else
  IO_LOG_DIR="${SCPCALL}_scale_init"
fi

mkdir -p $TMPOUT/${STIME}/bdy
if ((PNETCDF != 1)); then
  mkdir -p $TMPOUT/${STIME}/bdy/${MEM_BDY}
fi

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_init | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
        -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = $RESTART_OUTPUT," \
        -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${RESTART_OUT_BASENAME}\"," \
        -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE}\"," \
        -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPDAT_CONSTDB}/land/param.bucket.conf\"," \
        -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"${BASENAME_BOUNDARY}\"," \
        -e "/!--BASENAME_ORG--/a BASENAME_ORG = \"${BASENAME_ORG}\"," \
        -e "/!--FILETYPE_ORG--/a FILETYPE_ORG = \"${FILETYPE_ORG}\"," \
        -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${NUMBER_OF_FILES}," \
        -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${NUMBER_OF_TSTEPS}," \
        -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${NUMBER_OF_SKIP_TSTEPS}," \
        -e "/!--BOUNDARY_UPDATE_DT--/a BOUNDARY_UPDATE_DT = $BDYINT.D0," \
        -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${LATLON_CATALOGUE_FNAME}\"," \
        -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${OFFLINE_PARENT_BASENAME}\"," \
        -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
        -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
    > $TMPDIR/init.conf

if [ -e "$TMPDAT/conf/config.nml.grads_boundary" ]; then
  cat $TMPDAT/conf/config.nml.grads_boundary | \
      sed -e "s#--DIR--#${TMPDIR}#g" \
      > $TMPDIR/gradsbdy.conf
fi

#===============================================================================

exit 0
