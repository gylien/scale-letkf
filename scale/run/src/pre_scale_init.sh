#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE boundary creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main
. src/func_datetime.sh

if (($# < 12)); then
  cat >&2 << EOF

[pre_scale_init.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK MEM_NP TOPO LANDUSE BDYORG STIME FCSTLEN MKINIT MEM TMPDIR EXECDIR DATADIR [STARTFRAME]

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  BDYORG   Path of the source boundary files
           SCALE history: XXX
           WRF: Basename of WRF files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  MKINIT   Make initial condition as well?
            0: No
            1: Yes
  MEM      Name of the ensemble member
  TMPDIR   Temporary directory to run scale-les_init
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files
  STARTFRAME

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
BDYORG="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift
DATADIR="$1"; shift
STARTFRAME="${1:-1}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

TMPSUBDIR=$(basename "$(cd "$TMPDIR" && pwd)")

ln -fs ${TOPO}*.nc $TMPDIR
ln -fs ${LANDUSE}*.nc $TMPDIR

#for q in $(seq $MEM_NP); do
#  sfx=$(printf $SCALE_SFX $((q-1)))
##  if [ -e "$TOPO$sfx" ]; then
#    ln -fs $TOPO$sfx $TMPDIR/topo$sfx
##  fi
##  if [ -e "$LANDUSE$sfx" ]; then
#    ln -fs $LANDUSE$sfx $TMPDIR/landuse$sfx
##  fi
#done

RESTART_OUTPUT='.false.'
if ((MKINIT == 1 || (OCEAN_INPUT == 1 && OCEAN_FORMAT == 99))); then
  RESTART_OUTPUT='.true.'
fi

if ((BDY_FORMAT == 1)); then
  if [ ! -s "${BDYORG}.pe000000.nc" ]; then
    echo "[Error] $0: Cannot find source boundary file '${BDYORG}.pe000000.nc'."
    exit 1
  fi
  ln -fs ${BDYORG}*.nc $TMPDIR
#  NUMBER_OF_FILES=$(((FCSTLEN-1)/BDYINT+2))
#  NUMBER_OF_FILES=$(((FCSTLEN-1)/BDYINT+1+STARTFRAME))
  NUMBER_OF_FILES=1
  NUMBER_OF_TSTEPS=$(((FCSTLEN-1)/BDYINT+1+STARTFRAME))
  NUMBER_OF_SKIP_TSTEPS=$((STARTFRAME-1))

  BASENAME_ORG="${TMPSUBDIR}\/history"
  FILETYPE_ORG='SCALE-LES'
  USE_NESTING='.true.'
  OFFLINE='.true.'
elif ((BDY_FORMAT == 2)); then
  i=0
  time=$STIME
  etime_bdy=$(datetime $STIME $((FCSTLEN+BDYINT)) s)
#  tmp_etime_bdy=$(datetime $STIME $((BDYINT+BDYINT)) s)  # T. Honda (may be not necessary?)
#  if (( etime_bdy < tmp_etime_bdy )); then               #
#    etime_bdy=${tmp_etime_bdy}                           #
#  fi                                                     #
  while ((time < etime_bdy)); do
    if [ -s "${BDYORG}_${time}" ]; then
      ln -fs "${BDYORG}_${time}" $TMPDIR/wrfout_$(printf %05d $i)
    else
      echo "[Error] $0: Cannot find source boundary file '${BDYORG}_${time}'."
      exit 1
    fi
    i=$((i+1))
    time=$(datetime $time $BDYINT s)
  done
  NUMBER_OF_FILES=$i
  NUMBER_OF_TSTEPS=1
  NUMBER_OF_SKIP_TSTEPS=0

  BASENAME_ORG="${TMPSUBDIR}\/wrfout"
  FILETYPE_ORG='WRF-ARW'
  USE_NESTING='.false.'
  OFFLINE='.true.'
else
  echo "[Error] $0: Unsupport boundary file types" >&2
  exit 1
fi

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_init | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"${TMPSUBDIR}\/init_LOG\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = $RESTART_OUTPUT," \
        -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${TMPSUBDIR}\/init\"," \
        -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"${TMPSUBDIR}\/topo\"," \
        -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${TMPSUBDIR}\/landuse\"," \
        -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"${TMPSUBDIR}\/boundary\"," \
        -e "/!--BASENAME_ORG--/a BASENAME_ORG = \"${BASENAME_ORG}\"," \
        -e "/!--FILETYPE_ORG--/a FILETYPE_ORG = \"${FILETYPE_ORG}\"," \
        -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = $NUMBER_OF_FILES," \
        -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = $NUMBER_OF_TSTEPS," \
        -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = $NUMBER_OF_SKIP_TSTEPS," \
        -e "/!--BOUNDARY_UPDATE_DT--/a BOUNDARY_UPDATE_DT = $BDYINT.D0," \
        -e "/!--USE_NESTING--/a USE_NESTING = $USE_NESTING," \
        -e "/!--OFFLINE--/a OFFLINE = $OFFLINE," \
    > $TMPDIR/init.conf

#===============================================================================

exit 0
