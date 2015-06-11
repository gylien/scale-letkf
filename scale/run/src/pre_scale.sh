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

Usage: $0 MYRANK MEM_NP INIT OCEAN BDY TOPO LANDUSE STIME FCSTLEN FCSTINT HISTINT TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  INIT     Basename of SCALE initial files
  OCEAN    Basename of SCALE initial ocean files
  BDY      Basename of SCALE boundary files
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval of restart files (second)
  HISTINT  Output interval of history files (second)
  TMPDIR   Temporary directory to run the model
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
INIT="$1"; shift
OCEAN="$1"; shift
BDY="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
HISTINT="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift
DATADIR="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

#ln -fs $EXECDIR/scale-les $TMPDIR

#ln -fs $DATADIR/rad/PARAG.29 $TMPDIR
#ln -fs $DATADIR/rad/PARAPC.29 $TMPDIR
#ln -fs $DATADIR/rad/VARDATA.RM29 $TMPDIR
#ln -fs $DATADIR/rad/cira.nc $TMPDIR
#ln -fs $DATADIR/rad/MIPAS/day.atm $TMPDIR
#ln -fs $DATADIR/rad/MIPAS/equ.atm $TMPDIR
#ln -fs $DATADIR/rad/MIPAS/sum.atm $TMPDIR
#ln -fs $DATADIR/rad/MIPAS/win.atm $TMPDIR
#ln -fs $DATADIR/land/param.bucket.conf $TMPDIR

ln -fs ${INIT}*.nc $TMPDIR
if [ "$OCEAN" = '-' ]; then
  OCEAN=$INIT
else
  ln -fs ${OCEAN}*.nc $TMPDIR
fi
ln -fs ${BDY}*.nc $TMPDIR
ln -fs ${TOPO}*.nc $TMPDIR
ln -fs ${LANDUSE}*.nc $TMPDIR

###
### Given $mem_np, do exact loop, use exact process
###
#for q in $(seq $MEM_NP); do
#  sfx=$(printf $SCALE_SFX $((q-1)))
##  if [ -e "$INIT$sfx" ]; then
#    ln -fs $INIT$sfx $TMPDIR/init$sfx
##  fi
##  if [ -e "$BDY$sfx" ]; then
#    ln -fs $BDY$sfx $TMPDIR/boundary$sfx
##  fi
##  if [ -e "$TOPO$sfx" ]; then
#    ln -fs $TOPO$sfx $TMPDIR/topo$sfx
##  fi
##  if [ -e "$LANDUSE$sfx" ]; then
#    ln -fs $LANDUSE$sfx $TMPDIR/landuse$sfx
##  fi
#done

#===============================================================================

TMPSUBDIR=$(basename "$(cd "$TMPDIR" && pwd)")

cat $TMPDAT/conf/config.nml.scale | \
    sed -e "s/\[IO_LOG_BASENAME\]/ IO_LOG_BASENAME = \"${TMPSUBDIR}\/LOG\",/" \
        -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[TIME_DURATION\]/ TIME_DURATION = ${FCSTLEN}.D0,/" \
        -e "s/\[TIME_DT_ATMOS_RESTART\]/ TIME_DT_ATMOS_RESTART = ${FCSTINT}.D0,/" \
        -e "s/\[TIME_DT_OCEAN_RESTART\]/ TIME_DT_OCEAN_RESTART = ${FCSTINT}.D0,/" \
        -e "s/\[TIME_DT_LAND_RESTART\]/ TIME_DT_LAND_RESTART = ${FCSTINT}.D0,/" \
        -e "s/\[TIME_DT_URBAN_RESTART\]/ TIME_DT_URBAN_RESTART = ${FCSTINT}.D0,/" \
        -e "s/\[RESTART_IN_BASENAME\]/ RESTART_IN_BASENAME = \"${TMPSUBDIR}\/$(basename ${INIT})\",/" \
        -e "s/\[RESTART_OUT_BASENAME\]/ RESTART_OUT_BASENAME = \"${TMPSUBDIR}\/restart\",/" \
        -e "s/\[TOPO_IN_BASENAME\]/ TOPO_IN_BASENAME = \"${TMPSUBDIR}\/$(basename ${TOPO})\",/" \
        -e "s/\[LANDUSE_IN_BASENAME\]/ LANDUSE_IN_BASENAME = \"${TMPSUBDIR}\/$(basename ${LANDUSE})\",/" \
        -e "s/\[ATMOS_BOUNDARY_IN_BASENAME\]/ ATMOS_BOUNDARY_IN_BASENAME = \"${TMPSUBDIR}\/$(basename ${BDY})\",/" \
        -e "s/\[OCEAN_RESTART_IN_BASENAME\]/ OCEAN_RESTART_IN_BASENAME = \"${TMPSUBDIR}\/$(basename ${OCEAN})\",/" \
        -e "s/\[HISTORY_DEFAULT_BASENAME\]/ HISTORY_DEFAULT_BASENAME = \"${TMPSUBDIR}\/history\",/" \
        -e "s/\[HISTORY_DEFAULT_TINTERVAL\]/ HISTORY_DEFAULT_TINTERVAL = ${HISTINT}.D0,/" \
        -e "s/\[MONITOR_OUT_BASENAME\]/ MONITOR_OUT_BASENAME = \"${TMPSUBDIR}\/monitor\",/" \
    > $TMPDIR/run.conf

#===============================================================================

exit 0
