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

Usage: $0 MYRANK MEM_NP TOPO LANDUSE WRFOUT STIME FCSTLEN MKINIT MEM TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  WRFOUT   Basename of WRF files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  MKINIT   Make initial condition as well?
            0: No
            1: Yes
  MEM      Name of the ensemble member
  TMPDIR   Temporary directory to run scale-les_init
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
WRFOUT="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
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

ln -fs $EXECDIR/scale-les_init $TMPDIR

ln -fs $DATADIR/rad/PARAG.29 $TMPDIR
ln -fs $DATADIR/rad/PARAPC.29 $TMPDIR
ln -fs $DATADIR/rad/VARDATA.RM29 $TMPDIR
ln -fs $DATADIR/rad/cira.nc $TMPDIR
ln -fs $DATADIR/rad/MIPAS/day.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/equ.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/sum.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/win.atm $TMPDIR
ln -fs $DATADIR/land/param.bucket.conf $TMPDIR

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
if ((MKINIT == 1)); then
  RESTART_OUTPUT='.true.'
fi

i=0
time=$STIME
etime_bdy=$(datetime $STIME $((FCSTLEN+BDYINT)) s)
while ((time < etime_bdy)); do
  if ((BDY_ENS == 1)); then
    wrfoutfile="$TMPOUT/bdywrf/${MEM}/wrfout_${time}"
  else
    wrfoutfile="$TMPOUT/bdywrf/mean/wrfout_${time}"
  fi
  if [ ! -s "$wrfoutfile" ]; then
    echo "[Error] $0: Cannot find WRFOUT file '$wrfoutfile'."
    exit 1
  fi

  ln -fs $TMPOUT/bdywrf/${MEM}/wrfout_${time} $TMPDIR/wrfout_$(printf %05d $i)

  i=$((i+1))
  time=$(datetime $time $BDYINT s)
done
NUMBER_OF_FILES=$i

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_init | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[RESTART_OUTPUT\]/ RESTART_OUTPUT = $RESTART_OUTPUT,/" \
        -e "s/\[NUMBER_OF_FILES\]/ NUMBER_OF_FILES = $NUMBER_OF_FILES,/" \
        -e "s/\[BOUNDARY_UPDATE_DT\]/ BOUNDARY_UPDATE_DT = $BDYINT.D0,/" \
    > $TMPDIR/init.conf

#===============================================================================

exit 0
