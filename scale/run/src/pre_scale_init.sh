#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE boundary creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main
. src/func_datetime.sh

if (($# < 10)); then
  cat >&2 << EOF

[pre_scale_init.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK MEM_NP TOPO LANDUSE WRFOUT STIME FCSTLEN TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  WRFOUT   Basename of WRF files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  TMPDIR   Temporary directory to run scale-les_init
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"
MEM_NP="$2"
TOPO="$3"
LANDUSE="$4"
WRFOUT="$5"
STIME="$6"
FCSTLEN="$7"
TMPDIR="$8"
EXECDIR="$9"
DATADIR="${10}"

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

for q in $(seq $MEM_NP); do
  sfx=$(printf $SCALE_SFX $((q-1)))
#  if [ -e "$TOPO$sfx" ]; then
    ln -fs $TOPO$sfx $TMPDIR/topo$sfx
#  fi
#  if [ -e "$LANDUSE$sfx" ]; then
    ln -fs $LANDUSE$sfx $TMPDIR/landuse$sfx
#  fi
done

NUMBER_OF_FILES=$((FCSTLEN/ANLWRF_INT+2))

i=0
time=$STIME
for c in $(seq $NUMBER_OF_FILES); do
  wrfoutfile="${WRFOUT}_${time}"
  if [ -e "$wrfoutfile" ]; then
    ln -fs $wrfoutfile $TMPDIR/wrfout_$(printf %05d $i)
  else
    echo "[Error] $0: Cannot find WRFOUT file '$wrfoutfile'."
    exit 1
  fi
  i=$((i+1))
  time=$(datetime $time $ANLWRF_INT s)
done

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_init | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[NUMBER_OF_FILES\]/ NUMBER_OF_FILES = $NUMBER_OF_FILES,/" \
        -e "s/\[BOUNDARY_UPDATE_DT\]/ BOUNDARY_UPDATE_DT = $ANLWRF_INT.D0,/" \
    > $TMPDIR/init.conf

#===============================================================================

exit 0
