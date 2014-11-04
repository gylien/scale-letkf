#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE model run.
#  August  2014  created,  Guo-Yuan Lien
#  October 2014  modified, Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 12)); then
  cat >&2 << EOF

[scale_pre.sh] Prepare a temporary directory for SCALE model run,

Usage: $0 MYRANK MEM_NP INIT BDY TOPO LANDUSE STIME FCSTLEN FCSTINT TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  INIT     Basename of SCALE initial files
  BDY      Basename of SCALE boundary files
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval (second)
  TMPDIR   Temporary directory to run the model
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"
MEM_NP="$2"
INIT="$3"
BDY="$4"
TOPO="$5"
LANDUSE="$6"
STIME="$7"
FCSTLEN="$8"
FCSTINT="$9"
TMPDIR="${10}"
EXECDIR="${11}"
DATADIR="${12}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

ln -fs $EXECDIR/scale-les $TMPDIR

ln -fs $DATADIR/rad/PARAG.29 $TMPDIR
ln -fs $DATADIR/rad/PARAPC.29 $TMPDIR
ln -fs $DATADIR/rad/VARDATA.RM29 $TMPDIR
ln -fs $DATADIR/rad/cira.nc $TMPDIR
ln -fs $DATADIR/rad/MIPAS/day.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/equ.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/sum.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/win.atm $TMPDIR

###
### Given $mem_np, do exact loop, use exact process
###
for q in $(seq $MEM_NP); do
  sfx=$(printf $SCALE_SFX $((q-1)))
#  if [ -e "$INIT$sfx" ]; then
    ln -fs $INIT$sfx $TMPDIR/init$sfx
#  fi
#  if [ -e "$BDY$sfx" ]; then
    ln -fs $BDY$sfx $TMPDIR/boundary$sfx
#  fi
#  if [ -e "$TOPO$sfx" ]; then
    ln -fs $TOPO$sfx $TMPDIR/topo$sfx
#  fi
#  if [ -e "$LANDUSE$sfx" ]; then
    ln -fs $LANDUSE$sfx $TMPDIR/landuse$sfx
#  fi
done

#===============================================================================

cat $TMPDAT/conf/scale.conf | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[TIME_DURATION\]/ TIME_DURATION = ${FCSTLEN}.D0,/" \
        -e "s/\[HISTORY_DEFAULT_TINTERVAL\]/ HISTORY_DEFAULT_TINTERVAL = ${FCSTINT}.D0,/" \
    > $TMPDIR/run.conf

#===============================================================================

exit 0
