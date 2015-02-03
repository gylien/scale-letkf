#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE model run.
#  August  2014  created,  Guo-Yuan Lien
#  October 2014  modified, Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 13)); then
  cat >&2 << EOF

[pre_scale.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK MEM_NP INIT BDY TOPO LANDUSE STIME FCSTLEN FCSTINT HISTINT TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  INIT     Basename of SCALE initial files
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

ln -fs $EXECDIR/scale-les $TMPDIR

ln -fs $DATADIR/rad/PARAG.29 $TMPDIR
ln -fs $DATADIR/rad/PARAPC.29 $TMPDIR
ln -fs $DATADIR/rad/VARDATA.RM29 $TMPDIR
ln -fs $DATADIR/rad/cira.nc $TMPDIR
ln -fs $DATADIR/rad/MIPAS/day.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/equ.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/sum.atm $TMPDIR
ln -fs $DATADIR/rad/MIPAS/win.atm $TMPDIR
ln -fs $DATADIR/land/param.bucket.conf $TMPDIR

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

cat $TMPDAT/conf/config.nml.scale | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[TIME_DURATION\]/ TIME_DURATION = ${FCSTLEN}.D0,/" \
        -e "s/\[TIME_DT_ATMOS_RESTART\]/ TIME_DT_ATMOS_RESTART = ${FCSTLEN}.D0,/" \
        -e "s/\[HISTORY_DEFAULT_TINTERVAL\]/ HISTORY_DEFAULT_TINTERVAL = ${HISTINT}.D0,/" \
    > $TMPDIR/run.conf


#echo
#echo $TMPDIR
#echo '##### pre #####'
#ls -lL $TMPDIR
#echo



#===============================================================================

exit 0
