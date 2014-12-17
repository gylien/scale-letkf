#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 13)); then
  cat >&2 << EOF

[pre_obsope_node.sh] 

Usage: $0 MYRANK STIME ATIME TMPDIR EXECDIR OBSDIR MEM_NODES MEM_NP SLOT_START SLOT_END SLOT_BASE FCSTLEN FCSTINT

  MYRANK      My rank number (not used)
  STIME       Start time (format: YYYYMMDDHHMMSS)
  ATIME       Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR      Temporary directory to run the program
  EXECDIR     Directory of SCALE executable files
  OBSDIR      Directory of SCALE data files
  MEM_NODES   Number of nodes for a member
  MEM_NP      Number of processes for a member
  SLOT_START  Start observation timeslots
  SLOT_END    End observation timeslots
  SLOT_BASE   The base slot
  FCSTLEN     Forecast length (second)
  FCSTINT     Output interval (second)

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift 
OBSDIR="$1"; shift
MEM_NODES="$1"; shift
MEM_NP="$1"; shift
SLOT_START="$1"; shift
SLOT_END="$1"; shift
SLOT_BASE="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

ln -fs $EXECDIR/obsope $TMPDIR

ln -fs $OBSDIR/obs_${ATIME}.dat $TMPDIR/obs.dat

#===============================================================================

cat $TMPDAT/conf/config.nml.obsope | \
    sed -e "s/\[NNODES\]/ NNODES = $NNODES,/" \
        -e "s/\[PPN\]/ PPN = $PPN,/" \
        -e "s/\[MEM_NODES\]/ MEM_NODES = $MEM_NODES,/" \
        -e "s/\[MEM_NP\]/ MEM_NP = $MEM_NP,/" \
        -e "s/\[SLOT_START\]/ SLOT_START = $SLOT_START,/" \
        -e "s/\[SLOT_END\]/ SLOT_END = $SLOT_END,/" \
        -e "s/\[SLOT_BASE\]/ SLOT_BASE = $SLOT_BASE,/" \
        -e "s/\[SLOT_TINTERVAL\]/ SLOT_TINTERVAL = $LTIMESLOT.D0,/" \
    > $TMPDIR/obsope.conf

cat $TMPDAT/conf/config.nml.scale | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[TIME_DURATION\]/ TIME_DURATION = ${FCSTLEN}.D0,/" \
        -e "s/\[HISTORY_DEFAULT_TINTERVAL\]/ HISTORY_DEFAULT_TINTERVAL = ${FCSTINT}.D0,/" \
    >> $TMPDIR/obsope.conf

#===============================================================================

exit 0
