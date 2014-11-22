#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run.
#  November 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.all
. src/func_distribute.sh

if (($# < 12)); then
  cat >&2 << EOF

[pre_obsope.sh] Prepare a temporary directory for obsope

Usage: $0 MYRANK STIME FCSTLEN FCSTINT SLOT_BASE SLOT_NUM TMPDIR EXECDIR OBSDIR

  MYRANK     My rank number (not used)
  STIME      Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN    Forecast length (second)
  FCSTINT    Output interval (second)
  SLOT_BASE  The base slot
  SLOT_NUM   Number of observation timeslots
  TMPDIR     Temporary directory to run the model
  EXECDIR    Directory of SCALE executable files
  OBSDIR     Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
SLOT_BASE="$1"; shift
SLOT_NUM="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift 
OBSDIR="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

ATIME=$(datetime $STIME $LCYCLE s)

historybaselen=7

#initbaselen=4 ######

#-------------------------------------------------------------------------------

distribute_da_cycle - -

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

ln -fs $EXECDIR/obsope $TMPDIR

ln -fs $OBSDIR/obs_${ATIME}.dat $TMPDIR/obs.dat

for m in $(seq $MEMBER); do
  for ifile in $(cd $TMPOUT/${ATIME}/gues/${name_m[$m]} ; ls history*.nc); do
    ln -fs $TMPOUT/${ATIME}/gues/${name_m[$m]}/history*.nc $TMPDIR/hist.${name_m[$m]}${ifile:$historybaselen}
  done
#  for ifile in $(cd $TMPOUT/${ATIME}/gues/${name_m[$m]} ; ls init*.nc); do
#    ln -fs $TMPOUT/${ATIME}/gues/${name_m[$m]}/init*.nc $TMPDIR/init.${name_m[$m]}${ifile:$initbaselen}
#  done
done

#===============================================================================

cat $TMPDAT/conf/obsope.conf | \
    sed -e "s/\[NNODES\]/ NNODES = $NNODES,/" \
        -e "s/\[PPN\]/ PPN = $PPN,/" \
        -e "s/\[MEM_NODES\]/ MEM_NODES = $mem_nodes,/" \
        -e "s/\[MEM_NP\]/ MEM_NP = $mem_np,/" \
        -e "s/\[SLOT_NUM\]/ SLOT_NUM = $SLOT_NUM,/" \
        -e "s/\[SLOT_BASE\]/ SLOT_BASE = $SLOT_BASE,/" \
        -e "s/\[SLOT_TINTERVAL\]/ SLOT_TINTERVAL = $LTIMESLOT.D0,/" \
    > $TMPDIR/obsope.conf

cat $TMPDAT/conf/scale.conf | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
        -e "s/\[TIME_DURATION\]/ TIME_DURATION = ${FCSTLEN}.D0,/" \
        -e "s/\[HISTORY_DEFAULT_TINTERVAL\]/ HISTORY_DEFAULT_TINTERVAL = ${FCSTINT}.D0,/" \
    >> $TMPDIR/obsope.conf

#===============================================================================

exit 0
