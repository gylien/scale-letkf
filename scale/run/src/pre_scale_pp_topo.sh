#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE topo creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 5)); then
  cat >&2 << EOF

[scale_pre.sh] Prepare a temporary directory for SCALE model run,

Usage: $0 MYRANK STIME TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  TMPDIR   Temporary directory to run scale-les_pp
  EXECDIR  Directory of SCALE executable files
  DATADIR  Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"
STIME="$2"
TMPDIR="$3"
EXECDIR="$4"
DATADIR="$5"

#if ((MACHINE_TYPE == 10)); then
#  TMPDIR="$(pwd)/$TMPDIR"
#  EXECDIR="$(pwd)/$EXECDIR"
#  DATADIR="$(pwd)/$DATADIR"
#fi

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

ln -fs $EXECDIR/scale-les_pp $TMPDIR

ln -fs $DATADIR/topo/DEM50M/Products $TMPDIR/input

#===============================================================================


echo 777
pwd
ls -l
ls -l $TMPDAT
ls -l $TMPDAT/conf


cat $TMPDAT/conf/scale_pp_topo.conf | \
    sed -e "s/\[TIME_STARTDATE\]/ TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,/" \
    > $TMPDIR/pp.conf

#===============================================================================

exit 0
