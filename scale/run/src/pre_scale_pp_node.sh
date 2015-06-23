#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of scale run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 7)); then
  cat >&2 << EOF

[pre_scale_pp_node.sh] 

Usage: $0 MYRANK MEM_NODES MEM_NP TMPDIR EXECDIR DATADIR

  MYRANK     My rank number (not used)
  MEM_NODES  Number of nodes for a member
  MEM_NP     Number of processes per member
  TMPDIR     Temporary directory to run the model
  EXECDIR    Directory of SCALE executable files
  DATADIR    Directory of SCALE data files

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NODES="$1"; shift
MEM_NP="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift
DATADIR="$1"; shift
MEMBER_RUN="$1"

#===============================================================================

#mkdir -p $TMPDIR
#rm -fr $TMPDIR/*

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

#===============================================================================

cat $TMPDAT/conf/config.nml.letkf | \
    sed -e "s/\[MEMBER\]/ MEMBER = $MEMBER,/" \
        -e "s/\[MEMBER_RUN\]/ MEMBER_RUN = $MEMBER_RUN,/" \
        -e "s/\[SLOT_START\]/ SLOT_START = 1,/" \
        -e "s/\[SLOT_END\]/ SLOT_END = 1,/" \
        -e "s/\[SLOT_BASE\]/ SLOT_BASE = 1,/" \
        -e "s/\[SLOT_TINTERVAL\]/ SLOT_TINTERVAL = $LTIMESLOT.D0,/" \
        -e "s/\[NNODES\]/ NNODES = $NNODES,/" \
        -e "s/\[PPN\]/ PPN = $PPN,/" \
        -e "s/\[MEM_NODES\]/ MEM_NODES = $MEM_NODES,/" \
        -e "s/\[MEM_NP\]/ MEM_NP = $MEM_NP,/" \
    > $TMPDIR/scale-les_pp_ens.conf

#===============================================================================

exit 0
