#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 3)); then
  cat >&2 << EOF

[init_all_node.sh] 

Usage: $0 MYRANK DATADIR RUNDIR

  MYRANK      My rank number (not used)
  DATADIR     Directory of data files
  RUNDIR      Directory of runtime files

EOF
  exit 1
fi

MYRANK="$1"; shift
#DATADIR="$1"; shift
#RUNDIR="$1"

DATADIR="$TMPDAT"
RUNDIR="$TMPRUN"

#===============================================================================

mkdir -p $RUNDIR/scale_pp
ln -fs $DATADIR/exec/scale-les_pp_ens $RUNDIR/scale_pp
ln -fs $DATADIR/rad/PARAG.29 $RUNDIR/scale_pp
ln -fs $DATADIR/rad/PARAPC.29 $RUNDIR/scale_pp
ln -fs $DATADIR/rad/VARDATA.RM29 $RUNDIR/scale_pp
ln -fs $DATADIR/rad/cira.nc $RUNDIR/scale_pp
ln -fs $DATADIR/rad/MIPAS/day.atm $RUNDIR/scale_pp
ln -fs $DATADIR/rad/MIPAS/equ.atm $RUNDIR/scale_pp
ln -fs $DATADIR/rad/MIPAS/sum.atm $RUNDIR/scale_pp
ln -fs $DATADIR/rad/MIPAS/win.atm $RUNDIR/scale_pp
ln -fs $DATADIR/land/param.bucket.conf $RUNDIR/scale_pp
if [ "$TOPO_FORMAT" != 'prep' ]; then
  rm -f $RUNDIR/scale_pp/input_topo
  ln -fs $DATADIR/topo/${TOPO_FORMAT}/Products $RUNDIR/scale_pp/input_topo
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  rm -f $RUNDIR/scale_pp/input_landuse
  ln -fs $DATADIR/landuse/${LANDUSE_FORMAT}/Products $RUNDIR/scale_pp/input_landuse
fi

mkdir -p $RUNDIR/scale_init
ln -fs $DATADIR/exec/scale-les_init_ens $RUNDIR/scale_init
ln -fs $DATADIR/rad/PARAG.29 $RUNDIR/scale_init
ln -fs $DATADIR/rad/PARAPC.29 $RUNDIR/scale_init
ln -fs $DATADIR/rad/VARDATA.RM29 $RUNDIR/scale_init
ln -fs $DATADIR/rad/cira.nc $RUNDIR/scale_init
ln -fs $DATADIR/rad/MIPAS/day.atm $RUNDIR/scale_init
ln -fs $DATADIR/rad/MIPAS/equ.atm $RUNDIR/scale_init
ln -fs $DATADIR/rad/MIPAS/sum.atm $RUNDIR/scale_init
ln -fs $DATADIR/rad/MIPAS/win.atm $RUNDIR/scale_init
ln -fs $DATADIR/land/param.bucket.conf $RUNDIR/scale_init

mkdir -p $RUNDIR/scale
ln -fs $DATADIR/exec/scale-les_ens $RUNDIR/scale
ln -fs $DATADIR/rad/PARAG.29 $RUNDIR/scale
ln -fs $DATADIR/rad/PARAPC.29 $RUNDIR/scale
ln -fs $DATADIR/rad/VARDATA.RM29 $RUNDIR/scale
ln -fs $DATADIR/rad/cira.nc $RUNDIR/scale
ln -fs $DATADIR/rad/MIPAS/day.atm $RUNDIR/scale
ln -fs $DATADIR/rad/MIPAS/equ.atm $RUNDIR/scale
ln -fs $DATADIR/rad/MIPAS/sum.atm $RUNDIR/scale
ln -fs $DATADIR/rad/MIPAS/win.atm $RUNDIR/scale
ln -fs $DATADIR/land/param.bucket.conf $RUNDIR/scale

mkdir -p $RUNDIR/obsope
ln -fs $DATADIR/exec/obsope $RUNDIR/obsope

mkdir -p $RUNDIR/letkf
ln -fs $DATADIR/exec/letkf $RUNDIR/letkf

#===============================================================================

exit 0
