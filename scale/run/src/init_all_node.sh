#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 2)); then
  cat >&2 << EOF

[init_all_node.sh] 

Usage: $0 MYRANK SCPCALL

  MYRANK   My rank number (not used)
  SCPCALL  Called from which script? (fcst/cycle)

EOF
  exit 1
fi

MYRANK="$1"; shift
SCPCALL="$1"

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
if ((BDY_FORMAT == 1)); then
  if ((DATA_BDY_TMPLOC == 1)); then
    ln -fs $TMPDAT/bdyscale/latlon_domain_catalogue.txt $RUNDIR/scale_init/latlon_domain_catalogue.txt
  elif ((DATA_BDY_TMPLOC == 2)); then
    ln -fs $TMPOUT/bdyscale/latlon_domain_catalogue.txt $RUNDIR/scale_init/latlon_domain_catalogue.txt
  fi
fi

mkdir -p $RUNDIR/scale
cp -f $DATADIR/exec/scale-les_ens $RUNDIR/scale
#ln -fs $DATADIR/exec/scale-les_ens $RUNDIR/scale
#ln -fs $DATADIR/rad/PARAG.29 $RUNDIR/scale
#ln -fs $DATADIR/rad/PARAPC.29 $RUNDIR/scale
#ln -fs $DATADIR/rad/VARDATA.RM29 $RUNDIR/scale
#ln -fs $DATADIR/rad/cira.nc $RUNDIR/scale
#ln -fs $DATADIR/rad/MIPAS/day.atm $RUNDIR/scale
#ln -fs $DATADIR/rad/MIPAS/equ.atm $RUNDIR/scale
#ln -fs $DATADIR/rad/MIPAS/sum.atm $RUNDIR/scale
#ln -fs $DATADIR/rad/MIPAS/win.atm $RUNDIR/scale
#ln -fs $DATADIR/land/param.bucket.conf $RUNDIR/scale

if [ "$SCPCALL" == 'cycle' ]; then
  mkdir -p $RUNDIR/obsope
  cp -f $DATADIR/exec/obsope $RUNDIR/obsope
#  ln -fs $DATADIR/exec/obsope $RUNDIR/obsope

  #-- H08 --
  if [ -e "$DATADIR/rttov/rtcoef_himawari_8_ahi.dat" ]; then
    cp -f $DATADIR/rttov/rtcoef_himawari_8_ahi.dat $RUNDIR/obsope
  fi
  if [ -e "$DATADIR/rttov/sccldcoef_himawari_8_ahi.dat" ]; then
    cp -f $DATADIR/rttov/sccldcoef_himawari_8_ahi.dat $RUNDIR/obsope
  fi

  mkdir -p $RUNDIR/letkf
  cp -f $DATADIR/exec/letkf $RUNDIR/letkf
#  ln -fs $DATADIR/exec/letkf $RUNDIR/letkf
fi

#===============================================================================

exit 0
