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

#-------------------------------------------------------------------------------

. config.${SCPCALL}
. src/func_datetime.sh
. src/func_${SCPCALL}.sh

setting

#===============================================================================

mkdir -p $TMPRUN/scale_pp
ln -fs $TMPDAT/exec/scale-les_pp_ens $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/PARAG.29 $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/PARAPC.29 $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/VARDATA.RM29 $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/cira.nc $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/MIPAS/day.atm $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/MIPAS/equ.atm $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/MIPAS/sum.atm $TMPRUN/scale_pp
ln -fs $TMPDAT/rad/MIPAS/win.atm $TMPRUN/scale_pp
ln -fs $TMPDAT/land/param.bucket.conf $TMPRUN/scale_pp
if [ "$TOPO_FORMAT" != 'prep' ]; then
  rm -f $TMPRUN/scale_pp/input_topo
  ln -fs $TMPDAT/topo/${TOPO_FORMAT}/Products $TMPRUN/scale_pp/input_topo
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  rm -f $TMPRUN/scale_pp/input_landuse
  ln -fs $TMPDAT/landuse/${LANDUSE_FORMAT}/Products $TMPRUN/scale_pp/input_landuse
fi

mkdir -p $TMPRUN/scale_init
ln -fs $TMPDAT/exec/scale-les_init_ens $TMPRUN/scale_init
ln -fs $TMPDAT/rad/PARAG.29 $TMPRUN/scale_init
ln -fs $TMPDAT/rad/PARAPC.29 $TMPRUN/scale_init
ln -fs $TMPDAT/rad/VARDATA.RM29 $TMPRUN/scale_init
ln -fs $TMPDAT/rad/cira.nc $TMPRUN/scale_init
ln -fs $TMPDAT/rad/MIPAS/day.atm $TMPRUN/scale_init
ln -fs $TMPDAT/rad/MIPAS/equ.atm $TMPRUN/scale_init
ln -fs $TMPDAT/rad/MIPAS/sum.atm $TMPRUN/scale_init
ln -fs $TMPDAT/rad/MIPAS/win.atm $TMPRUN/scale_init
ln -fs $TMPDAT/land/param.bucket.conf $TMPRUN/scale_init
if ((BDY_FORMAT == 1)); then
  if ((DATA_BDY_TMPLOC == 1)); then
    ln -fs $TMPDAT/bdyscale/latlon_domain_catalogue.txt $TMPRUN/scale_init/latlon_domain_catalogue.txt
  elif ((DATA_BDY_TMPLOC == 2)); then
    ln -fs $TMPOUT/bdyscale/latlon_domain_catalogue.txt $TMPRUN/scale_init/latlon_domain_catalogue.txt
  fi
fi

mkdir -p $TMPRUN/scale
cp -f $TMPDAT/exec/scale-les_ens $TMPRUN/scale

if [ "$SCPCALL" = 'cycle' ]; then
  mkdir -p $TMPRUN/obsope
  cp -f $TMPDAT/exec/obsope $TMPRUN/obsope
  #-- H08 --
  if [ -e "$TMPDAT/rttov/rtcoef_himawari_8_ahi.dat" ]; then
    cp -f $TMPDAT/rttov/rtcoef_himawari_8_ahi.dat $TMPRUN/obsope
  fi
  if [ -e "$TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat" ]; then
    cp -f $TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat $TMPRUN/obsope
  fi

  mkdir -p $TMPRUN/letkf
  cp -f $TMPDAT/exec/letkf $TMPRUN/letkf
fi

if [ "$SCPCALL" = 'cycle' ]; then
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  while ((time <= ETIME)); do
    mkdir -p $TMPOUT/${time}/log/scale_pp
    mkdir -p $TMPOUT/${time}/log/scale_init
    mkdir -p $TMPOUT/${time}/log/scale
    mkdir -p $TMPOUT/${atime}/log/obsope
    mkdir -p $TMPOUT/${atime}/log/letkf

    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
elif [ "$SCPCALL" = 'fcst' ]; then   ###### -- not verified !! ######
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do
    mkdir -p $TMPOUT/${time}/log/scale_pp
    mkdir -p $TMPOUT/${time}/log/scale_init
    mkdir -p $TMPOUT/${time}/log/scale

    time=$(datetime $time $((lcycles * CYCLE)) s)
  done
fi

#===============================================================================

exit 0
