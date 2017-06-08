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
SCPCALL="$1"; shift

#-------------------------------------------------------------------------------

. config.${SCPCALL}
. src/func_datetime.sh
. src/func_${SCPCALL}.sh

setting

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-1" >&2
fi

#===============================================================================

mkdir -p $TMPRUN/scale_pp
cp -f $TMPDAT/exec/scale-rm_pp_ens $TMPRUN/scale_pp

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-2" >&2
fi

mkdir -p $TMPRUN/scale_init
cp -f $TMPDAT/exec/scale-rm_init_ens $TMPRUN/scale_init

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-3" >&2
fi

mkdir -p $TMPRUN/scale
cp -f $TMPDAT/exec/scale-rm_ens $TMPRUN/scale

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-4" >&2
fi

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
  if [ -e "$TMPDAT/rttov/rtcoef_himawari_8_ahi.bin" ]; then
    cp -f $TMPDAT/rttov/rtcoef_himawari_8_ahi.bin $TMPRUN/obsope
  fi
  if [ -e "$TMPDAT/rttov/sccldcoef_himawari_8_ahi.bin" ]; then
    cp -f $TMPDAT/rttov/sccldcoef_himawari_8_ahi.bin $TMPRUN/obsope
  fi

  mkdir -p $TMPRUN/letkf
  cp -f $TMPDAT/exec/letkf $TMPRUN/letkf
  #-- H08 --
  if [ -e "$TMPDAT/rttov/rtcoef_himawari_8_ahi.dat" ]; then
    cp -f $TMPDAT/rttov/rtcoef_himawari_8_ahi.dat $TMPRUN/letkf
  fi
  if [ -e "$TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat" ]; then
    cp -f $TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat $TMPRUN/letkf
  fi
  if [ -e "$TMPDAT/rttov/rtcoef_himawari_8_ahi.bin" ]; then
    cp -f $TMPDAT/rttov/rtcoef_himawari_8_ahi.bin $TMPRUN/letkf
  fi
  if [ -e "$TMPDAT/rttov/sccldcoef_himawari_8_ahi.bin" ]; then
    cp -f $TMPDAT/rttov/sccldcoef_himawari_8_ahi.bin $TMPRUN/letkf
  fi
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-5" >&2
fi

mkdir -p $TMPOUT/const/topo
if ((LANDUSE_UPDATE != 1)); then
  mkdir -p $TMPOUT/const/landuse
fi

if [ "$SCPCALL" = 'cycle' ]; then
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  while ((time <= ETIME)); do
    if ((LANDUSE_UPDATE == 1)); then
      mkdir -p $TMPOUT/${time}/landuse
    fi
    mkdir -p $TMPOUT/${time}/log/scale_pp
    mkdir -p $TMPOUT/${time}/log/scale_init
    mkdir -p $TMPOUT/${time}/log/scale
    mkdir -p $TMPOUT/${atime}/log/obsope
    mkdir -p $TMPOUT/${atime}/log/letkf

    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
elif [ "$SCPCALL" = 'fcst' ]; then
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do
    if ((LANDUSE_UPDATE == 1)); then
      mkdir -p $TMPOUT/${time}/landuse
    fi
    mkdir -p $TMPOUT/${time}/log/${SCPCALL}_scale_pp
    mkdir -p $TMPOUT/${time}/log/${SCPCALL}_scale_init
    mkdir -p $TMPOUT/${time}/log/${SCPCALL}_scale

    time=$(datetime $time $lcycles s)
  done
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 5-6" >&2
fi

#===============================================================================

exit 0
