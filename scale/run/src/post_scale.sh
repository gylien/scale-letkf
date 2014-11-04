#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 7)); then
  cat >&2 << EOF

[scale_pre.sh] Prepare a temporary directory for SCALE model run,

Usage: $0 MYRANK MEM_NP STIME FCSTLEN FCSTINT TMPDIR

  MYRANK   My rank number (not used)
  MEM_NP   Number of processes per member
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM      Name of the ensemble member
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval (second)
  TMPDIR   Temporary directory to run the model

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM_NP="$1"; shift
STIME="$1"; shift
MEN="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
TMPDIR="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

#mkdir -p $TMPDIR

if ((FOUT_OPT <= 2)); then
  mkdir -p $TMPOUT/fcst/${STIME}/${MEM}
  mv -f $TMPDIR/history*.nc $TMPOUT/fcst/${STIME}/${MEM}/${STIME}


#for c in `seq $CYCLES`; do
#  STIMEgrads=$(datetimegrads ${STIME[$c]})
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    if [ "$FOUT_OPT" -le 1 ]; then
#      mkdir -p $OUTDIR/fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}
#    fi
#    mkdir -p $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 x > \
#                     $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhx.ctl
#    if [ "$FOUT_OPT" -le 2 ]; then
#      mkdir -p $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}
#      $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 p > \
#                       $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhp.ctl
#    fi

#    fh=0
#    while [ "$fh" -le "$FCSTLEN" ]; do
#      fhhh=`printf '%03d' $fh`
#      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
#      mkdir -p $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
#      cd $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
#      ln -fs ../../../fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
#      if [ ! -s 'yyyymmddhhx.ctl' ]; then
#        $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
#                         yyyymmddhhx.ctl
#      fi
#      if [ "$FOUT_OPT" -le 2 ]; then
#        mkdir -p $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
#        cd $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
#        ln -fs ../../../fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
#        if [ ! -s 'yyyymmddhhp.ctl' ]; then
#          $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                           yyyymmddhhp.ctl
#        fi
#      fi
#    fh=$((fh + FCSTOUT))
#    done
#  done
#done

##-------------------------------------------------------------------------------
#if [ "$SHAREDISK" = '0' ]; then
##-------------------------------------------------------------------------------

#cd $TMPMPI
#mkdir -p $tmpstageout
#rm -f $tmpstageout/*

#for c in `seq $CYCLES`; do
#  for m in `seq $fmember`; do
#    mt=$(((c-1) * fmember + m))
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig" >> $tmpstageout/out.${node_m[$mt]}
#    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#    fh=0
#    while [ "$fh" -le "$FCSTLEN" ]; do
#      fhhh=`printf '%03d' $fh`
#      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
#      if [ "$FOUT_OPT" -le 1 ]; then
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$mt]}
#        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      if [ "$FOUT_OPT" -le 2 ]; then
#        echo "mv|fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
#      fi
#      echo "mv|verfo1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#      echo "mv|verfa2/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
#    fh=$((fh + FCSTOUT))
#    done
#  done
#done

#===============================================================================

exit 0
