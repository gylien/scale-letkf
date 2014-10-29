#!/bin/bash
#===============================================================================
#
#  Initialization functions
#  October 2014,    Guo-Yuan Lien
#
#  *Require source 'config.all' first.
#
#===============================================================================

function create_outdir {
#-------------------------------------------------------------------------------
# Create necessary subdirectories and files in the output directory.
#
# Usage: create_outdir
#
#-------------------------------------------------------------------------------

#if (($# < 1)); then
#  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
#  exit 1
#fi

#if [ -z "$1" ]; then
#  local INIT_TIME=''
#else
#  local INIT_TIME=$(datetime $1)
#fi

#===============================================================================

if [ -e "$OUTDIR" ]; then
  if [ -d "$OUTDIR" ]; then
    read -p "The output directory '$OUTDIR' exists. Do you want to overwrite? y/[n] " ans
    [ "$ans" != 'Y' ] && [ "$ans" != 'y' ] && exit 1
  else
    echo "[Error] $FUNCNAME: \$OUTDIR '$OUTDIR' exists and not a directory. Stop." >&2
    exit 1
  fi
fi

#===============================================================================

mkdir -p $OUTDIR
mkdir -p $OUTDIR/gues
#mkdir -p $OUTDIR/guesg
#mkdir -p $OUTDIR/guesgp
mkdir -p $OUTDIR/anal
#mkdir -p $OUTDIR/analg
#mkdir -p $OUTDIR/analgp
mkdir -p $OUTDIR/infl
mkdir -p $OUTDIR/log
mkdir -p $OUTDIR/obs

mkdir -p $OUTDIR/gues/mean
mkdir -p $OUTDIR/gues/sprd
#mkdir -p $OUTDIR/guesg/mean
#mkdir -p $OUTDIR/guesg/sprd
#mkdir -p $OUTDIR/guesgp/mean

mkdir -p $OUTDIR/anal/mean
mkdir -p $OUTDIR/anal/sprd
#mkdir -p $OUTDIR/analg/mean
#mkdir -p $OUTDIR/analg/sprd
#mkdir -p $OUTDIR/analgp/mean

#if [ "$INIT_TIME" != '' ]; then
#  STIMEgrads=$(datetimegrads $(datetime $INIT_TIME $LCYCLE h))
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                   $OUTDIR/guesg/mean/yyyymmddhh.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                   $OUTDIR/guesg/sprd/yyyymmddhh.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                   $OUTDIR/guesgp/mean/yyyymmddhhp.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                   $OUTDIR/analg/mean/yyyymmddhh.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                   $OUTDIR/analg/sprd/yyyymmddhh.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                   $OUTDIR/analgp/mean/yyyymmddhhp.ctl
#  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                   $OUTDIR/infl/yyyymmddhh.ctl
#fi

#-------------------------------------------------------------------------------

for m in $(seq $MEMBER); do
  mem=$(printf $MEMBER_FMT $m)
  mkdir -p $OUTDIR/gues/$mem
#  mkdir -p $OUTDIR/guesg/$mem
#  mkdir -p $OUTDIR/guesgp/$mem
  mkdir -p $OUTDIR/anal/$mem
#  mkdir -p $OUTDIR/analg/$mem
#  mkdir -p $OUTDIR/analgp/$mem

#  if [ "$INIT_TIME" != '' ]; then
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
#                     $OUTDIR/guesg/$mem/yyyymmddhhx.ctl
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                     $OUTDIR/guesgp/$mem/yyyymmddhhp.ctl
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
#                     $OUTDIR/analg/$mem/yyyymmddhh.ctl
#    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
#                     $OUTDIR/analgp/$mem/yyyymmddhhp.ctl
#  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================
