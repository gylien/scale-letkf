#!/bin/bash

# Load functions
. src/func_datetime.sh

EXP=2000m_InSnd_LT_SN14_Mac_0605
OUTPUT=/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/$EXP

STIME="20000101004000"
ETIME="20000101010000"

FINT=30 # forecast history interval(s)

ID1="4001"
NAME1=ref

ID2="4002"
NAME2=vr

mkdir -p $OUTPUT/$STIME/fcst/mean/obs

it=0
TIME=$STIME
while [ "$TIME" != "$ETIME" ]
do
  it=$((it + 1))
  itt=`printf "%03d" $it`

  TIME=$(datetime $TIME ${FINT} s)
  ORG1=$OUTPUT/$STIME/fcst/mean/obs_i${STIME}_mean_${ID1}_t${itt}.dat
  NEW1=$OUTPUT/$STIME/fcst/mean/obs/${NAME1}_${TIME}.dat

  ORG2=$OUTPUT/$STIME/fcst/mean/obs_i${STIME}_mean_${ID2}_t${itt}.dat
  NEW2=$OUTPUT/$STIME/fcst/mean/obs/${NAME2}_${TIME}.dat

  ln -sf $ORG1 $NEW1
  ln -sf $ORG2 $NEW2


done

