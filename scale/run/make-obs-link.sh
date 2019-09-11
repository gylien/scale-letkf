#!/bin/bash

# Load functions
. src/func_datetime.sh

EXP=2000m_WK1982_LT_SN14_NATURE
EXP=2000m_WK1982_LT_SN14_NATURE_1MIN
OUTPUT=/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/$EXP

OBSSIM_RADAR_LON=180
OBSSIM_RADAR_LAT=180

OHEAD=radar3d

STIME="20010101010000"


mkdir -p $OUTPUT/$STIME/fcst/mean/obs

FCSTLEN=3600
FCSTOUT=30

FCSTLEN=1800
FCSTOUT=60
FT=0

while [ $FT -le $FCSTLEN ] 
do
  CFT=$(printf %06d $FT)
  CTIME=$(datetime $STIME ${FT} s)
  #echo $FT" "$CFT

  ORG=$OUTPUT/$STIME/fcst/mean/${OHEAD}_i${STIME}_mean_${OBSSIM_RADAR_LON}_${OBSSIM_RADAR_LAT}_t${CFT}.dat
  NEW=$OUTPUT/$STIME/fcst/mean/obs/${OHEAD}_${CTIME}.dat

  echo $ORG" "$NEW
  ln -sf ${ORG}  ${NEW}

  FT=$((FT  + FCSTOUT))
done

exit



it=1
TIME=$STIME
while [ "$TIME" != "$ETIME" ]
do
  itt=`printf "%03d" $it`


  ORG1=$OUTPUT/$STIME/fcst/mean/obs_i${STIME}_mean_${ID1}_t${itt}.dat
  NEW1=$OUTPUT/$STIME/fcst/mean/obs/${NAME1}_${TIME}.dat

  ORG2=$OUTPUT/$STIME/fcst/mean/obs_i${STIME}_mean_${ID2}_t${itt}.dat
  NEW2=$OUTPUT/$STIME/fcst/mean/obs/${NAME2}_${TIME}.dat
 

  ln -sf $ORG1 $NEW1
  ln -sf $ORG2 $NEW2

  it=$((it + 1))
  TIME=$(datetime $TIME ${FINT} s)


done

