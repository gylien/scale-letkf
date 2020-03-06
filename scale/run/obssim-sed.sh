#/bin/bash

ORG=obssim-radar-lt-ofp.sh

TMEM_S=1
# Total member
TMEM=320
#TMEM=20

DMEM=10

rm -f ???[0-9]_obssim.sh


for SMEM in `seq ${TMEM_S} ${DMEM} ${TMEM}`; do

  if (( SMEM >= TMEM )); then
    exit
  fi

  EMEM=$((SMEM + DMEM - 1))
  if (( EMEM > TMEM )); then
    EMEM=$TMEM
  fi
  SMEM0=`printf %04d $SMEM`
  echo $SMEM" "$EMEM" "$SMEM0

  cat ${ORG} | \
  sed -e "/#--SMEM--/a SMEM=$SMEM" \
      -e "/#--EMEM--/a EMEM=$EMEM" \
  > ${SMEM0}_obssim.sh

  bash ${SMEM0}_obssim.sh

done


