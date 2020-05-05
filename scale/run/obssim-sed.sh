#/bin/bash

ORG=obssim-radar-lt-ofp.sh

TMEM_S=0
TMEM_S=21
# Total member
TMEM=320

DMEM=3
DMEM=5

#TMEM_S=111
#TMEM=215

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


