#!/bin/bash
#===============================================================================
#
#  Create stage-in and stage-out scripts for the K-computer
#
#===============================================================================

. config.all
#. src/func_util.sh

STAGING_DIR="$1"

#-------------------------------------------------------------------------------
# stage-in: Files in TMPDAT directory

if [ -s "$STAGING_DIR/stagein.dat" ]; then
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if ((USE_RANKDIR == 1)); then
          if ((TMPDAT_MODE == 3)); then
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPDAT_STG}/${destin}\""
          else
            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPDAT_STG}/${destin}\""
          fi
        else
          echo "#PJM --stgin-dir \"${source} ${TMPDAT_STG}/${destin}\""
        fi
      else
        if ((USE_RANKDIR == 1)); then
          if ((TMPDAT_MODE == 3)); then
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPDAT_STG}/${destin}\""
          else
            echo "#PJM --stgin \"rank=0 ${source} 0:.${TMPDAT_STG}/${destin}\""
          fi
        else
          echo "#PJM --stgin \"${source} ${TMPDAT_STG}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.dat"
fi

#-------------------------------------------------------------------------------
# stage-in: Files in TMPOUT directory

i=0
while [ -s "$STAGING_DIR/stagein.out.$((i+1))" ]; do
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if ((USE_RANKDIR == 1)); then
        echo "#PJM --stgin \"rank=${i} ${source} ${i}:${TMPOUT_STG}/${destin}\""
      else
        echo "#PJM --stgin \"${source} ${TMPOUT_STG}/${destin}\""
      fi
    fi
  done < "$STAGING_DIR/stagein.out.$((i+1))"
  i=$((i+1))
done

#-------------------------------------------------------------------------------
# stage-in: nodefiles

if ((USE_RANKDIR == 1)); then
#  echo "#PJM --stgin-dir \"rank=* $TMPS/node %r:./node\""
  echo "#PJM --stgin-dir \"rank=0 $TMPS/node 0:./node\""
else
  echo "#PJM --stgin-dir \"$TMPS/node ./node\""
fi

#-------------------------------------------------------------------------------
# stage-in: scripts

if ((USE_RANKDIR == 1)); then
  echo "#PJM --stgin \"rank=* $TMPS/config.all %r:./config.all\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.fcst %r:./config.fcst\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/fcst.sh %r:./fcst.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/src/* %r:./src/\""
else
  echo "#PJM --stgin \"$TMPS/config.all ./config.all\""
  echo "#PJM --stgin \"$SCRP_DIR/config.fcst ./config.fcst\""
  echo "#PJM --stgin \"$SCRP_DIR/fcst.sh ./fcst.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/src/* ./src/\""
fi

#===============================================================================

exit 0
