#!/bin/bash
#===============================================================================
#
#  Create stage-in and stage-out scripts for the K-computer
#
#===============================================================================

. config.main

STAGING_DIR="$1"; shift
PROGNAME="$1"

#===============================================================================
# stage-in: Files in TMPDAT directory

if [ -s "$STAGING_DIR/stagein.dat" ]; then
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if ((USE_RANKDIR == 1)); then
          if ((TMPDAT_MODE == 3)); then
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPDAT_STG}/${destin} recursive=10\""
          else
            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPDAT_STG}/${destin} recursive=10\""
          fi
        else
          echo "#PJM --stgin-dir \"${source} ${TMPDAT_STG}/${destin} recursive=10\""
        fi
      else
        if ((USE_RANKDIR == 1)); then
          if ((TMPDAT_MODE == 3)); then
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPDAT_STG}/${destin}\""
          else
            echo "#PJM --stgin \"rank=0 ${source} 0:${TMPDAT_STG}/${destin}\""
          fi
        else
          echo "#PJM --stgin \"${source} ${TMPDAT_STG}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.dat" # | sort | uniq
fi

#-------------------------------------------------------------------------------
# stage-in: Files in TMPOUT directory

if [ -s "$STAGING_DIR/stagein.out" ]; then
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if ((USE_RANKDIR == 1)); then
          if ((TMPOUT_MODE == 3)); then
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPOUT_STG}/${destin} recursive=10\""
          else
            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPOUT_STG}/${destin} recursive=10\""
          fi
        else
          echo "#PJM --stgin-dir \"${source} ${TMPOUT_STG}/${destin} recursive=10\""
        fi
      else
        if ((USE_RANKDIR == 1)); then
          if ((TMPOUT_MODE == 3)); then
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPOUT_STG}/${destin}\""
          else
            echo "#PJM --stgin \"rank=0 ${source} 0:${TMPOUT_STG}/${destin}\""
          fi
        else
          echo "#PJM --stgin \"${source} ${TMPOUT_STG}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.out" # | sort | uniq
fi

#-------------------

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
  done < "$STAGING_DIR/stagein.out.$((i+1))" # | sort | uniq
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

mkdir -p $LOGDIR
touch $LOGDIR/${PROGNAME}.err

if ((USE_RANKDIR == 1)); then
  echo "#PJM --stgin \"rank=* $TMPS/config.main %r:./config.main\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.rc %r:./config.rc\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.${PROGNAME} %r:./config.${PROGNAME}\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}.sh %r:./${PROGNAME}.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}_step.sh %r:./${PROGNAME}_step.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/src/* %r:./src/\""
  echo "#PJM --stgin \"rank=0 $LOGDIR/${PROGNAME}.err 0:./log/${PROGNAME}.err\""
else
  echo "#PJM --stgin \"$TMPS/config.main ./config.main\""
  echo "#PJM --stgin \"$SCRP_DIR/config.rc ./config.rc\""
  echo "#PJM --stgin \"$SCRP_DIR/config.${PROGNAME} ./config.${PROGNAME}\""
  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}.sh ./${PROGNAME}.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}_step.sh ./${PROGNAME}_step.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/src/* ./src/\""
  echo "#PJM --stgin \"$LOGDIR/${PROGNAME}.err ./log/${PROGNAME}.err\""
fi

#===============================================================================
# stage-out: Files in TMPOUT directory

i=0
while [ -s "$STAGING_DIR/stageout.out.$((i+1))" ]; do
  while read line; do
    destin="$(echo $line | cut -d '|' -s -f1)"
    source="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if ((USE_RANKDIR == 1)); then
        echo "#PJM --stgout \"rank=${i} ${i}:${TMPOUT_STG}/${source} ${destin}\""
      else
        echo "#PJM --stgout \"${TMPOUT_STG}/${source} ${destin}\""
      fi
    fi
  done < "$STAGING_DIR/stageout.out.$((i+1))" # | sort | uniq
  i=$((i+1))
done

#-------------------------------------------------------------------------------
# stage-out: standard log files

if ((USE_RANKDIR == 1)); then
  echo "#PJM --stgout \"rank=0 0:./log/* $LOGDIR/\""
else
  echo "#PJM --stgout \"./log/* $LOGDIR/\""
fi

#===============================================================================

exit 0
