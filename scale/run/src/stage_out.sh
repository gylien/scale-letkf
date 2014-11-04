#!/bin/bash
#===============================================================================
#
#  Built-in stage-out script
#
#===============================================================================

. config.all

MYRANK="$1"   # s: run on the server node (create directories)

if [ "$MYRANK" == '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------
# Files in TMPOUT directory

if ((TMPOUT_MODE >= 2)); then
  if [ -s "$STAGING_DIR/stageout.out.$((MYRANK+1))" ]; then
    while read line; do
      destin="$(echo $line | cut -d '|' -s -f1)"
      source="$(echo $line | cut -d '|' -s -f2)"
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        if [ "$MYRANK" == 's' ]; then
          mkdir -p "$(dirname ${destin})"
        else
          $SCP -r "${TMPOUT}/${source}" "${SCP_HOSTPREFIX}${destin}"
        fi
      fi
    done < "$STAGING_DIR/stagein.out.$((MYRANK+1))"
  fi
fi

#===============================================================================

exit 0
