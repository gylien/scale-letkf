#!/bin/bash
#===============================================================================
#
#  Built-in stage-out script
#
#===============================================================================

. config.main

MYRANK="$1"  # a: run on the server node, stage out all files
             # s: run on the server node, only create directories
LOOP="$2"    # Number of loop

if [ "$MYRANK" = '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------
# Files in TMPOUT directory

filelist=
if [ "$MYRANK" = 'a' ] || [ "$MYRANK" = 's' ]; then
  if ((ONLINE_STGOUT == 1)); then
    filelist=$(ls $STAGING_DIR/stageout.loop.${LOOP}.* 2> /dev/null)
  else
    filelist=$(ls $STAGING_DIR/stageout.out.* 2> /dev/null)
  fi
elif ((TMPOUT_MODE >= 2)); then
  if ((ONLINE_STGOUT == 1)); then
    filelist=$(ls $STAGING_DIR/stageout.loop.${LOOP}.$((MYRANK+1)) 2> /dev/null)
  else
    filelist=$(ls $STAGING_DIR/stageout.out.$((MYRANK+1)) 2> /dev/null)
  fi
fi

for ifile in $filelist; do
  while read line; do
    destin="$(echo $line | cut -d '|' -s -f1)"
    source="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ "$MYRANK" = 'a' ] || [ "$MYRANK" = 's' ]; then
        mkdir -p "$(dirname ${destin})"
      fi
      if [ "$MYRANK" = 'a' ] || [ "$MYRANK" != 's' ] && [ -e "${TMPOUT}/${source}" ]; then
        $SCP -r "${TMPOUT}/${source}" "${SCP_HOSTPREFIX}${destin}"
        if ((ONLINE_STGOUT == 1)); then
          flag="$(echo $line | cut -d '|' -s -f3)"
          if [ "$flag" = 'rm' ]; then
            rm -fr "${TMPOUT}/${source}"
          fi
        fi
      fi
    fi
  done < "$ifile"
done

#===============================================================================

exit 0
