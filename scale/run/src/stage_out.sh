#!/bin/bash
#===============================================================================
#
#  Built-in stage-out script
#
#===============================================================================

. config.all

MYRANK="$1"   # a: run on the server node, stage out all files
              # s: run on the server node, only create directories

if [ "$MYRANK" = '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------
# Files in TMPOUT directory

filelist=
if [ "$MYRANK" = 'a' ]; then
  filelist=$(ls $STAGING_DIR/stageout.out.* 2> /dev/null)
elif ((TMPOUT_MODE >= 2)); then
  filelist=$(ls $STAGING_DIR/stageout.out.$((MYRANK+1)) 2> /dev/null)
fi

for ifile in $filelist; do
  while read line; do
    destin="$(echo $line | cut -d '|' -s -f1)"
    source="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ] && [ -e "${TMPOUT}/${source}" ]; then
      if [ "$MYRANK" = 'a' ] || [ "$MYRANK" = 's' ]; then
        mkdir -p "$(dirname ${destin})"
      fi
      if [ "$MYRANK" = 'a' ] || [ "$MYRANK" != 's' ]; then
        $SCP -r "${TMPOUT}/${source}" "${SCP_HOSTPREFIX}${destin}"
      fi
    fi
  done < "$ifile"
done

#===============================================================================

exit 0
