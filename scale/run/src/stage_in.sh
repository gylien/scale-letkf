#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script
#
#===============================================================================

. config.all
. src/func_util.sh

MYRANK="$1"

if [ "$MYRANK" == '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------

if (((TMPDAT_MODE <= 2 && MYRANK == 0) || (TMPDAT_MODE == 3))); then
  safe_init_tmpdir $TMPDAT
fi
if (((TMPRUN_MODE <= 2 && MYRANK == 0) || (TMPRUN_MODE == 3))); then
  safe_init_tmpdir $TMPRUN
fi
if (((TMPOUT_MODE <= 2 && MYRANK == 0) || (TMPOUT_MODE == 3))); then
  safe_init_tmpdir $TMPOUT
fi

#-------------------------------------------------------------------------------
# Files in TMPDAT directory

if (((TMPDAT_MODE == 2 && MYRANK == 0) || (TMPDAT_MODE == 3))); then

  if [ -s "$STAGING_DIR/stagein.dat" ]; then
    cd $TMPDAT
    while read line; do
      source="$(echo $line | cut -d '|' -s -f1)"
      destin="$(echo $line | cut -d '|' -s -f2)"
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        mkdir -p "$(dirname $destin)"
        $SCP -r "${SCP_HOSTPREFIX}${source}" "$destin"
      fi
    done < "$STAGING_DIR/stagein.dat"
  fi

fi

#-------------------------------------------------------------------------------
# Files in TMPOUT directory

if ((TMPOUT_MODE >= 2)); then
  if [ -s "$STAGING_DIR/stagein.out.$((MYRANK+1))" ]; then
    cd $TMPOUT
    while read line; do
      source="$(echo $line | cut -d '|' -s -f1)"
      destin="$(echo $line | cut -d '|' -s -f2)"
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        mkdir -p "$(dirname $destin)"
        $SCP -r "${SCP_HOSTPREFIX}${source}" "$destin"
      fi
    done < "$STAGING_DIR/stagein.out.$((MYRANK+1))"
  fi
fi

#===============================================================================

exit 0
