#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run on the computing-node side)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 2)); then
  cat >&2 << EOF

Usage: $MYNAME MYRANK STGDIR [TYPE]

   MYRANK   (number): Rank of this node
            '-'     : Determinte the rank of this node automatically
   STGDIR   Directory where the files are staged into
   TYPE     'share': Stage-in to a shared directory (default)

EOF
  exit 1
fi

MYRANK="$1"; shift
STGDIR="$1"; shift
TYPE="${1:-share}"

# This may not work...
if [[ "$MYRANK" == '-' ]]; then
  MYRANK=$MPI_DRANK
fi

#-------------------------------------------------------------------------------
# Check the validity of $STGDIR for safety and clean $STGDIR

if ((MYRANK == 0)) || [[ "$TYPE" == 'local' ]] ; then
  mkdir -p $STGDIR || exit $?

  if [[ ! -d "$STGDIR" ]]; then
    echo "[Error] $MYNAME: '$STGDIR' is not a directory." >&2
    exit 1
  fi
  if [[ ! -O "$STGDIR" ]]; then
    echo "[Error] $MYNAME: '$STGDIR' is not owned by you." >&2
    exit 1
  fi

  rm -fr $STGDIR/* || exit $?
fi

#===============================================================================

exit 0
