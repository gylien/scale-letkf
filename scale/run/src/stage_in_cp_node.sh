#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run parallelly on all computing nodes)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 4)); then
  cat >&2 << EOF

Usage: $MYNAME MYRANK NRANKS STGLIST STGDIR [TYPE THREAD]

   MYRANK   (number): Rank of this node
            '-'     : Determinte the rank of this node automatically
   NRANKS   Total number of nodes
   STGLIST  File of the stage-in list
   STGDIR   Directory where the files are staged into
   TYPE     'share': Stage-in to a shared directory (default)
            'local': Stage-in to local directories
   THREAD   Number of parallel copying threads (default: 1)

EOF
  exit 1
fi

MYRANK="$1"; shift
NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
TYPE="${1:-share}"; shift
THREAD="${1:-1}"

# This may not work...
if [[ "$MYRANK" == '-' ]]; then
  MYRANK=$MPI_DRANK
fi

#-------------------------------------------------------------------------------

function stage_in_cp_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
  local destinstg="${STGDIR}/${destin}"
  if [[ -z "$destin" ]]; then
    : # do nothing
  elif [[ -z "$source" ]]; then
    if [[ "$destinstg" != */ ]]; then
      if (mkdir -p "$(dirname "$destinstg")"); then
        touch "$destinstg" # create an empty file
      fi
    else
      mkdir -p "$destinstg" # create an empty directory
    fi
  elif [[ "$source" != /* ]]; then
    echo "$MYNAME: source '$source' is not an absolute path" >&2
  elif [[ "$source" != */ && "$destinstg" != */ ]]; then # files
    if [[ -d "$destinstg" ]]; then
      echo "$MYNAME: source '$source' is a regular file, but destination '$destinstg' exists and is a directory" >&2
    elif (mkdir -p "$(dirname "$destinstg")"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        cp -L "$source" "$destinstg" &
      else
        cp -L "$source" "$destinstg"
      fi
    fi
  elif [[ "$source" == */ && "$destinstg" == */ ]]; then # directories
    if (mkdir -p "$destinstg"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        cp -rL "$source"* "$destinstg" > /dev/null 2>&1 &
      else
        cp -rL "$source"* "$destinstg" > /dev/null 2>&1
      fi
    fi
  else
    echo "$MYNAME: source '$source' and destination '$destinstg' need to be the same type" >&2
  fi
}

#-------------------------------------------------------------------------------
# Check the validity of $STGDIR for safety

mkdir -p $STGDIR || exit $?

if [[ ! -d "$STGDIR" ]]; then
  echo "[Error] $MYNAME: '$STGDIR' is not a directory." >&2
  exit 1
fi
if [[ ! -O "$STGDIR" ]]; then
  echo "[Error] $MYNAME: '$STGDIR' is not owned by you." >&2
  exit 1
fi

#-------------------------------------------------------------------------------
# Stage-in

if [[ -s "$STGLIST" ]]; then
  i=0
  while read line; do
    if [[ "$TYPE" == 'local' ]] || ((i % NRANKS == MYRANK)); then
      stage_in_cp_sub
    fi
    i=$((i+1))
  done < "$STGLIST" # | sort | uniq
fi

if [[ -s "$STGLIST.$((MYRANK+1))" ]]; then
  while read line; do
    stage_in_cp_sub
  done < "$STGLIST.$((MYRANK+1))" # | sort | uniq
fi

if ((THREAD > 1)); then
  wait
fi

#===============================================================================

exit 0
