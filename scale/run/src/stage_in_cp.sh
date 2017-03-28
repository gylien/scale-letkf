#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run on a single node)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 3)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST STGDIR [THREAD]

   NRANKS   Total number of nodes
   STGLIST  File of the stage-in list
   STGDIR   Directory where the files are staged into
   THREAD   Number of parallel copying threads (default: 1)

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
THREAD="${1:-1}"

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
  while read line; do
    stage_in_cp_sub
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    while read line; do
      stage_in_cp_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

if ((THREAD > 1)); then
  wait
fi

#===============================================================================

exit 0
