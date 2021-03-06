#!/bin/bash
#===============================================================================
#
#  Built-in stage-out script (run on a single node)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 3)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST STGDIR [THREAD STEP]

   NRANKS   Total number of nodes
   STGLIST  File of the stage-out list
   STGDIR   Directory where the files are staged out from
   THREAD   Number of parallel copying threads (default: 1)
   STEP     Step ID with which the files are processed
            'a': process all steps (default)

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
THREAD="${1:-1}"; shift
STEP="${1:-a}"

#-------------------------------------------------------------------------------

function stage_out_cp_sub () {
  local destin="$(echo $line | cut -d '|' -s -f1)"
  local source="$(echo $line | cut -d '|' -s -f2)"
  local sourcestg="${STGDIR}/${source}"
  local istep="$(echo $line | cut -d '|' -s -f3)"
  local flag=
  if [[ "$STEP" != 'a' ]]; then
    local flag="$(echo $line | cut -d '|' -s -f4)"
  fi
  if [[ "$STEP" != 'a' && "$istep" -ne "$STEP" ]]; then
    : # do nothing
  elif [[ -z "$source" || -z "$destin" ]]; then
    : # do nothing
  elif [[ "$destin" != /* ]]; then
    echo "$MYNAME: destination '$destin' is not an absolute path" >&2
#  elif [[ ! -e "$sourcestg" ]]; then
#    echo "$MYNAME: source '$sourcestg' does not exists" >&2
  elif [[ "$sourcestg" != */ && "$destin" != */ ]]; then # files
    if [[ -d "$destin" ]]; then
      echo "$MYNAME: source '$sourcestg' is a regular file, but destination '$destin' exists and is a directory" >&2
    elif (mkdir -p "$(dirname "$destin")"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        ( cp -L "$sourcestg" "$destin" && [ "$flag" = 'r' ] && rm -f "$sourcestg" ) &
      else
        cp -L "$sourcestg" "$destin" && [ "$flag" = 'r' ] && rm -f "$sourcestg"
      fi
    fi
  elif [[ "$sourcestg" == */ && "$destin" == */ ]]; then # directories
    if (mkdir -p "$destin"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        ( cp -rL "$sourcestg"* "$destin" > /dev/null 2>&1 && [ "$flag" = 'r' ] && rm -fr "$sourcestg" ) &
      else
        cp -rL "$sourcestg"* "$destin" > /dev/null 2>&1 && [ "$flag" = 'r' ] && rm -fr "$sourcestg"
      fi
    fi
  else
    echo "$MYNAME: source '$sourcestg' and destination '$destin' need to be the same type" >&2
  fi
}

#-------------------------------------------------------------------------------
# Check the validity of $STGDIR for safety

if [[ ! -d "$STGDIR" ]]; then
  echo "[Error] $MYNAME: '$STGDIR' is not a directory." >&2
  exit 1
fi
if [[ ! -O "$STGDIR" ]]; then
  echo "[Error] $MYNAME: '$STGDIR' is not owned by you." >&2
  exit 1
fi

#-------------------------------------------------------------------------------
# Stage-out

if [[ -s "$STGLIST" ]]; then
  while read line; do
    stage_out_cp_sub
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    while read line; do
      stage_out_cp_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

if ((THREAD > 1)); then
  wait
fi

#===============================================================================

exit 0
