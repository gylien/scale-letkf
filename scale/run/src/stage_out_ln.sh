#!/bin/bash
#===============================================================================
#
#  Built-in stage-out script (run on a single node)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 3)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST STGDIR

   NRANKS   Total number of nodes
   STGLIST  File of the stage-out list
   STGDIR   Directory where the files are staged out from

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"

#-------------------------------------------------------------------------------

function stage_out_ln_sub () {
  local destin="$(echo $line | cut -d '|' -s -f1)"
  local source="$(echo $line | cut -d '|' -s -f2)"
  local sourcestg="${STGDIR}/${source}"
  if [[ -z "$source" || -z "$destin" ]]; then
    : # do nothing
  elif [[ "$destin" != /* ]]; then
    echo "$MYNAME: destination '$destin' is not an absolute path" >&2
  elif [[ -e "$sourcestg" || -L "$sourcestg" ]]; then
    echo "$MYNAME: cannot use links for stage-out when source '$sourcestg' exists" >&2
  elif [[ -L "$(dirname "$sourcestg")" ]]; then
    echo "$MYNAME: should not create a link inside a link of directory '$(dirname "$sourcestg")'" >&2
  elif [[ "$sourcestg" != */ && "$destin" != */ ]]; then # files
    if [[ -d "$destin" ]]; then
      echo "$MYNAME: source '$sourcestg' is a regular file, but destination '$destin' exists and is a directory" >&2
    elif (mkdir -p "$(dirname "$destin")") && (mkdir -p "$(dirname "$sourcestg")"); then
      rm -f "$destin"
      ln -s "$destin" "$sourcestg"
    fi
  elif [[ "$sourcestg" == */ && "$destin" == */ ]]; then # directories
    if (mkdir -p "$destin") && (mkdir -p "$(dirname "$sourcestg")"); then
      rm -fr "$destin"*
      ln -s "${destin%/}" "${sourcestg%/}" # remove possible trailing slashes
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
    stage_out_ln_sub
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    while read line; do
      stage_out_ln_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

#===============================================================================

exit 0
