#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run on a single node)
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 3)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST STGDIR [TYPE]

   NRANKS   Total number of nodes
   STGLIST  File of the stage-in list
   STGDIR   Directory where the files are staged into
   TYPE     ''
            'local_rankdir': Stage-in to shared directory, but under rank directories

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
TYPE="$1"

#-------------------------------------------------------------------------------

function stage_in_ln_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
  if [[ "$TYPE" == 'local_rankdir' ]]; then
#    if ((n == 0)); then
#      local destinstg="${STGDIR}/all/${destin}"
#    else
      local destinstg="${STGDIR}/$((n-1))/${destin}"
#    fi
  else
    local destinstg="${STGDIR}/${destin}"
  fi
  if [[ -z "$destin" ]]; then
    : # do nothing
  elif [[ -z "$source" ]]; then
    if [[ "$destinstg" != */ ]]; then
      echo "$MYNAME: should not create an empty file '$destinstg'" >&2
    else
      mkdir -p "$destinstg" # create an empty directory
    fi
  elif [[ "$source" != /* ]]; then
    echo "$MYNAME: source '$source' is not an absolute path" >&2
  elif [[ ! -e "$source" ]]; then
    echo "$MYNAME: source '$source' does not exists" >&2
  elif [[ -e "$destinstg" ]]; then
    echo "$MYNAME: destination '$destinstg' exists" >&2
  elif [[ -L "$(dirname "$destinstg")" ]]; then
    echo "$MYNAME: should not create a link inside a link of directory '$(dirname "$destinstg")'" >&2
  elif [[ ("$source" != */ && "$destinstg" != */) || ("$source" == */ && "$destinstg" == */) ]]; then
    if (mkdir -p "$(dirname "$destinstg")"); then
      ln -s "${source%/}" "${destinstg%/}" # remove possible trailing slashes
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

#n=0
if [[ "$TYPE" != 'local_rankdir' ]]; then
if [[ -s "$STGLIST" ]]; then
  while read line; do
    stage_in_ln_sub
  done < "$STGLIST" # | sort | uniq
fi
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    while read line; do
      stage_in_ln_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

#===============================================================================

exit 0
