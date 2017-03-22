#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run on the server side)
#
#===============================================================================

if (($# < 3)); then
  cat >&2 << EOF

Usage: $0 NRANKS STGLIST STGDIR [CLEAN]

   NRANKS   Total number of nodes
   STGLIST  File of the stage-in list
   STGDIR   Directory where the files are staged into
   CLEAN    Clean the STGDIR before copying files?
            0: No
            1: Yes (default)

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
CLEAN="${1:-1}"

#-------------------------------------------------------------------------------

function stage_in_ln_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
  if [ -z "$destin" ]; then
    : # do nothing
  elif [ -z "$source" ]; then
    if [ "${destin: -1}" = '/' ]; then
      mkdir -p "${STGDIR}/${destin}" # create an empty directory
    else
      echo "$(basename $0): should not create an empty file '${STGDIR}/${destin}'" >&2
    fi
  elif [ "${source:0:1}" != '/' ]; then
    echo "$(basename $0): source '${source}' is not an absolute path" >&2
  elif [ ! -e "$source" ]; then
    echo "$(basename $0): source '${source}' does not exists" >&2
  elif [ -e "${STGDIR}/${destin}" ]; then
    echo "$(basename $0): destination '${STGDIR}/${destin}' exists" >&2
  elif [ -L "$(dirname ${STGDIR}/${destin})" ]; then
    echo "$(basename $0): should not create a link inside a link of directory '$(dirname ${STGDIR}/${destin})'" >&2
  elif (mkdir -p "$(dirname ${STGDIR}/${destin})"); then
    ln -s "$(dirname $source)/$(basename $source)" "${STGDIR}/$(dirname $destin)/$(basename $destin)" # remove possible trailing slashes
  fi
}

#-------------------------------------------------------------------------------
# Check the validity of $STGDIR for safety

mkdir -p $STGDIR || exit $?

if [ ! -d "$STGDIR" ]; then
  echo "[Error] $0: '$STGDIR' is not a directory." >&2
  exit 1
fi
if [ ! -O "$STGDIR" ]; then
  echo "[Error] $0: '$STGDIR' is not owned by you." >&2
  exit 1
fi

#-------------------------------------------------------------------------------
# Clean $STGDIR

if ((CLEAN == 1)); then
  rm -fr $STGDIR/* || exit $?
fi

#-------------------------------------------------------------------------------
# Stage-in

if [ -s "$STGLIST" ]; then
  while read line; do
    stage_in_ln_sub
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [ -s "$STGLIST.$n" ]; then
    while read line; do
      stage_in_ln_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

#===============================================================================

exit 0
