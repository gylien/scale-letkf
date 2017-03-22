#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script (run on the server side)
#
#===============================================================================

if (($# < 3)); then
  cat >&2 << EOF

Usage: $0 NRANKS STGLIST STGDIR [CLEAN THREAD]

   NRANKS   Total number of nodes
   STGLIST  File of the stage-in list
   STGDIR   Directory where the files are staged into
   CLEAN    Clean the STGDIR before copying files?
            0: No
            1: Yes (default)
   THREAD   Number of parallel copying threads (default: 1)

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
STGDIR="$1"; shift
CLEAN="${1:-1}"; shift
THREAD="${1:-1}"

#-------------------------------------------------------------------------------

function stage_in_cp_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
  if [ -z "$destin" ]; then
    : # do nothing
  elif [ -z "$source" ]; then
    if [ "${destin: -1}" != '/' ]; then
      if (mkdir -p "$(dirname ${STGDIR}/${destin})"); then
        touch "${STGDIR}/${destin}" # create an empty file
      fi
    else
      mkdir -p "${STGDIR}/${destin}" # create an empty directory
    fi
  elif [ "${source:0:1}" != '/' ]; then
    echo "$(basename $0): source '${source}' is not an absolute path" >&2
  elif [ "${source: -1}" != '/' ] && [ "${destin: -1}" != '/' ]; then # files
    if [ -d "${STGDIR}/${destin}" ]; then
      echo "$(basename $0): source '${source}' is a regular file, but destination '${STGDIR}/${destin}' exists and is a directory" >&2
    elif (mkdir -p "$(dirname ${STGDIR}/${destin})"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        cp -L "${source}" "${STGDIR}/${destin}" &
      else
        cp -L "${source}" "${STGDIR}/${destin}"
      fi
    fi
  elif [ "${source: -1}" = '/' ] && [ "${destin: -1}" = '/' ]; then # directories
    if (mkdir -p "${STGDIR}/${destin}"); then
      if ((THREAD > 1)); then
        while (($(jobs -p | wc -l) >= THREAD)); do
          sleep 1s
        done
        cp -rL "${source}"* "${STGDIR}/${destin}" &
      else
        cp -rL "${source}"* "${STGDIR}/${destin}"
      fi
    fi
  else
    echo "$(basename $0): source '${source}' and destination '${STGDIR}/${destin}' need to be the same type" >&2
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
    stage_in_cp_sub
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [ -s "$STGLIST.$n" ]; then
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
