#!/bin/bash
#===============================================================================
#
#  Create stage-out scripts for the K-computer
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 2)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST [USE_RANKDIR TYPE]

   NRANKS       Total number of nodes
   STGLIST      File of the stage-out list
   USE_RANKDIR  Whether enable the rank-directory?
                0: No (default)
                1: Yes
   TYPE         'share': Stage-in to a shared directory (default)
                'local': Stage-in to local directories

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
USE_RANKDIR="${1:-0}"; shift
TYPE="${1:-share}"

if ((USE_RANKDIR == 0)) && [[ "$TYPE" = 'local' ]]; then
  echo "[Error] $MYNAME: When the rank-directory is not enabled (\$USE_RANKDIR = 0), \$TYPE cannot be 'local'" >&2
  exit 1
fi
max_depth=10 # maximum depth for directory stage-out

#-------------------------------------------------------------------------------

function stage_out_K_sub () {
  local destin="$(echo $line | cut -d '|' -s -f1)"
  local source="$(echo $line | cut -d '|' -s -f2)"
  if [[ -z "$source" || -z "$destin" ]]; then
    : # do nothing
  elif [[ "$destin" != /* ]]; then
    echo "$MYNAME: destination '$destin' is not an absolute path" >&2
    exit 1
  elif [[ "$source" != */ && "$destin" != */ ]]; then # files
    if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
      echo "#PJM --stgout \"rank=$((i % NRANKS)) $((i % NRANKS)):../${source} $destin\""
    elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
      if ((n == 0)); then
        echo "#PJM --stgout \"rank=* %r:./${source} $destin\""
      else
        echo "#PJM --stgout \"rank=${i} ${i}:./${source} $destin\""
      fi
    else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
      echo "#PJM --stgout \"./${source} $destin\""
    fi
  elif [[ "$source" == */ && "$destin" == */ ]]; then # directories
    if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
      echo "#PJM --stgout-dir \"rank=$((i % NRANKS)) $((i % NRANKS)):../${source%/} ${destin%/} recursive=${max_depth}\"" # remove trailing slashes
    elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
      if ((n == 0)); then
        echo "#PJM --stgout-dir \"rank=* %r:./${source%/} ${destin%/} recursive=${max_depth}\""
      else
        echo "#PJM --stgout-dir \"rank=${i} ${i}:./${source%/} ${destin%/} recursive=${max_depth}\""
      fi
    else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
      echo "#PJM --stgout-dir \"./${source%/} ${destin%/} recursive=${max_depth}\""
    fi
  else
    echo "$MYNAME: source '$source' and destination '$destin' need to be the same type" >&2
    exit 1
  fi
}

#-------------------------------------------------------------------------------
# Stage-out

if [[ -s "$STGLIST" ]]; then
  n=0
  i=0
  while read line; do
    stage_out_K_sub
    i=$((i+1))
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    i=$((n-1))
    while read line; do
      stage_out_K_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

#===============================================================================

exit 0
