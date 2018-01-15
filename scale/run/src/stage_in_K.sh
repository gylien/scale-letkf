#!/bin/bash
#===============================================================================
#
#  Create stage-in scripts for the K-computer
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 2)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST [USE_RANKDIR TYPE TMPS]

   NRANKS       Total number of nodes
   STGLIST      File of the stage-in list
   USE_RANKDIR  Whether enable the rank-directory?
                0: No (default)
                1: Yes
   TYPE         'share': Stage-in to a shared directory (default)
                'local': Stage-in to local directories
   TMPS         Temporary directory for '.empty' file

EOF
  exit 1
fi

NRANKS="$1"; shift
STGLIST="$1"; shift
USE_RANKDIR="${1:-0}"; shift
TYPE="${1:-share}"; shift
TMPS="${1:-$(dirname $STGLIST)}"

if ((USE_RANKDIR == 0)) && [[ "$TYPE" = 'local' ]]; then
  echo "[Error] $MYNAME: When the rank-directory is not enabled (\$USE_RANKDIR = 0), \$TYPE cannot be 'local'" >&2
  exit 1
fi
max_depth=10 # maximum depth for directory stage-in

#-------------------------------------------------------------------------------

function stage_in_K_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
  local range="$(echo $line | cut -d '|' -s -f3)"
  local rank_s
  local rank_e
  if [[ -z "$destin" ]]; then
    : # do nothing
  else
    if [[ -z "$source" ]]; then
      if [[ "$destin" != */ ]]; then # stage-in an empty file
        source="$emptyfile"
      else # stage-in an empty file to create a directory
        source="$emptyfile"
        destin="${destin}."
      fi
    fi
    if [[ "$source" != /* && "$source" != "$emptyfile" ]]; then
      echo "$MYNAME: source '$source' is not an absolute path" >&2
      exit 1
    elif [[ "$source" != */ && "$destin" != */ ]]; then # files
      if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
        echo "#PJM --stgin \"rank=$((i % NRANKS)) $source $((i % NRANKS)):../${destin}\""
      elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
        if ((n == 0)); then
          echo "#PJM --stgin \"rank=* $source %r:./${destin}\""
        elif [[ -n "$range" ]]; then
          if [[ "$range" =~ ^[0-9]+-[0-9]+$ ]]; then
            rank_s="$(echo $range | cut -d '-' -s -f1)"
            rank_e="$(echo $range | cut -d '-' -s -f2)"
            echo "#PJM --stgin \"rank=$((rank_s-1))-$((rank_e-1)) $source %r:./${destin}\""
          fi
        else
          echo "#PJM --stgin \"rank=${i} $source ${i}:./${destin}\""
        fi
      else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
        echo "#PJM --stgin \"$source ./${destin}\""
      fi
    elif [[ "$source" == */ && "$destin" == */ ]]; then # directories
      if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
        echo "#PJM --stgin-dir \"rank=$((i % NRANKS)) ${source%/} $((i % NRANKS)):../${destin%/} recursive=${max_depth}\"" # remove trailing slashes
      elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
        if ((n == 0)); then
          echo "#PJM --stgin-dir \"rank=* ${source%/} %r:./${destin%/} recursive=${max_depth}\""
        elif [[ -n "$range" ]]; then
          if [[ "$range" =~ ^[0-9]+-[0-9]+$ ]]; then
            rank_s="$(echo $range | cut -d '-' -s -f1)"
            rank_e="$(echo $range | cut -d '-' -s -f2)"
            echo "#PJM --stgin-dir \"rank=$((rank_s-1))-$((rank_e-1)) ${source%/} ${i}:./${destin%/} recursive=${max_depth}\""
          fi
        else
          echo "#PJM --stgin-dir \"rank=${i} ${source%/} ${i}:./${destin%/} recursive=${max_depth}\""
        fi
      else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
        echo "#PJM --stgin-dir \"${source%/} ./${destin%/} recursive=${max_depth}\""
      fi
    else
      echo "$MYNAME: source '$source' and destination '$destin' need to be the same type" >&2
      exit 1
    fi
  fi
}

#-------------------------------------------------------------------------------
# Stage-in

emptyfile="$TMPS/.empty"
touch $emptyfile

if [[ -s "$STGLIST" ]]; then
  n=0
  i=0
  while read line; do
    stage_in_K_sub
    i=$((i+1))
  done < "$STGLIST" # | sort | uniq
fi

for n in $(seq $NRANKS); do
  if [[ -s "$STGLIST.$n" ]]; then
    i=$((n-1))
    while read line; do
      stage_in_K_sub
    done < "$STGLIST.$n" # | sort | uniq
  fi
done

#===============================================================================

exit 0
