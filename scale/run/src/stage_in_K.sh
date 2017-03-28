#!/bin/bash
#===============================================================================
#
#  Create stage-in scripts for the K-computer
#
#===============================================================================

MYNAME=$(basename $0)

if (($# < 2)); then
  cat >&2 << EOF

Usage: $MYNAME NRANKS STGLIST [USE_RANKDIR TYPE]

   NRANKS       Total number of nodes
   STGLIST      File of the stage-in list
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
max_depth=10 # maximum depth for directory stage-in

#-------------------------------------------------------------------------------

function stage_in_K_sub () {
  local source="$(echo $line | cut -d '|' -s -f1)"
  local destin="$(echo $line | cut -d '|' -s -f2)"
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
    elif [[ "$source" != */ && "$destin" != */ ]]; then # files
      if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
        echo "#PJM --stgin \"rank=$((i % NRANKS)) $source $((i % NRANKS)):../${destin}\""
      elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
        if ((n == 0)); then
          echo "#PJM --stgin \"rank=* $source %r:./${destin}\""
        else
          echo "#PJM --stgin \"rank=${i} $source ${i}:./${destin}\""
        fi
      else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
        echo "#PJM --stgin \"$source ./${destin}\""
      fi
    elif [[ "$source" == */ && "$destin" == */ ]]; then # directories
      # remove trailing slashes
      if ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'share' ]]; then
        echo "#PJM --stgin-dir \"rank=$((i % NRANKS)) $(dirname "$source")/$(basename "$source") $((i % NRANKS)):../$(dirname "$destin")/$(basename "$destin") recursive=${max_depth}\""
      elif ((USE_RANKDIR == 1)) && [[ "$TYPE" == 'local' ]]; then
        if ((n == 0)); then
          echo "#PJM --stgin-dir \"rank=* $(dirname "$source")/$(basename "$source") %r:./$(dirname "$destin")/$(basename "$destin") recursive=${max_depth}\""
        else
          echo "#PJM --stgin-dir \"rank=${i} $(dirname "$source")/$(basename "$source") ${i}:./$(dirname "$destin")/$(basename "$destin") recursive=${max_depth}\""
        fi
      else # ((USE_RANKDIR == 0)) && [[ "$TYPE" == 'share' ]]
        echo "#PJM --stgin-dir \"$(dirname "$source")/$(basename "$source") ./$(dirname "$destin")/$(basename "$destin") recursive=${max_depth}\""
      fi
    else
      echo "$MYNAME: source '$source' and destination '$destin' need to be the same type" >&2
    fi
  fi
}

#-------------------------------------------------------------------------------
# Stage-in

emptyfile="$(dirname $STGLIST)/.empty"
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

##-------------------------------------------------------------------------------
## stage-in: nodefiles

#if [ "$STG_TYPE" = 'K_rankdir' ]; then
#  echo "#PJM --stgin-dir \"rank=* $TMPS/node %r:./node\""
##  echo "#PJM --stgin-dir \"rank=0 $TMPS/node 0:./node\""
#else
#  echo "#PJM --stgin-dir \"$TMPS/node ./node\""
#fi

##-------------------------------------------------------------------------------
## stage-in: scripts

#if [ "$STG_TYPE" = 'K_rankdir' ]; then
#  echo "#PJM --stgin \"rank=* $TMPS/config.main %r:./config.main\""
#  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.rc %r:./config.rc\""
#  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.${PROGNAME} %r:./config.${PROGNAME}\""
#  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}.sh %r:./${PROGNAME}.sh\""
#  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}_step.sh %r:./${PROGNAME}_step.sh\""
#  echo "#PJM --stgin \"rank=* $SCRP_DIR/src/* %r:./src/\""
#else
#  echo "#PJM --stgin \"$TMPS/config.main ./config.main\""
#  echo "#PJM --stgin \"$SCRP_DIR/config.rc ./config.rc\""
#  echo "#PJM --stgin \"$SCRP_DIR/config.${PROGNAME} ./config.${PROGNAME}\""
#  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}.sh ./${PROGNAME}.sh\""
#  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}_step.sh ./${PROGNAME}_step.sh\""
#  echo "#PJM --stgin \"$SCRP_DIR/src/* ./src/\""
#fi

#===============================================================================

exit 0
