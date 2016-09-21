#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 2)); then
  cat >&2 << EOF

[final_all_node.sh] 

Usage: $0 MYRANK TNO SCPCALL [FCST_CYCLE]

  MYRANK   My rank number (not used)
  SCPCALL  Called from which script? (fcst/cycle)
  FCST_CYCLE

EOF
  exit 1
fi

MYRANK="$1"; shift
TNO="$1"; shift
SCPCALL="$1"; shift
FCST_CYCLE="${1:-CYCLE}"

#-------------------------------------------------------------------------------

. config.${SCPCALL}
. src/func_datetime.sh
. src/func_${SCPCALL}.sh

setting

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 8-1" >&2
fi

#===============================================================================

if [ "$SCPCALL" = 'cycle' ]; then
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  while ((time <= ETIME)); do
    mv -f $TMPOUT/${time}/log $TMPOUT/${time}/log_${TNO}
    mv -f $TMPOUT/${atime}/log $TMPOUT/${atime}/log_${TNO}
    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
elif [ "$SCPCALL" = 'fcst' ]; then
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do
    for c in $(seq $FCST_CYCLE); do
      time2=$(datetime $time $((lcycles * (c-1))) s)
      if ((time2 <= ETIME)); then
        mv -f $TMPOUT/${time2}/log $TMPOUT/${time2}/log_${TNO}
      fi
    done
    time=$(datetime $time $((lcycles * FCST_CYCLE)) s)
  done
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ### 8-2" >&2
fi

#===============================================================================

exit 0
