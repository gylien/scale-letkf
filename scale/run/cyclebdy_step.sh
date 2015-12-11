#!/bin/bash
#===============================================================================
#
#  Run one step of ensemble forecasts
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cyclebdy_step.sh [STEPFUNC MYRANK LOOP ITER]
#
#  Use settings:
#    config.main
#    config.cyclebdy
#    config.nml.scale_pp_topo
#    config.nml.scale_pp_landuse
#    config.nml.scale_init
#    config.nml.scale
#
#===============================================================================

cd "$(dirname "$0")"
#myname=$(basename "$0")
#myname1=${myname%.*}

#===============================================================================
# Configuration

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

. config.main
res=$? && ((res != 0)) && exit $res
. config.cyclebdy
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_cyclebdy.sh

#-------------------------------------------------------------------------------

setting

STEPFUNC="${1}"; shift
MYRANK="${1}"; shift
#TIME="${1}"; shift
LOOP="${1}"; shift
ITER="${1:-0}"

#===============================================================================
# Determine the distibution schemes

declare -a node
declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

distribute_fcst "$MEMBERS" $CYCLE machinefile - $NODEFILE_DIR/distr

#===============================================================================
# Run one step

loop=$LOOP

time=$(datetime $STIME $((lcycles * CYCLE * (LOOP-1))) s)
rcycle=0
for c in $(seq $CYCLE); do
  time2=$(datetime $time $((lcycles * (c-1))) s)
  if ((time2 <= ETIME)); then
    stimes[$c]=$time2
    stimesfmt[$c]="$(datetime_fmt ${stimes[$c]})"
    rcycle=$c  # The "real" number of cycles
  else
    stimes[$c]=
    stimesfmt[$c]=
  fi
done
if ((rcycle == 0)); then
  echo "[Error] Wrong loop number, LOOP = $LOOP" >&2
  exit 1
fi

iter=$ITER
if ((ITER == 0)); then
  its=1
  ite=$nitmax
else
  its=$iter
  ite=$iter
fi

#-------------------------------------------------------------------------------

#echo $STEPFUNC $MYRANK $TIME $LOOP 1>&2

$STEPFUNC
res=$? && ((res != 0)) && exit $res

#echo $STEPFUNC $MYRANK $TIME $LOOP ... done 1>&2

#===============================================================================

exit 0
