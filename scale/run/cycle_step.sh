#!/bin/bash
#===============================================================================
#
#  Run one step of data assimilation cycles.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_step.sh [STEPFUNC TIME LOOP]
#
#  Use settings:
#    config.main
#    config.cycle
#    config.nml.scale_pp_topo
#    config.nml.scale_pp_landuse
#    config.nml.scale_init
#    config.nml.scale
#    config.nml.obsope
#    config.nml.letkf
#
#===============================================================================

cd "$(dirname "$0")"
myname=$(basename "$0")
myname1=${myname%.*}

#===============================================================================
# Configuration

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

. config.main
res=$? && ((res != 0)) && exit $res
. config.cycle
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_cycle.sh

#-------------------------------------------------------------------------------

setting

STEPFUNC="${1}"; shift
MYRANK="${1}"; shift
TIME="${1}"; shift
LOOP="${1}"

#===============================================================================
# Determine the distibution schemes

declare -a procs
declare -a mem2node
declare -a node
declare -a name_m
declare -a node_m

distribute_da_cycle machinefile -

#===============================================================================
# Run one step

time=$TIME
loop=$LOOP

atime=$(datetime $time $LCYCLE s)
timefmt="$(datetime_fmt ${time})"
obstime $time

#-------------------------------------------------------------------------------

$STEPFUNC
res=$? && ((res != 0)) && exit $res

#===============================================================================

exit 0
