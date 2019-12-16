#!/bin/bash
#===============================================================================
#
#  Run one step of ensemble forecasts
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_step.sh [STEPFUNC MYRANK LOOP ITER]
#
#  Use settings:
#    config.main
#    config.fcst
#    config.nml.scale_pp
#    config.nml.scale_init
#    config.nml.scale
#    config.nml.scale_user
#    config.nml.grads_boundary
#
#===============================================================================

cd "$(dirname "$0")"
#myname=$(basename "$0")

#===============================================================================
# Configuration

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

. config.main
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_pp.sh

#-------------------------------------------------------------------------------

setting

STEPFUNC="${1}"; shift
MYRANK="${1}"; shift
#TIME="${1}"; shift
LOOP="${1}"; shift
ITER="${1:-0}"

######
if (( MYRANK == 0 )); then
  echo "[$(datetime_now)] pp_step: $@: start"
fi
######

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


distribute_fcst "$MEMBERS" 1 - $NODEFILE_DIR

#===============================================================================
# Run one step

#iter=1
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

######
if (( MYRANK == 0 )); then
  echo "[$(datetime_now)] pp_step: $@: end"
fi
######

exit 0
