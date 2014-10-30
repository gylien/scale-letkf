#!/bin/bash
#===============================================================================
#
#  Common utilities (using built-in 'datetime' program)
#  August 2014, Guo-Yuan Lien
#
#  *Require source 'config.all' first.
#
#===============================================================================

function safe_init_tmpdir {
#-------------------------------------------------------------------------------
# Safely initialize a temporary directory
#
# Usage: safe_init_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local DIRNAME="$1"



echo "###### $DIRNAME ######"



#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi

mkdir -p $DIRNAME
(($? != 0)) && exit $?

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

if ((MACHINE_TYPE != 10 || PREP != 0)); then 
  rm -fr $DIRNAME/*
  (($? != 0)) && exit $?
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function safe_rm_tmpdir {
#-------------------------------------------------------------------------------
# Safely remove a temporary directory
#
# Usage: safe_rm_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi
if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME
(($? != 0)) && exit $?

#-------------------------------------------------------------------------------
}

#===============================================================================

function mpirunf {
#-------------------------------------------------------------------------------
# Submit a MPI job according to nodefile
#
# Usage: mpiexec_nodefile NODEFILE RUNDIR PROG [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   RUNDIR    Working directory
#             '-': the current directory
#   PROG      Program
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#-------------------------------------------------------------------------------

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local RUNDIR="$1"; shift
local PROG="$1"; shift
local ARGS="$@"

#-------------------------------------------------------------------------------

if ((MACHINE_TYPE == 1)); then

  local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  if [ "$RUNDIR" == '-' ]; then
    $MPIRUN $HOSTLIST 1 $PROG $ARGS
  else
    $MPIRUN -d $RUNDIR $HOSTLIST 1 $PROG $ARGS
  fi

elif ((MACHINE_TYPE == 10)); then

  local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"

  if [ "$RUNDIR" == '-' ]; then
    mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $PROG $ARGS
  else
    ( cd $RUNDIR && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $PROG $ARGS )
  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function pdbash {
#-------------------------------------------------------------------------------
# Submit bash parallel scripts according to nodefile, only one process in each node
#
# Usage: pdbash NODEFILE PROC_OPT SCRIPT [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   PROC_OPT  Options of using processes
#             'all':  run the script in all processes listed in $NODEFILE
#             'alln': run the script in all nodes list in $NODEFILE, one process per node
#             'one':  run the script only in the first process and node in $NODEFILE
#   SCRIPT    Script (the working directory is set to $SCRP_DIR)
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#
#-------------------------------------------------------------------------------

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PROC_OPT="$1"; shift
local SCRIPT="$1"; shift
local ARGS="$@"

#-------------------------------------------------------------------------------

if ((MACHINE_TYPE == 1)); then

  if [ "$PROC_OPT" == 'all' ]; then
    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  elif [ "$PROC_OPT" == 'alln' ]; then
    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE} | sort | uniq)
  elif [ "$PROC_OPT" == 'one' ]; then
    local HOSTLIST=$(head -n 1 ${NODEFILE_DIR}/${NODEFILE})
  else
    exit 1
  fi
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  if [ -f "$TMPDAT/exec/pdbash" ]; then
    $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $TMPDAT/exec/pdbash $SCRIPT $ARGS
  else
    $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $COMMON_DIR/pdbash $SCRIPT $ARGS
  fi
#  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 bash $SCRIPT - $ARGS

elif ((MACHINE_TYPE == 10)); then


  ls -l .

  echo ${NODEFILE_DIR}
  ls ${NODEFILE_DIR}



  if [ "$PROC_OPT" == 'all' ]; then
    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"
  elif [ "$PROC_OPT" == 'alln' ]; then
    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}_tmp"
    cat ${NODEFILE_DIR}/${NODEFILE} | sort | uniq > $vcoordfile
  elif [ "$PROC_OPT" == 'one' ]; then
    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}_tmp"
    head -n 1 ${NODEFILE_DIR}/${NODEFILE} > $vcoordfile
  else
    exit 1
  fi

  if [ -f "$TMPDAT/exec/pdbash" ]; then
    ( cd $SCRP_DIR && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $TMPDAT/exec/pdbash $SCRIPT $ARGS )
  else
    ( cd $SCRP_DIR && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $COMMON_DIR/pdbash $SCRIPT $ARGS )
  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

#function parse_steps {
##-------------------------------------------------------------------------------
## Parse the input steps to be executed
##
## Usage: parse_steps NSTEPS [STEP]
##
##   NSTEPS  Number of steps
##   STEP    Steps of the script to be executed
##           'all': Run all steps
##           '2-3': Run steps 2 to 3
##           '2-' : Run steps after 2
##           (default: all)
##
## Return variables:
##   run_step[1...$NSTEPS]  Array of run/no run flags
##-------------------------------------------------------------------------------

#if (($# < 1)); then
#  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
#  exit 1
#fi

#local NSTEPS=$1
#local STEP="${2:-all}"

##-------------------------------------------------------------------------------

#local p1
#local p2
#if [ "$STEP" = 'all' ]; then
#  p1=1
#  p2=$NSTEPS
#elif [ $(echo $STEP | grep '-') ]; then
#  p1=$(echo $STEP | cut -d '-' -s -f1)
#  if [ -z "$p1" ]; then
#    p1=1
#  fi
#  p2=$(echo $STEP | cut -d '-' -s -f2)
#  if [ -z "$p2" ]; then
#    p2=$NSTEPS
#  fi
#else
#  p1=$STEP
#  p2=$STEP
#fi

#local p
#for p in $(seq $NSTEPS); do
#  if ((p >= p1)) && ((p <= p2)); then
#    run_step[$p]=1
#  else
#    run_step[$p]=0
#  fi
#done

##-------------------------------------------------------------------------------
#}

##===============================================================================

#function step_timing {
##-------------------------------------------------------------------------------
## Print the information and timing of a step
##
## Usage: step_timing ISTEP
##
##   ISTEP  Current step
##
## Other input variables:
##   $stepname[...]  Names of steps
##-------------------------------------------------------------------------------

#if (($# < 1)); then
#  echo "[Error] $FUNCNAME: Insufficient arguments." 1>&2
#  exit 1
#fi

#local ISTEP="$1"

##-------------------------------------------------------------------------------

#echo "[$(datetime_now)] ${stepname[$ISTEP]}" 1>&2
#echo
#printf " %2d. %-55s\n" $ISTEP "${stepname[$ISTEP]}"

##-------------------------------------------------------------------------------
#}

##===============================================================================
