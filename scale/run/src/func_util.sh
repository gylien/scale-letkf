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

local DIRNAME="$1"



echo "###### $DIRNAME ######"



#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." >&2
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

rm -fr $DIRNAME/*
(($? != 0)) && exit $?

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

local DIRNAME="$1"



echo "!!!!!! $DIRNAME !!!!!!"



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
# Usage: mpirunf NODEFILE RUNDIR PROG [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   RUNDIR    Working directory
#             -: the current directory
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

elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11)); then

#echo 21
  local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"

#echo 22
#echo $vcoordfile
#echo "mpirunf $NODEFILE $RUNDIR $PROG $ARGS"

  if [ "$RUNDIR" == '-' ]; then
    mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $PROG $ARGS
  else
    ( cd $RUNDIR && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $PROG $ARGS )
  fi

#echo 23

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
#             all:  run the script in all processes listed in $NODEFILE
#             alln: run the script in all nodes list in $NODEFILE, one process per node
#             one:  run the script only in the first process and node in $NODEFILE
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

elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11)); then

#echo 11
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

#echo 12
#echo "======"
#echo "pdbash $NODEFILE $PROC_OPT $SCRIPT $ARGS"
#echo $vcoordfile
#cat $vcoordfile
#echo "======"

  ( cd $SCRP_DIR && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $TMPDAT/exec/pdbash $SCRIPT $ARGS )

#echo 13

fi

#-------------------------------------------------------------------------------
}

#===============================================================================
