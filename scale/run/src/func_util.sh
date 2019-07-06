#!/bin/bash
#===============================================================================
#
#  Common utilities (using built-in 'datetime' program)
#  August 2014, Guo-Yuan Lien
#
#  *Require source 'config.main' first.
#
#===============================================================================

safe_init_tmpdir () {
#-------------------------------------------------------------------------------
# Safely initialize a temporary directory
#
# Usage: safe_init_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi

mkdir -p $DIRNAME || exit $?

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME/* || exit $?

#-------------------------------------------------------------------------------
}

#===============================================================================

safe_rm_tmpdir () {
#-------------------------------------------------------------------------------
# Safely remove a temporary directory
#
# Usage: safe_rm_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi
if [ ! -e "$DIRNAME" ]; then
  return 0
fi
if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME || exit $?

#-------------------------------------------------------------------------------
}

#===============================================================================

rev_path () {
#-------------------------------------------------------------------------------
# Compose the reverse path of a path
#
# Usage: rev_path PATH
#
#   PATH  The forward path
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local path="$1"

#-------------------------------------------------------------------------------

local rpath='.'
local base
while [ "$path" != '.' ]; do
  base=$(basename $path)
  res=$? && ((res != 0)) && exit $res
  path=$(dirname $path)
  if [ "$base" = '..' ]; then
    if [ -d "$path" ]; then
      rpath="$rpath/$(basename $(cd $path && pwd))"
    else
      echo "[Error] $FUNCNAME: Error in reverse path search." 1>&2
      exit 1
    fi
  elif [ "$base" != '.' ]; then
    rpath="$rpath/.."
  fi
done
if [ ${rpath:0:2} = './' ]; then
  echo ${rpath:2}
else
  echo $rpath
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

mpirunf () {
#-------------------------------------------------------------------------------
# Submit a MPI job according to nodefile
#
# Usage: mpirunf NODEFILE PROG [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   PROG      Program
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#-------------------------------------------------------------------------------

if (($# < 2)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PROG="$1"; shift
local CONF="$1"; shift
local STDOUT="$1"; shift
local ARGS="$@"

progbase=$(basename $PROG)
progdir=$(dirname $PROG)

#-------------------------------------------------------------------------------

if [ "$MPI_TYPE" = 'sgimpt' ]; then

  local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  $MPIRUN $HOSTLIST 1 $PROG $CONF $STDOUT $ARGS
#  $MPIRUN $HOSTLIST 1 omplace -nt ${THREADS} $PROG $CONF $STDOUT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN $HOSTLIST 1 $PROG $CONF $STDOUT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'openmpi' ]; then

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

  $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} $PROG $CONF $STDOUT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} $PROG $CONF $STDOUT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'impi' ]; then

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

  $MPIRUN -n $NNP -machinefile ${NODEFILE_DIR}/${NODEFILE} $PROG $CONF $STDOUT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -n $NNP -machinefile ${NODEFILE_DIR}/${NODEFILE} $PROG $CONF $STDOUT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'K' ]; then

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

  if [ "$PRESET" = 'K_rankdir' ]; then
    mpiexec -n $NNP -of-proc $STDOUT $PROG $CONF '' $ARGS
    res=$?
    if ((res != 0)); then
      echo "[Error] mpiexec -n $NNP -of-proc $STDOUT $PROG $CONF '' $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi
  else
    mpiexec -n $NNP -vcoordfile "${NODEFILE_DIR}/${NODEFILE}" -of-proc $STDOUT $PROG $CONF '' $ARGS
    res=$?
    if ((res != 0)); then 
      echo "[Error] mpiexec -n $NNP -vcoordfile \"${NODEFILE_DIR}/${NODEFILE}\" -of-proc $STDOUT $PROG $CONF '' $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi
  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

pdbash () {
#-------------------------------------------------------------------------------
# Submit bash parallel scripts according to nodefile
#
# Usage: pdbash NODEFILE PROC_OPT SCRIPT [ARGS]
#
#   NODEFILE  Name of nodefile (omit the directory $NODEFILE_DIR)
#   PROC_OPT  Options of using processes
#             all:  run the script in all processes listed in $NODEFILE
#             one:  run the script only in the first process and node in $NODEFILE
#   SCRIPT    Script (the working directory is set to $SCRP_DIR)
#   ARGS      Arguments passed into the program
#
# Other input variables:
#   $NODEFILE_DIR  Directory of nodefiles
#-------------------------------------------------------------------------------

if (($# < 2)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"; shift
local PROC_OPT="$1"; shift
local SCRIPT="$1"; shift
local ARGS="$@"

if ((RUN_LEVEL <= 2)); then
  pdbash_exec="$COMMON_DIR/pdbash"
else
  pdbash_exec="$TMPDAT/exec/pdbash"
fi
if [ ! -x "$pdbash_exec" ]; then
  echo "[Error] $FUNCNAME: Cannot find 'pdbash' program." >&2
  exit 1
fi

if [ "$PROC_OPT" != 'all' ] && [ "$PROC_OPT" != 'one' ]; then
  echo "[Error] $FUNCNAME: \$PROC_OPT needs to be {all|one}." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

if [ "$MPI_TYPE" = 'sgimpt' ]; then

  if [ "$PROC_OPT" == 'all' ]; then
    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  elif [ "$PROC_OPT" == 'one' ]; then
    local HOSTLIST=$(head -n 1 ${NODEFILE_DIR}/${NODEFILE})
  fi
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS
#  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 bash $SCRIPT - $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'openmpi' ]; then

  if [ "$PROC_OPT" == 'all' ]; then
    NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)
  elif [ "$PROC_OPT" == 'one' ]; then
    NNP=1
  fi

  $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -np $NNP -hostfile ${NODEFILE_DIR}/${NODEFILE} -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'impi' ]; then

  if [ "$PROC_OPT" == 'all' ]; then
    NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)
  elif [ "$PROC_OPT" == 'one' ]; then
    NNP=1
  fi

  $MPIRUN -n $NNP -machinefile ${NODEFILE_DIR}/${NODEFILE} -gwdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS
  res=$?
  if ((res != 0)); then
    echo "[Error] $MPIRUN -n $NNP -machinefile ${NODEFILE_DIR}/${NODEFILE} -gwdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS" >&2
    echo "        Exit code: $res" >&2
    exit $res
  fi

elif [ "$MPI_TYPE" = 'K' ]; then

  if [ "$PROC_OPT" == 'one' ]; then
    mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS
    res=$?
    if ((res != 0)); then
      echo "[Error] mpiexec -n 1 $pdbash_exec $SCRIPT $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi
  else
    mpiexec $pdbash_exec $SCRIPT $ARGS
    res=$?
    if ((res != 0)); then
      echo "[Error] mpiexec $pdbash_exec $SCRIPT $ARGS" >&2
      echo "        Exit code: $res" >&2
      exit $res
    fi
  fi

fi

#-------------------------------------------------------------------------------
}

#===============================================================================

pdrun () {
#-------------------------------------------------------------------------------
# Return if it is the case to run parallel scripts, according to nodefile
#
# Usage: pdrun GROUP OPT
#
#   GROUP   Group of processes
#           all:     all processes
#           (group): process group #(group)
#   OPT     Options of the ways to pick up processes
#           all:  run the script in all processes in the group
#           alln: run the script in all nodes in the group, one process per node (default)
#           one:  run the script only in the first process in the group
#
# Other input variables:
#   MYRANK  The rank of the current process
#
# Exit code:
#   0: This process is used
#   1: This process is not used
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local GROUP="$1"; shift
local OPT="${1:-alln}"

#-------------------------------------------------------------------------------

local mynode=${proc2node[$((MYRANK+1))]}
if [ -z "$mynode" ]; then
  exit 1
fi

local res=1
local n

if [ "$GROUP" = 'all' ]; then

  if [ "$OPT" = 'all' ]; then
    exit 0
  elif [ "$OPT" = 'alln' ]; then
    res=0
    for n in $(seq $MYRANK); do
      if ((${proc2node[$n]} == mynode)); then
        res=1
        break
      fi
    done
  elif [ "$OPT" = 'one' ]; then
    if ((MYRANK == 0)); then
      exit 0
    fi
  fi

elif ((GROUP <= parallel_mems)); then

  local mygroup=${proc2group[$((MYRANK+1))]}
  local mygrprank=${proc2grpproc[$((MYRANK+1))]}

  if ((mygroup = GROUP)); then
    if [ "$OPT" = 'all' ]; then
      exit 0
    elif [ "$OPT" = 'alln' ]; then
      res=0
      for n in $(seq $((mygrprank-1))); do
        if ((${mem2node[$(((GROUP-1)*mem_np+n))]} == mynode)); then
          res=1
          break
        fi
      done
    elif [ "$OPT" = 'one' ]; then
      if ((${mem2node[$(((GROUP-1)*mem_np+1))]} == mynode)); then
        res=0
        for n in $(seq $((mygrprank-1))); do
          if ((${mem2node[$(((GROUP-1)*mem_np+n))]} == mynode)); then
            res=1
            break
          fi
        done
      fi
    fi
  fi

fi

exit $res

#-------------------------------------------------------------------------------
}

#===============================================================================

stage_in () {
#-------------------------------------------------------------------------------
# Stage-in files into the runtime temporary directories based on the staging lists
#
# Usage: stage_in [RUN_ON]
#
#   RUN_ON  Run on which side?
#           'server': run on the server side
#           'node':   run on the computing-node side (using 'pdbash') (default)
#
# Other input variables:
#   $SCRP_DIR
#   $TMP
#   $TMPL
#   $STAGING_DIR
#   $STGINLIST_LINK
#   $STGINLIST_SHARE
#   $STGINLIST_LOCAL
#   $STGOUTLIST_LINK
#   $NNODES
#   $STAGE_THREAD
#-------------------------------------------------------------------------------

local RUN_ON="${1:-node}"

#-------------------------------------------------------------------------------

if ((DISK_MODE == 1)); then
  if [ -s "${STAGING_DIR}/${STGOUTLIST_LINK}.1" ] || [ -s "${STAGING_DIR}/${STGOUTLIST_LINK}" ]; then
    local errmsg=$(bash $SCRP_DIR/src/stage_out_ln.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_LINK} $TMP 2>&1) # code same for both server and computing-node sides
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
      return 1
    fi
  fi
fi

if [ -s "${STAGING_DIR}/${STGINLIST_LINK}.1" ] || [ -s "${STAGING_DIR}/${STGINLIST_LINK}" ]; then
#  safe_init_tmpdir $TMP || return $?
  local errmsg=$(bash $SCRP_DIR/src/stage_in_ln.sh $NNODES ${STAGING_DIR}/${STGINLIST_LINK} $TMP 2>&1) # code same for both server and computing-node sides
  if [ -n "$errmsg" ]; then
    echo "$errmsg" >&2
    return 1
  fi
fi

if [ -s "${STAGING_DIR}/${STGINLIST_SHARE}.1" ] || [ -s "${STAGING_DIR}/${STGINLIST_SHARE}" ]; then
#  safe_init_tmpdir $TMP || return $?
  if [ "$RUN_ON" = 'server' ]; then
    local errmsg=$(bash $SCRP_DIR/src/stage_in_cp.sh $NNODES ${STAGING_DIR}/${STGINLIST_SHARE} $TMP $STAGE_THREAD 2>&1)
  else
    local errmsg=$(pdbash node all $SCRP_DIR/src/stage_in_cp_node.sh $NNODES ${STAGING_DIR}/${STGINLIST_SHARE} $TMP share $STAGE_THREAD 2>&1)
  fi
  if [ -n "$errmsg" ]; then
    echo "$errmsg" >&2
    return 1
  fi
fi

if [ "$RUN_ON" = 'node' ]; then # stage-in to local directories can only be done on the computing-node side
  if [ -s "${STAGING_DIR}/${STGINLIST_LOCAL}.1" ] || [ -s "${STAGING_DIR}/${STGINLIST_LOCAL}" ]; then
    pdbash node all $SCRP_DIR/src/stage_in_init_stgdir_node.sh $TMPL local || return $?
    local errmsg=$(pdbash node all $SCRP_DIR/src/stage_in_cp_node.sh $NNODES ${STAGING_DIR}/${STGINLIST_LOCAL} $TMPL local $STAGE_THREAD 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
      return 1
    fi
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

stage_out () {
#-------------------------------------------------------------------------------
# Stage-out files from the runtime temporary directories based on the staging lists
#
# Usage: stage_out [RUN_ON STEP]
#
#   RUN_ON  Run on which side?
#           'server': run on the server side
#           'node':   run on the computing-node side (using 'pdbash') (default)
#   STEP    Step ID with which the files are processed
#           'a': process all steps (default)
#
# Other input variables:
#   $SCRP_DIR
#   $TMP
#   $TMPL
#   $STAGING_DIR
#   $STGOUTLIST_SHARE
#   $STGOUTLIST_LOCAL
#   $NNODES
#   $STAGE_THREAD
#-------------------------------------------------------------------------------

local RUN_ON="${1:-node}"; shift
local STEP="${1:-a}"

#-------------------------------------------------------------------------------

if [ -s "${STAGING_DIR}/${STGOUTLIST_SHARE}.1" ] || [ -s "${STAGING_DIR}/${STGOUTLIST_SHARE}" ]; then
  if [ "$RUN_ON" = 'server' ]; then
    errmsg=$(bash $SCRP_DIR/src/stage_out_cp.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_SHARE} $TMP $STAGE_THREAD $STEP 2>&1)
  else
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_SHARE} $TMP share $STAGE_THREAD $STEP 2>&1)
  fi
  if [ -n "$errmsg" ]; then
    echo "$errmsg" >&2
#    return 1
  fi
fi

if [ "$RUN_ON" = 'node' ]; then # stage-out from local directories can only be done on the computing-node side
  if [ -s "${STAGING_DIR}/${STGOUTLIST_LOCAL}.1" ] || [ -s "${STAGING_DIR}/${STGOUTLIST_LOCAL}" ]; then
    errmsg=$(pdbash node all $SCRP_DIR/src/stage_out_cp_node.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_LOCAL} $TMPL local $STAGE_THREAD $STEP 2>&1)
    if [ -n "$errmsg" ]; then
      echo "$errmsg" >&2
#      return 1
    fi
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

stage_K_inout () {
#-------------------------------------------------------------------------------
# Print stage-in/out scripts for K-computer jobs based on the staging lists
#
# Usage: stage_out [USE_RANKDIR]
#
#   USE_RANKDIR  Whether enable the rank-directory?
#                0: No (default)
#                1: Yes
#
# Other input variables:
#   $SCRP_DIR
#   $STAGING_DIR
#   $STGINLIST_SHARE
#   $STGINLIST_LOCAL
#   $STGOUTLIST_SHARE
#   $STGOUTLIST_LOCAL
#   $NNODES
#   $jobscrp
#-------------------------------------------------------------------------------

USE_RANKDIR="${1:-0}"

#-------------------------------------------------------------------------------

if [ -s "${STAGING_DIR}/${STGINLIST_SHARE}.1" ] || [ -s "${STAGING_DIR}/${STGINLIST_SHARE}" ]; then
  bash $SCRP_DIR/src/stage_in_K.sh $NNODES ${STAGING_DIR}/${STGINLIST_SHARE} $USE_RANKDIR share $TMPS 1>> $jobscrp || exit $?
fi
if [ -s "${STAGING_DIR}/${STGINLIST_LOCAL}.1" ] || [ -s "${STAGING_DIR}/${STGINLIST_LOCAL}" ]; then
  bash $SCRP_DIR/src/stage_in_K.sh $NNODES ${STAGING_DIR}/${STGINLIST_LOCAL} $USE_RANKDIR local $TMPS 1>> $jobscrp || exit $?
fi
if [ -s "${STAGING_DIR}/${STGOUTLIST_SHARE}.1" ] || [ -s "${STAGING_DIR}/${STGOUTLIST_SHARE}" ]; then
  bash $SCRP_DIR/src/stage_out_K.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_SHARE} $USE_RANKDIR share 1>> $jobscrp || exit $?
fi
if [ -s "${STAGING_DIR}/${STGOUTLIST_LOCAL}.1" ] || [ -s "${STAGING_DIR}/${STGOUTLIST_LOCAL}" ]; then
  bash $SCRP_DIR/src/stage_out_K.sh $NNODES ${STAGING_DIR}/${STGOUTLIST_LOCAL} $USE_RANKDIR local 1>> $jobscrp || exit $?
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

bdy_setting () {
#-------------------------------------------------------------------------------
# Calculate scale_init namelist settings for boundary files
#
# Usage: bdy_setting TIME FCSTLEN PARENT_LCYCLE [PARENT_FOUT] [PARENT_REF_TIME] [SINGLE_FILE]
#
#   TIME
#   FCSTLEN
#   PARENT_LCYCLE
#   PARENT_FOUT
#   PARENT_REF_TIME
#   SINGLE_FILE
#
# Return variables:
#   $nbdy
#   $ntsteps
#   $ntsteps_skip
#   $bdy_times[1...$nbdy]
#   $bdy_start_time
#   $parent_start_time
#
#  *Require source 'func_datetime' first.
#-------------------------------------------------------------------------------

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local TIME=$(datetime $1); shift
local FCSTLEN=$1; shift
local PARENT_LCYCLE=$1; shift
local PARENT_FOUT=${1:-$PARENT_LCYCLE}; shift
local PARENT_REF_TIME=${1:-$TIME}; shift
local SINGLE_FILE=${1:-0}

PARENT_REF_TIME=$(datetime $PARENT_REF_TIME)

#-------------------------------------------------------------------------------
# compute $parent_start_time based on $PARENT_REF_TIME and $PARENT_LCYCLE

parent_start_time=$PARENT_REF_TIME
local parent_start_time_prev=$parent_start_time
while ((parent_start_time <= TIME)); do
  parent_start_time_prev=$parent_start_time
  parent_start_time=$(datetime $parent_start_time $PARENT_LCYCLE s)
done
parent_start_time=$parent_start_time_prev

while ((parent_start_time > TIME)); do
  parent_start_time=$(datetime $parent_start_time -${PARENT_LCYCLE} s)
done

#-------------------------------------------------------------------------------
# compute $bdy_start_time, $ntsteps_skip, and $ntsteps_total based on $parent_start_time and $PARENT_FOUT
# (assume $bdy_start_time <= $TIME)

ntsteps_skip=-1
bdy_start_time=$parent_start_time
while ((bdy_start_time <= TIME)); do
  bdy_start_time_prev=$bdy_start_time
  bdy_start_time=$(datetime $bdy_start_time $PARENT_FOUT s)
  ntsteps_skip=$((ntsteps_skip+1))
done
bdy_start_time=$bdy_start_time_prev

local ntsteps_total=$(((FCSTLEN-1)/PARENT_FOUT+2 + ntsteps_skip))

if ((bdy_start_time != TIME)); then
  if (($(datetime $bdy_start_time $(((ntsteps_total-ntsteps_skip-1)*PARENT_FOUT)) s) < $(datetime $TIME $FCSTLEN s))); then
    ntsteps_total=$((ntsteps_total+1))
  fi
fi

#-------------------------------------------------------------------------------
# compute $ntsteps

if ((PARENT_LCYCLE % PARENT_FOUT != 0)); then
  echo "[Error] $FUNCNAME: $PARENT_LCYCLE needs to be an exact multiple of $PARENT_FOUT." >&2
  exit 1
fi

if ((SINGLE_FILE == 1)); then
  ntsteps=$ntsteps_total
else
  ntsteps=$((PARENT_LCYCLE / PARENT_FOUT))
fi

#-------------------------------------------------------------------------------
# compute $nbdy and $bdy_times[1...$nbdy]

nbdy=1
bdy_times[1]=$parent_start_time
while ((ntsteps_total > ntsteps)); do
  nbdy=$((nbdy+1))
  bdy_times[$nbdy]=$(datetime ${bdy_times[$((nbdy-1))]} $PARENT_LCYCLE s)
  ntsteps_total=$((ntsteps_total-ntsteps))
done

if ((nbdy == 1)); then
  ntsteps=$ntsteps_total
fi

#echo "\$nbdy              = $nbdy" >&2
#echo "\$ntsteps           = $ntsteps" >&2
#echo "\$ntsteps_skip      = $ntsteps_skip" >&2
#echo "\$ntsteps_total     = $ntsteps_total" >&2
#echo "\$bdy_start_time    = $bdy_start_time" >&2
#echo "\$parent_start_time = $parent_start_time" >&2

#-------------------------------------------------------------------------------
}

#===============================================================================

job_submit_PJM () {
#-------------------------------------------------------------------------------
# Submit a PJM job.
#
# Usage: job_submit_PJM
#
#   JOBSCRP  Job script
#
# Return variables:
#   $jobid  Job ID monitered
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBSCRP="$1"

local rundir=$(dirname $JOBSCRP)
local scrpname=$(basename $JOBSCRP)

#-------------------------------------------------------------------------------

if [ "$PRESET" = 'K_micro' ] ; then 
  res=$(cd $rundir && pjsub $scrpname -g $(id -ng)s 2>&1)
else
  res=$(cd $rundir && pjsub $scrpname 2>&1)
fi
echo $res

if [ -z "$(echo $res | grep 'ERR')" ]; then
  jobid=$(echo $res | grep 'submitted' | cut -d ' ' -f 6)
  if [ -z "$jobid" ]; then
    echo "[Error] $FUNCNAME: Error found when submitting a job." >&2
    exit 1
  fi
else
  echo "[Error] $FUNCNAME: Error found when submitting a job." >&2
  exit 1
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

job_end_check_PJM () {
#-------------------------------------------------------------------------------
# Check if a PJM job has ended.
#
# Usage: job_end_check_PJM JOBID
#
#   JOBID  Job ID monitored
#
# Return variables:
#   $jobstat    Job status
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBID="$1"

#-------------------------------------------------------------------------------

local res=0
local tmp
while true; do
  tmp=$(pjstat ${JOBID} | tail -n 1)
  if [ -z "$(echo $tmp | grep ${JOBID})" ]; then
    break
  fi
  sleep 10s
done

tmp=$(pjstat -H ${JOBID} | tail -n 1)
if [ -z "$(echo $tmp | grep ${JOBID})" ]; then
  echo "[Error] $FUNCNAME: Cannot find PJM job ${JOBID}." >&2
  return 99
else
  jobstat=$(echo $tmp | cut -d ' ' -f3)
  if [ "$jobstat" = 'REJECT' ] || [ "$jobstat" = 'CANCEL' ]; then
    res=98
  elif [ "$jobstat" = 'ERROR' ]; then
    res=97
  fi  
fi

if ((res != 0)); then
  echo "[Error] $FUNCNAME: PJM job $JOBID ended with errors." >&2
  echo "        status      = $jobstat" >&2
  return $res
fi
return 0

#-------------------------------------------------------------------------------
}

#===============================================================================

job_end_check_PJM_K () {
#-------------------------------------------------------------------------------
# Check if a PJM job has ended (specialized for the K computer).
#
# Usage: job_end_check_PJM_K JOBID
#
#   JOBID  Job ID monitored
#
# Return variables:
#   $jobstat    Job status
#   $jobec      Job exit code
#   $jobreason  Job exit reason
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBID="$1"

#-------------------------------------------------------------------------------

local res=0
local tmp
while true; do
  tmp=$(pjstat -H day=5 --choose ST,EC,REASON ${JOBID} | tail -n 1)
  if [ -z "$tmp" ]; then
    echo "[Error] $FUNCNAME: Cannot find PJM job ${JOBID}." >&2
    return 99
  fi

  jobstat=$(echo $tmp | cut -d ' ' -f1)
  jobec=$(echo $tmp | cut -d ' ' -f2)
  jobreason=$(echo $tmp | cut -d ' ' -f3-)

  if [ "$jobstat" = 'RJT' ] || [ "$jobstat" = 'CCL' ]; then
    res=98
    break
  elif [ "$jobstat" = 'EXT' ]; then
    if [ "$jobreason" != '-' ]; then
      res=97
    elif ((jobec != 0)); then
      res=$jobec
    fi
    break
  fi
  sleep 30s
done

if ((res != 0)); then
  echo "[Error] $FUNCNAME: PJM job $JOBID ended with errors." >&2
  echo "        status      = $jobstat" >&2
  echo "        exit code   = $jobec" >&2
  echo "        exit reason = $jobreason" >&2
  return $res
fi
return 0

#-------------------------------------------------------------------------------
}

#===============================================================================

job_submit_torque () {
#-------------------------------------------------------------------------------
# Submit a torque job.
#
# Usage: job_submit_torque
#
#   JOBSCRP  Job script
#
# Return variables:
#   $jobid  Job ID monitered
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBSCRP="$1"

local rundir=$(dirname $JOBSCRP)
local scrpname=$(basename $JOBSCRP)

#-------------------------------------------------------------------------------

res=$(cd $rundir && qsub $scrpname 2>&1)
jobid=$(echo $res | cut -d '.' -f1)

if ! [[ "$jobid" =~ ^[0-9]+$ ]] ; then
  jobid=
  echo "[Error] $FUNCNAME: Error found when submitting a job." >&2
  exit 1
fi

echo "qsub Job $jobid submitted."

#-------------------------------------------------------------------------------
}

#===============================================================================

job_end_check_torque () {
#-------------------------------------------------------------------------------
# Check if a torque job has ended.
#
# Usage: job_end_check_torque JOBID
#
#   JOBID  Job ID monitored
#
# * Do not support exit code yet
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBID="$1"

#-------------------------------------------------------------------------------

local res=0
local tmp
while true; do
  tmp=$(qstat ${JOBID} 2> /dev/null)
  if (($? != 0)); then
    break
  fi
  sleep 5s
done

return 0

#-------------------------------------------------------------------------------
}

#===============================================================================

backup_exp_setting () {
#-------------------------------------------------------------------------------
# Backup experimental settings
#
# Usage: backup_exp_setting JOBNAME JOB_DIR JOB_ID JOB_LOG_PREFIX JOB_LOG_TYPES [JOB_LOG_TYPE_WAIT]
#
#   JOBNAME
#   JOB_DIR
#   JOB_ID
#   JOB_LOG_PREFIX
#   JOB_LOG_TYPES
#   JOB_LOG_TYPE_WAIT
#
# Other input variables:
#   $OUTDIR
#   $SCRP_DIR
#   $SCALEDIR
#   $STIME
#-------------------------------------------------------------------------------

if (($# < 5)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBNAME="$1"; shift
local JOB_DIR="$1"; shift
local JOB_ID="$1"; shift
local JOB_LOG_PREFIX="$1"; shift
local JOB_LOG_TYPES="$1"; shift
local JOB_LOG_TYPE_WAIT="$1"

#-------------------------------------------------------------------------------

if [ -n "$JOB_LOG_TYPE_WAIT" ]; then
  n=0
  nmax=30
  while [ ! -s "$JOB_DIR/${JOB_LOG_PREFIX}.${JOB_LOG_TYPE_WAIT}${JOB_ID}" ] && ((n < nmax)); do
    n=$((n+1))
    sleep 2s
  done
fi

mkdir -p $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}
cp -f $SCRP_DIR/config.main $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}
cp -f $SCRP_DIR/config.${JOBNAME} $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}
cp -f $SCRP_DIR/config.nml.* $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}
cp -f $JOB_DIR/${JOBNAME}_job.sh $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}
for p in ${JOB_LOG_TYPES}; do
  if [ -f "$JOB_DIR/${JOB_LOG_PREFIX}.${p}${JOB_ID}" ]; then
    cp -f $JOB_DIR/${JOB_LOG_PREFIX}.${p}${JOB_ID} $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}/job.${p}
  fi
done

( cd $SCRP_DIR && git log -1 --format="SCALE-LETKF version %h (%ai)" > $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}/version )
( cd $SCALEDIR/scale-rm && git log -1 --format="SCALE       version %h (%ai)" >> $OUTDIR/exp/${JOB_ID}_${JOBNAME}_${STIME}/version )

return 0

#-------------------------------------------------------------------------------
}

#===============================================================================

scale_filename_sfx () {
#-------------------------------------------------------------------------------
# Return the suffix of SCALE file names for either split-file NetCDF or PnetCDF formats
#
# Usage: scale_filename_sfx [PE]
#
#   PE  Process number
#
# Other input variables:
#   $PNETCDF
#-------------------------------------------------------------------------------

local PE="${1:-0}"

#-------------------------------------------------------------------------------

if ((PNETCDF == 1)); then
  echo '.nc'
else
  printf $SCALE_SFX $PE
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

scale_filename_bdy_sfx () {
#-------------------------------------------------------------------------------
# Return the suffix of SCALE boundary file names (offline nesting) 
# for either split-file NetCDF or PnetCDF formats
#
# Usage: scale_filename_bdy_sfx [PE]
#
#   PE  Process number
#
# Other input variables:
#   $PNETCDF_BDY_SCALE
#-------------------------------------------------------------------------------

local PE="${1:-0}"

#-------------------------------------------------------------------------------

if ((PNETCDF_BDY_SCALE == 1)); then
  echo '.nc'
else
  printf $SCALE_SFX $PE
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
