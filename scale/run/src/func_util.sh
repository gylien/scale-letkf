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



#echo "###### $DIRNAME ######"



#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi

mkdir -p $DIRNAME
res=$? && ((res != 0)) && exit $res

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME/*
res=$? && ((res != 0)) && exit $res

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



#echo "!!!!!! $DIRNAME !!!!!!"



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

rm -fr $DIRNAME
res=$? && ((res != 0)) && exit $res

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
local ARGS="$@"

progbase=$(basename $PROG)
progdir=$(dirname $PROG)

#-------------------------------------------------------------------------------

if ((MACHINE_TYPE == 1)); then

  local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  $MPIRUN -d $progdir $HOSTLIST 1 ./$progbase $ARGS
#  $MPIRUN -d $progdir $HOSTLIST 1 omplace -nt ${THREADS} ./$progbase $ARGS

elif ((MACHINE_TYPE == 2)); then

#  local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
#  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

  NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)

#echo "$MPIRUN -np $NNP -wdir $progdir ./$progbase $ARGS"
  $MPIRUN -np $NNP -wdir $progdir ./$progbase $ARGS

elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11 || MACHINE_TYPE == 12)); then

  local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"

  if ((USE_RANKDIR == 1)); then

#pwd 1>&2
#mpiexec /work/system/bin/msh "/bin/ls -lL $progdir" 1>&2
#echo "mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./${progdir}/${progbase} $ARGS" 1>&2

    mpiexec -n $(cat $vcoordfile | wc -l) ./${progdir}/${progbase} $ARGS
#    mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./${progdir}/${progbase} $ARGS
#echo "fipp -C -d Fprofd_${progbase} -Icall,hwm mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./${progdir}/${progbase} $ARGS"
#    fipp -C -d Fprofd_${progbase} -Icall,hwm mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./${progdir}/${progbase} $ARGS
#    mpiexec -n $(cat $vcoordfile | wc -l) -of-proc std-file -vcoordfile $vcoordfile ./${progdir}/${progbase} $ARGS
#    mpiexec -n $(cat $vcoordfile | wc -l) -of-proc std-file -vcoordfile $vcoordfile -mca mpi_deadlock_timeout 60 -mca mpi_deadlock_timeout_delay 20 ./${progdir}/${progbase} $ARGS
    res=$?
    if ((res != 0)); then
      echo "[Error] mpiexec -n $(cat $vcoordfile | wc -l) ./${progdir}/${progbase} $ARGS" >&2
      echo "        Exit code: $res" >&2
    fi

  else

    ( cd $progdir && mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./$progbase $ARGS )
#echo "( cd $progdir && fipp -C -d Fprofd_${progbase} -Icall,hwm mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./$progbase $ARGS )"
#    ( cd $progdir && fipp -C -d Fprofd_${progbase} -Icall,hwm mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./$progbase $ARGS )
#    ( cd $progdir && mpiexec -n $(cat $vcoordfile | wc -l) -of-proc std-file -vcoordfile $vcoordfile ./$progbase $ARGS )
#    ( cd $progdir && mpiexec -n $(cat $vcoordfile | wc -l) -of-proc std-file -vcoordfile $vcoordfile -mca mpi_deadlock_timeout 60 -mca mpi_deadlock_timeout_delay 20 ./$progbase $ARGS )
    res=$?
    if ((res != 0)); then 
      echo "[Error] mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile ./$progbase $ARGS" >&2
      echo "        Exit code: $res" >&2
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
###             alln: run the script in all nodes list in $NODEFILE, one process per node
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

if [ -x "$TMPDAT/exec/pdbash" ]; then
  pdbash_exec="$TMPDAT/exec/pdbash"
elif [ -x "$COMMON_DIR/pdbash" ]; then
  pdbash_exec="$COMMON_DIR/pdbash"
else
  echo "[Error] $FUNCNAME: Cannot find 'pdbash' program." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

if ((MACHINE_TYPE == 1)); then

  if [ "$PROC_OPT" == 'all' ]; then
    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
###  elif [ "$PROC_OPT" == 'alln' ]; then
###    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE} | sort | uniq)
  elif [ "$PROC_OPT" == 'one' ]; then
    local HOSTLIST=$(head -n 1 ${NODEFILE_DIR}/${NODEFILE})
  else
    exit 1
  fi
  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

#echo "  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS"

  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 $pdbash_exec $SCRIPT $ARGS
#  $MPIRUN -d $SCRP_DIR $HOSTLIST 1 bash $SCRIPT - $ARGS

elif ((MACHINE_TYPE == 2)); then

  if [ "$PROC_OPT" == 'all' ]; then
#    local HOSTLIST=$(cat ${NODEFILE_DIR}/${NODEFILE})
    NNP=$(cat ${NODEFILE_DIR}/${NODEFILE} | wc -l)
  elif [ "$PROC_OPT" == 'one' ]; then
#    local HOSTLIST=$(head -n 1 ${NODEFILE_DIR}/${NODEFILE})
    NNP=1
  else
    exit 1
  fi
#  HOSTLIST=$(echo $HOSTLIST | sed 's/  */,/g')

#echo "$MPIRUN -np $NNP -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS"
  $MPIRUN -np $NNP -wdir $SCRP_DIR $pdbash_exec $SCRIPT $ARGS

elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11 || MACHINE_TYPE == 12)); then

#echo 11
  if [ "$PROC_OPT" == 'all' ]; then
    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"
###  elif [ "$PROC_OPT" == 'alln' ]; then
###    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}_tmp"
###    cat ${NODEFILE_DIR}/${NODEFILE} | sort | uniq > $vcoordfile
  elif [ "$PROC_OPT" == 'one' ]; then
    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}"
###    local vcoordfile="${NODEFILE_DIR}/${NODEFILE}_tmp"
###    head -n 1 ${NODEFILE_DIR}/${NODEFILE} > $vcoordfile
  else
    exit 1
  fi

#echo 12
#echo "======"
#echo "pdbash $NODEFILE $PROC_OPT $SCRIPT $ARGS"
#echo $vcoordfile
#cat $vcoordfile
#echo "======"



  if ((USE_RANKDIR == 1)); then


#    pdbash_exec="./dat/exec/pdbash"


#pwd 1>&2
#ls -l .. 1>&2
#ls -l 1>&2
#ls -l src 1>&2
#mpiexec /work/system/bin/msh "/bin/ls -l dat/exec"
#echo "mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS" 1>&2
#cat $vcoordfile 1>&2


#    mpiexec -n $(cat $vcoordfile | wc -l) -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS
#    mpiexec -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS
#    mpiexec $pdbash_exec $SCRIPT $ARGS

    if [ "$PROC_OPT" == 'one' ]; then
      mpiexec -n 1 -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS
    else
      mpiexec -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS
    fi

  else



    if [ "$PROC_OPT" == 'one' ]; then
      ( cd $SCRP_DIR && mpiexec -n 1 -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS )
    else
      ( cd $SCRP_DIR && mpiexec -vcoordfile $vcoordfile $pdbash_exec $SCRIPT $ARGS )
    fi


  fi

#echo 13

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

history_files_for_bdy () {
#-------------------------------------------------------------------------------
# Find the corresponding history files for preparing boundary files
#
# Usage: history_files_for_bdy
#
#   TIME
#   FCSTLEN
#   PARENT_LCYCLE
#   PARENT_FOUT
#   PARENT_REF_TIME
#   ONEFILE
#
# Return variables:
#   $nfiles
#   $ntsteps
#   $ntsteps_skip
#   $history_times[1...$nfiles]
#
#  *Require source 'func_datetime' first.
#-------------------------------------------------------------------------------

if (($# < 5)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local TIME=$1; shift
local FCSTLEN=$1; shift
local PARENT_LCYCLE=$1; shift
local PARENT_FOUT=$1; shift
local PARENT_REF_TIME=$1; shift
local ONEFILE=${1:-0}

#-------------------------------------------------------------------------------

local parent_time_start=$PARENT_REF_TIME
local parent_time_start_prev=$parent_time_start
while ((parent_time_start <= TIME)); do
  parent_time_start_prev=$parent_time_start
  parent_time_start=$(datetime $parent_time_start $PARENT_LCYCLE s)
done
parent_time_start=$parent_time_start_prev

while ((parent_time_start > TIME)); do
  parent_time_start=$(datetime $parent_time_start -${PARENT_LCYCLE} s)
done

ntsteps_skip=0
local itime=$parent_time_start
while ((itime < TIME)); do
  ntsteps_skip=$((ntsteps_skip+1))
  itime=$(datetime $itime ${PARENT_FOUT} s)
done
if ((itime > TIME)); then
  echo "[Error] $FUNCNAME: Cannot not find the requested timeframe ($TIME) in history files." >&2
  exit 1
fi

local ntsteps_total=$(((FCSTLEN-1)/PARENT_FOUT+2 + ntsteps_skip))

if ((ONEFILE == 1)); then
  ntsteps=$ntsteps_total
  nfiles=1
  history_times[1]=$(datetime $parent_time_start $PARENT_LCYCLE s)
else
  ntsteps=$((PARENT_LCYCLE / PARENT_FOUT))
  nfiles=1
  history_times[1]=$(datetime $parent_time_start $PARENT_LCYCLE s)
  while ((ntsteps_total > ntsteps)); do
    nfiles=$((nfiles+1))
    history_times[$nfiles]=$(datetime ${history_times[$((nfiles-1))]} $PARENT_LCYCLE s)
    ntsteps_total=$((ntsteps_total-ntsteps))
  done
fi

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

res=$(cd $rundir && pjsub $scrpname 2>&1)
echo $res

if [ -z "$(echo $res | grep '\[ERR.\]')" ]; then
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
# Check if a K-computer job has ended.
#
# Usage: job_end_check_PJM JOBID
#
#   JOBID  Job ID monitored
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBID="$1"

#-------------------------------------------------------------------------------

while true; do
  jobnum=$(pjstat $JOBID | sed -n '2p' | awk '{print $10}')
  if [[ "$jobnum" =~ ^[0-9]+$ ]]; then
    if ((jobnum == 0)); then
      break
    fi
  fi
  sleep 5s
done

#-------------------------------------------------------------------------------
}

#===============================================================================
