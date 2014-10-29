#!/bin/bash
#===============================================================================
#
#  Adaptively distribute members on nodes
#  April   2013,          Guo-Yuan Lien
#  August  2014, modified Guo-Yuan Lien
#
#  *Require source 'config.all' first.
#
#===============================================================================

function set_mem_np {
#-------------------------------------------------------------------------------
# Set up numbers of nodes and processes for a member.
#
# Usage: set_mem_np [MEM] [MIN_NP] [MAX_NP]
#
#   MEM     Number of members
#           (default: $MEMBER)
#   MIN_NP  Minimum number of processes for a member
#           (default: 1)
#   MAX_NP  Maximum number of processes for a member
#           (default: infinity)
#
# Other input variables:
#   $MEMBER  Ensemble size (if $MEM is not given)
#   $NNODES  Number of total nodes
#   $PPN     Number of processes per node
#
# Return variables:
#   $mem_nodes        Number of nodes for a member
#   $mem_np           Number of processes for a member
#   $repeat_mems      Number of members that uses one round of nodes
#   $parallel_mems    Number of members that can run in parallel
#   $parallel_cycles  Number of parallel cycles needed to finish all members
#-------------------------------------------------------------------------------

local MEM=${1:-$MEMBER}
local MIN_NP=${2:-none}
local MAX_NP=${3:-none}

#-------------------------------------------------------------------------------
# Determine optimal numbers of nodes ($mem_nodes) and processes ($mem_np)
# based on the number of members ($MEM) and total nodes ($NNODES) and processes per node ($PPN).
# The 'optimal' here means minimum $mem_np but still can occupy all avaiable processes at once.

if ((NNODES >= MEM)); then
  mem_nodes=$((NNODES/MEM))
  mem_np=$((PPN*mem_nodes))
else
  mem_nodes=1
  mempn=$(((MEM-1)/NNODES+1))
  if ((mempn > PPN)); then
    mem_np=1
  else
    mem_np=$((PPN/mempn))
  fi
fi

#-------------------------------------------------------------------------------
# Limited to the minimum processes for a member ($MIN_NP)
#        and the maximum processes for a member ($MAX_NP),
# and require it occupy full nodes if multiple nodes are used.

if [ "$MIN_NP" != 'none' ] && ((mem_np < MIN_NP)); then
  mem_np=$MIN_NP
fi
if [ "$MAX_NP" != 'none' ] && ((mem_np > MAX_NP)); then
  mem_np=$MAX_NP
fi
if ((mem_np > PPN)); then
  mem_nodes=$(((mem_np-1)/PPN+1))
  if ((mem_nodes > NNODES)); then
    echo "[Error] Total number of nodes is insufficient." >&2
    exit 1
  fi
  mem_np=$((PPN*mem_nodes))
  repeat_mems=$((NNODES/mem_nodes))
  parallel_mems=$repeat_mems
else
  mem_nodes=1
  mempn=$((PPN/mem_np))
  if ((mempn == 1)); then
    mem_np=$PPN
  fi
  repeat_mems=$NNODES
  parallel_mems=$((repeat_mems*mempn))
fi
parallel_cycles=$(((MEM-1)/parallel_mems+1))

#-------------------------------------------------------------------------------
}

#===============================================================================

function set_mem2proc {
#-------------------------------------------------------------------------------
# Set up the relation from members to nodes and processes
#
# Usage: set_mem2proc [MEM]
#
#   MEM  Number of members
#        (default: $MEMBER)
#
# Input variables:
#   $MEMBER             Ensemble size (if $MEM is not given)
#   $NNODES             Number of total nodes
#   $PPN                Number of processes per node
#   $totalnp            Total number of processes
#   $mem_nodes          Number of nodes for a member
#   $mem_np             Number of processes for a member
#   $node[1...$NNODES]  Name of nodes
#
# Return variables:
#   $procs[1...$totalnp]           Sequence of (total) processes
#   $mem2proc[1...($MEM*$mem_np)]  Relation from members to nodes and processes (pseudo 2-D array)
#   $node_m[1...$MEM]              Name of node(s) for each member
#-------------------------------------------------------------------------------

local MEM=${1:-$MEMBER}

#-------------------------------------------------------------------------------

local m=0
local i
local n
local nn
local p

while ((m < MEM)); do
  for i in $(seq $PPN); do
    for n in $(seq $NNODES); do
      m=$((m+1))
      if ((mem_nodes == 1)) && ((m <= MEM)); then
        for p in $(seq $mem_np); do
          mem2proc[$(((m-1)*mem_np+p))]=$n
          node_m[$m]="${node[$n]}*$mem_np"
        done
      fi
      if ((m <= totalnp)); then
        procs[$m]=$n
      fi
    done
  done
done

#-------------------------------------------------------------------------------

if ((mem_nodes > 1)); then
  n=0
  node_m[$m]=''
  for m in $(seq $MEM); do
    for nn in $(seq $mem_nodes); do
      for p in $(seq $((PPN*(nn-1)+1)) $((PPN*nn))); do
        mem2proc[$(((m-1)*mem_np+p))]=$((n+nn))
      done
      node_m[$m]="${node_m[$m]}${node[$((n+nn))]}*$PPN "
    done
    n=$((n+mem_nodes))
    if ((n+mem_nodes > NNODES)); then
      n=0
    fi
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function distribute_da_cycle {
#-------------------------------------------------------------------------------
# Distribute members on nodes for DA cycling run.
#
# Usage: distribute_da_cycle [NODEFILE NODEFILEDIR]
#
#   NODEFILE     The pre-determined nodefile (required when $MACHINE_TYPE = 1)
#   NODEFILEDIR  Directory to output nodefiles
#                '-': No output (default)
#
# Other input variables:
#   $MEMBER        Ensemble size
#   $NNODES        Number of total nodes
#   $PPN           Number of processes per node
#   $MEMBER_FMT
#   $MIN_NP_SCALE
#   $MAX_NP_SCALE
#   
# Return variables:
#   $totalnp                              Total number of processes
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $procs[1...$totalnp]                  Sequence of (total) processes
#   $mem2proc[1...(($MEMBER+2)*$mem_np)]  Relation from members to nodes and processes (pseudo 2-D array)
#   $mmean                                Index of the ensemble mean ($MEMBER+1)
#   $msprd                                Index of the ensemble spread ($MEMBER+2)
#   $node[1...$nnodes]                    Name of each node
#   $name_m[1...$MEMBER+2]                Name of each member
#   $node_m[1...$MEMBER+2]                Name of node(s) for each member
#   $repeat_mems                          Number of members that uses one round of nodes
#   $parallel_mems                        Number of members that can run in parallel
#   $parallel_cycles                      Number of parallel cycles needed to finish all members
#
# Output files:
#   [$TMP/node/proc]       All processes
#   [$TMP/node/node]       One process per node
#   [$TMP/node/proc.MMMM]  Processes for each member
#-------------------------------------------------------------------------------

local NODEFILE=${1:-machinefile}
local NODEFILEDIR=${2:-'-'}

#-------------------------------------------------------------------------------
# Set up node names and member names

totalnp=$((PPN*NNODES))
mmean=$((MEMBER+1))
msprd=$((MEMBER+2))

local n
for n in $(seq $NNODES); do
  node[$n]="($((n-1)))"
done

local m
for m in $(seq $MEMBER); do
  name_m[$m]=$(printf $MEMBER_FMT $m)
done
name_m[$mmean]='mean'
name_m[$msprd]='sprd'

#-------------------------------------------------------------------------------
# If a nodefile is used, parse it and get the node names.

if ((MACHINE_TYPE == 1)); then
  read_nodefile_pbs "$NODEFILE"
fi

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $((MEMBER+1)) $MIN_NP_SCALE $MAX_NP_SCALE

set_mem2proc $((MEMBER+1))

local p
for p in $(seq $mem_np); do
  mem2proc[$(((msprd-1)*mem_np+p))]=${mem2proc[$(((mmean-1)*mem_np+p))]}
done
node_m[$msprd]=${node_m[$mmean]}

#-------------------------------------------------------------------------------
# Create nodefiles

if [ "$NODEFILEDIR" != '-' ] && [ -d "$NODEFILEDIR" ]; then
  local p
  for p in $(seq $totalnp); do
    echo ${node[${procs[$p]}]} >> $NODEFILEDIR/proc
  done
  for n in $(seq $NNODES); do
    echo ${node[$n]} >> $NODEFILEDIR/node
  done
  for m in $(seq $((MEMBER+1))); do
    for p in $(seq $mem_np); do
      echo ${node[${mem2proc[$(((m-1)*mem_np+p))]}]} >> $NODEFILEDIR/proc.${name_m[$m]}
    done
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function distribute_fcst {
#-------------------------------------------------------------------------------
# Distribute members on nodes for ensemble forecasts.
#
# Usage: distribute_fcst MEMBERS [CYCLE NODEFILE NODEFILEDIR]
#
#   MEMBERS      List of forecast members
#   CYCLE        Number of forecast cycles run in parallel
#                (default: 1)
#   NODEFILE     The pre-determined nodefile (required when $MACHINE_TYPE = 1)
#   NODEFILEDIR  Directory to output nodefiles
#                '-': No output (default)
#
# Other input variables:
#   $NNODES        Number of total nodes
#   $PPN           Number of processes per node
#   $MIN_NP_SCALE
#   $MAX_NP_SCALE
#   $MACHINE_TYPE
#   
# Return variables:
#   $fmember                              Number of forecast members
#   $fmembertot                           Total number of forecast numbers for all cycles
#   $totalnp                              Total number of processes
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $procs[1...$totalnp]                  Sequence of (total) processes
#   $mem2proc[1...($fmembertot*$mem_np)]  Relation from members to nodes and processes (pseudo 2-D array)
#   $node[1...$nnodes]                    Name of each node
#   $name_m[1...$fmembertot]              Name of each member
#   $node_m[1...$fmembertot]              Name of node(s) for each member
#   $repeat_mems                          Number of members that uses one round of nodes
#   $parallel_mems                        Number of members that can run in parallel
#   $parallel_cycles                      Number of parallel cycles needed to finish all members
#
# Output files:
#   [$TMP/node/proc]            All processes
#   [$TMP/node/node]            One process per node
#   [$TMP/node/proc.CCCC.MMMM]  Processes for each cycle and each member
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local MEMBERS="$1"
local CYCLE=${2:-1}
local NODEFILE=${3:-machinefile}
local NODEFILEDIR=${4:-'-'}

#-------------------------------------------------------------------------------
# Set up node names and member names, and also get the number of members

totalnp=$((PPN*NNODES))

local n
for n in $(seq $NNODES); do
  node[$n]="($((n-1)))"
done

fmember=0
for iname in $MEMBERS; do
  fmember=$((fmember+1))
  name_m[$fmember]=$iname
done

fmembertot=$((fmember * CYCLE))

#-------------------------------------------------------------------------------
# If a nodefile is used, parse it and get the node names.

if ((MACHINE_TYPE == 1)); then
  read_nodefile_pbs "$NODEFILE"
fi

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $fmembertot $MIN_NP_SCALE $MAX_NP_SCALE

set_mem2proc $fmembertot

#-------------------------------------------------------------------------------
# Create nodefiles

if [ "$NODEFILEDIR" != '-' ] && [ -d "$NODEFILEDIR" ]; then
  local p
  for p in $(seq $totalnp); do
    echo ${node[${procs[$p]}]} >> $NODEFILEDIR/proc
  done
  for n in $(seq $NNODES); do
    echo ${node[$n]} >> $NODEFILEDIR/node
  done
  for c in $(seq $CYCLE); do
    cf=$(printf $CYCLE_FMT $c)
    for m in $(seq $fmember); do
      mm=$(((c-1) * fmember + m))
      if ((mm != m)); then
        name_m[$mm]=${name_m[$m]}
      fi
      for p in $(seq $mem_np); do
        echo ${node[${mem2proc[$(((mm-1)*mem_np+p))]}]} >> $NODEFILEDIR/proc.${cf}.${name_m[$mm]}
      done
    done
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function read_nodefile_pbs {
#-------------------------------------------------------------------------------
# Parse the PBS-type nodefile.
# check if it is consistent to the $NNODES and $PPN settings and get the node names.
#
# Usage: read_nodefile_pbs NODEFILE
#
#   NODEFILE  PBS-type Nodefile for mpiexec
#
# Other input variables:
#   $NNODES  Number of total nodes
#   $PPN     Number of processes per node
#   
# Return variables:
#   $node[1...$nnodes]  Name of each node
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local NODEFILE="$1"

if [ ! -s "$NODEFILE" ]; then
  echo "[Error] $FUNCNAME: Can not find \$NODEFILE '$NODEFILE'" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

local nodelist=$(cat $NODEFILE | sort | uniq)
local inode
local ippn
local n=0
for inode in $nodelist; do
  n=$((n+1))
  node[$n]=$inode
  ippn=`cat $NODEFILE | grep $inode | wc -l`
  if ((ippn != PPN)); then
    echo "[Error] $FUNCNAME: Number of processes per node in \$NODEFILE" >&2
    echo "          is not consistent to the setting in 'configure.sh'" >&2
    exit 1
  fi
done
if ((n != NNODES)); then
  echo "[Error] $FUNCNAME: Number of nodes in \$NODEFILE" >&2
  echo "          is not consistent to the setting in 'configure.sh'" >&2
  exit 1
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
