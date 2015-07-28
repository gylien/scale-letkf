#!/bin/bash
#===============================================================================
#
#  Adaptively distribute members on nodes
#  April   2013,          Guo-Yuan Lien
#  August  2014, modified Guo-Yuan Lien
#
#  *Require source 'config.main' first.
#
#===============================================================================

set_mem_np () {
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
  local mempn=$(((MEM-1)/NNODES+1))
  if ((mempn > PPN)); then
    mem_np=1
  else
    mem_np=$((PPN/mempn))
  fi
fi

#-------------------------------------------------------------------------------
# Limited to the minimum processes for a member ($MIN_NP)
#        and the maximum processes for a member ($MAX_NP),

if [ "$MIN_NP" != 'none' ] && ((mem_np < MIN_NP)); then
  mem_np=$MIN_NP
fi
if [ "$MAX_NP" != 'none' ] && ((mem_np > MAX_NP)); then
  mem_np=$MAX_NP
fi

#-------------------------------------------------------------------------------
# Re-calculate $mem_nodes, and calculate
# the period that members run on repeated nodes ($repeat_mems), and
# the number of members that can be run parallelly ($parallel_mems).

mem_nodes=$(((mem_np-1)/PPN+1))
if ((mem_nodes > NNODES)); then
  echo "[Error] Total number of nodes is insufficient." >&2
  exit 1
fi
repeat_mems=$((NNODES/mem_nodes))
if ((mem_nodes == 1)); then
  parallel_mems=$((repeat_mems * (PPN/mem_np)))
else
  parallel_mems=$repeat_mems
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

set_mem2node () {
#-------------------------------------------------------------------------------
# Set up the relation from members to nodes and processes
#
# Usage: set_mem2node [MEM]
#
#   MEM  Number of members
#        (default: $MEMBER)
#
# Input variables:
#   $MEMBER             Ensemble size (if $MEM is not given)
#   $NNODES             Number of total nodes
#   $PPN                Number of processes per node
#   $mem_nodes          Number of nodes for a member
#   $mem_np             Number of processes for a member
#   $node[1...$NNODES]  Name of nodes
#
# Return variables:
#   $totalnp                       Total number of processes
#   $procs[1...$totalnp]           Sequence of (total) processes
#   $mem2node[1...($MEM*$mem_np)]  Relation from (members, processes) to nodes (pseudo 2-D array)
#   $node_m[1...$MEM]              Name of node(s) for each member
#-------------------------------------------------------------------------------

local MEM=${1:-$MEMBER}

#-------------------------------------------------------------------------------

local tppn=$((mem_np / mem_nodes))
local tmod=$((mem_np % mem_nodes))

local m
local nn
local q
local qn
local n=0
local procs_add=1
totalnp=0

for m in $(seq $MEM); do

  ###### SHORT node list description
  if ((mem_nodes == 1)); then
    node_m[$m]="${node[$((n+1))]}*$tppn"
  elif ((tmod == 0)); then
    node_m[$m]="[${node[$((n+1))]}-${node[$((n+mem_nodes))]}]*$tppn"
  else
    if ((tmod == 1)); then
      node_m[$m]="${node[$((n+1))]}*$((tppn+1))"
    else
      node_m[$m]="[${node[$((n+1))]}-${node[$((n+tmod))]}]*$((tppn+1))"
    fi
    if (($((mem_nodes - tmod)) == 1)); then
      node_m[$m]="${node_m[$m]} ${node[$((n+mem_nodes))]}*$tppn"
    else
      node_m[$m]="${node_m[$m]} [${node[$((n+tmod+1))]}-${node[$((n+mem_nodes))]}]*$tppn"
    fi
  fi
  ######

#  node_m[$m]=''
  qn=0
  for nn in $(seq $mem_nodes); do
    if ((nn <= tmod)); then
      tppnt=$((tppn+1))
    else
      tppnt=$tppn
    fi
    for q in $(seq $((qn+1)) $((qn+tppnt))); do
      mem2node[$(((m-1)*mem_np+q))]=$((n+nn))
    done
    qn=$((qn+tppnt))
#    node_m[$m]="${node_m[$m]}${node[$((n+nn))]}*$tppnt "

    if ((procs_add == 1)); then
      for p in $(seq $((totalnp+1)) $((totalnp+PPN))); do
        procs[$p]=$((n+nn))
      done
      totalnp=$((totalnp+PPN))
    fi
  done

  n=$((n+mem_nodes))
  if ((n+mem_nodes > NNODES)); then
    n=0
    procs_add=0
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

distribute_da_cycle () {
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
#   $SCALE_NP
#   
# Return variables:
#   $totalnp                              Total number of processes
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $procs[1...$totalnp]                  Sequence of (total) processes
#   $mem2node[1...(($MEMBER+2)*$mem_np)]  Relation from members to nodes and processes (pseudo 2-D array)
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

mmean=$((MEMBER+1))
msprd=$((MEMBER+2))

if ((MACHINE_TYPE == 1)); then
  read_nodefile_pbs "$NODEFILE"
elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11 || MACHINE_TYPE == 12)); then
  local n
  local p
  for n in $(seq $NNODES_real); do
    for p in $(seq $PPN_real); do
      node[$(((n-1)*PPN_real+p))]="($((n-1)))"
    done
  done
else
  echo "[Error] Unsupported \$MACHINE_TYPE." >&2
  exit 1
fi

local m
for m in $(seq $MEMBER); do
  name_m[$m]=$(printf $MEMBER_FMT $m)
done
name_m[$mmean]='mean'
name_m[$msprd]='sprd'

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $((MEMBER+1)) $SCALE_NP $SCALE_NP

set_mem2node $((MEMBER+1))

local p
for p in $(seq $mem_np); do
  mem2node[$(((msprd-1)*mem_np+p))]=${mem2node[$(((mmean-1)*mem_np+p))]}
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
      echo ${node[${mem2node[$(((m-1)*mem_np+p))]}]} >> $NODEFILEDIR/proc.${name_m[$m]}
    done
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

distribute_fcst () {
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
#   $NNODES_real   XXXXXX
#   $PPN_real      XXXXXX
#   $SCALE_NP
#   $MACHINE_TYPE
#   
# Return variables:
#   $totalnp                              Total number of processes
#   $fmember                              Number of forecast members
#   $fmembertot                           Total number of forecast numbers for all cycles
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $procs[1...$totalnp]                  Sequence of (total) processes
#   $mem2node[1...($fmembertot*$mem_np)]  Relation from (members, processes) to nodes (pseudo 2-D array)
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

if ((MACHINE_TYPE == 1)); then
  read_nodefile_pbs "$NODEFILE"
elif ((MACHINE_TYPE == 10 || MACHINE_TYPE == 11 || MACHINE_TYPE == 12)); then
  local n
  local p
  for n in $(seq $NNODES_real); do
    for p in $(seq $PPN_real); do
      node[$(((n-1)*PPN_real+p))]="($((n-1)))"
    done
  done
else
  echo "[Error] Unsupported \$MACHINE_TYPE." >&2
  exit 1
fi

fmember=0
for iname in $MEMBERS; do
  fmember=$((fmember+1))
  name_m[$fmember]=$iname
done

for c in $(seq 2 $CYCLE); do
  for m in $(seq $fmember); do
    name_m[$(((c-1)*fmember+m))]=${name_m[$m]}
  done
done

fmembertot=$((fmember * CYCLE))

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $fmembertot $SCALE_NP $SCALE_NP

set_mem2node $fmembertot

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
      for p in $(seq $mem_np); do
        echo ${node[${mem2node[$(((mm-1)*mem_np+p))]}]} >> $NODEFILEDIR/proc.${cf}.${name_m[$mm]}
      done
    done
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

read_nodefile_pbs () {
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
