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
#   $mem_nodes  Number of nodes for a member
#   $mem_np     Number of processes for a member
#   $totalnp    Total number of processes
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

#-------------------------------------------------------------------------------

totalnp=$((NNODES*PPN))

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
#   $n_mem                         Number of members that use one round of nodes
#   $n_mempn                       Number of members that run in parallel in a node
#   $repeat_mems                   (= $n_mem)
#   $parallel_mems                 Number of members that run in parallel
#   $nitmax                        Number of parallel cycles needed to finish all members
#   $mem2node[1...($MEM*$mem_np)]  Relation from (members, m_processes) to nodes (pseudo 2-D array)
#   $mem2proc[1...($MEM*$mem_np)]  Relation from (members, m_processes) to processes (pseudo 2-D array)
#   $proc2node[1...$totalnp]       Relation from processes to nodes
#   $proc2group[1...$totalnp]      Relation from processes to groups
#   $proc2grpproc[1...$totalnp]    Relation from processes to m_processes
#   $node_m[1...$MEM]              Short node list description of each member
#-------------------------------------------------------------------------------

local MEM=${1:-$MEMBER}

#-------------------------------------------------------------------------------

local ns=0
local n
local p
for n in $(seq $NNODES); do
  for p in $(seq $((ns+1)) $((ns+PPN))); do
    proc2node[$p]=$n
  done
  ns=$((ns+PPN))
done

if ((mem_nodes > 1)); then
  n_mem=$((NNODES / mem_nodes))
  n_mempn=1
else
  n_mem=$NNODES
  n_mempn=$((PPN / mem_np))
fi
repeat_mems=$n_mem
parallel_mems=$((n_mem * n_mempn))
nitmax=$(((MEM - 1) / (n_mem * n_mempn) + 1))
local tppn=$((mem_np / mem_nodes))
local tmod=$((mem_np % mem_nodes))

local m
local nn
local q
local qs

m=1
for it in $(seq $nitmax); do
  for i in $(seq 0 $((n_mempn-1))); do
    n=0
    for j in $(seq 0 $((n_mem-1))); do
      if ((m > MEM && it > 1)); then break; fi

      qs=0
      for nn in $(seq 0 $((mem_nodes-1))); do
        if ((nn < tmod)); then
          tppnt=$((tppn+1))
        else
          tppnt=$tppn
        fi
        for q in $(seq 0 $((tppnt-1))); do
          ip=$(((n+nn)*PPN + i*mem_np + q))
          if ((m <= MEM)); then
            mem2node[$(((m-1)*mem_np+qs+1))]=$((n+nn+1))
            mem2proc[$(((m-1)*mem_np+qs+1))]=$((ip+1))
            if ((it == 1)); then
              proc2group[$((ip+1))]=$m
              proc2grpproc[$((ip+1))]=$((qs+1))
            fi
          fi
          qs=$((qs+1))
        done
      done

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

      m=$((m+1))
      n=$((n+mem_nodes))
    done
    if ((m > MEM && it > 1)); then break; fi
  done
  if ((m > MEM && it > 1)); then break; fi
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
#   $MEMBER      Ensemble size
#   $NNODES      Number of total nodes
#   $PPN         Number of processes per node
#   $MEMBER_FMT
#   $SCALE_NP
#   
# Return variables:
#   $node[1...$nnodes]                    Name of nodes
#
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $totalnp                              Total number of processes
#
#   $mmean                                Index of the ensemble mean ($MEMBER+1)
#   $msprd                                Index of the ensemble spread ($MEMBER+2)
#   $node_m[1...$MEMBER+2]                Short node list description of each member
#   $name_m[1...$MEMBER+2]                Name of members
#
#   $n_mem                                Number of members that use one round of nodes
#   $n_mempn                              Number of members that run in parallel in a node
#   $repeat_mems                          (= $n_mem)
#   $parallel_mems                        Number of members that run in parallel
#   $nitmax                               Number of parallel cycles needed to finish all members
#   $mem2node[1...(($MEMBER+2)*$mem_np)]  Relation from (members, m_processes) to nodes (pseudo 2-D array)
#   $mem2proc[1...(($MEMBER+2)*$mem_np)]  Relation from (members, m_processes) to processes (pseudo 2-D array)
#   $proc2node[1...$totalnp]              Relation from processes to nodes
#   $proc2group[1...$totalnp]             Relation from processes to groups
#   $proc2grpproc[1...$totalnp]           Relation from processes to m_processes
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
  mem2proc[$(((msprd-1)*mem_np+p))]=${mem2proc[$(((mmean-1)*mem_np+p))]}
done
node_m[$msprd]=${node_m[$mmean]}

#-------------------------------------------------------------------------------
# Create nodefiles

if [ "$NODEFILEDIR" != '-' ] && [ -d "$NODEFILEDIR" ]; then
  local p
  for p in $(seq $totalnp); do  
    echo ${node[${proc2node[$p]}]} >> $NODEFILEDIR/proc
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
#   $node[1...$nnodes]                    Name of nodes
#   $fmember                              Number of forecast members
#   $fmembertot                           Total number of forecast numbers for all cycles
#
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $totalnp                              Total number of processes
#
#   $node_m[1...$fmembertot]              Short node list description of each member
#   $name_m[1...$fmembertot]              Name of members
#
#   $n_mem                                Number of members that use one round of nodes
#   $n_mempn                              Number of members that run in parallel in a node
#   $repeat_mems                          (= $n_mem)
#   $parallel_mems                        Number of members that run in parallel
#   $nitmax                               Number of parallel cycles needed to finish all members
#   $mem2node[1...($fmembertot*$mem_np)]  Relation from (members, m_processes) to nodes (pseudo 2-D array)
#   $mem2proc[1...($fmembertot*$mem_np)]  Relation from (members, m_processes) to processes (pseudo 2-D array)
#   $proc2node[1...$totalnp]              Relation from processes to nodes
#   $proc2group[1...$totalnp]             Relation from processes to groups
#   $proc2grpproc[1...$totalnp]           Relation from processes to m_processes
#
# Output files:
#   [$TMP/node/proc]            All processes
#   [$TMP/node/node]            One process per node
#   [$TMP/node/proc.CCCC.MMMM]  Processes for each member
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
    echo ${node[${proc2node[$p]}]} >> $NODEFILEDIR/proc
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
