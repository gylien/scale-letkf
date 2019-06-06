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
#   $MEMBER        Ensemble size (if $MEM is not given)
#   $NNODES_APPAR  Apparent number of total nodes
#   $PPN_APPAR     Apparent number of processes per node
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
# based on the number of members ($MEM) and total nodes ($NNODES_APPAR) and processes per node ($PPN_APPAR).
# The 'optimal' here means minimum $mem_np but still can occupy all avaiable processes at once.

if ((NNODES_APPAR >= MEM)); then
  mem_nodes=$((NNODES_APPAR/MEM))
  mem_np=$((PPN_APPAR*mem_nodes))
else
  mem_nodes=1
  local mempn=$(((MEM-1)/NNODES_APPAR+1))
  if ((mempn > PPN_APPAR)); then
    mem_np=1
  else
    mem_np=$((PPN_APPAR/mempn))
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

mem_nodes=$(((mem_np-1)/PPN_APPAR+1))
if ((mem_nodes > NNODES_APPAR)); then
  echo "[Error] Total number of nodes is insufficient." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

totalnp=$((NNODES_APPAR*PPN_APPAR))

#-------------------------------------------------------------------------------
}

#===============================================================================

set_mem2node () {
#-------------------------------------------------------------------------------
# Set up the relation from members to nodes and processes
#
# Usage: set_mem2node [MEM USE_CACHE SAVE_CACHE]
#
#   MEM         Number of members
#               (default: $MEMBER)
#   USE_CACHE   Use 'distr' file cache?
#               0: No (default)
#               1: Yes
#   SAVE_CACHE  Save 'distr' file cache?
#               0: No
#               1: Yes (default)
#
# Input variables:
#   $MEMBER                   Ensemble size (if $MEM is not given)
#   $NNODES_APPAR             Apparent number of total nodes
#   $PPN_APPAR                Apparent number of processes per node
#   $mem_nodes                Number of nodes for a member
#   $mem_np                   Number of processes for a member
#   $node[1...$NNODES_APPAR]  Name of nodes
#   $NODEFILEDIR
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

local MEM=${1:-$MEMBER}; shift
local USE_CACHE=${1:-0}; shift
local SAVE_CACHE=${1:-1}

#-------------------------------------------------------------------------------

local ns=0
local n
local p

if ((USE_CACHE == 0)); then
  for n in $(seq $NNODES_APPAR); do
    for p in $(seq $((ns+1)) $((ns+PPN_APPAR))); do
      proc2node[$p]=$n
      if ((SAVE_CACHE == 1)); then
        echo "proc2node[$p]=$n" >> $NODEFILEDIR/distr
      fi
    done
    ns=$((ns+PPN_APPAR))
  done
fi # ((USE_CACHE == 0))

if ((mem_nodes > 1)); then
  n_mem=$((NNODES_APPAR / mem_nodes))
  n_mempn=1
else
  n_mem=$NNODES_APPAR
  n_mempn=$((PPN_APPAR / mem_np))
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

if ((USE_CACHE == 0)); then
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
            ip=$(((n+nn)*PPN_APPAR + i*mem_np + q))
            if ((m <= MEM)); then
              mem2node[$(((m-1)*mem_np+qs+1))]=$((n+nn+1))
              mem2proc[$(((m-1)*mem_np+qs+1))]=$((ip+1))
              if ((SAVE_CACHE == 1)); then
                echo "mem2node[$(((m-1)*mem_np+qs+1))]=$((n+nn+1))" >> $NODEFILEDIR/distr
                echo "mem2proc[$(((m-1)*mem_np+qs+1))]=$((ip+1))" >> $NODEFILEDIR/distr
              fi
              if ((it == 1)); then
                proc2group[$((ip+1))]=$m
                proc2grpproc[$((ip+1))]=$((qs+1))
                if ((SAVE_CACHE == 1)); then
                  echo "proc2group[$((ip+1))]=$m" >> $NODEFILEDIR/distr
                  echo "proc2grpproc[$((ip+1))]=$((qs+1))" >> $NODEFILEDIR/distr
                fi
              fi
            fi
            qs=$((qs+1))
          done
        done

        ###### SHORT node list description
        if ((mem_nodes == 1)); then
          node_m[$m]="${node[$((n+1))]}*$tppn"
          if ((SAVE_CACHE == 1)); then
            node_m_out[$m]="\${node[$((n+1))]}*$tppn"
          fi
        elif ((tmod == 0)); then
          node_m[$m]="[${node[$((n+1))]}-${node[$((n+mem_nodes))]}]*$tppn"
          if ((SAVE_CACHE == 1)); then
            node_m_out[$m]="[\${node[$((n+1))]}-\${node[$((n+mem_nodes))]}]*$tppn"
          fi
        else
          if ((tmod == 1)); then
            node_m[$m]="${node[$((n+1))]}*$((tppn+1))"
            if ((SAVE_CACHE == 1)); then
              node_m_out[$m]="\${node[$((n+1))]}*$((tppn+1))"
            fi
          else
            node_m[$m]="[${node[$((n+1))]}-${node[$((n+tmod))]}]*$((tppn+1))"
            if ((SAVE_CACHE == 1)); then
              node_m_out[$m]="[\${node[$((n+1))]}-\${node[$((n+tmod))]}]*$((tppn+1))"
            fi
          fi
          if (($((mem_nodes - tmod)) == 1)); then
            node_m[$m]="${node_m[$m]} ${node[$((n+mem_nodes))]}*$tppn"
            if ((SAVE_CACHE == 1)); then
              node_m_out[$m]="${node_m_out[$m]} \${node[$((n+mem_nodes))]}*$tppn"
            fi
          else
            node_m[$m]="${node_m[$m]} [${node[$((n+tmod+1))]}-${node[$((n+mem_nodes))]}]*$tppn"
            if ((SAVE_CACHE == 1)); then
              node_m_out[$m]="${node_m_out[$m]} [\${node[$((n+tmod+1))]}-\${node[$((n+mem_nodes))]}]*$tppn"
            fi
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

  if ((SAVE_CACHE == 1)); then
    for m in $(seq $MEM); do
      echo "node_m[$m]=\"${node_m_out[$m]}\"" >> $NODEFILEDIR/distr
    done
  fi
fi # ((USE_CACHE == 0))


#-------------------------------------------------------------------------------
}

#===============================================================================

distribute_da_cycle () {
#-------------------------------------------------------------------------------
# Distribute members on nodes for DA cycling run.
#
# Usage: distribute_da_cycle [NODELIST NODEFILEDIR SAVE_CACHE MEMBERS]
#
#   NODELIST     List of node names
#                '-':    The node names has not been determined yet (default)
#                '(0)':  Ordered sequence '(*)' starting from 0
#                (other) A given nodelist file
#   NODEFILEDIR  Directory to output nodefiles
#                '-': No output (default)
#   SAVE_CACHE   Save 'distr' file cache?
#                0: No
#                1: Yes (default)
#   MEMBERS      List of forecast members
#                'all': All sequential numbers (default)
#
# Other input variables:
#   $MEMBER        Ensemble size
#   $DET_RUN       Whether the deterministic run is enabled
#   $NNODES        Number of total nodes
#   $PPN           Number of processes per node
#   $NNODES_APPAR  Apparent number of total nodes
#   $PPN_APPAR     Apparent number of processes per node
#   $MEMBER_FMT
#   $SCALE_NP_TOT
#   
# Return variables:
#   $node[1...$nnodes]                    Name of nodes
#
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $totalnp                              Total number of processes
#
#   $mtot                                 Number of members + ensemble mean (+ deterministic run)
#   $mmean                                Index of the ensemble mean ($MEMBER+1)
#   $mmdet                                Index of the deterministic run ($MEMBER+2)
#   $node_m[1...$mtot]                    Short node list description of each member
#   $name_m[1...$mtot]                    Name of members
#
#   $n_mem                                Number of members that use one round of nodes
#   $n_mempn                              Number of members that run in parallel in a node
#   $repeat_mems                          (= $n_mem)
#   $parallel_mems                        Number of members that run in parallel
#   $nitmax                               Number of parallel cycles needed to finish all members
#   $mem2node[1...($mtot*$mem_np)]        Relation from (members, m_processes) to nodes (pseudo 2-D array)
#   $mem2proc[1...($mtot*$mem_np)]        Relation from (members, m_processes) to processes (pseudo 2-D array)
#   $proc2node[1...$totalnp]              Relation from processes to nodes
#   $proc2group[1...$totalnp]             Relation from processes to groups
#   $proc2grpproc[1...$totalnp]           Relation from processes to m_processes
#
# Output files:
#   [$NODEFILEDIR/proc]       All processes
#   [$NODEFILEDIR/node]       One process per node
#   [$NODEFILEDIR/distr]      File cache for shell array variables
#-------------------------------------------------------------------------------

local NODELIST=${1:--}; shift
local NODEFILEDIR=${1:--}; shift
local SAVE_CACHE=${1:-1}; shift
local MEMBERS="${1:-all}"

if [ "$NODEFILEDIR" != '-' ] && [ ! -d "$NODEFILEDIR" ]; then
  echo "[Error] $FUNCNAME: \$NODEFILEDIR is given but is not an existing directory: '$NODEFILEDIR'" >&2
  exit 1
fi

if [ -s "${NODEFILEDIR}/proc" ] && [ -s "${NODEFILEDIR}/node" ]; then
#  echo "[INFO] $FUNCNAME: Skip creating 'proc' and 'node' files because they already exist." >&2
  NODELIST='-'
fi

if [ "$NODEFILEDIR" = '-' ]; then
  SAVE_CACHE=0
fi

local use_cache=0
if [ -s "${NODEFILEDIR}/distr" ]; then
  use_cache=1
  SAVE_CACHE=0
fi

#if ((SAVE_CACHE == 1)); then
#  echo "[INFO] $FUNCNAME: Save 'distr' file cache." >&2
#fi

#-------------------------------------------------------------------------------
# Set up node names

if [ "$NODELIST" = '-' ]; then
  : # do nothing
elif [ "$NODELIST" = '(0)' ]; then
  local n
  local p
  local appar_npn=$((NNODES_APPAR/NNODES))
  if ((use_cache == 0)); then
    for n in $(seq $NNODES); do
      for p in $(seq $appar_npn); do
        node[$(((n-1)*appar_npn+p))]="($((n-1)))"
        if ((SAVE_CACHE == 1)); then
          echo "node[$(((n-1)*appar_npn+p))]=\"($((n-1)))\"" >> $NODEFILEDIR/distr
        fi
      done
    done
  fi # ((use_cache == 0))
else
  read_nodefile_pbs "$NODELIST"
fi

#-------------------------------------------------------------------------------
# Set up member names

if [ "$MEMBERS" = 'all' ]; then
  local m
  for m in $(seq $MEMBER); do
    name_m[$m]=$(printf $MEMBER_FMT $m)
  done
else
  local m=0
  for iname in $MEMBERS; do
    m=$((m+1))
    name_m[$m]=$iname
  done
  if ((m != MEMBER)); then
    echo "[Error] Number of members (\$MEMBERS) is not equal to \$MEMBER." >&2
    exit 1
  fi
fi

mtot=$((MEMBER+1))
mmean=$((MEMBER+1))
name_m[$mmean]='mean'
if ((DET_RUN == 1)); then
  mtot=$((MEMBER+2))
  mmdet=$((MEMBER+2))
  name_m[$mmdet]='mdet'
fi

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $mtot $SCALE_NP_TOT $SCALE_NP_TOT

set_mem2node $mtot $use_cache $SAVE_CACHE

if ((use_cache == 1)); then
#  echo "[INFO] $FUNCNAME: Use 'distr' file cache." >&2
  . ${NODEFILEDIR}/distr
fi

#-------------------------------------------------------------------------------
# Create nodefiles

if [ "$NODELIST" != '-' ] && [ "$NODEFILEDIR" != '-' ]; then
#  echo "[INFO] $FUNCNAME: Save 'proc', 'node' files." >&2
  local p
  for p in $(seq $totalnp); do  
    echo ${node[${proc2node[$p]}]} >> $NODEFILEDIR/proc
  done
  for n in $(seq $NNODES_APPAR); do
    echo ${node[$n]} >> $NODEFILEDIR/node
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

distribute_da_cycle_set () {
#-------------------------------------------------------------------------------
# Distribute members on nodes for DA cycling run.
#
# Usage: distribute_da_cycle_set [NODELIST NODEFILEDIR SAVE_CACHE]
#
#   NODELIST     List of node names
#                '-':    The node names has not been determined yet (default)
#                '(0)':  Ordered sequence '(*)' starting from 0
#                (other) A given nodelist file
#   NODEFILEDIR  Directory to output nodefiles
#                '-': No output (default)
#   SAVE_CACHE   Save 'distr' file cache?
#                0: No
#                1: Yes (default)
#
# Other input variables:
#   $MEMBER        Ensemble size
#   $DET_RUN       Whether the deterministic run is enabled
#   $NNODES        Number of total nodes
#   $PPN           Number of processes per node
#   $NNODES_APPAR  Apparent number of total nodes
#   $PPN_APPAR     Apparent number of processes per node
#   $MEMBER_FMT
#   $SCALE_NP_TOT
#   
# Return variables:
#   $node[1...$nnodes]                    Name of nodes
#
#   $mem_nodes                            Number of nodes for a member
#   $mem_np                               Number of processes for a member
#   $totalnp                              Total number of processes
#
#   $mtot                                 Number of members + ensemble mean (+ deterministic run)
#   $mmean                                Index of the ensemble mean ($MEMBER+1)
#   $mmdet                                Index of the deterministic run ($MEMBER+2)
#   $node_m[1...$mtot]                    Short node list description of each member
#   $name_m[1...$mtot]                    Name of members
#
#   $n_mem                                Number of members that use one round of nodes
#   $n_mempn                              Number of members that run in parallel in a node
#   $repeat_mems                          (= $n_mem)
#   $parallel_mems                        Number of members that run in parallel
#   $nitmax                               Number of parallel cycles needed to finish all members
#   $mem2node[1...($mtot*$mem_np)]        Relation from (members, m_processes) to nodes (pseudo 2-D array)
#   $mem2proc[1...($mtot*$mem_np)]        Relation from (members, m_processes) to processes (pseudo 2-D array)
#   $proc2node[1...$totalnp]              Relation from processes to nodes
#   $proc2group[1...$totalnp]             Relation from processes to groups
#   $proc2grpproc[1...$totalnp]           Relation from processes to m_processes
#
# Output files:
#   [$NODEFILEDIR/proc]       All processes
#   [$NODEFILEDIR/node]       One process per node
#   [$NODEFILEDIR/distr]      File cache for shell array variables
#-------------------------------------------------------------------------------

local NODELIST=${1:--}; shift
local NODEFILEDIR=${1:--}; shift
local SAVE_CACHE=${1:-1}

if [ "$NODEFILEDIR" != '-' ] && [ ! -d "$NODEFILEDIR" ]; then
  echo "[Error] $FUNCNAME: \$NODEFILEDIR is given but is not an existing directory: '$NODEFILEDIR'" >&2
  exit 1
fi

if [ -s "${NODEFILEDIR}/proc" ] && [ -s "${NODEFILEDIR}/node" ]; then
#  echo "[INFO] $FUNCNAME: Skip creating 'proc' and 'node' files because they already exist." >&2
  NODELIST='-'
fi

if [ "$NODEFILEDIR" = '-' ]; then
  SAVE_CACHE=0
fi

local use_cache=0
if [ -s "${NODEFILEDIR}/distr" ]; then
  use_cache=1
  SAVE_CACHE=0
fi

#if ((SAVE_CACHE == 1)); then
#  echo "[INFO] $FUNCNAME: Save 'distr' file cache." >&2
#fi

#-------------------------------------------------------------------------------
# Set up node names

if [ "$NODELIST" = '-' ]; then
  : # do nothing
elif [ "$NODELIST" = '(0)' ]; then
  local n
  local p
  local appar_npn=$((NNODES_APPAR/NNODES))
######
  local s
  for s in $(seq 2); do
######
  if ((use_cache == 0)); then
    for n in $(seq $NNODES); do
      for p in $(seq $appar_npn); do
        node[$(((s-1)*NNODES*appar_npn+(n-1)*appar_npn+p))]="($(((s-1)*NNODES+n-1)))"
        if ((SAVE_CACHE == 1)); then
          echo "node[$(((s-1)*NNODES*appar_npn+(n-1)*appar_npn+p))]=\"($(((s-1)*NNODES+n-1)))\"" >> $NODEFILEDIR/distr
        fi
      done
    done
  fi # ((use_cache == 0))
######
  done
######
else
  read_nodefile_pbs "$NODELIST"
fi

#-------------------------------------------------------------------------------
# Set up member names

local m
for m in $(seq $MEMBER); do
  name_m[$m]=$(printf $MEMBER_FMT $m)
done

mtot=$((MEMBER+1))
mmean=$((MEMBER+1))
name_m[$mmean]='mean'
if ((DET_RUN == 1)); then
  mtot=$((MEMBER+2))
  mmdet=$((MEMBER+2))
  name_m[$mmdet]='mdet'
fi

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $mtot $SCALE_NP_TOT $SCALE_NP_TOT

set_mem2node $mtot $use_cache $SAVE_CACHE

if ((use_cache == 1)); then
#  echo "[INFO] $FUNCNAME: Use 'distr' file cache." >&2
  . ${NODEFILEDIR}/distr
fi

#-------------------------------------------------------------------------------
# Create nodefiles

######
for s in $(seq 2); do
######
if [ "$NODELIST" != '-' ] && [ "$NODEFILEDIR" != '-' ]; then
#  echo "[INFO] $FUNCNAME: Save 'proc', 'node' files." >&2
  local p
  for p in $(seq $totalnp); do  
    if ((s == 1)); then ###
      echo ${node[${proc2node[$p]}]} >> $NODEFILEDIR/proc
    fi
    echo ${node[$(((s-1)*NNODES+${proc2node[$p]}))]} >> $NODEFILEDIR/set${s}.proc
  done
  for n in $(seq $NNODES_APPAR); do
    if ((s == 1)); then ###
      echo ${node[$n]} >> $NODEFILEDIR/node
    fi
    echo ${node[$(((s-1)*NNODES+$n))]} >> $NODEFILEDIR/set${s}.node
  done
fi
######
done
######

#-------------------------------------------------------------------------------
}

#===============================================================================

distribute_fcst () {
#-------------------------------------------------------------------------------
# Distribute members on nodes for ensemble forecasts.
#
# Usage: distribute_fcst MEMBERS [CYCLE NODELIST NODEFILEDIR SAVE_CACHE]
#
#   MEMBERS      List of forecast members
#   CYCLE        Number of forecast cycles run in parallel
#                (default: 1)
#   NODELIST     List of node names
#                '-':    The node names has not been determined yet (default)
#                '(0)':  Ordered sequence '(*)' starting from 0
#                (other) A given nodelist file
#   NODEFILEDIR  Directory to output nodefiles
#                '-': No output (default)
#   SAVE_CACHE   Save 'distr' file cache?
#                0: No
#                1: Yes (default)
#
# Other input variables:
#   $NNODES        Number of total nodes
#   $PPN           Number of processes per node
#   $NNODES_APPAR  Apparent number of total nodes
#   $PPN_APPAR     Apparent number of processes per node
#   $SCALE_NP_TOT
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
#   $cycle_auto                           Automatically determined $CYCLE value
#
# Output files:
#   [$NODEFILEDIR/proc]            All processes
#   [$NODEFILEDIR/node]            One process per node
#   [$NODEFILEDIR/distr]           File cache for shell array variables
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local MEMBERS="$1"; shift
local CYCLE=${1:-0}; shift
local NODELIST=${1:--}; shift
local NODEFILEDIR=${1:--}; shift
local SAVE_CACHE=${1:-1}

if [ "$NODEFILEDIR" != '-' ] && [ ! -d "$NODEFILEDIR" ]; then
  echo "[Error] $FUNCNAME: \$NODEFILEDIR is given but is not an existing directory: '$NODEFILEDIR'" >&2
  exit 1
fi

if [ -s "${NODEFILEDIR}/proc" ] && [ -s "${NODEFILEDIR}/node" ]; then
#  echo "[INFO] $FUNCNAME: Skip creating 'proc' and 'node' files because they already exist." >&2
  NODELIST='-'
fi

if [ "$NODEFILEDIR" = '-' ]; then
  SAVE_CACHE=0
fi

local use_cache=0
if [ -s "${NODEFILEDIR}/distr" ]; then
  use_cache=1
  SAVE_CACHE=0
fi

#if ((SAVE_CACHE == 1)); then
#  echo "[INFO] $FUNCNAME: Save 'distr' file cache." >&2
#fi

#-------------------------------------------------------------------------------
# Set up node names

if [ "$NODELIST" = '-' ]; then
  : # do nothing
elif [ "$NODELIST" = '(0)' ]; then
  local n
  local p
  if ((use_cache == 0)); then
    for n in $(seq $NNODES); do
      for p in $(seq $PPN); do
        node[$(((n-1)*PPN+p))]="($((n-1)))"
        if ((SAVE_CACHE == 1)); then
          echo "node[$(((n-1)*PPN+p))]=\"($((n-1)))\"" >> $NODEFILEDIR/distr
        fi
      done
    done
  fi # ((use_cache == 0))
else
  read_nodefile_pbs "$NODELIST"
fi

#-------------------------------------------------------------------------------
# Set up member names and also get the number of members

fmember=0
for iname in $MEMBERS; do
  fmember=$((fmember+1))
  name_m[$fmember]=$iname
done

if ((CYCLE == 0)); then
  set_mem_np $fmember $SCALE_NP_TOT $SCALE_NP_TOT
  set_mem2node $fmember 0 0
  CYCLE=$((parallel_mems / fmember))
  if ((CYCLE < 1)); then
    CYCLE=1
  fi
  cycle_auto=$CYCLE
fi

for c in $(seq 2 $CYCLE); do
  for m in $(seq $fmember); do
    name_m[$(((c-1)*fmember+m))]=${name_m[$m]}
  done
done

fmembertot=$((fmember * CYCLE))

#-------------------------------------------------------------------------------
# Set up the distribution of members on nodes

set_mem_np $fmembertot $SCALE_NP_TOT $SCALE_NP_TOT

set_mem2node $fmembertot $use_cache $SAVE_CACHE

if ((use_cache == 1)); then
#  echo "[INFO] $FUNCNAME: Use 'distr' file cache." >&2
  . ${NODEFILEDIR}/distr
fi

#-------------------------------------------------------------------------------
# Create nodefiles

if [ "$NODELIST" != '-' ] && [ "$NODEFILEDIR" != '-' ]; then
#  echo "[INFO] $FUNCNAME: Save 'proc', 'node' files." >&2
  local p
  for p in $(seq $totalnp); do
    echo ${node[${proc2node[$p]}]} >> $NODEFILEDIR/proc
  done
  for n in $(seq $NNODES_APPAR); do
    echo ${node[$n]} >> $NODEFILEDIR/node
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

read_nodefile_pbs () {
#-------------------------------------------------------------------------------
# Parse the PBS-type nodefile.
# check if it is consistent to the $NNODES_APPAR and $PPN_APPAR settings and get the node names.
#
# Usage: read_nodefile_pbs NODEFILE
#
#   NODEFILE  PBS-type Nodefile for mpiexec
#
# Other input variables:
#   $NNODES_APPAR  Apparent number of total nodes
#   $PPN_APPAR     Apparent number of processes per node
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
  ippn=`cat $NODEFILE | grep -Fx $inode | wc -l`
  if ((ippn != PPN_APPAR)); then
    echo "[Error] $FUNCNAME: Number of processes per node in \$NODEFILE" >&2
    echo "          is not consistent to the setting in 'config.main'" >&2
    exit 1
  fi
done
if ((n != NNODES_APPAR)); then
  echo "[Error] $FUNCNAME: Number of nodes in \$NODEFILE" >&2
  echo "          is not consistent to the setting in 'config.main'" >&2
  exit 1
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
