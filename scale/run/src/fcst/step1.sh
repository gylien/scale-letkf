#!/bin/bash
#===============================================================================
#
#  SCALE forecasts - Initialization
#  August 2014, Guo-Yuan Lien
#
#===============================================================================

for c in $(seq $CYCLE); do
  if ((c == 1)); then
    stimef[$c]=$(datetime $STIME)
  else
    stimef[$c]=$(datetime ${stimef[$((c-1))]} $((LCYCLE*CYCLE_SKIP)) s)
  fi
  stimefr[$c]="$(datetime_fmt ${stimef[$c]})"
done

#if [ "$IF_EFSO" = '1' ]; then
#  FCSTLEN=$EFSOFLEN
#  FCSTOUT=$EFSOFOUT
#fi

#-------------------------------------------------------------------------------

mkdir -p $TMP1
tmpnode="$TMP1/node"
tmpstagein="$TMP1/stagein"
tmpstageout="$TMP1/stageout"

# Files that only a copy is needed
tmpprog="$STMP1/prog"
tmpdata="$STMP1/data"

# Files that are different in each node
if ((LTMPUSE == 0)); then
  tmpscale="$STMP1/scale"
  tmpverify="$STMP1/verify"
  tmpout="$STMP1/out"
  if ((STMPHOST == 0)); then
    ln -s $OUTDIR $tmpout
  fi
elif ((LTMPUSE == 1)); then
  tmpscale="$LTMP2/scale"
  tmpverify="$LTMP2/verify"
  tmpout="$LTMP2/out"
#  if [ "$IF_VERF" = '1' ]; then
#    echo "[Warning] $0: \$OBS, \$ANLGRDP, \$ANLGRDP2 need to be in a shared disk that can be seen from all computing nodes." 1>&2
#  fi
else
  echo "[Error] $0: Unsupported \$LTMPUSE setting." 1>&2
  exit 1
fi

echo "[$(datetime_now)] Start fcst.sh $@" 1>&2

#-------------------------------------------------------------------------------

if [ "$MEMBERS" = 'all' ] || [ "$MEMBERS" = 'mems' ]; then
  if [ "$MEMBERS" = 'all' ]; then
    MEMBERS="mean "
  else
    MEMBERS=''
  fi
  for m in $(seq $MEMBER); do
    MEMBERS="${MEMBERS}$(printf $MEMBER_FMT $m) "
  done
else
  memberstmp=''
  for m in $MEMBERS; do
    if [ "$m" = 'mean' ] || [ "$m" = 'sprd' ]; then
      memberstmp="$memberstmp$m "
    else
      memberstmp="$memberstmp$(printf $MEMBER_FMT $((10#$m))) "
      (($? != 0)) && exit 1
    fi
  done
  MEMBERS="$memberstmp"
fi

#-------------------------------------------------------------------------------

#===============================================================================
# Print job information

echo
echo " +----------------------------------------------------------------+"
echo " |                        SCALE-Forecasts                         |"
echo " +----------------------------------------------------------------+"
for s in $(seq $nsteps); do
  printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
done
echo " +----------------------------------------------------------------+"
echo
echo "  Number of cycles:         $CYCLE"
echo "  Forecast start time:"
for c in $(seq $CYCLE); do
  printf "    Cycle %-5s %s\n" "$c:" "${stimefr[$c]}"
done
echo
echo "  Forecast length:          $FCSTLEN s"
echo "  Output interval:          $FCSTOUT s"
echo

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m
if ((primary_script == 1)); then
  distribute_fcst "$MEMBERS" $CYCLE $MTYPE 1 1 machinefile
else
  distribute_fcst "$MEMBERS" $CYCLE $MTYPE 0 1 machinefile
fi

echo
echo "===================================================================="

#===============================================================================
