#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#  August  2014, modified from GFS-LETKF, Guo-Yuan Lien
#  October 2014, modified,                Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP PREP]
#    (see variable explanation in 'config.fcst')
#
#  Use settings:
#    config.all
#    config.fcst (optional)
#
#===============================================================================

cd "$(dirname "$0")"
myname=$(basename "$0")
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.all
(($? != 0)) && exit $?

if [ -f config.$myname1 ]; then
  . config.$myname1
fi

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh

. src/func_$myname1.sh

#-------------------------------------------------------------------------------

nsteps=5
declare -a stepname
declare -a stepfunc

stepname[1]='Prepare topo files'
stepfunc[1]='topo'
stepname[2]='Prepare landuse files'
stepfunc[2]='landuse'
stepname[3]='Prepare boundary files'
stepfunc[3]='boundary'
stepname[4]='Perturb boundaries'
stepfunc[4]='pertbdy'
stepname[5]='Run ensemble forecasts'
stepfunc[5]='ensfcst'
stepname[6]='Run verification'
stepfunc[6]='verf'

#-------------------------------------------------------------------------------

USAGE="
[$myname] Run ensemble forecasts and (optional) verifications.

Configuration files:
  config.all
  config.$myname1 (optional)

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP PREP]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  MEMBERS     List of forecast members ('mean' for ensemble mean)
               all:     Run all members including ensemble mean (default)
               mems:    Run all members but not including ensemble mean
               '2 4 6': Run members 2, 4, 6
  CYCLE       Number of forecast cycles run in parallel
               (default: 1)
  CYCLE_SKIP  Run forecasts every ? cycles
               (default: 1)
  IF_VERF     Run verification? [Not finished!]
               0: No (default)
               1: Yes
              * to run the verification, a shared disk storing observations
                and reference model analyses needs to be used
  IF_EFSO     Use EFSO forecast length and output interval? [Not finished!]
               0: No (default)
               1: Yes
  ISTEP       The initial step in the first cycle from which this script starts
               (default: the first step)
  FSTEP       The final step in the last cycle by which this script ends
               (default: the last step)
  PREP        Is this a preparation run?
               1: Yes (default when \$MACHINE_TYPE = 10)
               0: No  (default otherwise)
"

#-------------------------------------------------------------------------------

if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
  echo "$USAGE"
  exit 0
fi
if (($# < 1)) && [ -z "$STIME" ]; then
  echo "$USAGE" >&2
  exit 1
fi

STIME=$(datetime ${1:-$STIME})
ETIME=$(datetime ${2:-$ETIME})
MEMBERS=${3:-$MEMBERS}
CYCLE=${4:-$CYCLE}
CYCLE_SKIP=${5:-$CYCLE_SKIP}
IF_VERF=${6:-$IF_VERF}
IF_EFSO=${7:-$IF_EFSO}
ISTEP=${8:-$ISTEP}
FSTEP=${9:-$FSTEP}
PREP=${10:-$PREP}

ETIME=${ETIME:-$STIME}
MEMBERS="${MEMBERS:-all}"
CYCLE=${CYCLE:-1}
CYCLE_SKIP=${CYCLE_SKIP:-1}
IF_VERF=${IF_VERF:-0}
IF_EFSO=${IF_EFSO:-0}
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-$nsteps}
if ((MACHINE_TYPE == 10)); then
  PREP=${PREP:-1}
else
  PREP=${PREP:-0}
fi

#-------------------------------------------------------------------------------

mkdir -p $LOGDIR

if ((MACHINE_TYPE == 10 && PREP == 0)); then
  exec 2>> $TMP/runlog/${myname1}.err
else
  sleep 0.01s
####exec 2>> $LOGDIR/${myname1}.err
fi
echo "[$(datetime_now)] Start $myname" >&2

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT FOUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '                      %-10s = %s\n' $vname ${!vname} >&2
done

#===============================================================================
# More configuration

safe_init_tmpdir $TMP

#-------------------------------------------------------------------------------

# more variable difinition???

CYCLE_FMT='%04d'


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

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m

#if ((MACHINE_TYPE == 10)); then # K-computer: create temporary nodefiles in local disk, then stage in
#  safe_init_tmpdir $TMPS/node
#  distribute_fcst "$MEMBERS" $CYCLE machinefile $TMPS/node
#else


if ((MACHINE_TYPE == 10)); then

  if ((PREP == 1)); then
    NODEFILE_DIR="$TMPS/node"
    safe_init_tmpdir $NODEFILE_DIR
    distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR
  else
    distribute_fcst "$MEMBERS" $CYCLE - -
  fi

else

  safe_init_tmpdir $NODEFILE_DIR
  distribute_fcst "$MEMBERS" $CYCLE machinefile $NODEFILE_DIR

fi


#fi

#===============================================================================

if ((MACHINE_TYPE == 10)); then

  if ((PREP == 1)); then

    STAGING_DIR="$TMPS/staging"

    init

    ## create K job script

    ## submit job

    ## wait for job to finish
    
  fi  

else
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  init

####
#safe_init_tmpdir $TMPDAT
#safe_init_tmpdir $TMPOUT
#safe_init_tmpdir $TMPRUN
####

  pdbash node all $SCRP_DIR/src/stage_in.sh
fi

#===============================================================================

if ((PREP == 0)); then

#===============================================================================
# Run cycles of forecasts


declare -a stimes
declare -a stimesfmt
lcycles=$((LCYCLE * CYCLE_SKIP))
s_flag=1
e_flag=0
time=$STIME

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

  for c in $(seq $CYCLE); do
    time2=$(datetime $time $((lcycles * (c-1))) s)
    if (($(datetime $time2 $lcycles s) > ETIME)); then
      e_flag=1
    fi

    if ((time2 <= ETIME)); then
      stimes[$c]=$time2
      stimesfmt[$c]="$(datetime_fmt ${stimes[$c]})"
      rcycle=$c  # The "real" number of cycles
    else
      stimes[$c]=
      stimesfmt[$c]=
    fi
  done

#-------------------------------------------------------------------------------
# Write the header of the log file

if ((MACHINE_TYPE == 10)); then
  exec > $TMP/runlog/${myname1}_${stimes[1]}.log
else
  sleep 0.01s
####  exec > $LOGDIR/${myname1}_${stimes[1]}.log
fi

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                        SCALE-Forecasts                         |"
  echo " +----------------------------------------------------------------+"
  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then
      printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
    fi
  done
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Number of cycles:         $rcycle"
  echo "  Forecast start time:"
  for c in $(seq $rcycle); do
    printf "    Cycle %-5s %s\n" "$c:" "${stimesfmt[$c]}"
  done
  echo
  echo "  Forecast length:          $FCSTLEN s"
  echo "  Output interval:          $FCSTOUT s"
  echo
  echo "  Nodes used:               $NNODES"
#  if ((MTYPE == 1)); then
    for n in `seq $NNODES`; do
      echo "    ${node[$n]}"
    done
#  fi
  echo
  echo "  Processes per node:       $PPN"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Number of members:        $fmember"
  for c in $(seq $rcycle); do
    echo "    Cycle $c:"
    for m in $(seq $fmember); do
      mm=$(((c-1) * fmember + m))
      echo "      ${name_m[$m]}: ${node_m[$mm]}"
    done
  done
  echo
  echo "===================================================================="

#-------------------------------------------------------------------------------
# Call functions to run the job

  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

      echo "[$(datetime_now)] ${stepname[$s]}" >&2
      echo
      printf " %2d. %-55s\n" $s "${stepname[$s]}"

      ${stepfunc[$s]}

      echo
      echo "===================================================================="

    fi
  done

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo
  echo " +----------------------------------------------------------------+"
  echo " |             SCALE-Forecasts successfully completed             |"
  echo " +----------------------------------------------------------------+"
  echo

#-------------------------------------------------------------------------------

  time=$(datetime $time $((lcycles * CYCLE)) s)
  s_flag=0

#-------------------------------------------------------------------------------
done
#-------------------------------------------------------------------------------

#===============================================================================

fi # ((PREP == 0))

#===============================================================================
# Finalization

#safe_rm_tmpdir $TMP
#safe_rm_tmpdir $TMPS

echo "[$(datetime_now)] Finish fcst.sh $@" >&2

#===============================================================================

exit 0
