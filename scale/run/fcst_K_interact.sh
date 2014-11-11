#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a K-computer job script (interact) and run it.
#
#  November 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_K_interact.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]
#
#===============================================================================

cd "$(dirname "$0")"
myname1='fcst'

#===============================================================================
# Configuration

. config.all
(($? != 0)) && exit $?
. config.$myname1
(($? != 0)) && exit $?

. src/func_datetime.sh
. src/func_$myname1.sh

#-------------------------------------------------------------------------------

if ((TMPDAT_MODE >= 3 || TMPRUN_MODE >= 3 || TMPOUT_MODE >= 3)); then
  echo "[Error] $0: When using 'interact' resource group," >&2
  echo "        \$TMPDAT_MODE, \$TMPRUN_MODE, \$TMPOUT_MODE all need to be 1 or 2." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}"

jobscrp="${myname1}_job.sh"

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $(basename $0) $@"
echo

for vname in DIR OUTDIR ANLWRF OBS OBSNCEP MEMBER NNODES PPN \
             FCSTLEN FCSTOUT EFSOFLEN EFSOFOUT OUT_OPT \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP; do
  printf '  %-10s = %s\n' $vname "${!vname}"
done

echo

#===============================================================================
# Creat a job script

echo "[$(datetime_now)] Create a job script '$jobscrp'"

cat > $jobscrp << EOF
#!/bin/sh
##PJM -g ra000015
#PJM -N ${myname1}_${SYSNAME}
#PJM -s
#PJM --rsc-list "node=${NNODES_real}"
#PJM --rsc-list "elapse=${TIME_LIMIT}"
#PJM --mpi "shape=${NNODES_real}"
#PJM --mpi "proc=$NNODES"
#PJM --mpi assign-online-node

. /work/system/Env_base
export OMP_NUM_THREADS=${THREADS}
export PARALLEL=${THREADS}

cd $(pwd)

./${myname1}.sh $STIME "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP"
EOF

echo "[$(datetime_now)] Entering the interactive mode..."
echo

pjsub --interact --sparam wait-time=3600 $jobscrp

#===============================================================================

echo
echo "[$(datetime_now)] Finish $(basename $0) $@"

exit 0
