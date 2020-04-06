#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in a K-computer job script and run it.
#
#  October 2014, created,                 Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_K.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='fcst'

RCSGRP=${RCSGRP:-"regular-flat"}
GNAME=${GNAME:-`id -ng`}

#===============================================================================
# Configuration

. config.main || exit $?
. config.${job} || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_${job}.sh || exit $?

#STAGING_DIR="$TMPSL/staging"
#NODEFILE_DIR="$TMPS/node"

#-------------------------------------------------------------------------------

if [ `echo $PRESET | cut -c 1` = 'K' ] && ((USE_TMP_LINK == 1)); then
  echo "[Error] $0: Wrong disk mode for K computer staged jobs." >&2
  exit 1
fi

if [ "$PRESET" = 'K_rankdir' ] && ((PNETCDF == 1)); then
  echo "[Error] When PNETCDF is enabled, 'K_rankdir' preset cannot be used." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------

statfile=${myname%.*}.stat.${PARENT_REF_TIME}.${STIME}

echo "prep" > $statfile
rm temp.lock

echo "[$(datetime_now)] Start $myname $@"

setting "$@" || exit $?

if [ "$CONF_MODE" = 'static' ]; then
  . src/func_common_static.sh || exit $?
  . src/func_${job}_static.sh || exit $?
fi

echo
print_setting || exit $?
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"
 
safe_init_tmpdir $TMPS || exit $?

#===============================================================================
# Determine the distribution schemes

echo "[$(datetime_now)] Determine the distribution schemes"

declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

safe_init_tmpdir $NODEFILE_DIR || exit $?
#distribute_fcst "$MEMBERS" $CYCLE "$NODELIST_TYPE" $NODEFILE_DIR || exit $?
distribute_fcst "$MEMBERS" $CYCLE - $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

#cp -L $SCRP_DIR/config.main $TMPS/config.main

cat $SCRP_DIR/config.main | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMPS/config.main
echo "NODEFILE_DIR=\"\$TMPROOT/node\"" >> $TMPS/config.main
echo "RUN_LEVEL=4" >> $TMPS/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMPS/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
if [ "$CONF_MODE" = 'static' ]; then
  staging_list_static || exit $?
  config_file_list $TMPS/config || exit $?
else
  staging_list || exit $?
fi

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

#${TMPS}/config.main|config.main
#
cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${SCRP_DIR}/config.rc|config.rc
${SCRP_DIR}/${job}.sh|${job}.sh
${SCRP_DIR}/src/|src/
EOF

if [ "$DISK_MODE" != "1" ] ; then
cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${SCRP_DIR}/config.${job}|config.${job}
${NODEFILE_DIR}/|node/
EOF
else
cp ${SCRP_DIR}/config.${job} $TMP/config.${job}
fi

if [ "$CONF_MODE" != 'static' ]; then
  echo "${SCRP_DIR}/${job}_step.sh|${job}_step.sh" >> ${STAGING_DIR}/${STGINLIST}
fi

#===============================================================================
# Staging in

stage_in server || exit $?


#===============================================================================
# Creat a job script

jobscrp="$TMP/${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

NPIN=`expr 255 / \( $PPN \) + 1`

cat > $jobscrp << EOF
#!/bin/sh -l
#PJM -L rscgrp=${RCSGRP}
#PJM -L node=${NNODES}
#PJM -L elapse=${TIME_LIMIT}
#PJM --mpi "proc=$((NNODES*PPN))"
#PJM --omp "thread=${THREADS}"
#PJM -g ${GNAME}


rm -f machinefile 
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
 for ippn in \$(seq $PPN); do
   echo "\$inode" >> machinefile
 done
done

module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5

ulimit -s unlimited

#export OMP_STACKSIZE=128m
export OMP_NUM_THREADS=1

export I_MPI_PIN_PROCESSER_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_PJOLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi

export I_MPI_PERHOST=${PPN}
export I_MPI_PIN_DOMAIN=${NPIN}

export KMP_HW_SUBSET=1t

export HFI_NO_CPUAFFINITY=1
unset KMP_AFFINITY

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" "$TIME_LIMIT" || exit \$?
EOF

#===============================================================================
# Run the job

echo "[$(datetime_now)] Run ${job} job on PJM"
echo

job_submit_PJM $jobscrp

jobid=$(grep 'pjsub Job' fcst_ofp.log.${PARENT_REF_TIME}.${STIME} | grep -v ++ | grep -v + | cut -d ' ' -f6)
echo "submit" $jobid > $statfile
echo

job_end_check_PJM $jobid
res=$?

[ $res = 77 ] && exit $res

echo "plot" > $statfile


#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_job.sh 'o e'

###if [ "$CONF_MODE" = 'static' ]; then
###  config_file_save $TMPS/config || exit $?
###fi

archive_log

#if ((CLEAR_TMP == 1)); then
#  safe_rm_tmpdir $TMPS
#  safe_rm_tmpdir $TMPSL
#fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
