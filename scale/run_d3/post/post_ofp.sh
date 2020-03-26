#!/bin/bash

#
# SNO is executed by a bulk job
#

mydir=`dirname $0`
POSTDIR=`cd $mydir ; pwd`

cd $POSTDIR/..

. config.main || exit $?
. src/func_util.sh || exit $?


PARENT_REF_TIME=$1
STIME=$2
FCSTLEN=$3

time_int=600 ### plot interval

ntime=`expr $FCSTLEN \/ $time_int`

RUNDIR="${TMP}/../sno_grads_d3_ref_${PARENT_REF_TIME}_${STIME}"

# Which domain do you want to convert?
DOM=3 ### fixed

OUTDIR="${OUTPUT}/${EXP3}/ref_${PARENT_REF_TIME}"

cd $POSTDIR

PPN=64 # Process per node

# Output file (X & Y process number) for each member
NP_OFILE_X=1
NP_OFILE_Y=1

# Do not edit!
NP_OFILE=$((${NP_OFILE_X} * ${NP_OFILE_Y})) # Output file (process number) for each member

#SNO_MEMBERS=${MEMBER}
SNO_MEMBERS=1
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean mdet" # All members + mean + mdet
#SNO_MEM_L=$(seq -f %04g ${SNO_MEMBERS})" mean" # All members + mean
SNO_MEM_L="mean"

# Total SNO processes  
NP_TOTAL=$((${SNO_MEMBERS} * ${NP_OFILE}))

# Convert variables (Other variables will NOT be included in converted files)
#VARS='"U", "V", "W", "T", "RH" , "QV",  "QC", "QI", "QR", "QS", "QG","QHYD", "PREC", "MSLP", "SFC_TEMP","U10", "V10", "T2"'
VARS='"T","PREC", "MSLP", "SFC_TEMP","U10", "V10", "T2"'

#VARS='"Uprs", "Vprs", "Gprs", "Tprs" , "QVprs", "QHYDprs", "PREC", "MSLP", "SFC_TEMP","U10", "V10", "T2"'

TOPO=0 # Process topography file? # 1: Yes, 0: No
if (( TOPO > 0 )) ; then
  VARS='"TOPO"'
  SNO_MEM_L="mean"
fi

# Get total SNO NODEs
if (( NP_TOTAL < PPN )) ; then
  SNO_NODE=1
else
  SNO_NODE=$((${NP_TOTAL} / ${PPN}))
fi

###############################

# Path for SNO binary
SNOBIN_ORG=${TOPDIR}/scale_develop/bin/sno
SNOBIN=${RUNDIR}/sno
if [ ! -e ${SNOBIN_ORG} ] ; then
  echo "No SNO binary!"
  exit
fi

###############################


rm -rf ${RUNDIR}

mkdir -p ${RUNDIR}/conf
mkdir -p ${RUNDIR}/log

# copy binary 
cp ${SNOBIN_ORG} ${SNOBIN}

conf_bulk="${RUNDIR}/conf/bulk_sno.conf"

cnt=0
for mem in  ${SNO_MEM_L} # member loop
do 
  cnt=$((${cnt} + 1))
  echo $mem

  mkdir -p ${RUNDIR}/grads/${mem}

  SNO_OUTPUT=${OUTDIR}/${STIME}/fcst_sno_np$(printf %05d ${NP_OFILE})/${mem}
  SNO_BASENAME_IN="${OUTDIR}/${STIME}/fcst/${mem}/history"
  SNO_BASENAME_OUT="${SNO_OUTPUT}/history"

  if (( TOPO > 0 )) ; then
    SNO_OUTPUT=${OUTDIR}/const/topo_sno_np$(printf %05d ${NP_OFILE})
    SNO_BASENAME_IN="${OUTDIR}/const/topo/topo"
    SNO_BASENAME_OUT="${SNO_OUTPUT}/topo"
  fi

  if [ ! -e ${SNO_OUTPUT} ] ; then
    mkdir -p ${SNO_OUTPUT}
  fi

  conf="${RUNDIR}/conf/sno_${mem}.conf"

cat << EOF >> $conf
&PARAM_IO
 IO_LOG_BASENAME = "log/LOG_${mem}",
 IO_LOG_ALLNODE = .false.,
 IO_LOG_SUPPRESS = .false.,
 IO_LOG_NML_SUPPRESS = .false.,
/

&PARAM_SNO
 basename_in  = "${SNO_BASENAME_IN}",
 basename_out  = "${SNO_BASENAME_OUT}",
 vars         = ${VARS},
 nprocs_x_out = ${NP_OFILE_X},
 nprocs_y_out = ${NP_OFILE_Y},
 output_grads = .true.,
 output_gradsctl = .false.,
 output_single = .true.,
 dirpath_out="grads/${mem}",
/


EOF

  ln -s ${conf} ${conf_bulk}.${cnt}

done # member loop



jobsh="${RUNDIR}/job_sno.sh"

cat << EOF >> $jobsh
#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L node=${SNO_NODE}
#PJM -L elapse="00:30:00"
#PJM --mpi proc=${NP_TOTAL}
#PJM --omp thread=1
#PJM -g $(echo $(id -ng))

rm -f machinefile
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
 for ippn in \$(seq ${PPN}); do
    echo "\$inode" >> machinefile
  done
done
module load hdf5
module load netcdf
module load netcdf-fortran

export FORT_FMT_RECL=400

export HFI_NO_CPUAFFINITY=1
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_preferred,,
export I_MPI_FABRICS_LIST=tmi
unset KMP_AFFINITY

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=${NPIN}
export I_MPI_PERHOST=${PPN}
export KMP_HW_SUBSET=1t


ulimit -s unlimited


echo "[\$(date "+%Y/%m/%d %H:%M:%S")] Start SNO"
mpirun -np ${NP_OFILE} ${SNOBIN} ${conf_bulk}.\${PJM_BULKNUM}
echo "[\$(date "+%Y/%m/%d %H:%M:%S")] End SNO"

EOF

echo 'convert...'
cd ${RUNDIR}
pjsub --bulk --sparam 1-${cnt} job_sno.sh 
cd $POSTDIR

jobid=$(grep 'pjsub Job' post.log.${PARENT_REF_TIME}.${STIME} | cut -d ' ' -f6)
job_end_check_PJM $jobid
res=$?

STIMEf="${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:0:2}"
export LANG=en_us
tstamp=`date -d "$STIMEf" +%H:%mZ%d%b%Y | tr '[a-z]' '[A-Z]'`

echo 'merge...'
for mem in  ${SNO_MEM_L} ;do # member loop
cp $POSTDIR/ctl/*.ctl $RUNDIR/grads/${mem}/
cd $RUNDIR/grads/$mem
sed -i -e "s/<--DATE-->/$tstamp/g" *.ctl
sed -i -e "s/<--NT-->/$ntime/g" *.ctl
grads -bcl "$POSTDIR/merge.gs ${ntime}" 
mkdir -p $OUTDIR/$STIME/fcstgp/$mem
cp history.ctl $OUTDIR/$STIME/fcstgp/$mem
cp history.grd $OUTDIR/$STIME/fcstgp/$mem
done


cd $POSTDIR/plot_grads
[ -f plot.lock ] && echo 'wait...'
while [ -f plot.lock ] ;do
 sleep 60
done

echo 'plot...'
echo $STIME > plot.lock 
 [ ! -z "`ls out/`" ] && rm out/*
 grads -bcl "plot_driver_d3_10min.gs $OUTDIR/$STIME/fcstgp/$mem/history.ctl 1 25 1" &> plot.log 
 mkdir -p $OUTDIR/$STIME/fcstgpi/$mem 
 mv out/*.png $OUTDIR/$STIME/fcstgpi/$mem/
 rm plot.lock

echo 'done.'

cd $POSTDIR


