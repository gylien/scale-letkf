#!/bin/bash
#===============================================================================

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------
PLACE=Saitama
DX=$1
STIME=$2

NCYCLE=120
WTIME_L="01:30:00"
NMEM=50
DACYCLE_RUN_FCST_TIME=1800
#MAX_DACYCLE_RUN_FCST=$NCYCLE
MAX_DACYCLE_RUN_FCST=110
NUM_DACYCLE_FCST_MEM=10

intv_sec=`expr \( $NCYCLE - 1 \) \* 30`
STIME_in="${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}:${STIME:10:2}:${STIME:12:2}"
ETIME=`date -d "${intv_sec} second ${STIME_in}" +'%Y%m%d%H%M%S'`

[ "$DX" == "1km" ]            && NNODES=`expr \( $NMEM + 2 + $NUM_DACYCLE_FCST_MEM \) ` ### 1km / 64domain
[ "$DX" == "500m_verysmall" ] && NNODES=`expr \( $NMEM + 2 + $NUM_DACYCLE_FCST_MEM \) ` ### 500m / 64domain
[ "$DX" == "500m_small" ]     && NNODES=`expr \( \( $NMEM + 2 + $NUM_DACYCLE_FCST_MEM \) \* 4 \) ` ### 500m / 256domain
[ "$DX" == "500m" ]           && NNODES=`expr \( \( $NMEM + 2 + $NUM_DACYCLE_FCST_MEM \) \* 16 \) ` ### 500m / 1024domain

#-------------------------------------------------------------------------------

CONFIG="${PLACE}/d4_${DX}"
if [ ! -d config/${CONFIG} ];then
 echo "EXPTYPE "$EXPTYPE" not supported !"
 exit
fi

cp config/${CONFIG}/config.* .

cat config.main.ofp | \
    sed -e "s/<MEMBER>/${NMEM}/g" | \
    sed -e "s/<NNODES>/${NNODES}/g"| \
    sed -e "s#<INDIR>#\${OUTDIR}#g" | \
    sed -e "s#<BGDIR>#\${DIR}/run/bgdata#g"  \
    > config.main

cat config/${CONFIG}/config.cycle | \
    sed -e "s/<STIME>/${STIME}/g" | \
    sed -e "s/<ETIME>/${ETIME}/g" | \
    sed -e "s/<WTIME_L>/${WTIME_L}/g" | \
    sed -e "s/<MEMBER>/${NMEM}/g" | \
    sed -e "s/<MAX_DACYCLE_RUN_FCST>/${MAX_DACYCLE_RUN_FCST}/g" | \
    sed -e "s/<NUM_DACYCLE_FCST_MEM>/${NUM_DACYCLE_FCST_MEM}/g" | \
    sed -e "s/<DACYCLE_RUN_FCST_TIME>/${DACYCLE_RUN_FCST_TIME}/g"  \
    > config.cycle


. config.main || exit $?

#-------------------------------------------------------------------------------
### prepare latest init and boundary files
##./prep.sh next 
 ./prep.sh $STIME 

#-------------------------------------------------------------------------------

./cycle_ofp.sh > cycle_ofp.log 2>&1 || exit $?

intv_sec_h=`expr \( $NCYCLE \) \* 30`
ETIMEh=`date -d "${intv_sec_h} second ${STIME_in}" +'%Y%m%d%H%M%S'`

./move_restart.sh $STIME $ETIMEh

./store_images.sh dacycle_${DX}_${STIME} $OUTDIR/${STIME}/dafcst &> log_store_images & 
#
#-------------------------------------------------------------------------------

  jobname="cycle_${SYSNAME}"
  jobid=$(grep 'pjsub Job' cycle_ofp.log | cut -d ' ' -f6)
  logdir="$OUTDIR/exp/${jobid}_cycle_${STIME}"
  stdout="$logdir/job.o"
  stderr="$logdir/job.e"
  jobinfo="$logdir/job.i"

#-------------------------------------------------------------------------------

rm -f cycle_job.sh
rm -f ${jobname}.?${jobid}

mkdir -p exp
rm -f exp/*
ln -s $OUTDIR/exp/${jobid}_cycle_${STIME} exp

#-------------------------------------------------------------------------------
