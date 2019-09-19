#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 12)); then
  cat >&2 << EOF

[pre_letkf_node.sh]

Usage: $0 MYRANK STIME ATIME TMPDIR OBSDIR MEM_NODES MEM_NP SLOT_START SLOT_END SLOT_BASE TOPO OBSOUT_OPT [ADAPTINFL SPRD_OUT RTPS_INFL_OUT NOBS_OUT]

  MYRANK      My rank number (not used)
  STIME
  ATIME       Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR      Temporary directory to run the program
  OBSDIR      Directory of SCALE data files
  MEM_NODES   Number of nodes for a member
  MEM_NP      Number of processes for a member
  SLOT_START  Start observation timeslots
  SLOT_END    End observation timeslots
  SLOT_BASE   The base slot
  TOPO        Basename of SCALE topography files
  OBSOUT_OPT
  ADAPTINFL
  SPRD_OUT
  RTPS_INFL_OUT
  NOBS_OUT

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
OBSDIR="$1"; shift
MEM_NODES="$1"; shift
MEM_NP="$1"; shift
SLOT_START="$1"; shift
SLOT_END="$1"; shift
SLOT_BASE="$1"; shift
TOPO="$1"; shift
OBSOUT_OPT="$1"; shift
ADAPTINFL="${1:-0}"; shift
SPRD_OUT="${1:-1}"; shift
RTPS_INFL_OUT="${1:-0}"; shift
NOBS_OUT="${1:-0}"

#===============================================================================

FILE_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
fi

OBS_IN_NAME_LIST=
for iobs in $(seq $OBSNUM); do
  if [ "${OBSNAME[$iobs]}" != '' ]; then
    OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'$OBSDIR/${OBSNAME[$iobs]}_${ATIME}.${OBSFOOT[$iobs]}', "
  fi

  rm -f $OBSDIR/${OBSNAME[$iobs]}_${STIME}.${OBSFOOT[$iobs]}

done

OBSDA_RUN_LIST=
for iobs in $(seq $OBSNUM); do
  if [ -n "${OBSOPE_SEPARATE[$iobs]}" ] && ((${OBSOPE_SEPARATE[$iobs]} == 1)); then
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.false., "
  else
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.true., "
  fi
done

DET_RUN_TF='.false.'
if ((DET_RUN == 1)); then
  DET_RUN_TF='.true.'
fi
OBSDA_IN='.false.'
if ((OBSOPE_RUN == 1)); then
  OBSDA_IN='.true.'
fi
OBSDEP_OUT_TF='.false.'
if ((OBSOUT_OPT <= 3)); then
  OBSDEP_OUT_TF='.true.'
fi
SPRD_OUT_TF='.true.'
if ((SPRD_OUT == 0)); then
  SPRD_OUT_TF='.false.'
fi
INFL_MUL_ADAPTIVE='.false.'
if ((ADAPTINFL == 1)); then
  INFL_MUL_ADAPTIVE='.true.'
fi
RTPS_INFL_OUT_TF='.false.'
if ((RTPS_INFL_OUT == 1)); then
  RTPS_INFL_OUT_TF='.true.'
fi
NOBS_OUT_TF='.false.'
if ((NOBS_OUT == 1)); then
  NOBS_OUT_TF='.true.'
fi

if ((PNETCDF == 1)); then
  HISTORY_IN_BASENAME="${TMPOUT}/${STIME}/hist/@@@@.history"
  GUES_IN_BASENAME="${TMPOUT}/${ATIME}/anal/@@@@.init"
  GUES_MEAN_INOUT_BASENAME="${TMPOUT}/${ATIME}/gues/mean.init"
  GUES_SPRD_OUT_BASENAME="${TMPOUT}/${ATIME}/gues/sprd.init"
  ANAL_OUT_BASENAME="${TMPOUT}/${ATIME}/anal/@@@@.init"
  INFL_MUL_IN_BASENAME="${TMPOUT}/${ATIME}/diag/infl"
  INFL_MUL_OUT_BASENAME="${TMPOUT}/${ATIME}/diag/infl"
  INFL_ADD_IN_BASENAME="${TMPOUT}/const/addi/@@@@.init"
  RELAX_SPREAD_OUT_BASENAME="${TMPOUT}/${ATIME}/diag/rtps"
  NOBS_OUT_BASENAME="${TMPOUT}/${ATIME}/diag/nobs"
else
  HISTORY_IN_BASENAME="${OUTDIR}/${STIME}/hist/@@@@/history"
  GUES_IN_BASENAME="${OUTDIR}/${ATIME}/anal/@@@@/init"
  GUES_MEAN_INOUT_BASENAME="${OUTDIR}/${ATIME}/gues/mean/init"
  GUES_SPRD_OUT_BASENAME="${OUTDIR}/${ATIME}/gues/sprd/init"
  ANAL_OUT_BASENAME="${OUTDIR}/${ATIME}/anal/@@@@/init"
  INFL_MUL_IN_BASENAME="${OUTDIR}/${ATIME}/diag/infl/init"
  INFL_MUL_OUT_BASENAME="${OUTDIR}/${ATIME}/diag/infl/init"
  INFL_ADD_IN_BASENAME="${OUTDIR}/const/addi/@@@@/init"
  RELAX_SPREAD_OUT_BASENAME="${OUTDIR}/${ATIME}/diag/rtps/init"
  NOBS_OUT_BASENAME="${OUTDIR}/${ATIME}/diag/nobs/init"
fi

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

cat $TMPDAT/conf/config.nml.ensmodel | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
        -e "/!--DET_RUN--/a DET_RUN = ${DET_RUN_TF}," \
        -e "/!--PPN--/a PPN = $PPN_APPAR," \
        -e "/!--MEM_NODES--/a MEM_NODES = $MEM_NODES," \
        -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = 1," \
        -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $MEM_NP," \
    > $TMPDIR/letkf.conf

cat $TMPDAT/conf/config.nml.letkf | \
    sed -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
        -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
        -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
        -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = \"${HISTORY_IN_BASENAME}\"," \
        -e "/!--SLOT_START--/a SLOT_START = $SLOT_START," \
        -e "/!--SLOT_END--/a SLOT_END = $SLOT_END," \
        -e "/!--SLOT_BASE--/a SLOT_BASE = $SLOT_BASE," \
        -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = $LTIMESLOT.D0," \
        -e "/!--OBSDA_IN--/a OBSDA_IN = $OBSDA_IN," \
        -e "/!--OBSDA_IN_BASENAME--/a OBSDA_IN_BASENAME = \"${TMPOUT}/${ATIME}/obsgues/@@@@/obsda.ext\"," \
        -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"${GUES_IN_BASENAME}\"," \
        -e "/!--GUES_MEAN_INOUT_BASENAME--/a GUES_MEAN_INOUT_BASENAME = \"${GUES_MEAN_INOUT_BASENAME}\"," \
        -e "/!--GUES_SPRD_OUT_BASENAME--/a GUES_SPRD_OUT_BASENAME = \"${GUES_SPRD_OUT_BASENAME}\"," \
        -e "/!--GUES_SPRD_OUT--/a GUES_SPRD_OUT = ${SPRD_OUT_TF}," \
        -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"${ANAL_OUT_BASENAME}\"," \
        -e "/!--ANAL_SPRD_OUT--/a ANAL_SPRD_OUT = ${SPRD_OUT_TF}," \
        -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--INFL_MUL_ADAPTIVE--/a INFL_MUL_ADAPTIVE = ${INFL_MUL_ADAPTIVE}," \
        -e "/!--INFL_MUL_IN_BASENAME--/a INFL_MUL_IN_BASENAME = \"${INFL_MUL_IN_BASENAME}\"," \
        -e "/!--INFL_MUL_OUT_BASENAME--/a INFL_MUL_OUT_BASENAME = \"${INFL_MUL_OUT_BASENAME}\"," \
        -e "/!--INFL_ADD_IN_BASENAME--/a INFL_ADD_IN_BASENAME = \"${INFL_ADD_IN_BASENAME}\"," \
        -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RTPS_INFL_OUT_TF}," \
        -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"${RELAX_SPREAD_OUT_BASENAME}\"," \
        -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
        -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"${NOBS_OUT_BASENAME}\"," \
        -e "/!--OBSDEP_OUT--/a OBSDEP_OUT = ${OBSDEP_OUT_TF}," \
        -e "/!--OBSDEP_OUT_BASENAME--/a OBSDEP_OUT_BASENAME = \"${TMPOUT}/${ATIME}/obs/obsdep\"," \
        -e "/!--H08_NOWDATE--/a H08_NOWDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--H08_RTTOV_COEF_PATH--/a H08_RTTOV_COEF_PATH = \"${TMPDIR}\"," \
        -e "/!--H08_OUTFILE_BASENAME--/a H08_OUTFILE_BASENAME = \"${TMPOUT}/${ATIME}/Him8/Him8_${ATIME}\"," \
    >> $TMPDIR/letkf.conf
#        -e "/!--H08_VBC_PATH--/a H08_VBC_PATH = \"${TMPOUT}/vbc\"," \

# Most of these parameters are not important for letkf
cat $TMPDAT/conf/config.nml.scale | \
    sed -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = ${IO_AGGREGATE}," \
        -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/PARAG.29\"," \
        -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/PARAPC.29\"," \
        -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/VARDATA.RM29\"," \
        -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPDAT_CONSTDB}/rad/cira.nc\"," \
        -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPDAT_CONSTDB}/rad/MIPAS\"," \
        -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
    >> $TMPDIR/letkf.conf


#if [ -f "${TMPOUT}/vbc/Him8_vbca_${STIME}.dat" ] ; then
#  cp ${TMPOUT}/vbc/Him8_vbca_${STIME}.dat ${TMPOUT}/vbc/Him8_vbcf.dat
#elif [ ! -d "${TMPOUT}/vbc" ] ; then
#  mkdir "${TMPOUT}/vbc"
#fi

if [ ! -d "${TMPOUT}/${ATIME}/Him8" ] ; then
  mkdir "${TMPOUT}/${ATIME}/Him8"
fi


#===============================================================================

exit 0
