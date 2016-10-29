#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 14)); then
  cat >&2 << EOF

[pre_letkf_node.sh]

Usage: $0 MYRANK STIME ATIME TMPDIR OBSDIR MEM_NODES MEM_NP SLOT_START SLOT_END SLOT_BASE TOPO ADAPTINFL RTPS_INFL_OUT NOBS_OUT [MEMBERSEQ]

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
  ADAPTINFL
  RTPS_INFL_OUT
  NOBS_OUT
  MEMBERSEQ

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
ADAPTINFL="$1"; shift
RTPS_INFL_OUT="$1"; shift
NOBS_OUT="$1"; shift
MEMBERSEQ=${1:-$MEMBER}

#===============================================================================

OBS_IN_NAME_LIST=
for iobs in $(seq $OBSNUM); do
  if [ "${OBSNAME[$iobs]}" != '' ]; then
#    ln -fs $OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat $TMPDIR/${OBSNAME[$iobs]}.dat
    OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'$OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat', "
  fi
done

OBSDA_RUN_LIST=
for iobs in $(seq $OBSNUM); do
  if [ -n "${OBSOPE_SEPARATE[$iobs]}" ] && ((${OBSOPE_SEPARATE[$iobs]} == 1)); then
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.false., "
  else
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.true., "
  fi
done

OBSDA_IN='.false.'
if ((OBSOPE_RUN == 1)); then
  OBSDA_IN='.true.'
fi
INFL_MUL_ADAPTIVE='.false.'
if ((ADAPTINFL == 1)); then
  INFL_MUL_ADAPTIVE='.true.'
fi
RELAX_SPREAD_OUT='.false.'
if ((RTPS_INFL_OUT == 1)); then
  RELAX_SPREAD_OUT='.true.'
fi
NOBS_OUT_TF='.false.'
if ((NOBS_OUT == 1)); then
  NOBS_OUT_TF='.true.'
fi

# --- Parameter estimation (Tomita 2008) ----

if ((MYRANK == 0 )) && [ "$PARAM_EST" == "T" ] ; then

  if [ -e ${TMPDAT}/param/EPARAM_TOMITA_ANAL${STIME}.txt ] ; then
    PARAM_FILE="${TMPDAT}/param/EPARAM_TOMITA_ANAL${STIME}.txt"
  elif [ -e ${TMPOUT}/${STIME}/log/letkf/EPARAM_TOMITA_ANAL${STIME}.txt ] ; then
    PARAM_FILE="${TMPOUT}/${STIME}/log/letkf/EPARAM_TOMITA_ANAL${STIME}.txt"
  else
    echo "No parameter input!! check! ""${TMPOUT}/${STIME}/log/letkf/EPARAM_TOMITA_ANAL${STIME}.txt"
    echo "No parameter input!! check! ""${TMPDAT}/param/EPARAM_TOMITA_ANAL${STIME}.txt"
    echo "No parameter input!! check! "$PARAM_FILE
    exit 1
  fi
 
  if [ ! -e ${TMPOUT}/param ] ; then
    mkdir -p ${TMPOUT}/param
  fi
 
  cp $PARAM_FILE $TMPDIR/EPARAM_TOMITA_GUES.txt

fi

# -------------------------------------------
# -- Cloud dependent obs err --

BB_LIST="07 08 09 10 11 12 13 14 15 16"

if [ ! -e ${TMPDAT}/Him8 ] ; then
  mkdir ${TMPDAT}/Him8 
fi

for BB in ${BB_LIST} ; do
  CA_FILE1_A="${TMPDAT}/Him8/Him8_ERR_CA_A_B${BB}_${STIME}.dat"
  CA_FILE2_A="Him8_ERR_CA_A_B${BB}.dat"
  CA_FILE1_B="${TMPDAT}/Him8/Him8_ERR_CA_A_B${BB}_${STIME}.dat"
  CA_FILE2_B="Him8_ERR_CA_A_B${BB}.dat"
  if [ -e ${CA_FILE1_A} ] && [ -e ${CA_FILE1_B} ] ; then
    cp $CA_FILE1_A ${TMPDIR}/$CA_FILE2_A
    cp $CA_FILE1_B ${TMPDIR}/$CA_FILE2_B
  fi
done

#===============================================================================

cat $TMPDAT/conf/config.nml.letkf | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBERSEQ," \
        -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
        -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
        -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
        -e "/!--OBSDA_OUT--/a OBSDA_OUT = .true." \
        -e "/!--OBSDA_OUT_BASENAME--/a OBSDA_OUT_BASENAME = \"${TMPOUT}/${ATIME}/obsgues/@@@@/obsda\"," \
        -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = '${TMPOUT}/${STIME}/hist/@@@@/history'," \
        -e "/!--SLOT_START--/a SLOT_START = $SLOT_START," \
        -e "/!--SLOT_END--/a SLOT_END = $SLOT_END," \
        -e "/!--SLOT_BASE--/a SLOT_BASE = $SLOT_BASE," \
        -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = $LTIMESLOT.D0," \
        -e "/!--OBSDA_IN--/a OBSDA_IN = $OBSDA_IN," \
        -e "/!--OBSDA_IN_BASENAME--/a OBSDA_IN_BASENAME = \"${TMPOUT}/${ATIME}/obsgues/@@@@/obsda.ext\"," \
        -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"${TMPOUT}/${ATIME}/gues/@@@@/init\"," \
        -e "/!--GUES_OUT_MEAN_BASENAME--/a GUES_OUT_MEAN_BASENAME = \"${TMPOUT}/${ATIME}/gues/mean/init\"," \
        -e "/!--GUES_OUT_SPRD_BASENAME--/a GUES_OUT_SPRD_BASENAME = \"${TMPOUT}/${ATIME}/gues/sprd/init\"," \
        -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"${TMPOUT}/${ATIME}/anal/@@@@/init\"," \
        -e "/!--ANAL_OUT_MEAN_BASENAME--/a ANAL_OUT_MEAN_BASENAME = \"${TMPOUT}/${ATIME}/anal/mean/init\"," \
        -e "/!--ANAL_OUT_SPRD_BASENAME--/a ANAL_OUT_SPRD_BASENAME = \"${TMPOUT}/${ATIME}/anal/sprd/init\"," \
        -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--INFL_MUL_ADAPTIVE--/a INFL_MUL_ADAPTIVE = ${INFL_MUL_ADAPTIVE}," \
        -e "/!--INFL_MUL_IN_BASENAME--/a INFL_MUL_IN_BASENAME = \"${TMPOUT}/${ATIME}/diag/infl/init\"," \
        -e "/!--INFL_MUL_OUT_BASENAME--/a INFL_MUL_OUT_BASENAME = \"${TMPOUT}/${ATIME}/diag/infl/init\"," \
        -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RELAX_SPREAD_OUT}," \
        -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"${TMPOUT}/${ATIME}/diag/rtps/init\"," \
        -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
        -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"${TMPOUT}/${ATIME}/diag/nobs/init\"," \
        -e "/!--NNODES--/a NNODES = $NNODES," \
        -e "/!--PPN--/a PPN = $PPN," \
        -e "/!--MEM_NODES--/a MEM_NODES = $MEM_NODES," \
        -e "/!--MEM_NP--/a MEM_NP = $MEM_NP," \
    > $TMPDIR/letkf.conf

# These parameters are not important for letkf
cat $TMPDAT/conf/config.nml.scale >> $TMPDIR/letkf.conf

#===============================================================================

exit 0
