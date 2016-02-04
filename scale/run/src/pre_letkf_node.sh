#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of LETKF run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 11)); then
  cat >&2 << EOF

[pre_letkf_node.sh]

Usage: $0 MYRANK ATIME TMPDIR EXECDIR OBSDIR MEM_NODES MEM_NP SLOT_START SLOT_END SLOT_BASE TOPO MEMBERSEQ

  MYRANK      My rank number (not used)
  ATIME       Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR      Temporary directory to run the program
  EXECDIR     Directory of SCALE executable files
  OBSDIR      Directory of SCALE data files
  MEM_NODES   Number of nodes for a member
  MEM_NP      Number of processes for a member
  SLOT_START  Start observation timeslots
  SLOT_END    End observation timeslots
  SLOT_BASE   The base slot
  TOPO        Basename of SCALE topography files
  MEMBERSEQ

EOF
  exit 1
fi

MYRANK="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
EXECDIR="$1"; shift 
OBSDIR="$1"; shift
MEM_NODES="$1"; shift
MEM_NP="$1"; shift
SLOT_START="$1"; shift
SLOT_END="$1"; shift
SLOT_BASE="$1"; shift
TOPO="$1"; shift
MEMBERSEQ=${1:-$MEMBER}

#===============================================================================

#mkdir -p $TMPDIR
#rm -fr $TMPDIR/*

#ln -fs $EXECDIR/letkf $TMPDIR

OBS_IN_NAME_LIST=
for iobs in $(seq $OBSNUM); do
  if [ "${OBSNAME[$iobs]}" != '' ]; then
#    ln -fs $OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat $TMPDIR/${OBSNAME[$iobs]}.dat
    OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'$OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat', "
  fi
done

#===============================================================================

cat $TMPDAT/conf/config.nml.letkf | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBERSEQ," \
        -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
        -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
        -e "/!--OBSDA_IN_BASENAME--/a OBSDA_IN_BASENAME = \"${TMPOUT}/${ATIME}/obsgues/@@@@/obsda\"," \
        -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"${TMPOUT}/${ATIME}/gues/@@@@/init\"," \
        -e "/!--GUES_OUT_MEAN_BASENAME--/a GUES_OUT_MEAN_BASENAME = \"${TMPOUT}/${ATIME}/gues/mean/init\"," \
        -e "/!--GUES_OUT_SPRD_BASENAME--/a GUES_OUT_SPRD_BASENAME = \"${TMPOUT}/${ATIME}/gues/sprd/init\"," \
        -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"${TMPOUT}/${ATIME}/anal/@@@@/init\"," \
        -e "/!--ANAL_OUT_MEAN_BASENAME--/a ANAL_OUT_MEAN_BASENAME = \"${TMPOUT}/${ATIME}/anal/mean/init\"," \
        -e "/!--ANAL_OUT_SPRD_BASENAME--/a ANAL_OUT_SPRD_BASENAME = \"${TMPOUT}/${ATIME}/anal/sprd/init\"," \
        -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--SLOT_START--/a SLOT_START = $SLOT_START," \
        -e "/!--SLOT_END--/a SLOT_END = $SLOT_END," \
        -e "/!--SLOT_BASE--/a SLOT_BASE = $SLOT_BASE," \
        -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = $LTIMESLOT.D0," \
        -e "/!--NNODES--/a NNODES = $NNODES," \
        -e "/!--PPN--/a PPN = $PPN," \
        -e "/!--MEM_NODES--/a MEM_NODES = $MEM_NODES," \
        -e "/!--MEM_NP--/a MEM_NP = $MEM_NP," \
    > $TMPDIR/letkf.conf

# These parameters are not important for obsope
cat $TMPDAT/conf/config.nml.scale | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"LOG\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = 2014, 1, 1, 0, 0, 0," \
        -e "/!--TIME_DURATION--/a TIME_DURATION = $LTIMESLOT.D0," \
        -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = $LTIMESLOT.D0," \
        -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = $LTIMESLOT.D0," \
        -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = $LTIMESLOT.D0," \
        -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = .D0," \
        -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"init\"," \
        -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"restart\"," \
        -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo\"," \
        -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse\"," \
        -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"boundary\"," \
        -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"init_ocean\"," \
        -e "/!--HISTORY_DEFAULT_BASENAME--/a HISTORY_DEFAULT_BASENAME = \"history\"," \
        -e "/!--HISTORY_DEFAULT_TINTERVAL--/a HISTORY_DEFAULT_TINTERVAL = $LTIMESLOT.D0," \
        -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"monitor\"," \
    >> $TMPDIR/letkf.conf

mkdir -p $TMPOUT/${ATIME}/log/letkf

#===============================================================================

exit 0
