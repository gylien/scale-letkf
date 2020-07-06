#!/bin/bash
#===============================================================================
#
#  Common functions for 'cycle/fcst' jobs
#
#===============================================================================

staging_list_common_static () {
#-------------------------------------------------------------------------------
# Usage: staging_list_common_static JOBTYPE
#
#   JOBTYPE  Job type (cycle/fcst)
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local JOBTYPE="$1"

#-------------------------------------------------------------------------------
# common variables

for d in $(seq $DOMNUM); do
  if ((PNETCDF == 1)); then
    mem_np_[$d]=1
  else
    mem_np_[$d]=${SCALE_NP[$d]}
  fi
done



## H08
#if [ "$JOBTYPE" = 'cycle' ]; then
#  if [ -e "${RTTOV_COEF}" ] && [ -e "${RTTOV_SCCOEF}" ]; then
#    cat >> ${STAGING_DIR}/${STGINLIST_CONSTDB} << EOF
#${RTTOV_COEF}|dat/rttov/rtcoef_himawari_8_ahi.dat
#${RTTOV_SCCOEF}|dat/rttov/sccldcoef_himawari_8_ahi.dat
#EOF
#  fi
#fi

if [ "$TOPO_FORMAT" != 'prep' ]; then
  mkdir -p $TMP/dat/topo/${TOPO_FORMAT}
  ln -sf ${DATADIR}/topo/${TOPO_FORMAT}/Products $TMP/dat/topo/${TOPO_FORMAT}/Products
  #echo "${DATADIR}/topo/${TOPO_FORMAT}/Products/|dat/topo/${TOPO_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  mkdir -p $TMP/dat/landuse/${LANDUSE_FORMAT}
  ln -sf ${DATADIR}/landuse/${LANDUSE_FORMAT}/Products $TMP/dat/landuse/${LANDUSE_FORMAT}/Products
  #echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products/|dat/landuse/${LANDUSE_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

#-------------------------------------------------------------------------------
# observations

if [ "$JOBTYPE" = 'cycle' ]; then
  mkdir -p $TMP/obs

  time=$(datetime $STIME $LCYCLE s)
  while ((time <= $(datetime $ETIME $LCYCLE s))); do
    for iobs in $(seq $OBSNUM); do
      if [ "${OBSNAME[$iobs]}" != '' ] && [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
        #echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs.${OBSNAME[$iobs]}_${time}.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
        ln -sf ${OBS}/${OBSNAME[$iobs]}_${time}.dat $TMP/obs/${OBSNAME[$iobs]}_${time}.dat
      fi
    done
    time=$(datetime $time $LCYCLE s)
  done
fi

#-------------------------------------------------------------------------------
# create empty directories

#cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#|mean/
#|log/
#EOF
#
#if [ "$JOBTYPE" = 'cycle' ]; then
#  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
#|sprd/
#EOF
#fi

#-------------------------------------------------------------------------------
# time-invariant outputs

#-------------------
# stage-in
#-------------------

# domain catalogue
#-------------------
#if ((BDY_FORMAT == 1)); then
#  if [ -s "$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt" ]; then
#    pathin="$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt"
#    path="latlon_domain_catalogue.bdy.txt"
#    echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#  else
#    echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
#    echo "        '$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt'" >&2
#    exit 1
#  fi
#fi

#-------------------
# stage-out
#-------------------

# domain catalogue
#-------------------
if ((LOG_OPT <= 3)); then
  for d in $(seq $DOMNUM); do
    path="latlon_domain_catalogue.d$(printf $DOMAIN_FMT $d).txt"
    pathout="${OUTDIR[$d]}/const/log/latlon_domain_catalogue.txt"
    echo "${pathout}|${path}|1" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+1))]}
  done
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

config_file_scale_launcher () {
#-------------------------------------------------------------------------------
# Generate the launcher configuration files for scale_pp/scale_init/scale
#
# Usage: config_file_scale_launcher MODEL_NAME CONF_NAME
#
#   JOBTYPE     Job type (cycle/fcst)
#   MODEL_NAME  (scale-rm_pp/scale-rm_init/scale-rm)
#   CONF_NAME   (pp/init/run)
#   MEMBER_RUN  Number of members needed to run
#
# Other input variables:
#   $nitmax
#   $time
#   $SCRP_DIR
#   $MEMBER
#   $mtot
#   $PPN_APPAR
#   $mem_nodes
#   $DOMNUM
#   $PRC_DOMAINS_LIST
#   $STAGING_DIR
#-------------------------------------------------------------------------------

local JOBTYPE="$1"; shift
local MODEL_NAME="$1"; shift
local CONF_NAME="$1"; shift
local MEMBER_RUN="$1"

#-------------------------------------------------------------------------------

local it
local conf_file

if [ "$JOBTYPE" = 'cycle' ]; then
  CONF_FILES_SEQNUM='.false.'
else
  CONF_FILES_SEQNUM='.true.'
fi

DET_RUN_TF='.false.'
if ((DET_RUN == 1)); then
  DET_RUN_TF='.true.'
fi

for it in $(seq $nitmax); do
  if ((nitmax == 1)); then
    conf_file="${MODEL_NAME}_${time}.conf"
  else
    conf_file="${MODEL_NAME}_${time}_${it}.conf"
  fi
  echo "  $conf_file"
  cat $SCRP_DIR/config.nml.ensmodel | \
      sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
          -e "/!--MEMBER_RUN--/a MEMBER_RUN = $MEMBER_RUN," \
          -e "/!--MEMBER_ITER--/a MEMBER_ITER = $it," \
          -e "/!--CONF_FILES--/a CONF_FILES = \"${CONF_NAME}.d<domain>_${time}.conf\"," \
          -e "/!--CONF_FILES_SEQNUM--/a CONF_FILES_SEQNUM = $CONF_FILES_SEQNUM," \
          -e "/!--DET_RUN--/a DET_RUN = $DET_RUN_TF," \
          -e "/!--PPN--/a PPN = $PPN_APPAR," \
          -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
          -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
          -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
      > $TMP/${conf_file}
#  if ((stage_config == 1)); then
#    echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
#  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

config_file_save () {
#-------------------------------------------------------------------------------
# Save the runtime configuration files in $OUTDIR
#
# Usage: config_file_save [CONFIG_DIR]
#
#   CONFIG_DIR  Temporary directory of configuration files
#               '-': Use $TMPROOT
#
# Other input variables:
#   $TMPROOT
#-------------------------------------------------------------------------------

local CONFIG_DIR="${1:--}"

if [ "$CONFIG_DIR" = '-' ]; then
  CONFIG_DIR="$TMPROOT"
fi

#-------------------------------------------------------------------------------

mkdir -p ${OUTDIR[1]}/config

cp -fr $CONFIG_DIR/*.conf ${OUTDIR[1]}/config/

#-------------------------------------------------------------------------------
}

#===============================================================================
