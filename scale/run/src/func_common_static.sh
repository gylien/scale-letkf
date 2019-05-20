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

#-------------------------------------------------------------------------------
# executable files

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${COMMON_DIR}/pdbash|pdbash
${COMMON_DIR}/datetime|datetime
${ENSMODEL_DIR}/scale-rm_pp_ens|scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|scale-rm_ens
EOF

if [ "$JOBTYPE" = 'cycle' ]; then
  if ((DACYCLE == 1)); then
    cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${ENSMODEL_DIR}/../dacycle/dacycle|dacycle
EOF
  fi
  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${OBSUTIL_DIR}/obsope|obsope
${LETKF_DIR}/letkf|letkf
EOF
fi

#-------------------------------------------------------------------------------
# database

cat >> ${STAGING_DIR}/${STGINLIST_CONSTDB} << EOF
${SCALEDIR}/scale-rm/test/data/rad/cira.nc|dat/rad/cira.nc
${SCALEDIR}/scale-rm/test/data/rad/PARAG.29|dat/rad/PARAG.29
${SCALEDIR}/scale-rm/test/data/rad/PARAPC.29|dat/rad/PARAPC.29
${SCALEDIR}/scale-rm/test/data/rad/rad_o3_profs.txt|dat/rad/rad_o3_profs.txt
${SCALEDIR}/scale-rm/test/data/rad/VARDATA.RM29|dat/rad/VARDATA.RM29
${SCALEDIR}/scale-rm/test/data/rad/MIPAS/|dat/rad/MIPAS/
${SCALEDIR}/scale-rm/test/data/land/|dat/land/
EOF

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
  echo "${DATADIR}/topo/${TOPO_FORMAT}/Products/|dat/topo/${TOPO_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products/|dat/landuse/${LANDUSE_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

#-------------------------------------------------------------------------------
# observations

if [ "$JOBTYPE" = 'cycle' ] && ((OBS_USE_JITDT != 1)); then
  time=$(datetime $STIME $LCYCLE s)
  while ((time <= $(datetime $ETIME $LCYCLE s))); do
    for iobs in $(seq $OBSNUM); do
      if [ "${OBSNAME[$iobs]}" != '' ]; then
        if [ "${OBS_FORMAT[$iobs]}" = 'PAWR_TOSHIBA' ]; then
          if [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.10000000.dat ] && \
             [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.20000000.dat ] && \
             [ -e ${OBS}/${OBSNAME[$iobs]}_${time}_pawr_qcf.dat ]; then
            if ((DACYCLE == 1)); then
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.10000000.dat|obs.${OBSNAME[$iobs]}.10000000_${time:0:8}-${time:8:6}.000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.20000000.dat|obs.${OBSNAME[$iobs]}.20000000_${time:0:8}-${time:8:6}.000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
              echo "${OBS}/${OBSNAME[$iobs]}_${time}_pawr_qcf.dat|obs.${OBSNAME[$iobs]}_pawr_qcf_${time:0:8}-${time:8:6}.000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
            else
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.10000000.dat|obs.${OBSNAME[$iobs]}_${time}.10000000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.20000000.dat|obs.${OBSNAME[$iobs]}_${time}.20000000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
              echo "${OBS}/${OBSNAME[$iobs]}_${time}_pawr_qcf.dat|obs.${OBSNAME[$iobs]}_${time}_pawr_qcf.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
            fi
          fi
        else
          if [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
            if ((DACYCLE == 1)); then
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs.${OBSNAME[$iobs]}_${time:0:8}-${time:8:6}.000.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
            else
              echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs.${OBSNAME[$iobs]}_${time}.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
            fi
          fi
        fi
      fi
    done
    time=$(datetime $time $LCYCLE s)
  done
fi

#-------------------------------------------------------------------------------
# create empty directories

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
|mean/
|log/
EOF

if [ "$JOBTYPE" = 'cycle' ]; then
  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
|sprd/
EOF

  if ((MAKEINIT == 1)); then
    echo "|mdet/" >> ${STAGING_DIR}/${STGINLIST}
    for m in $(seq $MEMBER); do
      echo "|${name_m[$m]}/" >> ${STAGING_DIR}/${STGINLIST}
    done
  fi
fi

#-------------------------------------------------------------------------------
# time-invariant outputs

#-------------------
# stage-in
#-------------------

# domain catalogue
#-------------------
if ((BDY_FORMAT == 1)); then
  if [ -s "$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt" ]; then
    pathin="$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt"
    path="latlon_domain_catalogue.bdy.txt"
    echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
  else
    echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
    echo "        '$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt'" >&2
    exit 1
  fi
fi

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
# Generate the launcher configuration files
#
# Usage: config_file_scale_launcher JOBTYPE CONF_NAME CONF_FILES MEMBER_RUN
#
#   JOBTYPE     Job type (cycle/fcst)
#   CONF_NAME   Name of this configuration file
#   CONF_FILES  Name pattern of configuration files for each member/domain
#   MEMBER_RUN  Number of members needed to run
#
# Other input variables:
#   $SCRP_DIR
#   $MEMBER
#   $PPN_APPAR
#   $mem_nodes
#   $DOMNUM
#   $PRC_DOMAINS_LIST
#-------------------------------------------------------------------------------

local JOBTYPE="$1"; shift
local CONF_NAME="$1"; shift
local CONF_FILES="$1"; shift
local MEMBER_RUN="$1"

#-------------------------------------------------------------------------------

local MEMBER_TOT
local ENS_WITH_MEAN_TF
local ENS_WITH_MDET_TF
local MDET_CYCLED_TF

if [ "$JOBTYPE" = 'cycle' ] || [ "$JOBTYPE" = 'letkf' ]; then
  MEMBER_TOT=$MEMBER
  ENS_WITH_MEAN_TF='.true.'
  ENS_WITH_MDET_TF='.false.'
  if ((DET_RUN == 1)); then
    ENS_WITH_MDET_TF='.true.'
  fi
elif [ "$JOBTYPE" = 'fcst' ]; then
  MEMBER_TOT=$((MEMBER+1))
  if ((DET_RUN == 1)); then
    MEMBER_TOT=$((MEMBER+2))
  fi
  ENS_WITH_MEAN_TF='.false.'
  ENS_WITH_MDET_TF='.false.'
fi

if ((DET_RUN_CYCLED == 1)); then
  MDET_CYCLED_TF='.true.'
else
  MDET_CYCLED_TF='.false.'
fi

cat $SCRP_DIR/config.nml.ensmodel | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBER_TOT," \
        -e "/!--ENS_WITH_MEAN--/a ENS_WITH_MEAN = $ENS_WITH_MEAN_TF," \
        -e "/!--ENS_WITH_MDET--/a ENS_WITH_MDET = $ENS_WITH_MDET_TF," \
        -e "/!--MEMBER_RUN--/a MEMBER_RUN = $MEMBER_RUN," \
        -e "/!--MDET_CYCLED--/a MDET_CYCLED = ${MDET_CYCLED_TF}," \
        -e "/!--CONF_FILES--/a CONF_FILES = \"${CONF_FILES}\"," \
        -e "/!--PPN--/a PPN = $PPN_APPAR," \
        -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
        -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
        -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
    > $CONFIG_DIR/${CONF_NAME}

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

cp -fr $CONFIG_DIR/* ${OUTDIR[1]}/config

#-------------------------------------------------------------------------------
}

#===============================================================================
