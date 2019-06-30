#!/bin/bash
#
#===============================================================================
#  Steps of 'cycle.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=5
stepname[1]='Run SCALE pp'
stepexecdir[1]="$TMPRUN/scale_pp"
stepexecname[1]="scale-rm_pp_ens"
stepname[2]='Run SCALE init'
stepexecdir[2]="$TMPRUN/scale_init"
stepexecname[2]="scale-rm_init_ens"
stepname[3]='Run ensemble forecasts'
stepexecdir[3]="$TMPRUN/scale"
stepexecname[3]="scale-rm_ens"
stepname[4]='Run observation operator'
stepexecdir[4]="$TMPRUN/obsope"
stepexecname[4]="obsope"
stepname[5]='Run LETKF'
stepexecdir[5]="$TMPRUN/letkf"
stepexecname[5]="letkf"

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run data assimilation cycles.

Configuration files:
  config.main
  config.cycle

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME ISTEP FSTEP CONF_MODE TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  ISTEP       The initial step in the first cycle from which this script starts
               (default: the first step)
  FSTEP       The final step in the last cycle by which this script ends
               (default: the last step)
  CONF_MODE   Mode of creating runtime configuration files: 'dynamic' or 'static'
               (default: 'dynamic')
  TIME_LIMIT  Requested time limit (only used when using a job scheduler)
               (default: 30 minutes)
"

#if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
#  echo "$USAGE" >&2
#  exit 1
#fi

#-------------------------------------------------------------------------------
# set parameters from command line

STIME=${1:-$STIME}; shift
ETIME=${1:-$ETIME}; shift
ISTEP=${1:-$ISTEP}; shift
FSTEP=${1:-$FSTEP}; shift
CONF_MODE=${1:-$CONF_MODE}; shift
TIME_LIMIT="${1:-$TIME_LIMIT}"

#if [ -z "$STIME" ]; then
#  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
#  echo "$USAGE" >&2
#  exit 1
#fi

#-------------------------------------------------------------------------------
# assign default values to and standardize the parameters

STIME=$(datetime $STIME)
ETIME=$(datetime ${ETIME:-$STIME})
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-$nsteps}
CONF_MODE=${CONF_MODE:-"dynamic"}
TIME_LIMIT=${TIME_LIMIT:-"0:30:00"}

#-------------------------------------------------------------------------------
# error detection

#if ((MACHINE_TYPE == 10 && ONLINE_STGOUT != 0)); then
#  echo "[Error] $myname: When \$MACHINE_TYPE = 10, \$ONLINE_STGOUT needs to be 0." >&2
#  exit 1
#fi

if [ "$CONF_MODE" != 'static' ] && ((DOMNUM > 1)); then
  echo "[Error] Online nesting with multiple domains is only allowed when the static-config mode is used (\$CONF_MODE = 'static')." 1>&2
  exit 1
fi

if ((RUN_LEVEL == 0)); then
  if ((ENABLE_PARAM_USER == 1)) && [ ! -e "$SCRP_DIR/config.nml.scale_user" ]; then
    echo "[Error] $myname: When \$ENABLE_PARAM_USER = 1, 'config.nml.scale_user' file is required." >&2
    exit 1
  fi
  if ((BDY_FORMAT == 4)) && [ ! -e "$SCRP_DIR/config.nml.grads_boundary" ]; then
    echo "[Error] $myname: When \$BDY_FORMAT = 4, 'config.nml.grads_boundary' file is required." >&2
    exit 1
  fi

  if ((MAKEINIT == 1)); then
    if [ -d "${OUTDIR}/${STIME}/anal" ]; then
      if [ -n "$(ls ${OUTDIR}/${STIME}/anal 2> /dev/null)" ]; then
        echo "[Error] $myname: Initial ensemble is to be generated (\$MAKEINIT = 1) at \"${OUTDIR}/${STIME}/anal/\", but existing data are found there;" >&2
        echo "        Set \$MAKEINIT = 0 or remove \"${OUTDIR}/${STIME}/anal/*\" before running this job." >&2
        exit 1
      fi
    fi
  fi
fi

#... more detections...

#-------------------------------------------------------------------------------
# common variables

OUT_CYCLE_SKIP=${OUT_CYCLE_SKIP:-1}

CYCLEFLEN=$WINDOW_E     # Model forecast length in a cycle (second)
if [ -z "$FCSTOUT" ] || ((FCSTOUT >= LTIMESLOT)); then
  CYCLEFOUT=$LTIMESLOT  # Model forecast output interval (second)
elif ((LTIMESLOT % FCSTOUT == 0)); then
  CYCLEFOUT=$FCSTOUT
else
  echo "[Error] If \$FCSTOUT < \$LTIMESLOT, \$LTIMESLOT needs to be an exact multiple of \$FCSTOUT" >&2
  exit 1
fi

if ((BDY_FORMAT >= 1)); then
  if ((BDYCYCLE_INT % BDYINT != 0)); then
    echo "[Error] \$BDYCYCLE_INT needs to be an exact multiple of \$BDYINT" >&2
    exit 1
  fi
  BDY_STARTFRAME_MAX=$((BDYCYCLE_INT / BDYINT))
  if [ -z "$PARENT_REF_TIME" ]; then
    PARENT_REF_TIME=$STIME
    for bdy_startframe in $(seq $BDY_STARTFRAME_MAX); do
      if ((BDY_FORMAT == 1)) && [ -s "$DATA_BDY_SCALE/${PARENT_REF_TIME}/hist/${BDY_MEAN}/history${SCALE_SFX_0}" ]; then
        break
      elif ((BDY_FORMAT == 2 && BDY_ROTATING == 1)) && [ -s "$DATA_BDY_WRF/${PARENT_REF_TIME}/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}" ]; then
        break
      elif ((BDY_FORMAT == 2 && BDY_ROTATING != 1)) && [ -s "$DATA_BDY_WRF/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}" ]; then
        break
      elif ((BDY_FORMAT == 4 && BDY_ROTATING == 1)) && [ -s "$DATA_BDY_GRADS/${PARENT_REF_TIME}/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd" ]; then
        break
      elif ((BDY_FORMAT == 4 && BDY_ROTATING != 1)) && [ -s "$DATA_BDY_GRADS/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd" ]; then
        break
      elif ((bdy_startframe == BDY_STARTFRAME_MAX)); then
        echo "[Error] Cannot find boundary files." >&2
        exit 1
      fi
      PARENT_REF_TIME=$(datetime $PARENT_REF_TIME -${BDYINT} s)
    done
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

print_setting () {
#-------------------------------------------------------------------------------

for vname in DIR DOMAIN @INDIR @OUTDIR @DATA_TOPO DATA_TOPO_BDY_SCALE @DATA_LANDUSE DATA_BDY_SCALE \
             DATA_BDY_SCALE_PREP DATA_BDY_WRF DATA_BDY_NICAM OBS OBSNCEP DET_RUN TOPO_FORMAT \
             LANDUSE_FORMAT LANDUSE_UPDATE BDY_FORMAT BDY_ENS BDYINT BDYCYCLE_INT PARENT_REF_TIME \
             ENABLE_PARAM_USER OCEAN_INPUT OCEAN_FORMAT LAND_INPUT LAND_FORMAT OBSNUM WINDOW_S WINDOW_E \
             LCYCLE LTIMESLOT MEMBER NNODES NNODES_APPAR PPN PPN_APPAR THREADS @SCALE_NP \
             STIME ETIME ISTEP FSTEP CONF_MODE FCSTOUT MAKEINIT OUT_OPT TOPOOUT_OPT \
             LANDUSEOUT_OPT BDYOUT_OPT OBSOUT_OPT LOG_OPT LOG_TYPE; do
  if [ "${vname:0:1}" = '@' ]; then
    for d in $(seq $DOMNUM); do
      vname_d="${vname:1}[$d]"
      printf '  %-20s = %s\n' "$vname_d" "${!vname_d}"
    done
  else
    printf '  %-20s = %s\n' $vname "${!vname}"
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

staging_list () {
#-------------------------------------------------------------------------------
# common variables

if ((PNETCDF == 1)); then
  local mem_np_=1
else
  local mem_np_=$mem_np
fi

#-------------------------------------------------------------------------------
# TMPDAT

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${COMMON_DIR}/pdbash|${DAT_SUBDIR}/exec/pdbash
${SCRP_DIR}/config.nml.scale_pp|${DAT_SUBDIR}/conf/config.nml.scale_pp
${SCRP_DIR}/config.nml.scale_init|${DAT_SUBDIR}/conf/config.nml.scale_init
${SCRP_DIR}/config.nml.scale|${DAT_SUBDIR}/conf/config.nml.scale
${SCRP_DIR}/config.nml.ensmodel|${DAT_SUBDIR}/conf/config.nml.ensmodel
${SCRP_DIR}/config.nml.letkf|${DAT_SUBDIR}/conf/config.nml.letkf
EOF

cat >> ${STAGING_DIR}/${STGINLIST_CONSTDB} << EOF
${SCALEDIR}/scale-rm/test/data/rad/cira.nc|${DAT_SUBDIR}/rad/cira.nc
${SCALEDIR}/scale-rm/test/data/rad/PARAG.29|${DAT_SUBDIR}/rad/PARAG.29
${SCALEDIR}/scale-rm/test/data/rad/PARAPC.29|${DAT_SUBDIR}/rad/PARAPC.29
${SCALEDIR}/scale-rm/test/data/rad/rad_o3_profs.txt|${DAT_SUBDIR}/rad/rad_o3_profs.txt
${SCALEDIR}/scale-rm/test/data/rad/VARDATA.RM29|${DAT_SUBDIR}/rad/VARDATA.RM29
${SCALEDIR}/scale-rm/test/data/rad/MIPAS/|${DAT_SUBDIR}/rad/MIPAS/
${SCALEDIR}/scale-rm/test/data/land/|${DAT_SUBDIR}/land/
EOF

if [ -e "${SCRP_DIR}/config.nml.scale_user" ]; then
  echo "${SCRP_DIR}/config.nml.scale_user|${DAT_SUBDIR}/conf/config.nml.scale_user" >> ${STAGING_DIR}/${STGINLIST}
fi
if [ -e "${SCRP_DIR}/config.nml.obsope" ]; then
  echo "${SCRP_DIR}/config.nml.obsope|${DAT_SUBDIR}/conf/config.nml.obsope" >> ${STAGING_DIR}/${STGINLIST}
fi
if [ -e "${SCRP_DIR}/config.nml.grads_boundary" ]; then
  echo "${SCRP_DIR}/config.nml.grads_boundary|${DAT_SUBDIR}/conf/config.nml.grads_boundary" >> ${STAGING_DIR}/${STGINLIST}
fi

if [ "$TOPO_FORMAT" != 'prep' ]; then
  echo "${DATADIR}/topo/${TOPO_FORMAT}/Products/|${DAT_SUBDIR}/topo/${TOPO_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products/|${DAT_SUBDIR}/landuse/${LANDUSE_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

time=$(datetime $STIME $LCYCLE s)
while ((time <= $(datetime $ETIME $LCYCLE s))); do
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ] && [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
      echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|${DAT_SUBDIR}/obs/${OBSNAME[$iobs]}_${time}.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
    fi
  done
  time=$(datetime $time $LCYCLE s)
done

if [ "$PRESET" = 'K' ] || [ "$PRESET" = 'K_rankdir' ]; then
  echo "${COMMON_DIR}/datetime|${DAT_SUBDIR}/exec/datetime" >> ${STAGING_DIR}/${STGINLIST}
fi

#-------------------------------------------------------------------------------
# TMPRUN

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${ENSMODEL_DIR}/scale-rm_pp_ens|${RUN_SUBDIR}/scale_pp/scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|${RUN_SUBDIR}/scale_init/scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|${RUN_SUBDIR}/scale/scale-rm_ens
${OBSUTIL_DIR}/obsope|${RUN_SUBDIR}/obsope/obsope
${LETKF_DIR}/letkf|${RUN_SUBDIR}/letkf/letkf
EOF

# H08
#-------------------
if [ -e "${RTTOV_COEF}" ] && [ -e "${RTTOV_SCCOEF}" ]; then
  cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${RTTOV_COEF}|${RUN_SUBDIR}/obsope/rtcoef_himawari_8_ahi.dat
${RTTOV_COEF}|${RUN_SUBDIR}/letkf/rtcoef_himawari_8_ahi.dat
${RTTOV_SCCOEF}|${RUN_SUBDIR}/obsope/sccldcoef_himawari_8_ahi.dat
${RTTOV_SCCOEF}|${RUN_SUBDIR}/letkf/sccldcoef_himawari_8_ahi.dat
EOF
fi

#-------------------------------------------------------------------------------
# TMPOUT

# empty directories
#-------------------

echo "|${OUT_SUBDIR}/const/log/" >> ${STAGING_DIR}/${STGINLIST}
if ((PNETCDF != 1)); then
  echo "|${OUT_SUBDIR}/const/topo/" >> ${STAGING_DIR}/${STGINLIST}
  if ((LANDUSE_UPDATE != 1)); then
    echo "|${OUT_SUBDIR}/const/landuse/" >> ${STAGING_DIR}/${STGINLIST}
  fi
fi

#-------------------

time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))

  # empty directories
  #-------------------

  if ((PNETCDF != 1)); then
    if ((LANDUSE_UPDATE == 1)); then
      echo "|${OUT_SUBDIR}/${time}/landuse/" >> ${STAGING_DIR}/${STGINLIST}
    fi
  fi

  echo "|${OUT_SUBDIR}/${time}/log/scale_pp/" >> ${STAGING_DIR}/${STGINLIST}
  echo "|${OUT_SUBDIR}/${time}/log/scale_init/" >> ${STAGING_DIR}/${STGINLIST}
  echo "|${OUT_SUBDIR}/${time}/log/scale/" >> ${STAGING_DIR}/${STGINLIST}
  echo "|${OUT_SUBDIR}/${atime}/log/obsope/" >> ${STAGING_DIR}/${STGINLIST}
  echo "|${OUT_SUBDIR}/${atime}/log/letkf/" >> ${STAGING_DIR}/${STGINLIST}
  echo "|${OUT_SUBDIR}/${atime}/obs/" >> ${STAGING_DIR}/${STGINLIST}

  #-------------------
  # stage-in
  #-------------------

  # anal
  #-------------------
  if ((loop == 1 && MAKEINIT != 1)); then
    for m in $(seq $mtot); do
      for q in $(seq $mem_np_); do
        path="${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
        echo "${INDIR}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  # anal_ocean
  #-------------------
  if ((OCEAN_INPUT == 1)) && ((OCEAN_FORMAT == 0)); then
    for m in $(seq $mtot); do
      for q in $(seq $mem_np_); do
        path="${time}/anal/${name_m[$m]}${CONNECTOR}init_ocean$(scale_filename_sfx $((q-1)))"
        echo "${INDIR}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  # anal_land
  #-------------------
  if ((LAND_INPUT == 1)) && ((LAND_FORMAT == 0)); then
    for m in $(seq $mtot); do
      for q in $(seq $mem_np_); do
        path="${time}/anal/${name_m[$m]}${CONNECTOR}init_land$(scale_filename_sfx $((q-1)))"
        echo "${INDIR}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  # topo
  #-------------------
  if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
    if ((DISK_MODE == 3)); then
      for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
        for q in $(seq $mem_np_); do
          path="const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
          echo "${DATA_TOPO}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    else
      for q in $(seq $mem_np_); do
        path="const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
        echo "${DATA_TOPO}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}
      done
    fi
  fi

  # topo (bdy_scale)
  #-------------------
  if ((loop == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
    if ((PNETCDF_BDY_SCALE == 1)); then
      pathin="${DATA_TOPO_BDY_SCALE}.nc"
      path="bdytopo/const/topo.nc"
    else
      pathin="${DATA_TOPO_BDY_SCALE}/"
      path="bdytopo/const/"
    fi
    echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
  fi

  # landuse
  #-------------------
  if ((loop == 1 || LANDUSE_UPDATE == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
    if ((DISK_MODE == 3)); then
      for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
        for q in $(seq $mem_np_); do
          if ((LANDUSE_UPDATE == 1)); then
            path="${time}/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
          else
            path="const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
          fi
          echo "${DATA_LANDUSE}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    else
      for q in $(seq $mem_np_); do
        if ((LANDUSE_UPDATE == 1)); then
          path="${time}/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
        else
          path="const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
        fi
        echo "${DATA_LANDUSE}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}
      done
    fi
  fi

  # bdy (prepared)
  #-------------------
  if ((BDY_FORMAT == 0)); then
    if ((BDY_ENS == 0)); then
      if ((DISK_MODE == 3)); then
        for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
          for q in $(seq $mem_np_); do
            pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
            path="${time}/bdy/mean${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
            if ((USE_INIT_FROM_BDY == 1)); then
              pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
              path="${time}/bdy/mean${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
              echo "${pathin}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
            fi
          done
        done
      else
        for q in $(seq $mem_np_); do
          pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          path="${time}/bdy/mean${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}
          if ((USE_INIT_FROM_BDY == 1)); then
            pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
            path="${time}/bdy/mean${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}
          fi
        done
      fi
    elif ((BDY_ENS == 1)); then
      for m in $(seq $mtot); do
        for q in $(seq $mem_np_); do
          path="${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          echo "${DATA_BDY_SCALE_PREP}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
          if ((USE_INIT_FROM_BDY == 1)); then
            path="${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
            echo "${DATA_BDY_SCALE_PREP}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
          fi
        done
      done
    fi
  fi

  # additive inflation
  #-------------------
  if ((loop == 1 && ADDINFL == 1)); then
    for m in $(seq $MEMBER); do
      for q in $(seq $mem_np_); do
        path="const/addi/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
        echo "${DATA_ADDINFL}/${path}|${OUT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  #-------------------
  # stage-out
  #-------------------

  # anal (initial time)
  #-------------------
  if ((loop == 1 && MAKEINIT == 1)); then
    path="${time}/anal/"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
  fi

  # topo
  #-------------------
  if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
    if ((PNETCDF == 1)); then
      path="const/topo.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      path="const/topo/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi

  # landuse
  #-------------------
  if ((loop == 1 || LANDUSE_UPDATE == 1)) && ((LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
    if ((PNETCDF == 1)); then
      if ((LANDUSE_UPDATE == 1)); then
        path="${time}/landuse.nc"
      else
        path="const/landuse.nc"
      fi
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      if ((LANDUSE_UPDATE == 1)); then
        path="${time}/landuse/"
      else
        path="const/landuse/"
      fi
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi

  # bdy
  #-------------------
  if ((BDY_FORMAT != 0)); then
    if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
      path="${time}/bdy/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    elif ((BDYOUT_OPT <= 2)); then
      if ((PNETCDF == 1)); then
        path="${time}/bdy/mean.boundary.nc"
#        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
        if ((USE_INIT_FROM_BDY == 1)); then
          path="${time}/bdy/mean.init_bdy.nc"
#          echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
          echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
        fi
      else
        path="${time}/bdy/mean/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        if ((DET_RUN == 1)); then
          path="${time}/bdy/mdet/"
          echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        fi
      fi
    fi
  fi

  # hist
  #-------------------
  if ((OUT_OPT <= 1)); then
    path="${time}/hist/"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
  elif ((OUT_OPT <= 2)); then
    if ((PNETCDF == 1)); then
      path="${time}/hist/mean.history.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      path="${time}/hist/mean/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      if ((DET_RUN == 1)); then
        path="${time}/hist/mdet/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      fi
    fi
  fi

  # gues
  #-------------------
  if ((OUT_OPT <= 3)); then
    path="${atime}/gues/"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
  elif ((OUT_OPT <= 6)); then
    if ((PNETCDF == 1)); then
      path="${atime}/gues/mean.init.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      if ((SPRD_OUT == 1)); then
        path="${atime}/gues/sprd.init.nc"
#        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      fi
      if ((DET_RUN == 1)); then
        path="${atime}/gues/mdet.init.nc"
#        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      fi
    else
      path="${atime}/gues/mean/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      if ((SPRD_OUT == 1)); then
        path="${atime}/gues/sprd/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      fi
      if ((DET_RUN == 1)); then
        path="${atime}/gues/mdet/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      fi
    fi
  fi

  # anal
  #-------------------
  if ((OUT_OPT <= 4 || (OUT_OPT <= 5 && loop % OUT_CYCLE_SKIP == 0) || atime > ETIME)); then
    path="${atime}/anal/"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
  elif ((OUT_OPT <= 7)); then
    if ((PNETCDF == 1)); then
      path="${atime}/anal/mean.init.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      if ((SPRD_OUT == 1)); then
        path="${atime}/anal/sprd.init.nc"
#        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      fi
      if ((DET_RUN == 1)); then
        path="${atime}/anal/mdet.init.nc"
#        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
      fi
    else
      path="${atime}/anal/mean/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      if ((SPRD_OUT == 1)); then
        path="${atime}/anal/sprd/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      fi
      if ((DET_RUN == 1)); then
        path="${atime}/anal/mdet/"
        echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      fi
    fi
  fi

  ### anal_ocean [mean]

  # diag
  #-------------------
  if ((ADAPTINFL == 1)); then
    if ((PNETCDF == 1)); then
      path="${atime}/diag/infl.init.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      path="${atime}/diag/infl/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi
  if ((RTPS_INFL_OUT == 1)); then
    if ((PNETCDF == 1)); then
      path="${atime}/diag/rtps.init.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      path="${atime}/diag/rtps/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi
  if ((NOBS_OUT == 1)); then
    if ((PNETCDF == 1)); then
      path="${atime}/diag/nobs.init.nc"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
    else
      path="${atime}/diag/nobs/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi

  # obs
  #-------------------
  if ((OBSOUT_OPT <= 3)); then
    path="${atime}/obs/"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
  fi

  # log
  #-------------------
  if [ "$MPI_TYPE" = 'K' ]; then
    log_zeros='0'
  else
    log_zeros="$PROCESS_FMT_0"
  fi

  if ((loop == 1 && LOG_OPT <= 3)); then
    path="const/log/latlon_domain_catalogue.txt"
    echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
  fi

  if ((LOG_OPT <= 2)); then
    if ((LOG_TYPE == 1)); then
#      path="${time}/log/scale_pp/0001_pp.conf"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${time}/log/scale_pp/0001_LOG${SCALE_SFX_NONC_0}"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${time}/log/scale_pp/NOUT.${log_zeros}"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${time}/log/scale_init/0001_init.conf"
      path="${time}/log/scale_init/mean_init.conf"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${time}/log/scale_init/0001_gradsbdy.conf"
      path="${time}/log/scale_init/mean_gradsbdy.conf"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
      path="${time}/log/scale_init/0001_LOG${SCALE_SFX_NONC_0}"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
      if ((BDY_ENS == 1)); then
#          path="${time}/log/scale_init/NOUT-1.${log_zeros}"
	  path="${time}/log/scale_init/NOUT-1-${log_zeros}"
      else
        path="${time}/log/scale_init/NOUT.${log_zeros}"
      fi
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
    else
      path="${time}/log/scale_pp/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      path="${time}/log/scale_init/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi
  if ((LOG_OPT <= 3)); then
    if ((LOG_TYPE == 1)); then
      path="${time}/log/scale/0001_run.conf"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
      path="${time}/log/scale/0001_LOG${SCALE_SFX_NONC_0}"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${time}/log/scale/NOUT-1.${log_zeros}"
      path="${time}/log/scale/NOUT-1-${log_zeros}"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
    else
      path="${time}/log/scale/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi
  if ((LOG_OPT <= 4)); then
    if ((LOG_TYPE == 1)); then
#      path="${atime}/log/obsope/obsope.conf"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${atime}/log/obsope/NOUT.${log_zeros}"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
      path="${atime}/log/letkf/letkf.conf"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#      path="${atime}/log/letkf/NOUT.${log_zeros}"
      path="${atime}/log/letkf/NOUT-${log_zeros}"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
    else
#      path="${atime}/log/obsope/"
#      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
      path="${atime}/log/letkf/"
      echo "${OUTDIR}/${path}|${OUT_SUBDIR}/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
    fi
  fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

#-------------------
# stage-in
#-------------------

# bdy
#-------------------
if ((BDY_FORMAT >= 1)); then
  if ((BDY_FORMAT == 1)); then
    if [ -s "$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt" ]; then
      pathin="$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt"
      path="bdyorg/latlon_domain_catalogue.txt"
      echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
    else
      echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
      echo "        '$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt'" >&2
      exit 1
    fi
  fi

  nbdy_all=0
  time=$STIME
  while ((time <= ETIME)); do
    bdy_setting $time $CYCLEFLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

    for ibdy in $(seq $nbdy); do
      time_bdy=${bdy_times[$ibdy]}

      bdy_processed=0
      for ibdy2 in $(seq $nbdy_all); do
        if ((${bdy_times_all[$ibdy2]} == $time_bdy)); then
          bdy_processed=1
          break
        fi
      done

      if ((bdy_processed == 0)); then
        nbdy_all=$((nbdy_all+1))
        bdy_times_all[${nbdy_all}]=$time_bdy
      fi

      if ((bdy_processed == 0 || BDY_ROTATING == 1)); then
        if ((BDY_FORMAT == 1)); then

          if ((BDY_ENS == 1)); then
            for m in $(seq $mtot); do
              mem=${name_m[$m]}
              if [ "$mem" = 'mean' ]; then
                mem="$BDY_MEAN"
              fi
              if ((PNETCDF_BDY_SCALE == 1)); then
                pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${mem}.history.nc"
                if ((BDY_ROTATING == 1)); then
                  path="bdyorg/${time_bdy}/${time_bdy}/${name_m[$m]}.history.nc"
                else
                  path="bdyorg/const/${time_bdy}/${name_m[$m]}.history.nc"
                fi
              else
                pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${mem}/"
                if ((BDY_ROTATING == 1)); then
                  path="bdyorg/${time_bdy}/${time_bdy}/${name_m[$m]}/"
                else
                  path="bdyorg/const/${time_bdy}/${name_m[$m]}/"
                fi
              fi
              echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            done
          else
            if ((PNETCDF_BDY_SCALE == 1)); then
              pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}.history.nc"
              if ((BDY_ROTATING == 1)); then
                path="bdyorg/${time_bdy}/${time_bdy}/mean.history.nc"
              else
                path="bdyorg/const/${time_bdy}/mean.history.nc"
              fi
            else             
              pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}/"
              if ((BDY_ROTATING == 1)); then
                path="bdyorg/${time_bdy}/mean/${time_bdy}/"
              else
                path="bdyorg/const/mean/${time_bdy}/"
              fi
            fi
            echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
          fi

        elif ((BDY_FORMAT == 2 || BDY_FORMAT == 4)); then

          if ((BDY_FORMAT == 2)); then
            data_bdy_i=$DATA_BDY_WRF
            filenum=1
            filename_prefix[1]='wrfout_'
            filename_suffix[1]=''
          elif ((BDY_FORMAT == 4)); then
            data_bdy_i=$DATA_BDY_GRADS
            filenum=3
            filename_prefix[1]='atm_'
            filename_suffix[1]='.grd'
            filename_prefix[2]='sfc_'
            filename_suffix[2]='.grd'
            filename_prefix[3]='land_'
            filename_suffix[3]='.grd'
          fi

          if ((BDY_ENS == 1)); then
            for m in $(seq $mtot); do
              for ifile in $(seq $filenum); do
                if ((BDY_ROTATING == 1)); then
                  pathin="$data_bdy_i/${time}/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                  path="bdyorg/${time}/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                else
                  pathin="$data_bdy_i/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                  path="bdyorg/const/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                fi
                echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
              done
            done
          else
            for ifile in $(seq $filenum); do
              if ((BDY_ROTATING == 1)); then
                pathin="$data_bdy_i/${time}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                path="bdyorg/${time}/mean/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              else
                pathin="$data_bdy_i/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                path="bdyorg/const/mean/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              fi
              echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            done
          fi

        fi
      fi # ((bdy_processed == 0 || BDY_ROTATING == 1))
    done

    time=$(datetime $time $LCYCLE s)
  done
fi # ((BDY_FORMAT >= 1))

#-------------------------------------------------------------------------------
}

#===============================================================================

enspp_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[1]}: Pre-processing script start" >&2
fi

if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
  echo "  ... skip this step (use prepared topo and landuse files)"
  exit 1
elif ((BDY_FORMAT == 0)); then
  echo "  ... skip this step (use prepared boundaries)"
  exit 1
elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
  echo "  ... skip this step (already done in the first cycle)"
  exit 1
fi

if ((BDY_FORMAT == 1)); then
  bdycatalogue=${TMPDAT_BDYDATA}/bdyorg/latlon_domain_catalogue.txt
  bdytopo=${TMPDAT_BDYDATA}/bdytopo/const/topo
fi

if ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= mtot ? repeat_mems : mtot))
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_pp_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale_pp $MEMBER_RUN $iter cycle
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[1]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= MEMBER_RUN)); then
      bash $SCRP_DIR/src/pre_scale_pp.sh $MYRANK $time ${name_m[$m]} \
           $TMPRUN/scale_pp/${name_m[$m]} $TMPDAT \
           cycle ${bdytopo} ${bdycatalogue}
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

enspp_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
  return 1
elif ((BDY_FORMAT == 0)); then
  return 1
fi

if ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= mtot ? repeat_mems : mtot))
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= MEMBER_RUN)); then
      bash $SCRP_DIR/src/post_scale_pp.sh $MYRANK $time \
           ${name_m[$m]} $TMPRUN/scale_pp/${name_m[$m]} $LOG_OPT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensinit_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[2]}: Pre-processing script start" >&2
fi

if ((BDY_FORMAT == 0)); then
  echo "  ... skip this step (use prepared boundaries)"
  exit 1
fi

bdy_setting $time $CYCLEFLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"
bdy_time_list=''
for ibdy in $(seq $nbdy); do
  bdy_time_list="${bdy_time_list}${bdy_times[$ibdy]} "
done

bdyorgf=${TMPDAT_BDYDATA}/bdyorg

if ((BDY_ENS == 1)); then
  MEMBER_RUN=$mtot
elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= mtot ? repeat_mems : mtot))
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

if ((LANDUSE_UPDATE == 1)); then
  time_l=${time}
else
  time_l='const'
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_init_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale_init $MEMBER_RUN $iter cycle
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[2]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= MEMBER_RUN)); then
      if ((BDY_ENS == 1)); then
        mem_bdy=${name_m[$m]}
      else
        mem_bdy='mean'
      fi

      bash $SCRP_DIR/src/pre_scale_init.sh $MYRANK \
           $TMPOUT/const/${CONNECTOR_TOPO}topo $TMPOUT/${time_l}/${CONNECTOR_LANDUSE}landuse \
           ${bdyorgf} $time $mkinit ${name_m[$m]} $mem_bdy \
           $TMPRUN/scale_init/${name_m[$m]} \
           "$bdy_time_list" $ntsteps $ntsteps_skip cycle
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensinit_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

if ((BDY_FORMAT == 0)); then
  return 1
fi

if ((BDY_ENS == 1)); then
  MEMBER_RUN=$mtot
elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= mtot ? repeat_mems : mtot))
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= MEMBER_RUN)); then
      if ((BDY_ENS == 1)); then
        mem_bdy=${name_m[$m]}
      else
        mem_bdy='mean'
      fi

      bash $SCRP_DIR/src/post_scale_init.sh $MYRANK $time \
           $mkinit $mem_bdy $TMPRUN/scale_init/${name_m[$m]} $LOG_OPT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Post-processing script (member) start" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[3]}: Pre-processing script start" >&2
fi

bdy_setting $time $CYCLEFLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

############
#if ((BDY_FORMAT == 1)); then
#  if ((DATA_BDY_TMPLOC == 1)); then
#    bdyorgf=$TMPDAT/bdyorg
#  elif ((DATA_BDY_TMPLOC == 2)); then
#    bdyorgf=$TMPOUT/bdyorg
#  fi
#  time_bdy=$(datetime $time $BDYCYCLE_INT s)
#  for bdy_startframe in $(seq $BDY_STARTFRAME_MAX); do
#    if [ -s "$bdyorgf/${time_bdy}/mean/history${SCALE_SFX_0}" ]; then
#      break
#    elif ((bdy_startframe == BDY_STARTFRAME_MAX)); then
#      echo "[Error] Cannot find boundary files from the SCALE history files." >&2
#      exit 1
#    fi
#    time_bdy=$(datetime $time_bdy -${BDYINT} s)
#  done
#  time_bdy=$(datetime $time_bdy -$BDYCYCLE_INT s)
#fi
############

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale $mtot $iter cycle
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

if ((LANDUSE_UPDATE == 1)); then
  time_l=${time}
else
  time_l='const'
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[3]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
#      if ((PERTURB_BDY == 1)); then
#        ...
#      fi

      if ((BDY_ENS == 1)); then
        mem_bdy=${name_m[$m]}
      else
        mem_bdy='mean'
      fi

      ocean_base='-'
      if ((OCEAN_INPUT == 1)); then
        if ((OCEAN_FORMAT == 0)); then
          ocean_base="$TMPOUT/${time}/anal/${mem_bdy}${CONNECTOR}init_ocean"
        elif ((OCEAN_FORMAT == 99 && mkinit != 1)); then
          ocean_base="$TMPOUT/${time}/bdy/${mem_bdy}${CONNECTOR}init_bdy"
        fi
      fi

      land_base='-'
      if ((LAND_INPUT == 1)); then
        if ((LAND_FORMAT == 0)); then
          land_base="$TMPOUT/${time}/anal/${mem_bdy}${CONNECTOR}init_land"
        elif ((LAND_FORMAT == 99 && mkinit != 1)); then
          land_base="$TMPOUT/${time}/bdy/${mem_bdy}${CONNECTOR}init_bdy"
        fi
      fi

      bdy_base="$TMPOUT/${time}/bdy/${mem_bdy}${CONNECTOR}boundary"

      bash $SCRP_DIR/src/pre_scale.sh $MYRANK ${name_m[$m]} \
           $TMPOUT/${time}/anal/${name_m[$m]}${CONNECTOR}init $ocean_base $land_base $bdy_base \
           $TMPOUT/const/${CONNECTOR_TOPO}topo $TMPOUT/${time_l}/${CONNECTOR_LANDUSE}landuse \
           $time $CYCLEFLEN $LCYCLE $CYCLEFOUT $TMPRUN/scale/${name_m[$m]} $OUT_OPT \
           cycle $bdy_start_time $SPRD_OUT $RTPS_INFL_OUT $NOBS_OUT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

DELETE_MEMBER=0
if ((OUT_OPT >= 5 && ((loop - 1) % OUT_CYCLE_SKIP != 0))); then
  DELETE_MEMBER=1
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
#      if ((PERTURB_BDY == 1)); then
#        ...
#      fi

      bash $SCRP_DIR/src/post_scale.sh $MYRANK $time \
           ${name_m[$m]} $CYCLEFLEN $TMPRUN/scale/${name_m[$m]} $LOG_OPT $OUT_OPT \
           cycle $DELETE_MEMBER $SPRD_OUT $RTPS_INFL_OUT $NOBS_OUT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

obsope_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[4]}: Pre-processing script start" >&2
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_obsope_node.sh $MYRANK \
       $time $atime $TMPRUN/obsope ${TMPDAT_OBS}/obs \
       $mem_nodes $mem_np $slot_s $slot_e $slot_b $MEMBER
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[4]}: Pre-processing script end" >&2
fi

for it in $(seq $nitmax); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[4]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
      bash $SCRP_DIR/src/pre_obsope.sh $MYRANK \
           $atime ${name_m[$m]}
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[4]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

obsope_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

for it in $(seq $nitmax); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[4]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
      bash $SCRP_DIR/src/post_obsope.sh $MYRANK \
           ${time} ${atime} ${name_m[$m]} $TMPRUN/obsope $LOG_OPT $OUT_OPT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[4]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

letkf_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((IO_ARB == 1)); then     ##
  if ((MYRANK == 0)); then   ##
    echo "[$(datetime_now)] ${time}: ${stepname[5]}: Wait for 360 seconds" >&2 ##
  fi                         ##
  sleep 360s                 ##
fi                           ##

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[5]}: Pre-processing script start" >&2
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_letkf_node.sh $MYRANK \
       $time $atime $TMPRUN/letkf ${TMPDAT_OBS}/obs \
       $mem_nodes $mem_np $slot_s $slot_e $slot_b $TMPOUT/const/${CONNECTOR_TOPO}topo $OBSOUT_OPT \
       $ADAPTINFL $SPRD_OUT $RTPS_INFL_OUT $NOBS_OUT \
       $MEMBER
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[5]}: Pre-processing script end" >&2
fi

for it in $(seq $nitmax); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[5]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
      bash $SCRP_DIR/src/pre_letkf.sh $MYRANK \
           $atime ${name_m[$m]} $OUT_OPT \
           $ADAPTINFL $SPRD_OUT $RTPS_INFL_OUT $NOBS_OUT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[5]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

letkf_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

for it in $(seq $nitmax); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[5]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= mtot)); then
      bash $SCRP_DIR/src/post_letkf.sh $MYRANK \
           ${atime} $TMPRUN/letkf $LOG_OPT
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[5]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

obstime () {
#-------------------------------------------------------------------------------
# Determine the observation time slots
#  *Require source 'func_datetime.sh'
#
# Usage: obstime TIME
#
#   TIME  Forecast start time
#
# Other input variables:
#   $LTIMESLOT
#   $WINDOW_S
#   $WINDOW_E
#   $LCYCLE
#
# Return variables:
#   $slot_s
#   $slot_e
#   $slot_b
#   $time_sl[1...$slot_e]
#   $timefmt_sl[1...$slot_e]
#-------------------------------------------------------------------------------

if (($# < 1)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

local TIME="$1"

#-------------------------------------------------------------------------------

local otime=$(datetime $TIME)               # HISTORY_OUTPUT_STEP0 = .true.,
#local otime=$(datetime $TIME $LTIMESLOT s)  # HISTORY_OUTPUT_STEP0 = .false.,
local otime_s=$(datetime $TIME $WINDOW_S s)
local otime_e=$(datetime $TIME $WINDOW_E s)
local otime_a=$(datetime $TIME $LCYCLE s)
local is=0
slot_s=0
while ((otime <= otime_e)); do
  is=$((is+1))
  time_sl[$is]=$otime
  timefmt_sl[$is]="$(datetime_fmt ${otime})"
  if ((slot_s == 0 && otime >= otime_s)); then
    slot_s=$is
  fi
  if ((otime == otime_a)); then
    slot_b=$is
  fi
otime=$(datetime $otime $LTIMESLOT s)
done
slot_e=$is

#-------------------------------------------------------------------------------
}

#===============================================================================

archive_log () {
#-------------------------------------------------------------------------------

if ((LOG_TYPE >= 3)); then
  time=$STIME
  atime=$(datetime $time $LCYCLE s)
  while ((time <= ETIME)); do
    for d in $(seq $DOMNUM); do
      if ((LOG_OPT <= 2)) && [ -d "${OUTDIR[$d]}/${time}/log/scale_pp" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_pp.tar scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_pp.tar.gz scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_pp.tar scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_pp.tar.gz scale_pp && rm -fr ${OUTDIR[$d]}/${time}/log/scale_pp
          fi
        fi
      fi

      if ((LOG_OPT <= 2)) && [ -d "${OUTDIR[$d]}/${time}/log/scale_init" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_init.tar scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_init.tar.gz scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale_init.tar scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale_init.tar.gz scale_init && rm -fr ${OUTDIR[$d]}/${time}/log/scale_init
          fi
        fi
      fi

      if ((LOG_OPT <= 3)) && [ -d "${OUTDIR[$d]}/${time}/log/scale" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale.tar scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale.tar.gz scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${time}/log -cf ${OUTDIR[$d]}/${time}/log/scale.tar scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${time}/log -czf ${OUTDIR[$d]}/${time}/log/scale.tar.gz scale && rm -fr ${OUTDIR[$d]}/${time}/log/scale
          fi
        fi
      fi

      if ((LOG_OPT <= 4)) && [ -d "${OUTDIR[$d]}/${atime}/log/obsope" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/obsope.tar obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/obsope.tar.gz obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/obsope.tar obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/obsope.tar.gz obsope && rm -fr ${OUTDIR[$d]}/${atime}/log/obsope
          fi
        fi
      fi

      if ((LOG_OPT <= 4)) && [ -d "${OUTDIR[$d]}/${atime}/log/letkf" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/letkf.tar letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/letkf.tar.gz letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -cf ${OUTDIR[$d]}/${atime}/log/letkf.tar letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf
          elif ((LOG_TYPE == 4)); then
            tar -C ${OUTDIR[$d]}/${atime}/log -czf ${OUTDIR[$d]}/${atime}/log/letkf.tar.gz letkf && rm -fr ${OUTDIR[$d]}/${atime}/log/letkf
          fi
        fi
      fi
    done # [ d in $(seq $DOMNUM) ]

    time=$(datetime $time $LCYCLE s)
    atime=$(datetime $time $LCYCLE s)
  done
  if ((TAR_THREAD > 1)); then
    wait
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
