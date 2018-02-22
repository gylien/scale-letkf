#!/bin/bash
#===============================================================================
#
#  Steps of 'cycle.sh'
#
#===============================================================================

staging_list_static () {
#-------------------------------------------------------------------------------
# Prepare all the staging list files
#
# Usage: staging_list_static
#
# Other input variables:
#   $STAGING_DIR
#   ...
#-------------------------------------------------------------------------------
# common variables

if ((PNETCDF == 1)); then
  local mem_np_=1
else
  local mem_np_=$mem_np
fi

#-------------------------------------------------------------------------------
# executable files

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${ENSMODEL_DIR}/scale-rm_pp_ens|scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|scale-rm_ens
${OBSUTIL_DIR}/obsope|obsope
${LETKF_DIR}/letkf|letkf
${COMMON_DIR}/pdbash|pdbash
${COMMON_DIR}/datetime|datetime
EOF

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
#  if [ -e "${RTTOV_COEF}" ] && [ -e "${RTTOV_SCCOEF}" ]; then
#    cat >> ${STAGING_DIR}/${STGINLIST_CONSTDB} << EOF
#${RTTOV_COEF}|dat/rttov/rtcoef_himawari_8_ahi.dat
#${RTTOV_SCCOEF}|dat/rttov/sccldcoef_himawari_8_ahi.dat
#EOF
#  fi

if [ "$TOPO_FORMAT" != 'prep' ]; then
  echo "${DATADIR}/topo/${TOPO_FORMAT}/Products/|dat/topo/${TOPO_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products/|dat/landuse/${LANDUSE_FORMAT}/Products/" >> ${STAGING_DIR}/${STGINLIST_CONSTDB}
fi

#-------------------------------------------------------------------------------
# observations

time=$(datetime $STIME $LCYCLE s)
while ((time <= $(datetime $ETIME $LCYCLE s))); do
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ] && [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
      echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs.${OBSNAME[$iobs]}_${time}.dat" >> ${STAGING_DIR}/${STGINLIST_OBS}
    fi
  done
  time=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
# create empty directories

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
|sprd/
|log/
EOF

#-------------------------------------------------------------------------------
# time-invariant outputs

#-------------------
# stage-out
#-------------------

# domain catalogue
#-------------------
if ((LOG_OPT <= 3)); then
  path="latlon_domain_catalogue.txt"
  pathout="${OUTDIR}/const/log/latlon_domain_catalogue.txt"
  echo "${pathout}|${path}|1" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[1]}
fi

#-------------------------------------------------------------------------------
# time-variant outputs

time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))

  #-------------------
  # stage-in
  #-------------------

  # anal
  #-------------------
  if ((loop == 1 && MAKEINIT != 1)); then
    for m in $(seq $mtot); do
      for q in $(seq $mem_np_); do
        pathin="${INDIR}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
        path="${name_m[$m]}/anal.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  # topo
  #-------------------
  if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
    if ((DISK_MODE == 3)); then
      for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
        for q in $(seq $mem_np_); do
          pathin="${DATA_TOPO}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
          path="topo.d01$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    else
      for q in $(seq $mem_np_); do
        pathin="${DATA_TOPO}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
        path="topo.d01$(scale_filename_sfx $((q-1)))"
        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
      done
    fi
  fi

#    # topo (bdy_scale)
#    #-------------------
#    if ((loop == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#      pathin="${DATA_TOPO_BDY_SCALE}.nc"
#      path="bdytopo.nc"
#      echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#    fi

  # landuse
  #-------------------
  if ((loop == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
    if ((DISK_MODE == 3)); then
      for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
        for q in $(seq $mem_np_); do
          pathin="${DATA_LANDUSE}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
          path="landuse.d01$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
      done
    else
      for q in $(seq $mem_np_); do
        pathin="${DATA_LANDUSE}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
        path="landuse.d01$(scale_filename_sfx $((q-1)))"
        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
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
            path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#            if ((USE_INIT_FROM_BDY == 1)); then
#              pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#              path="mean/analbdy.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#            fi
          done
        done
      else
        for q in $(seq $mem_np_); do
          pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
#          if ((USE_INIT_FROM_BDY == 1)); then
#            pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#            path="mean/analbdy.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
#          fi
        done
      fi
    elif ((BDY_ENS == 1)); then
      for m in $(seq $mtot); do
        for q in $(seq $mem_np_); do
          pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          path="${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#          if ((USE_INIT_FROM_BDY == 1)); then
#            pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#            path="${name_m[$m]}/analbdy.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#          fi
        done
      done
    fi
  fi

  #-------------------
  # stage-out
  #-------------------

#    # topo
#    #-------------------
#    if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#      path="topo.d01.nc"
#      pathout="${OUTDIR}/const/topo.nc"
##      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
#    fi

#    # landuse
#    #-------------------
#    if ((loop == 1 && LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
#      path="landuse.d01.nc"
#      pathout="${OUTDIR}/const/landuse.nc"
##      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
#    fi

#    # bdy
#    #-------------------
#    if ((loop == 1 && BDY_FORMAT != 0)); then
#      if ((BDY_ENS == 0)); then
#        path="mean/bdy.nc"
#        pathout="${OUTDIR}/${time}/bdy/mean.boundary.nc"
##        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}
#      elif ((BDY_ENS == 1)); then
#        for m in $(seq $mtot); do
#          path="${name_m[$m]}/bdy.nc"
#          pathout="${OUTDIR}/${time}/bdy/${name_m[$m]}.boundary.nc"
##          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+1))]}
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+1))]}
#        done
#      fi
#    fi

  # anal (initial time)
  #-------------------
  if ((loop == 1 && MAKEINIT == 1)); then
    for m in $(seq $mtot); do
      for q in $(seq $mem_np_); do
        path="${name_m[$m]}/anal.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
      done
    done
  fi

  # anal
  #-------------------
  if ((OUT_OPT <= 4 || (OUT_OPT <= 5 && loop % OUT_CYCLE_SKIP == 0) || atime > ETIME)); then
    mlist=$(seq $mtot)
  elif ((OUT_OPT <= 7)); then
    mlist="$mmean"
    if ((DET_RUN == 1)); then
      mlist="$mlist $mmdet"
    fi
  fi
  for m in $mlist; do
    for q in $(seq $mem_np_); do
      path="${name_m[$m]}/anal.d01_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
      pathout="${OUTDIR}/${atime}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
      if ((m == mmean && SPRD_OUT == 1)); then
        path="sprd/anal.d01_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR}/${atime}/anal/sprd${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
      fi
    done
  done

  # gues
  #-------------------
  if ((OUT_OPT <= 3)); then
    mlist=$(seq $mtot)
  elif ((OUT_OPT <= 6)); then
    mlist="$mmean"
    if ((DET_RUN == 1)); then
      mlist="$mlist $mmdet"
    fi
  fi
  for m in $mlist; do
    for q in $(seq $mem_np_); do
      path="${name_m[$m]}/gues.d01_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
      pathout="${OUTDIR}/${atime}/gues/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
      if ((m == mmean && SPRD_OUT == 1)); then
        path="sprd/gues.d01_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR}/${atime}/gues/sprd${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
      fi
    done
  done

  # hist
  #-------------------
  if ((OUT_OPT <= 1)); then
    mlist=$(seq $mtot)
  elif ((OUT_OPT <= 2)); then
    mlist="$mmean"
    if ((DET_RUN == 1)); then
      mlist="$mlist $mmdet"
    fi
  fi
  for m in $mlist; do
    for q in $(seq $mem_np_); do
      path="${name_m[$m]}/hist.d01_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
      pathout="${OUTDIR}/${time}/hist/${name_m[$m]}${CONNECTOR}history$(scale_filename_sfx $((q-1)))"
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
    done
  done

#    # diag
#    #-------------------
#    if ((RTPS_INFL_OUT == 1)); then
#      path="rtpsinfl.d01_$(datetime_scale $atime).nc"
#      pathout="${OUTDIR}/${atime}/diag/rtpsinfl.init.nc"
##      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+1))]}
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+1))]}
#    fi
#    if ((NOBS_OUT == 1)); then
#      path="nobs.d01_$(datetime_scale $atime).nc"
#      pathout="${OUTDIR}/${atime}/diag/nobs.init.nc"
##      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+1))]}
#      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+1))]}
#    fi

#    # obsgues
#    #-------------------
#    if ((OBSOUT_OPT <= 2)); then
#      for m in $(seq $mtot); do ###### either $mmean or $mmdet ? ######
#        path="${name_m[$m]}/obsgues.d01_${atime}.dat"
#        pathout="${OUTDIR}/${atime}/obsgues/${name_m[$m]}.obsda.dat"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+1))]}
#      done
#    fi

  # log
  #-------------------
  if [ "$MPI_TYPE" = 'K' ]; then
    log_nfmt='.%d'
  else
    log_nfmt="-${PROCESS_FMT}"
  fi

  if ((LOG_OPT <= 3)); then
    if ((LOG_TYPE == 1)); then
      mlist='1'
      plist='1'
    else
      mlist=$(seq $mtot)
      plist=$(seq $totalnp)
    fi
    for m in $mlist; do
      path="log/scale.${name_m[$m]}.LOG_${time}${SCALE_SFX_NONC_0}"
      pathout="${OUTDIR}/${time}/log/scale/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+1))]}
      path="log/scale.${name_m[$m]}.monitor_${time}${SCALE_SFX_NONC_0}"
      pathout="${OUTDIR}/${time}/log/scale/${name_m[$m]}_monitor${SCALE_SFX_NONC_0}"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+1))]}
    done
    for p in $plist; do
      if ((nitmax == 1)); then
        path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
        pathout="${OUTDIR}/${time}/log/scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
      else
        for it in $(seq $nitmax); do
          path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
          pathout="${OUTDIR}/${time}/log/scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
        done
      fi
    done
  fi

  if ((LOG_OPT <= 4)); then
    if ((LOG_TYPE == 1)); then
      plist='1'
    else
      plist=$(seq $totalnp)
    fi
    for p in $plist; do
      path="log/letkf.NOUT_${atime}$(printf -- "${log_nfmt}" $((p-1)))"
      pathout="${OUTDIR}/${atime}/log/letkf/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
    done
  fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
}

#===============================================================================

config_file_list () {
#-------------------------------------------------------------------------------
# Prepare all runtime configuration files
#
# Usage: config_file_list [CONFIG_DIR]
#
#   CONFIG_DIR  Temporary directory of configuration files to be staged to $TMPROOT
#               '-': Do not use a temporary directory and stage;
#                    output configuration files directly to $TMPROOT
#
# Other input variables:
#   $TMPROOT
#   $STAGING_DIR
#-------------------------------------------------------------------------------

local CONFIG_DIR="${1:--}"

local stage_config=1
if [ "$CONFIG_DIR" = '-' ]; then
  CONFIG_DIR="$TMPROOT"
  stage_config=0
fi

#-------------------------------------------------------------------------------

if ((PNETCDF == 1)); then
  local mem_np_=1
else
  local mem_np_=$mem_np
fi

#-------------------------------------------------------------------------------

echo
echo "Generate configration files..."

mkdir -p $CONFIG_DIR

FILE_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
fi

time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))

#  for s in $(seq $nsteps); do
#    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

#      if ((s == 1)); then
#        if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
#          continue
#        elif ((BDY_FORMAT == 0)); then
#          continue
#        elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
#          continue
#        fi
#      fi
#      if ((s == 2)); then
#        if ((BDY_FORMAT == 0)); then
#          continue
#        fi
#      fi
#      if ((s == 4)); then
#        if ((OBSOPE_RUN == 0)); then
#          continue
#        fi
#      fi

#    fi
#  done

  obstime $time

  #-----------------------------------------------------------------------------
  # scale (launcher)
  #-----------------------------------------------------------------------------

  for it in $(seq $nitmax); do
    if ((nitmax == 1)); then
      conf_file="scale-rm_ens_${time}.conf"
    else
      conf_file="scale-rm_ens_${time}_${it}.conf"
    fi
    echo "  $conf_file"
    cat $SCRP_DIR/config.nml.ensmodel | \
        sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
            -e "/!--MEMBER_RUN--/a MEMBER_RUN = $mtot," \
            -e "/!--MEMBER_ITER--/a MEMBER_ITER = $it," \
            -e "/!--CONF_FILES--/a CONF_FILES = \"@@@@/run_${time}.conf\"," \
            -e "/!--NNODES--/a NNODES = $NNODES_APPAR," \
            -e "/!--PPN--/a PPN = $PPN_APPAR," \
            -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
            -e "/!--MEM_NP--/a MEM_NP = $mem_np," \
        > $CONFIG_DIR/${conf_file}
    if ((stage_config == 1)); then
      echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
    fi
  done

  #-----------------------------------------------------------------------------
  # scale (each member)
  #-----------------------------------------------------------------------------

  for m in $(seq $mtot); do
    if [ "${name_m[$m]}" = 'mean' ]; then ###### using a variable for 'mean', 'mdet', 'sprd'
      RESTART_OUT_ADDITIONAL_COPIES=1
      RESTART_OUT_ADDITIONAL_BASENAME="\"mean/gues.d01\", "
      if ((SPRD_OUT == 1)); then
        RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+2))
        RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"sprd/anal.d01\", "
        RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"sprd/gues.d01\", "
      fi
#          if ((RTPS_INFL_OUT == 1)); then
#            RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
#            RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"rtpsinfl.d01\", "
#          fi
#          if ((NOBS_OUT == 1)); then
#            RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
#            RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"nobs.d01\", "
#          fi
    elif [ "${name_m[$m]}" = 'mdet' ]; then
      RESTART_OUT_ADDITIONAL_COPIES=1
      RESTART_OUT_ADDITIONAL_BASENAME="\"mdet/gues.d01\", "
    elif ((OUT_OPT <= 3)); then
      RESTART_OUT_ADDITIONAL_COPIES=1
      RESTART_OUT_ADDITIONAL_BASENAME="\"${name_m[$m]}/gues.d01\", "
    else
      RESTART_OUT_ADDITIONAL_COPIES=0
      RESTART_OUT_ADDITIONAL_BASENAME=
    fi
    if ((BDY_ENS == 1)); then
      mem_bdy=${name_m[$m]}
    else
      mem_bdy='mean'
    fi
    DOMAIN_CATALOGUE_OUTPUT=".false."
    if ((m == 1)); then
      DOMAIN_CATALOGUE_OUTPUT=".true."
    fi

    conf_file="${name_m[$m]}/run_${time}.conf"
    echo "  $conf_file"
    mkdir -p $CONFIG_DIR/${name_m[$m]}
    cat $SCRP_DIR/config.nml.scale | \
        sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale.${name_m[$m]}.LOG_${time}\"," \
            -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
            -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
            -e "/!--TIME_DURATION--/a TIME_DURATION = ${CYCLEFLEN}.D0," \
            -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${LCYCLE}.D0," \
            -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${LCYCLE}.D0," \
            -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${LCYCLE}.D0," \
            -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${LCYCLE}.D0," \
            -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${name_m[$m]}/anal.d01\"," \
            -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = .true.," \
            -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${name_m[$m]}/anal.d01\"," \
            -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo.d01\"," \
            -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse.d01\"," \
            -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${mem_bdy}/bdy_$(datetime_scale $time)\"," \
            -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
            -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0," \
            -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${name_m[$m]}/hist.d01_$(datetime_scale $time)\"," \
            -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${CYCLEFOUT}.D0," \
            -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"log/scale.${name_m[$m]}.monitor_${time}\"," \
            -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\"," \
            -e "/!--DOMAIN_CATALOGUE_FNAME--/a DOMAIN_CATALOGUE_FNAME = \"latlon_domain_catalogue.txt\"," \
            -e "/!--DOMAIN_CATALOGUE_OUTPUT--/a DOMAIN_CATALOGUE_OUTPUT = ${DOMAIN_CATALOGUE_OUTPUT}," \
            -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAG.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAPC.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/VARDATA.RM29\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/cira.nc\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPROOT_CONSTDB}/dat/rad/MIPAS\"," \
            -e "/!--TIME_END_RESTART_OUT--/a TIME_END_RESTART_OUT = .false.," \
            -e "/!--RESTART_OUT_ADDITIONAL_COPIES--/a RESTART_OUT_ADDITIONAL_COPIES = ${RESTART_OUT_ADDITIONAL_COPIES}," \
            -e "/!--RESTART_OUT_ADDITIONAL_BASENAME--/a RESTART_OUT_ADDITIONAL_BASENAME = ${RESTART_OUT_ADDITIONAL_BASENAME}" \
        > $CONFIG_DIR/${conf_file}
#    cat $SCRP_DIR/config.nml.scale_user | \
#        sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"XXXXXX\"," \
#        sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"XXXXXX\"," \
#        >> $CONFIG_DIR/${conf_file}
    if ((stage_config == 1)); then
      if ((DISK_MODE == 3)); then
        for q in $(seq $mem_np_); do
          echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
      else
        echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
      fi
    fi
  done

  #-----------------------------------------------------------------------------
  # letkf
  #-----------------------------------------------------------------------------

  OBS_IN_NAME_LIST=
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ]; then
      OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'${TMPROOT_OBS}/obs.${OBSNAME[$iobs]}_${atime}.dat', "
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

  DET_RUN_TF='.false.'
  if ((DET_RUN == 1)); then
    DET_RUN_TF='.true.'
  fi
  OBSDA_OUT='.false.'
  if ((OBSOUT_OPT <= 2)); then
    OBSDA_OUT='.true.'
  fi
  SPRD_OUT_TF='.true.'
  if ((SPRD_OUT == 0)); then
    SPRD_OUT_TF='.false.'
  fi
  RTPS_INFL_OUT_TF='.false.'
  if ((RTPS_INFL_OUT == 1)); then
    RTPS_INFL_OUT_TF='.true.'
  fi
  NOBS_OUT_TF='.false.'
  if ((NOBS_OUT == 1)); then
    NOBS_OUT_TF='.true.'
  fi

  conf_file="letkf_${atime}.conf"
  echo "  $conf_file"
  cat $SCRP_DIR/config.nml.letkf | \
      sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
          -e "/!--DET_RUN--/a DET_RUN = ${DET_RUN_TF}," \
          -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
          -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
          -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
          -e "/!--OBSDA_OUT--/a OBSDA_OUT = $OBSDA_OUT" \
          -e "/!--OBSDA_OUT_BASENAME--/a OBSDA_OUT_BASENAME = \"@@@@/obsgues.d01_${atime}\"," \
          -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = \"@@@@/hist.d01_$(datetime_scale $time)\"," \
          -e "/!--SLOT_START--/a SLOT_START = $slot_s," \
          -e "/!--SLOT_END--/a SLOT_END = $slot_e," \
          -e "/!--SLOT_BASE--/a SLOT_BASE = $slot_b," \
          -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = ${LTIMESLOT}.D0," \
          -e "/!--OBSDA_IN--/a OBSDA_IN = .false.," \
          -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"@@@@/anal.d01_$(datetime_scale $atime)\"," \
          -e "/!--GUES_MEAN_INOUT_BASENAME--/a GUES_MEAN_INOUT_BASENAME = \"mean/gues.d01_$(datetime_scale $atime)\"," \
          -e "/!--GUES_SPRD_OUT_BASENAME--/a GUES_SPRD_OUT_BASENAME = \"sprd/gues.d01_$(datetime_scale $atime)\"," \
          -e "/!--GUES_SPRD_OUT--/a GUES_SPRD_OUT = ${SPRD_OUT_TF}," \
          -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"@@@@/anal.d01_$(datetime_scale $atime)\"," \
          -e "/!--ANAL_SPRD_OUT--/a ANAL_SPRD_OUT = ${SPRD_OUT_TF}," \
          -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"topo.d01\"," \
          -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RTPS_INFL_OUT_TF}," \
          -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"rtpsinfl.d01_$(datetime_scale $atime).nc\"," \
          -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
          -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"nobs.d01_$(datetime_scale $atime).nc\"," \
          -e "/!--NNODES--/a NNODES = $NNODES_APPAR," \
          -e "/!--PPN--/a PPN = $PPN_APPAR," \
          -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
          -e "/!--MEM_NP--/a MEM_NP = $mem_np," \
          -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
          -e "/!--STIME--/a STIME = ${STIME:0:4}, ${STIME:4:2}, ${STIME:6:2}, ${STIME:8:2}, ${STIME:10:2}, ${STIME:12:2}, " \
          -e "/!--ETIME--/a ETIME = ${ETIME:0:4}, ${ETIME:4:2}, ${ETIME:6:2}, ${ETIME:8:2}, ${ETIME:10:2}, ${ETIME:12:2}, " \
          -e "/!--LCYCLE--/a LCYCLE = ${LCYCLE}.D0, " \
      > $CONFIG_DIR/${conf_file}
  # Most of these parameters are not important for letkf
  cat $SCRP_DIR/config.nml.scale | \
      sed -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
      >> $CONFIG_DIR/${conf_file}
  if ((stage_config == 1)); then
    echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
  fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

echo

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

if ((CLEAR_TMP == 1)); then
  local save_cmd='mv -f'
else
  local save_cmd='cp -f'
fi

#-------------------------------------------------------------------------------

time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))

  #-----------------------------------------------------------------------------
  # scale
  #-----------------------------------------------------------------------------

  if [ -d "${OUTDIR}/${time}/log/scale" ]; then
    # launcher
    for it in $(seq $nitmax); do
      if ((nitmax == 1)); then
        $save_cmd $CONFIG_DIR/scale-rm_ens_${time}.conf ${OUTDIR}/${time}/log/scale/scale-rm_ens.conf
      else
        $save_cmd $CONFIG_DIR/scale-rm_ens_${time}_${it}.conf ${OUTDIR}/${time}/log/scale/scale-rm_ens_${it}.conf
      fi
    done

    # each member
    for m in $(seq $mtot); do
      $save_cmd $CONFIG_DIR/${name_m[$m]}/run_${time}.conf ${OUTDIR}/${time}/log/scale/${name_m[$m]}_run.conf
    done
  fi

  #-----------------------------------------------------------------------------
  # letkf
  #-----------------------------------------------------------------------------

  if [ -d "${OUTDIR}/${atime}/log/letkf" ]; then
    $save_cmd $CONFIG_DIR/letkf_${atime}.conf ${OUTDIR}/${atime}/log/letkf/letkf.conf
  fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
}

#===============================================================================
