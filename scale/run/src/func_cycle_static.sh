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

declare -a mem_np_
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
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${INDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
          path="${name_m[$m]}/anal.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        done
      done
    done
  fi

  # topo
  #-------------------
  if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
    if ((DISK_MODE == 3)); then
      for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
            path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      done
    else
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
          path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
        done
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
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
            path="landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      done
    else
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
          path="landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
        done
      done
    fi
  fi

  # bdy (prepared)
  #-------------------
  if ((BDY_FORMAT == 0)); then
    if ((BDY_ENS == 0)); then
      if ((DISK_MODE == 3)); then
        for m in $(seq $((repeat_mems <= mtot ? repeat_mems : mtot))); do
          for q in $(seq ${mem_np_[1]}); do
            pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
            path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
          done
          if ((USE_INIT_FROM_BDY == 1)); then
            for d in $(seq $DOMNUM); do
              for q in $(seq ${mem_np_[$d]}); do
                pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
                path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              done
            done
          fi
        done
      else
        for q in $(seq ${mem_np_[1]}); do
          pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
        done
        if ((USE_INIT_FROM_BDY == 1)); then
          for d in $(seq $DOMNUM); do
            for q in $(seq ${mem_np_[$d]}); do
              pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
              path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
            done
          done
        fi
      fi
    elif ((BDY_ENS == 1)); then
      for m in $(seq $mtot); do
        for q in $(seq ${mem_np_[1]}); do
          pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
          path="${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
        done
        if ((USE_INIT_FROM_BDY == 1)); then
          for d in $(seq $DOMNUM); do
            for q in $(seq ${mem_np_[$d]}); do
              pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
              path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            done
          done
        fi
      done
    fi
  fi

  #-------------------
  # stage-out
  #-------------------

  # anal (initial time)
  #-------------------
  if ((loop == 1 && MAKEINIT == 1)); then
    for m in $(seq $mtot); do
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          pathout="${OUTDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        done
      done
    done
  fi

  # topo
  #-------------------
  if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
      done
    done
  fi

  # landuse
  #-------------------
  if ((loop == 1 && LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
      done
    done
  fi

  # bdy
  #-------------------
  if ((BDY_FORMAT != 0)); then
    if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
      for m in $(seq $mtot); do
        for q in $(seq ${mem_np_[1]}); do
          path="${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
          pathout="${OUTDIR[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
        done
        if ((USE_INIT_FROM_BDY == 1)); then
          for d in $(seq $DOMNUM); do
            for q in $(seq ${mem_np_[$d]}); do
              path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              pathout="${OUTDIR[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            done
          done
        fi
      done
    elif ((BDYOUT_OPT <= 2)); then
      for q in $(seq ${mem_np_[1]}); do
        path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$q]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+q))]}
      done
      if ((USE_INIT_FROM_BDY == 1)); then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      fi
    fi
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
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="${name_m[$m]}/anal.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[$d]}/${atime}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        if ((m == mmean && SPRD_OUT == 1)); then
          path="sprd/anal.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
          pathout="${OUTDIR[$d]}/${atime}/anal/sprd${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        fi
      done
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
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="${name_m[$m]}/gues.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[$d]}/${atime}/gues/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        if ((m == mmean && SPRD_OUT == 1)); then
          path="sprd/gues.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $atime)$(scale_filename_sfx $((q-1)))"
          pathout="${OUTDIR[$d]}/${atime}/gues/sprd${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        fi
      done
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
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="${name_m[$m]}/hist.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
        pathout="${OUTDIR[$d]}/${time}/hist/${name_m[$m]}${CONNECTOR}history$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
      done
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
      for d in $(seq $DOMNUM); do
        path="log/scale_init.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
        pathout="${OUTDIR[$d]}/${time}/log/scale_init/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
        path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
        pathout="${OUTDIR[$d]}/${time}/log/scale/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
        path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).monitor_${time}${SCALE_SFX_NONC_0}"
        pathout="${OUTDIR[$d]}/${time}/log/scale/${name_m[$m]}_monitor${SCALE_SFX_NONC_0}"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
      done
    done
    for p in $plist; do
      if ((nitmax == 1)); then
        path="log/scale-rm_init_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
        pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
        path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
        pathout="${OUTDIR[1]}/${time}/log/scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
      else
        for it in $(seq $nitmax); do
          path="log/scale-rm_init_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
          pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
          path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
          pathout="${OUTDIR[1]}/${time}/log/scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
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
      pathout="${OUTDIR[1]}/${atime}/log/letkf/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
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

config_file_scale_launcher () {
#-------------------------------------------------------------------------------
# Generate the launcher configuration files for scale_pp/scale_init/scale
#
# Usage: config_file_scale_launcher MODEL_NAME CONF_NAME
#
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

local MODEL_NAME="$1"; shift
local CONF_NAME="$1"; shift
local MEMBER_RUN="$1"

#-------------------------------------------------------------------------------

local it
local conf_file

for it in $(seq $nitmax); do
  if ((nitmax == 1)); then
    conf_file="${MODEL_NAME}_ens_${time}.conf"
  else
    conf_file="${MODEL_NAME}_ens_${time}_${it}.conf"
  fi
  echo "  $conf_file"
  cat $SCRP_DIR/config.nml.ensmodel | \
      sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
          -e "/!--MEMBER_RUN--/a MEMBER_RUN = $MEMBER_RUN," \
          -e "/!--MEMBER_ITER--/a MEMBER_ITER = $it," \
          -e "/!--CONF_FILES--/a CONF_FILES = \"<member>/${CONF_NAME}.d<domain>_${time}.conf\"," \
          -e "/!--PPN--/a PPN = $PPN_APPAR," \
          -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
          -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
          -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
      > $CONFIG_DIR/${conf_file}
  if ((stage_config == 1)); then
    echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
  fi
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
if ((PNETCDF_BDY_SCALE == 1)); then
  local mem_np_bdy_=1
else
  local mem_np_bdy_=$((DATA_BDY_SCALE_PRC_NUM_X*DATA_BDY_SCALE_PRC_NUM_Y))
fi

#-------------------------------------------------------------------------------

echo
echo "Generate configration files..."

mkdir -p $CONFIG_DIR

IO_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  IO_AGGREGATE=".true."
fi

PRC_DOMAINS_LIST=
for d in $(seq $DOMNUM); do
  PRC_DOMAINS_LIST="$PRC_DOMAINS_LIST${SCALE_NP[$d]}, "
done

time=$STIME
atime=$(datetime $time $LCYCLE s)
time_bdy_start_prev=0
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

  bdy_setting $time $CYCLEFLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

  if ((BDY_FORMAT != 0)); then

    #---------------------------------------------------------------------------
    # scale_init (launcher)
    #---------------------------------------------------------------------------

    if ((BDY_ENS == 1)); then
      config_file_scale_launcher scale-rm_init init $mtot
    elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
      config_file_scale_launcher scale-rm_init init 1
    else # local run directory: run multiple members as needed
      config_file_scale_launcher scale-rm_init init $((repeat_mems <= mtot ? repeat_mems : mtot))
    fi

    #---------------------------------------------------------------------------
    # scale_init (each member)
    #---------------------------------------------------------------------------

    if (((loop == 1 && MAKEINIT == 1) || USE_INIT_FROM_BDY == 1)); then
      RESTART_OUTPUT='.true.'
    else
      RESTART_OUTPUT='.false.'
    fi
    if (((loop == 1 && MAKEINIT == 1) && ${bdy_times[1]} != time)); then
      echo "[Error] $0: Unable to generate initial analyses (MAKEINIT) at this time" >&2
      echo "        that does not fit to any boundary data." >&2
      exit 1
    fi

    if ((BDY_ROTATING == 1 || ${bdy_times[1]} != time_bdy_start_prev)); then
      time_bdy_start_prev=${bdy_times[1]}
      nbdy_max=0
    fi
    if ((nbdy > nbdy_max)); then
      for ibdy in $(seq $((nbdy_max+1)) $nbdy); do
        time_bdy=${bdy_times[$ibdy]}

#        if ((BDY_FORMAT == 1)); then

          if ((BDY_ENS == 1)); then
            for m in $(seq $mtot); do
              if ((m == mmean)); then
                mem_bdy="$BDY_MEAN"
              else
                mem_bdy="${name_m[$m]}"
              fi
              for q in $(seq $mem_np_bdy_); do
                pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${mem_bdy}${CONNECTOR}history$(scale_filename_bdy_sfx $((q-1)))"
                path="${name_m[$m]}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
              done
            done
          else
            for q in $(seq $mem_np_bdy_); do
              pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}${CONNECTOR}history$(scale_filename_bdy_sfx $((q-1)))"
              path="mean/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            done
          fi

#        elif ((BDY_FORMAT == 2 || BDY_FORMAT == 4)); then

#          if ((BDY_FORMAT == 2)); then
#            data_bdy_i=$DATA_BDY_WRF
#            filenum=1
#            filename_prefix[1]='wrfout_'
#            filename_suffix[1]=''
#          elif ((BDY_FORMAT == 4)); then
#            data_bdy_i=$DATA_BDY_GRADS
#            filenum=3
#            filename_prefix[1]='atm_'
#            filename_suffix[1]='.grd'
#            filename_prefix[2]='sfc_'
#            filename_suffix[2]='.grd'
#            filename_prefix[3]='land_'
#            filename_suffix[3]='.grd'
#          fi

#          if ((BDY_ENS == 1)); then
#            for m in $(seq $mtot); do
#              for ifile in $(seq $filenum); do
#                if ((BDY_ROTATING == 1)); then
#                  pathin="$data_bdy_i/${time}/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
#                  path="bdyorg/${time}/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}"
#                else
#                  pathin="$data_bdy_i/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}"
#                  path="bdyorg/const/${name_m[$m]}/${filename_prefix[$ifile]}${time_bdy}"
#                fi
#                echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#              done
#            done
#          else
#            for ifile in $(seq $filenum); do
#              if ((BDY_ROTATING == 1)); then
#                pathin="$data_bdy_i/${time}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
#                path="bdyorg/${time}/mean/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
#              else
#                pathin="$data_bdy_i/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
#                path="bdyorg/const/mean/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
#              fi
#              echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#            done
#          fi

#        fi

  #      if ((BDY_ROTATING == 1)); then
  #        bdyorg_path="$(cd "${BDYORG}/${STIME}" && pwd)"
  #      else
  #        bdyorg_path="$(cd "${BDYORG}/const" && pwd)"
  #      fi
  #      if ((BDY_FORMAT == 2)); then
  #        if [ -s "${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}" ]; then
  #          ln -fs "${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}" $TMPDIR/bdydata${file_number}
  #        else
  #          echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/wrfout_${time_bdy}'."
  #          exit 1
  #        fi
  #      elif ((BDY_FORMAT == 4)); then
  #        if [ -s "${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd" ]; then
  #          ln -fs "${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd" $TMPDIR/bdyatm${file_number}.grd
  #        else
  #          echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/atm_${time_bdy}.grd'."
  #          exit 1
  #        fi
  #        if [ -s "${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd" ]; then
  #          ln -fs "${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd" $TMPDIR/bdysfc${file_number}.grd
  #        else
  #          echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/sfc_${time_bdy}.grd'."
  #          exit 1
  #        fi
  #        if [ -s "${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd" ]; then
  #          ln -fs "${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd" $TMPDIR/bdyland${file_number}.grd
  #        else
  #          echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${MEM_BDY}/land_${time_bdy}.grd'."
  #          exit 1
  #        fi
  #      fi
  #      i=$((i+1))

      done # [ ibdy in $(seq $((nbdy_max+1)) $nbdy) ]
      nbdy_max=$nbdy
    fi

    for m in $(seq $mtot); do
      if ((BDY_ENS == 1)); then
        mem_bdy=${name_m[$m]}
      else
        mem_bdy='mean'
      fi

      if ((BDY_FORMAT == 1)); then
        FILETYPE_ORG='SCALE-RM'
      elif ((BDY_FORMAT == 2)); then
        FILETYPE_ORG='WRF-ARW'
      elif ((BDY_FORMAT == 4)); then
        FILETYPE_ORG='GrADS'
      else
        echo "[Error] $0: Unsupport boundary file types." >&2
        exit 1
      fi
      if ((BDY_FORMAT == 4)); then
        BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/gradsbdy.conf"
      else
        if ((nbdy <= 1)); then
          BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d 0)"
        else
          BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)"
        fi
      fi

      for d in $(seq $DOMNUM); do
        dfmt=$(printf $DOMAIN_FMT $d)

        if ((d == 1)); then
          conf_file_src=$SCRP_DIR/config.nml.scale_init
        else
          conf_file_src=$SCRP_DIR/config.nml.scale_init.d$d
        fi
        mkdir -p $CONFIG_DIR/${name_m[$m]}
        conf="$(cat $conf_file_src | \
            sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale_init.${name_m[$m]}.d${dfmt}.LOG_${time}\"," \
                -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = ${IO_AGGREGATE}," \
                -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
                -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = $RESTART_OUTPUT," \
                -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${name_m[$m]}/init.d${dfmt}\"," \
                -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo.d${dfmt}\"," \
                -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse.d${dfmt}\"," \
                -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\",")"
        if ((BDY_FORMAT == 1)); then
          conf="$(echo "$conf" | \
              sed -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d 0)\"," \
                  -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
                  -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
                  -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${TMPROOT_BDYDATA}/latlon_domain_catalogue.bdy.txt\",")"
        fi
        conf="$(echo "$conf" | \
          sed -e "/!--BASENAME_ORG--/a BASENAME_ORG = \"${BASENAME_ORG}\"," \
              -e "/!--FILETYPE_ORG--/a FILETYPE_ORG = \"${FILETYPE_ORG}\"," \
              -e "/!--BOUNDARY_UPDATE_DT--/a BOUNDARY_UPDATE_DT = ${BDYINT}.D0,")"
        if ((d == 1)); then
          conf="$(echo "$conf" | \
            sed -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"${mem_bdy}/bdy_$(datetime_scale $time)\"," \
                -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${nbdy}," \
                -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${ntsteps}," \
                -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${ntsteps_skip},")"
        else
          #------ before SCALE 5.2
          conf="$(echo "$conf" | \
            sed -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"${mem_bdy}/bdy.d${dfmt}_$(datetime_scale $time)\"," \
                -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${nbdy}," \
                -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${ntsteps}," \
                -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${ntsteps_skip},")"
          #------ after SCALE 5.3
          #conf="$(echo "$conf" | \
          #  sed -e "/!--MAKE_BOUNDARY--/a MAKE_BOUNDARY = .false.,")"
          #------
        fi
        conf_file="${name_m[$m]}/init.d${dfmt}_${time}.conf"
        echo "  $conf_file"
        echo "$conf" > $CONFIG_DIR/${conf_file}

        if ((stage_config == 1)); then
          if ((DISK_MODE == 3)); then
            for q in $(seq $mem_np_); do
              echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
            done
          else
            echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
          fi
        fi

#        if ((BDY_FORMAT == 4)); then
#          conf_file="${mem_bdy}/gradsbdy.conf"
#          cat $SCRP_DIR/config.nml.grads_boundary | \
#              sed -e "s#--DIR--#${TMPROOT_BDYDATA}/${mem_bdy}/......#g" \
#              > $CONFIG_DIR/${conf_file}

#          if ((stage_config == 1)); then
#            echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#          fi
#        fi
      done # [ d in $(seq $DOMNUM) ]
    done # [ m in $(seq $mtot) ]

  fi # [ BDY_FORMAT != 0 ]

  #-----------------------------------------------------------------------------
  # scale (launcher)
  #-----------------------------------------------------------------------------

  config_file_scale_launcher scale-rm run $mtot

  #-----------------------------------------------------------------------------
  # scale (each member)
  #-----------------------------------------------------------------------------

  for m in $(seq $mtot); do
    DOMAIN_CATALOGUE_OUTPUT=".false."
    if ((m == 1)); then
      DOMAIN_CATALOGUE_OUTPUT=".true."
    fi

    for d in $(seq $DOMNUM); do
      dfmt=$(printf $DOMAIN_FMT $d)

      ONLINE_IAM_PARENT=".false."
      if ((d < DOMNUM)); then
        ONLINE_IAM_PARENT=".true."
      fi
      ONLINE_IAM_DAUGHTER=".false."
      if ((d > 1)); then
        ONLINE_IAM_DAUGHTER=".true."
      fi
      if ((loop == 1 && MAKEINIT == 1)); then
        RESTART_IN_BASENAME="${name_m[$m]}/init.d${dfmt}"
      else
        RESTART_IN_BASENAME="${name_m[$m]}/anal.d${dfmt}"
      fi
      if [ "${name_m[$m]}" = 'mean' ]; then ###### using a variable for 'mean', 'mdet', 'sprd'
        RESTART_OUT_ADDITIONAL_COPIES=1
        RESTART_OUT_ADDITIONAL_BASENAME="\"mean/gues.d$dfmt\", "
        if ((SPRD_OUT == 1)); then
          RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+2))
          RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"sprd/anal.d$dfmt\", "
          RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"sprd/gues.d$dfmt\", "
        fi
#        if ((RTPS_INFL_OUT == 1)); then
#          RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
#          RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"rtpsinfl.d$dfmt\", "
#        fi
#        if ((NOBS_OUT == 1)); then
#          RESTART_OUT_ADDITIONAL_COPIES=$((RESTART_OUT_ADDITIONAL_COPIES+1))
#          RESTART_OUT_ADDITIONAL_BASENAME="$RESTART_OUT_ADDITIONAL_BASENAME\"nobs.d$dfmt\", "
#        fi
      elif [ "${name_m[$m]}" = 'mdet' ]; then
        RESTART_OUT_ADDITIONAL_COPIES=1
        RESTART_OUT_ADDITIONAL_BASENAME="\"mdet/gues.d$dfmt\", "
      elif ((OUT_OPT <= 3)); then
        RESTART_OUT_ADDITIONAL_COPIES=1
        RESTART_OUT_ADDITIONAL_BASENAME="\"${name_m[$m]}/gues.d$dfmt\", "
      else
        RESTART_OUT_ADDITIONAL_COPIES=0
        RESTART_OUT_ADDITIONAL_BASENAME=
      fi

      if ((d == 1)); then
        conf_file_src=$SCRP_DIR/config.nml.scale
      else
        conf_file_src=$SCRP_DIR/config.nml.scale.d$d
      fi
      mkdir -p $CONFIG_DIR/${name_m[$m]}
      conf="$(cat $conf_file_src | \
          sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale.${name_m[$m]}.d${dfmt}.LOG_${time}\"," \
              -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = ${IO_AGGREGATE}," \
              -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
              -e "/!--TIME_DURATION--/a TIME_DURATION = ${CYCLEFLEN}.D0," \
              -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${LCYCLE}.D0," \
              -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${LCYCLE}.D0," \
              -e "/!--ONLINE_DOMAIN_NUM--/a ONLINE_DOMAIN_NUM = ${d}," \
              -e "/!--ONLINE_IAM_PARENT--/a ONLINE_IAM_PARENT = ${ONLINE_IAM_PARENT}," \
              -e "/!--ONLINE_IAM_DAUGHTER--/a ONLINE_IAM_DAUGHTER = ${ONLINE_IAM_DAUGHTER}," \
              -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME}\"," \
              -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = .true.," \
              -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${name_m[$m]}/anal.d${dfmt}\"," \
              -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo.d${dfmt}\"," \
              -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse.d${dfmt}\"," \
              -e "/!--HISTORY_DEFAULT_BASENAME--/a HISTORY_DEFAULT_BASENAME = \"${name_m[$m]}/hist.d${dfmt}_$(datetime_scale $time)\"," \
              -e "/!--HISTORY_DEFAULT_TINTERVAL--/a HISTORY_DEFAULT_TINTERVAL = ${CYCLEFOUT}.D0," \
              -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"log/scale.${name_m[$m]}.d${dfmt}.monitor_${time}\"," \
              -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\"," \
              -e "/!--DOMAIN_CATALOGUE_FNAME--/a DOMAIN_CATALOGUE_FNAME = \"latlon_domain_catalogue.d${dfmt}.txt\"," \
              -e "/!--DOMAIN_CATALOGUE_OUTPUT--/a DOMAIN_CATALOGUE_OUTPUT = ${DOMAIN_CATALOGUE_OUTPUT}," \
              -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAG.29\"," \
              -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAPC.29\"," \
              -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/VARDATA.RM29\"," \
              -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/cira.nc\"," \
              -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPROOT_CONSTDB}/dat/rad/MIPAS\"," \
              -e "/!--TIME_END_RESTART_OUT--/a TIME_END_RESTART_OUT = .false.," \
              -e "/!--RESTART_OUT_ADDITIONAL_COPIES--/a RESTART_OUT_ADDITIONAL_COPIES = ${RESTART_OUT_ADDITIONAL_COPIES}," \
              -e "/!--RESTART_OUT_ADDITIONAL_BASENAME--/a RESTART_OUT_ADDITIONAL_BASENAME = ${RESTART_OUT_ADDITIONAL_BASENAME}")"
      if ((d == 1)); then
        if ((BDY_ENS == 1)); then
          mem_bdy=${name_m[$m]}
        else
          mem_bdy='mean'
        fi
        conf="$(echo "$conf" | \
            sed -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${mem_bdy}/bdy_$(datetime_scale $time)\"," \
                -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = ${bdy_start_time:0:4}, ${bdy_start_time:4:2}, ${bdy_start_time:6:2}, ${bdy_start_time:8:2}, ${bdy_start_time:10:2}, ${bdy_start_time:12:2}," \
                -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0,")"
      fi
      conf_file="${name_m[$m]}/run.d${dfmt}_${time}.conf"
      echo "  $conf_file"
      echo "$conf" > $CONFIG_DIR/${conf_file}

      conf="$(cat $SCRP_DIR/config.nml.scale_user)"
      if ((OCEAN_INPUT == 1)); then
        if ((OCEAN_FORMAT == 99)); then
          conf="$(echo "$conf" | \
              sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${name_m[$m]}/init.d${dfmt}\",")"
        fi
      fi
      if ((LAND_INPUT == 1)); then
        if ((LAND_FORMAT == 99)); then
          conf="$(echo "$conf" | \
              sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${name_m[$m]}/init.d${dfmt}\",")"
        fi
      fi
      echo "$conf" >> $CONFIG_DIR/${conf_file}

      if ((stage_config == 1)); then
        if ((DISK_MODE == 3)); then
          for q in $(seq $mem_np_); do
            echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
          done
        else
          echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
        fi
      fi
    done # [ d in $(seq $DOMNUM) ]
  done # [ m in $(seq $mtot) ]

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

  for d in $(seq $DOMNUM); do
    dfmt=$(printf $DOMAIN_FMT $d)

    if ((d == 1)); then
      conf_file_src=$SCRP_DIR/config.nml.letkf
      conf_file_src2=$SCRP_DIR/config.nml.scale
      conf_file="letkf_${atime}.conf"
    else
      conf_file_src=$SCRP_DIR/config.nml.letkf.d$d
      conf_file_src2=$SCRP_DIR/config.nml.scale.d$d
      conf_file="letkf.d${dfmt}_${atime}.conf"
    fi
    echo "  $conf_file"
    cat $SCRP_DIR/config.nml.ensmodel | \
        sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
            -e "/!--CONF_FILES--/a CONF_FILES = \"letkf.d<domain>_${atime}.conf\"," \
            -e "/!--DET_RUN--/a DET_RUN = ${DET_RUN_TF}," \
            -e "/!--PPN--/a PPN = $PPN_APPAR," \
            -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
            -e "/!--NUM_DOMAIN--/a NUM_DOMAIN = $DOMNUM," \
            -e "/!--PRC_DOMAINS--/a PRC_DOMAINS = $PRC_DOMAINS_LIST" \
        > $CONFIG_DIR/${conf_file}

    cat $conf_file_src | \
        sed -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
            -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
            -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
            -e "/!--OBSDA_OUT--/a OBSDA_OUT = $OBSDA_OUT" \
            -e "/!--OBSDA_OUT_BASENAME--/a OBSDA_OUT_BASENAME = \"<member>/obsgues.d${dfmt}_${atime}\"," \
            -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = \"<member>/hist.d${dfmt}_$(datetime_scale $time)\"," \
            -e "/!--SLOT_START--/a SLOT_START = $slot_s," \
            -e "/!--SLOT_END--/a SLOT_END = $slot_e," \
            -e "/!--SLOT_BASE--/a SLOT_BASE = $slot_b," \
            -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = ${LTIMESLOT}.D0," \
            -e "/!--OBSDA_IN--/a OBSDA_IN = .false.," \
            -e "/!--GUES_IN_BASENAME--/a GUES_IN_BASENAME = \"<member>/anal.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--GUES_MEAN_INOUT_BASENAME--/a GUES_MEAN_INOUT_BASENAME = \"mean/gues.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--GUES_SPRD_OUT_BASENAME--/a GUES_SPRD_OUT_BASENAME = \"sprd/gues.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--GUES_SPRD_OUT--/a GUES_SPRD_OUT = ${SPRD_OUT_TF}," \
            -e "/!--ANAL_OUT_BASENAME--/a ANAL_OUT_BASENAME = \"<member>/anal.d${dfmt}_$(datetime_scale $atime)\"," \
            -e "/!--ANAL_SPRD_OUT--/a ANAL_SPRD_OUT = ${SPRD_OUT_TF}," \
            -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"topo.d${dfmt}\"," \
            -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RTPS_INFL_OUT_TF}," \
            -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"rtpsinfl.d${dfmt}_$(datetime_scale $atime).nc\"," \
            -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
            -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"nobs.d${dfmt}_$(datetime_scale $atime).nc\"," \
        >> $CONFIG_DIR/${conf_file}

    # Most of these parameters are not important for letkf
    cat $conf_file_src2 | \
        sed -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = ${IO_AGGREGATE}," \
        >> $CONFIG_DIR/${conf_file}

    if ((stage_config == 1)); then
      echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
    fi
  done # [ d in $(seq $DOMNUM) ]

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

mkdir -p ${OUTDIR[1]}/config

cp -fr $CONFIG_DIR/* ${OUTDIR[1]}/config

#-------------------------------------------------------------------------------
}

#===============================================================================
