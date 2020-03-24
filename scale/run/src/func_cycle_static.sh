#!/bin/bash
#===============================================================================
#
#  Functions for 'cycle' jobs
#
#===============================================================================

staging_list_static () {
#-------------------------------------------------------------------------------
# Prepare all the staging list files
#
# Usage: staging_list_static
#-------------------------------------------------------------------------------
# common section

declare -a mem_np_

staging_list_common_static cycle

mkdir -p ${TMPROOT}/topo
mkdir -p ${TMPROOT}/landuse

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
      mkdir -p ${TMPROOT}/${name_m[$m]}
      for d in $(seq $DOMNUM); do
        ln -sf ${INDIR[$d]}/${time}/anal/${name_m[$m]}/*.nc ${TMPROOT}/${name_m[$m]}/
#        for q in $(seq ${mem_np_[$d]}); do
#          pathin="${INDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#          path="${name_m[$m]}/anal.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
##          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
#          ln -sf ${pathin} ${TMPROOT}/${name_m[$m]}/anal.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))
#        done
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
            path="topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      done
    else
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
          path="topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
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
            path="landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      done
    else
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
          path="landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
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
          pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
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

  # additive inflation
  #-------------------
  if ((loop == 1 && ADDINFL == 1)); then
    for m in $(seq $MEMBER); do
      for d in $(seq $DOMNUM); do
        for q in $(seq ${mem_np_[$d]}); do
          pathin="${DATA_ADDINFL[$d]}/const/addi/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
          path="${name_m[$m]}/addi.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
          echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
        done
      done
    done
  fi

  #-------------------
  # stage-out
  #-------------------

  # anal (initial time)
  #-------------------
#  if ((loop == 1 && MAKEINIT == 1)); then
#    for m in $(seq $mtot); do
#      for d in $(seq $DOMNUM); do
#        for q in $(seq ${mem_np_[$d]}); do
#          path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#          pathout="${OUTDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
##          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
#          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
#        done
#      done
#    done
#  fi

  # topo
  #-------------------
  if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
    for d in $(seq $DOMNUM); do
      for q in $(seq ${mem_np_[$d]}); do
        path="topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
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
        path="landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
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
        pathout="${OUTDIR[1]}/${time}/bdy/mean${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
#        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$q]}
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+q))]}
      done
      if ((USE_INIT_FROM_BDY == 1)); then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/${time}/bdy/mean${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
          done
        done
      fi
    fi
  fi

  # anal
  #-------------------
  mlist=
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
  mlist=
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
  mlist=
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
  if ((LOG_TYPE == 1)); then
    mlist='1'
    plist='1'
  else
    mlist=$(seq $mtot)
    if ((BDY_ENS == 1)); then
      mlist_init=$(seq $mtot)
    elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
      mlist_init='1'
    else # local run directory: run multiple members as needed
      mlist_init=$(seq $((repeat_mems <= mtot ? repeat_mems : mtot)))
    fi
    plist=$(seq $totalnp)
  fi

  if ((BDY_FORMAT != 0 && LOG_OPT <= 2)); then
    for m in $mlist_init; do
      for d in $(seq $DOMNUM); do
        path="log/scale_init.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
        pathout="${OUTDIR[$d]}/${time}/log/scale_init/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
      done
    done
    for p in $plist; do
      if ((nitmax == 1)); then
        path="log/scale-rm_init_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
        pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
      else
        for it in $(seq $((BDY_ENS == 1 ? nitmax : 1))); do
          path="log/scale-rm_init_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
          pathout="${OUTDIR[1]}/${time}/log/scale_init/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
        done
      fi
    done
  fi
  if ((LOG_OPT <= 3)); then
    for m in $mlist; do
      for d in $(seq $DOMNUM); do
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
        path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
        pathout="${OUTDIR[1]}/${time}/log/scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
        echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
      else
        for it in $(seq $nitmax); do
          path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
          pathout="${OUTDIR[1]}/${time}/log/scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
          echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
        done
      fi
    done
  fi
  if ((LOG_OPT <= 4)); then
    for p in $plist; do
      path="log/letkf.NOUT_${atime}$(printf -- "${log_nfmt}" $((p-1)))"
      pathout="${OUTDIR[1]}/${atime}/log/letkf/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
      echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
    done
  fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done # [ ((time <= ETIME)) ]

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

FILE_AGGREGATE=".false"
if ((PNETCDF == 1)); then
  FILE_AGGREGATE=".true."
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
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" $mtot
    elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" 1
    else # local run directory: run multiple members as needed
      config_file_scale_launcher cycle scale-rm_init_ens "<member>/init" $((repeat_mems <= mtot ? repeat_mems : mtot))
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

        if ((BDY_FORMAT == 1)); then

          if ((BDY_ENS == 1)); then
            for m in $(seq $mtot); do
              if ((m == mmean)); then
                mem_bdy="$BDY_MEAN"
              else
                mem_bdy="${name_m[$m]}"
              fi
              for q in $(seq $mem_np_bdy_); do
                pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${mem_bdy}${CONNECTOR_BDY}history$(scale_filename_bdy_sfx $((q-1)))"
                path="${name_m[$m]}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
              done
            done
          else
            for q in $(seq $mem_np_bdy_); do
              pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}${CONNECTOR_BDY}history$(scale_filename_bdy_sfx $((q-1)))"
              path="mean/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            done
          fi

        elif ((BDY_FORMAT == 2 || BDY_FORMAT == 4)); then

          if ((BDY_FORMAT == 2)); then
            data_bdy_i="$DATA_BDY_WRF"
            filenum=1
            filename_prefix[1]='wrfout_'
            filename_suffix[1]=''
            filenamein_prefix[1]=''
            filenamein_suffix[1]=''
          elif ((BDY_FORMAT == 4)); then
            data_bdy_i="$DATA_BDY_GRADS"
            filenum=3
            filename_prefix[1]='atm_'
            filename_suffix[1]='.grd'
            filenamein_prefix[1]='atm_'
            filenamein_suffix[1]='.grd'
            filename_prefix[2]='sfc_'
            filename_suffix[2]='.grd'
            filenamein_prefix[2]='sfc_'
            filenamein_suffix[2]='.grd'
            filename_prefix[3]='land_'
            filename_suffix[3]='.grd'
            filenamein_prefix[3]='lnd_'
            filenamein_suffix[3]='.grd'
          fi

          if ((BDY_ENS == 1)); then
            for m in $(seq $mtot); do
              if ((m == mmean)); then
                mem_bdy="$BDY_MEAN"
              else
                mem_bdy="${name_m[$m]}"
              fi
              for ifile in $(seq $filenum); do
                if ((BDY_ROTATING == 1)); then
                  pathin="${data_bdy_i}/${time}/${mem_bdy}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                else
                  pathin="${data_bdy_i}/${mem_bdy}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                fi
                path="${name_m[$m]}/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
              done
            done
          else
            for ifile in $(seq $filenum); do
              if ((BDY_ROTATING == 1)); then
                pathin="${data_bdy_i}/${time}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              else
                pathin="${data_bdy_i}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
              fi
              path="mean/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            done
          fi

        fi # [ BDY_FORMAT == 2 || BDY_FORMAT == 4 ]
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
        LATLON_CATALOGUE_FNAME="${TMPROOT_BDYDATA}/latlon_domain_catalogue.bdy.txt"
      elif ((BDY_FORMAT == 2)); then
        FILETYPE_ORG='WRF-ARW'
        LATLON_CATALOGUE_FNAME=
      elif ((BDY_FORMAT == 4)); then
        FILETYPE_ORG='GrADS'
        LATLON_CATALOGUE_FNAME=
      else
        echo "[Error] $0: Unsupport boundary file types." >&2
        exit 1
      fi
      if ((BDY_FORMAT == 4)); then
        BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/gradsbdy.conf"
      else
        if ((nbdy <= 1)); then
          bdy_no_suffix="_$(printf %05d 0)"
        else
          bdy_no_suffix=
        fi
        BASENAME_ORG="${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}"
      fi

      for d in $(seq $DOMNUM); do
        dfmt=$(printf $DOMAIN_FMT $d)

        if ((d == 1)); then
          conf_file_src=$SCRP_DIR/config.nml.scale_init
        else
          conf_file_src=$SCRP_DIR/config.nml.scale_init.d$d
        fi
        conf="$(cat $conf_file_src | \
            sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale_init.${name_m[$m]}.d${dfmt}.LOG_${time}\"," \
                -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
                -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
                -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
                -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${mem_bdy}/init.d${dfmt}\"," \
                -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo/topo.d${dfmt}\"," \
                -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse/landuse.d${dfmt}\"," \
                -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\",")"
        if ((BDY_FORMAT == 1)); then
          conf="$(echo "$conf" | \
              sed -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d 0)\"," \
                  -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
                  -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
                  -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${LATLON_CATALOGUE_FNAME}\",")"
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
        mkdir -p $CONFIG_DIR/${name_m[$m]}
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
      done # [ d in $(seq $DOMNUM) ]

      if ((BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1))); then
        mkdir -p $CONFIG_DIR/${mem_bdy}
        conf_file="${mem_bdy}/gradsbdy.conf"
        echo "  $conf_file"
        if ((nbdy <= 1)); then
          bdy_no_suffix="_$(printf %05d 0)"
        else
          bdy_no_suffix=
        fi
        cat $SCRP_DIR/config.nml.grads_boundary | \
            sed -e "s#--DIR--/bdyatm#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_atm_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                -e "s#--DIR--/bdysfc#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_sfc_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                -e "s#--DIR--/bdyland#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_lnd_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
            > $CONFIG_DIR/${conf_file}

        if ((stage_config == 1)); then
          echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
        fi
      fi # [ BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1) ]
    done # [ m in $(seq $mtot) ]

  fi # [ BDY_FORMAT != 0 ]

  #-----------------------------------------------------------------------------
  # scale (launcher)
  #-----------------------------------------------------------------------------

  config_file_scale_launcher cycle scale-rm_ens "<member>/run" $mtot

  #-----------------------------------------------------------------------------
  # scale (each member)
  #-----------------------------------------------------------------------------

  for m in $(seq $mtot); do
    if ((BDY_ENS == 1)); then
      mem_bdy=${name_m[$m]}
    else
      mem_bdy='mean'
    fi
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

      if ((loop == 1 )); then
        RESTART_IN_POSTFIX_TIMELABEL_TF=".false."
      else
        RESTART_IN_POSTFIX_TIMELABEL_TF=".true."
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
      conf="$(cat $conf_file_src | \
          sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale.${name_m[$m]}.d${dfmt}.LOG_${time}\"," \
              -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
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
              -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = ${RESTART_IN_POSTFIX_TIMELABEL_TF}," \
              -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = .true.," \
              -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${name_m[$m]}/anal.d${dfmt}\"," \
              -e "/!--RESTART_OUT_POSTFIX_TIMELABEL--/a RESTART_OUT_POSTFIX_TIMELABEL = .true.," \
              -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo/topo.d${dfmt}\"," \
              -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse/landuse.d${dfmt}\"," \
              -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${name_m[$m]}/hist.d${dfmt}_$(datetime_scale $time)\"," \
              -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${CYCLEFOUT}.D0," \
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
        conf="$(echo "$conf" | \
            sed -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${mem_bdy}/bdy_$(datetime_scale $time)\"," \
                -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = ${bdy_start_time:0:4}, ${bdy_start_time:4:2}, ${bdy_start_time:6:2}, ${bdy_start_time:8:2}, ${bdy_start_time:10:2}, ${bdy_start_time:12:2}," \
                -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0,")"
      fi
      if [ ! -e "$SCRP_DIR/config.nml.scale_user" ]; then
        if ((OCEAN_INPUT == 1)); then
          if ((OCEAN_FORMAT == 99)); then
            conf="$(echo "$conf" | \
                sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${mem_bdy}/init.d${dfmt}_$(datetime_scale $time)\",")"
          fi
        fi
        if ((LAND_INPUT == 1)); then
          if ((LAND_FORMAT == 99)); then
            conf="$(echo "$conf" | \
                sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${mem_bdy}/init.d${dfmt}_$(datetime_scale $time)\",")"
          fi
        fi
      fi
      mkdir -p $CONFIG_DIR/${name_m[$m]}
      conf_file="${name_m[$m]}/run.d${dfmt}_${time}.conf"
      echo "  $conf_file"
      echo "$conf" > $CONFIG_DIR/${conf_file}

      if [ -e "$SCRP_DIR/config.nml.scale_user" ]; then
        conf="$(cat $SCRP_DIR/config.nml.scale_user)"
        if ((OCEAN_INPUT == 1)); then
          if ((OCEAN_FORMAT == 99)); then
            conf="$(echo "$conf" | \
                sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${mem_bdy}/init.d${dfmt}_$(datetime_scale $time)\",")"
          fi
        fi
        if ((LAND_INPUT == 1)); then
          if ((LAND_FORMAT == 99)); then
            conf="$(echo "$conf" | \
                sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"${mem_bdy}/init.d${dfmt}_$(datetime_scale $time)\",")"
          fi
        fi
        echo "$conf" >> $CONFIG_DIR/${conf_file}
      fi

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
            -e "/!--LETKF_TOPO_IN_BASENAME--/a LETKF_TOPO_IN_BASENAME = \"topo/topo.d${dfmt}\"," \
            -e "/!--INFL_ADD_IN_BASENAME--/a INFL_ADD_IN_BASENAME = \"<member>/addi.d${dfmt}\"," \
            -e "/!--RELAX_SPREAD_OUT--/a RELAX_SPREAD_OUT = ${RTPS_INFL_OUT_TF}," \
            -e "/!--RELAX_SPREAD_OUT_BASENAME--/a RELAX_SPREAD_OUT_BASENAME = \"rtpsinfl.d${dfmt}_$(datetime_scale $atime).nc\"," \
            -e "/!--NOBS_OUT--/a NOBS_OUT = ${NOBS_OUT_TF}," \
            -e "/!--NOBS_OUT_BASENAME--/a NOBS_OUT_BASENAME = \"nobs.d${dfmt}_$(datetime_scale $atime).nc\"," \
        >> $CONFIG_DIR/${conf_file}

    # Most of these parameters are not important for letkf
    cat $conf_file_src2 | \
        sed -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
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
