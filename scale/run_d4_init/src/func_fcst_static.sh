#!/bin/bash
#===============================================================================
#
#  Functions for 'fcst' jobs
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

staging_list_common_static fcst

#-------------------------------------------------------------------------------
# time-variant outputs

lcycles=$((LCYCLE * CYCLE_SKIP))
time_s=$STIME
loop=0
while ((time_s <= ETIME)); do
  loop=$((loop+1))

  for c in $(seq $CYCLE); do
    time=$(datetime $time_s $((lcycles * (c-1))) s)
    if ((time <= ETIME)); then
      ftime=$(datetime $time $FCSTLEN s)

      #-------------------
      # stage-in
      #-------------------

      # anal
      #-------------------
      if ((MAKEINIT != 1)); then
        for mm in $(seq $fmember); do
          m=$(((c-1) * fmember + mm))
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
          for mm in $(seq $fmember); do
            m=$(((c-1) * fmember + mm))
            for d in $(seq $DOMNUM); do
              for q in $(seq ${mem_np_[$d]}); do
                pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
                path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              done
            done
          done
        elif ((c == 1)); then
          for d in $(seq $DOMNUM); do
            for q in $(seq ${mem_np_[$d]}); do
              pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
              path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
            done
          done
        fi
      fi

#      # topo (bdy_scale)
#      #-------------------
#      if ((loop == 1 && c == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#        pathin="${DATA_TOPO_BDY_SCALE}.nc"
#        path="bdytopo.nc"
#        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#      fi

      # landuse
      #-------------------
      if ((loop == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
        if ((DISK_MODE == 3)); then
          for mm in $(seq $fmember); do
            m=$(((c-1) * fmember + mm))
            for d in $(seq $DOMNUM); do
              for q in $(seq ${mem_np_[$d]}); do
                pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
                path="landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              done
            done
          done
        elif ((c == 1)); then
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
            for mm in $(seq $fmember); do
              m=$(((c-1) * fmember + mm))
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
          for mm in $(seq $fmember); do
            m=$(((c-1) * fmember + mm))
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

      #-------------------
      # stage-out
      #-------------------

      # anal
      #-------------------
      if ((MAKEINIT == 1)); then
        for mm in $(seq $fmember); do
          m=$(((c-1) * fmember + mm))
          for d in $(seq $DOMNUM); do
            for q in $(seq ${mem_np_[$d]}); do
              path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              pathout="${OUTDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            done
          done
        done
      fi

      # topo
      #-------------------
      if ((loop == 1 && c == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            path="topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
          done
        done
      fi

      # landuse
      #-------------------
      if ((loop == 1 && c == 1 && LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            path="landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
          done
        done
      fi

      # bdy
      #-------------------
      if ((BDY_FORMAT != 0)); then
        if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
          for mm in $(seq $fmember); do
            m=$(((c-1) * fmember + mm))
            for q in $(seq ${mem_np_[1]}); do
              path="${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              pathout="${OUTDIR[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
            done
            if ((USE_INIT_FROM_BDY == 1)); then
              for d in $(seq $DOMNUM); do
                for q in $(seq ${mem_np_[$d]}); do
                  path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                  pathout="${OUTDIR[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#                  echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
                  echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
                done
              done
            fi
          done
        elif ((BDYOUT_OPT <= 2)); then
          for q in $(seq ${mem_np_[1]}); do
            path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[1]}/${time}/bdy/mean${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$q]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+q))]}
          done
          if ((USE_INIT_FROM_BDY == 1)); then
            for d in $(seq $DOMNUM); do
              for q in $(seq ${mem_np_[$d]}); do
                path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                pathout="${OUTDIR[$d]}/${time}/bdy/mean${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
#                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              done
            done
          fi
        fi
      fi

      # fcst
      #-------------------
      for mm in $(seq $fmember); do
        m=$(((c-1) * fmember + mm))
        for d in $(seq $DOMNUM); do
          for q in $(seq ${mem_np_[$d]}); do
            path="${name_m[$m]}/hist.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/${time}/fcst/${name_m[$m]}${CONNECTOR}history$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            if ((OUT_OPT <= 1)); then
              path="${name_m[$m]}/fcst.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)_$(datetime_scale $ftime)$(scale_filename_sfx $((q-1)))"
              pathout="${OUTDIR[$d]}/${time}/fcst/${name_m[$m]}${CONNECTOR}init_${ftime}$(scale_filename_sfx $((q-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
            fi
          done
        done
      done

      # log
      #-------------------
      if [ "$MPI_TYPE" = 'K' ]; then
        log_nfmt='.%d'
      else
        log_nfmt="-${PROCESS_FMT}"
      fi
      if ((LOG_TYPE == 1)); then
        mlist='1'  ####  the first of the next time??
        plist='1'
      else
        mlist=$(seq $fmember)
        if ((BDY_ENS == 1)); then
          mlist_init=$(seq $fmember)
        elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
          mlist_init='1'
        else # local run directory: run multiple members as needed
          mlist_init=$(seq $((repeat_mems <= fmember ? $repeat_mems : $fmember)))
        fi
        plist=$(seq $totalnp)
      fi

      if ((BDY_FORMAT != 0 && LOG_OPT <= 2)); then
        for mm in $mlist_init; do
          m=$(((c-1) * fmember + mm))
          for d in $(seq $DOMNUM); do
            path="log/scale_init.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale_init/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
          done
        done
        if ((c == 1)); then
          for p in $plist; do
            if ((nitmax == 1)); then
              path="log/scale-rm_init_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
              pathout="${OUTDIR[1]}/${time}/log/fcst_scale_init/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
            else
              for it in $(seq $nitmax); do
                path="log/scale-rm_init_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
                pathout="${OUTDIR[1]}/${time}/log/fcst_scale_init/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
              done
            fi
          done
        fi
      fi
      if ((LOG_OPT <= 3)); then
        for mm in $mlist; do
          m=$(((c-1) * fmember + mm))
          for d in $(seq $DOMNUM); do
            path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
            path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).monitor_${time}${SCALE_SFX_NONC_0}"
            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale/${name_m[$m]}_monitor${SCALE_SFX_NONC_0}"
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+1))]}
          done
        done
        if ((c == 1)); then
          for p in $plist; do
            if ((nitmax == 1)); then
              path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
              pathout="${OUTDIR[1]}/${time}/log/fcst_scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
            else
              for it in $(seq $nitmax); do
                path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
                pathout="${OUTDIR[1]}/${time}/log/fcst_scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
              done
            fi
          done
        fi
      fi

      #-------------------
    fi # [ ((time <= ETIME)) ]
  done # [ c in $(seq $CYCLE) ]

  time_s=$(datetime $time_s $((lcycles * CYCLE)) s)
done # [ ((time_s <= ETIME)) ]

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

lcycles=$((LCYCLE * CYCLE_SKIP))
time_s=$STIME
time_bdy_start_prev=0
loop=0
while ((time_s <= ETIME)); do
  loop=$((loop+1))

  rcycle=1
  for c in $(seq 2 $CYCLE); do
    if (($(datetime $time_s $((lcycles * (c-1))) s) <= ETIME)); then
      rcycle=$c  # The "real" number of cycles
    else
      break
    fi
  done

  if ((BDY_FORMAT != 0)); then

    if ((BDY_ENS == 1)); then
      m_run_onecycle=$fmember
    elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
      m_run_onecycle=1
    else # local run directory: run multiple members as needed
      m_run_onecycle=$((repeat_mems <= fmember ? $repeat_mems : $fmember))
    fi

    #---------------------------------------------------------------------------
    # scale_init (launcher)
    #---------------------------------------------------------------------------

    time=$time_s
    config_file_scale_launcher fcst fcst_scale-rm_init_ens "f<member>/init" $((m_run_onecycle*rcycle))

    #---------------------------------------------------------------------------
    # scale_init (each member)
    #---------------------------------------------------------------------------

    for c in $(seq $CYCLE); do
      time=$(datetime $time_s $((lcycles * (c-1))) s)
      if ((time <= ETIME)); then

        bdy_setting $time $FCSTLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

        if ((MAKEINIT == 1 || USE_INIT_FROM_BDY == 1)); then
          RESTART_OUTPUT='.true.'
        else
          RESTART_OUTPUT='.false.'
        fi
        if ((MAKEINIT == 1 && ${bdy_times[1]} != time)); then
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
                for m in $(seq $fmember); do
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
                for m in $(seq $fmember); do
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

        for mm in $(seq $m_run_onecycle); do
          m=$(((c-1) * m_run_onecycle + mm))
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
                    -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo.d${dfmt}\"," \
                    -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse.d${dfmt}\"," \
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
            mkdir -p $CONFIG_DIR/f$(printf $MEMBER_FMT $m)
            conf_file="f$(printf $MEMBER_FMT $m)/init.d${dfmt}_${time_s}.conf"
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
                sed -e "s#--FNAME_ATMOS--#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_atm_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                    -e "s#--FNAME_SFC--#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_sfc_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                    -e "s#--FNAME_LAND--#${TMPROOT_BDYDATA}/${mem_bdy}/bdyorg_lnd_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                > $CONFIG_DIR/${conf_file}

            if ((stage_config == 1)); then
              echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
            fi
          fi # [ BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1) ]
        done # [ mm in $(seq $m_run_onecycle) ]

      fi # [ ((time <= ETIME)) ]
    done # [ c in $(seq $CYCLE) ]

  fi # [ BDY_FORMAT != 0 ]

  #-----------------------------------------------------------------------------
  # scale (launcher)
  #-----------------------------------------------------------------------------

  time=$time_s
  config_file_scale_launcher fcst fcst_scale-rm_ens "f<member>/run" $((fmember*rcycle))

  #-----------------------------------------------------------------------------
  # scale (each member)
  #-----------------------------------------------------------------------------

  if ((OUT_OPT <= 1)); then
    RESTART_OUTPUT='.true.'
  else
    RESTART_OUTPUT='.false.'
  fi

  for c in $(seq $CYCLE); do
    time=$(datetime $time_s $((lcycles * (c-1))) s)
    if ((time <= ETIME)); then

      bdy_setting $time $FCSTLEN $BDYCYCLE_INT "$BDYINT" "$PARENT_REF_TIME" "$BDY_SINGLE_FILE"

      for mm in $(seq $fmember); do
        m=$(((c-1) * fmember + mm))
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
            RESTART_IN_BASENAME="${name_m[$m]}/init.d${dfmt}_$(datetime_scale $time)"
          else
            RESTART_IN_BASENAME="${name_m[$m]}/anal.d${dfmt}_$(datetime_scale $time)"
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
                  -e "/!--TIME_DURATION--/a TIME_DURATION = ${FCSTLEN}.D0," \
                  -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${FCSTLEN}.D0," \
                  -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${FCSTLEN}.D0," \
                  -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${FCSTLEN}.D0," \
                  -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${FCSTLEN}.D0," \
                  -e "/!--ONLINE_DOMAIN_NUM--/a ONLINE_DOMAIN_NUM = ${d}," \
                  -e "/!--ONLINE_IAM_PARENT--/a ONLINE_IAM_PARENT = ${ONLINE_IAM_PARENT}," \
                  -e "/!--ONLINE_IAM_DAUGHTER--/a ONLINE_IAM_DAUGHTER = ${ONLINE_IAM_DAUGHTER}," \
                  -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${RESTART_IN_BASENAME}\"," \
                  -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = .false.," \
                  -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
                  -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${name_m[$m]}/fcst.d${dfmt}_$(datetime_scale $time)\"," \
                  -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"topo.d${dfmt}\"," \
                  -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"landuse.d${dfmt}\"," \
                  -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${name_m[$m]}/hist.d${dfmt}_$(datetime_scale $time)\"," \
                  -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${FCSTOUT}.D0," \
                  -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"log/scale.${name_m[$m]}.d${dfmt}.monitor_${time}\"," \
                  -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\"," \
                  -e "/!--DOMAIN_CATALOGUE_FNAME--/a DOMAIN_CATALOGUE_FNAME = \"latlon_domain_catalogue.d${dfmt}.txt\"," \
                  -e "/!--DOMAIN_CATALOGUE_OUTPUT--/a DOMAIN_CATALOGUE_OUTPUT = ${DOMAIN_CATALOGUE_OUTPUT}," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAG.29\"," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAPC.29\"," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/VARDATA.RM29\"," \
                  -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/cira.nc\"," \
                  -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPROOT_CONSTDB}/dat/rad/MIPAS\",")"
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
          mkdir -p $CONFIG_DIR/f$(printf $MEMBER_FMT $m)
          conf_file="f$(printf $MEMBER_FMT $m)/run.d${dfmt}_${time_s}.conf"
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
      done # [ mm in $(seq $fmember) ]

    fi # [ ((time <= ETIME)) ]
  done # [ c in $(seq $CYCLE) ]

  #-------------------
  time_s=$(datetime $time_s $((lcycles * CYCLE)) s)
done # [ ((time_s <= ETIME)) ]

echo

#-------------------------------------------------------------------------------
}

#===============================================================================
