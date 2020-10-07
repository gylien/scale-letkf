#!ubin/bash
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


fmember=0
for iname in $MEMBERS; do
  fmember=$((fmember+1))
  name_m[$fmember]=$iname
  mkdir -p $TMP/$iname
  mkdir -p $TMP/bdy/$iname
done

totalnp=$((PPN*NNODES))
SCALE_NP_TOTAL=0
for d in `seq $DOMNUM`; do
  SCALE_NP_TOTAL=$((SCALE_NP_TOTAL+SCALE_NP[$d]))
done

mkdir -p $TMP/mean
mkdir -p $TMP/mdet


CYCLE=$((fmember*SCALE_NP_TOTAL/totalnp))
if (( CYCLE < 1 )) ; then
  CYCLE=1
fi

repeat_mems=$((fmember*SCALE_NP_TOTAL/totalnp))
nitmax=$((fmember*SCALE_NP_TOTAL/totalnp))
if (( nitmax < 1 )) ; then
  echo "Make sure nitmax " >&2
  exit 1
fi

mkdir -p ${TMPROOT}/topo
mkdir -p ${TMPROOT}/landuse
mkdir -p ${TMPROOT}/dat
mkdir -p ${TMPROOT}/log

#-------------------------------------------------------------------------------
# executable files

cp ${COMMON_DIR}/pdbash ${TMPROOT}/pdbash
cp ${COMMON_DIR}/datetime ${TMPROOT}/datetime
cp ${ENSMODEL_DIR}/scale-rm_pp_ens ${TMPROOT}/scale-rm_pp_ens
cp ${ENSMODEL_DIR}/scale-rm_init_ens ${TMPROOT}/scale-rm_init_ens
cp ${ENSMODEL_DIR}/scale-rm_ens ${TMPROOT}/scale-rm_ens

#-------------------------------------------------------------------------------
# database

cp -r ${SCALEDIR}/scale-rm/test/data/rad ${TMPROOT}/dat/rad
cp -r ${SCALEDIR}/scale-rm/test/data/land ${TMPROOT}/dat/land
cp -r ${SCALEDIR}/scale-rm/test/data/urban ${TMPROOT}/dat/urban

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
          mkdir -p ${TMPROOT}/${name_m[$m]}
          for d in $(seq $DOMNUM); do
            pathin="${INDIR[$d]}/${time}/anal/${name_m[$m]}"
            path="$TMPROOT/${name_m[$m]}_d$(printf $DOMAIN_FMT $d)"
          done
#            for q in $(seq ${SCALE_NP[$d]}); do
#              pathin="${INDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#              path="$TMPROOT/${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              ln -sf $pathin $path
#              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#            done
        done
      fi

#      # topo
#      #-------------------
#      if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
#        if ((DISK_MODE == 3)); then
#          for mm in $(seq $fmember); do
#            m=$(((c-1) * fmember + mm))
#            for d in $(seq $DOMNUM); do
#              for q in $(seq ${SCALE_NP[$d]}); do
#                pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#                path="${TMPROOT}/topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#                ln -sf $pathin $path
##                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#              done
#            done
#          done
#        elif ((c == 1)); then
#          for d in $(seq $DOMNUM); do
#            for q in $(seq ${SCALE_NP[$d]}); do
#              pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#              path="${TMPROOT}/topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#              ln -sf $pathin $path
##              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
#            done
#          done
#        fi
#      fi

#      # topo (bdy_scale)
#      #-------------------
#      if ((loop == 1 && c == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#        pathin="${DATA_TOPO_BDY_SCALE}.nc"
#        path="bdytopo.nc"
#        echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#      fi

#      # landuse
#      #-------------------
#      if ((loop == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
#        if ((DISK_MODE == 3)); then
#          for mm in $(seq $fmember); do
#            m=$(((c-1) * fmember + mm))
#            for d in $(seq $DOMNUM); do
#              for q in $(seq ${SCALE_NP[$d]}); do
#                pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#                path="${TMPROOT}/landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#                ln -sf $pathin $path
##                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#              done
#            done
#          done
#        elif ((c == 1)); then
#          for d in $(seq $DOMNUM); do
#            for q in $(seq ${SCALE_NP[$d]}); do
#              pathin="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#              path="${TMPROOT}/landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#              ln -sf $pathin $path
##              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
#            done
#          done
#        fi
#      fi

      # bdy (prepared)
      #-------------------
      if ((BDY_FORMAT == 0)); then
        if ((BDY_ENS == 0)); then
          if ((DISK_MODE == 3)); then
            for mm in $(seq $fmember); do
              m=$(((c-1) * fmember + mm))
              for q in $(seq ${SCALE_NP[1]}); do
                pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
                path="${TMPROOT}/mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                ln -sf $pathin $path
#                echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
              done
              if ((USE_INIT_FROM_BDY == 1)); then
                for d in $(seq $DOMNUM); do
                  for q in $(seq ${SCALE_NP[$d]}); do
                    pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
                    path="${TMPROOT}/mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                    ln -sf $pathin $path
#                    echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
                  done
                done
              fi
            done
          else
            for q in $(seq ${SCALE_NP[1]}); do
              pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
              path="${TMPROOT}/mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              ln -sf $pathin $path
#              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
            done
            if ((USE_INIT_FROM_BDY == 1)); then
              for d in $(seq $DOMNUM); do
                for q in $(seq ${SCALE_NP[$d]}); do
                  pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${BDY_MEAN}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
                  path="${TMPROOT}/mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                  ln -sf $pathin $path
#                  echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}
                done
              done
            fi
          fi
        elif ((BDY_ENS == 1)); then
          for mm in $(seq $fmember); do
            m=$(((c-1) * fmember + mm))
            for q in $(seq ${SCALE_NP[1]}); do
              pathin="${DATA_BDY_SCALE_PREP[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
              path="${TMPROOT}/${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              ln -sf $pathin $path
#              echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
            done
            if ((USE_INIT_FROM_BDY == 1)); then
              for d in $(seq $DOMNUM); do
                for q in $(seq ${SCALE_NP[$d]}); do
                  pathin="${DATA_BDY_SCALE_PREP[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
                  path="${TMPROOT}/${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
                  ln -sf $pathin $path
#                  echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
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
            for q in $(seq ${SCALE_NP[$d]}); do
              path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
              pathout="${OUTDIR[$d]}/${time}/anal/${name_m[$m]}${CONNECTOR}init$(scale_filename_sfx $((q-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
            done
          done
        done
      fi

      # topo
      #-------------------
      if ((loop == 1 && c == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${SCALE_NP[$d]}); do
            path="topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP[$d]}+q))]}
          done
        done
      fi

      # landuse
      #-------------------
      if ((loop == 1 && c == 1 && LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
        for d in $(seq $DOMNUM); do
          for q in $(seq ${SCALE_NP[$d]}); do
            path="llanduse/anduse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
            pathout="${OUTDIR[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$((${SCALE_NP_S[$d]}+q))]}
            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$((${SCALE_NP[$d]}+q))]}
          done
        done
      fi

#      # bdy
#      #-------------------
#      if ((BDY_FORMAT != 0)); then
#        if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
#          for mm in $(seq $fmember); do
#            m=$(((c-1) * fmember + mm))
#            for q in $(seq ${SCALE_NP[1]}); do
#              path="${name_m[$m]}/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#              pathout="${OUTDIR[1]}/${time}/bdy/${name_m[$m]}${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
##              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+q))]}
#            done
#            if ((USE_INIT_FROM_BDY == 1)); then
#              for d in $(seq $DOMNUM); do
#                for q in $(seq ${SCALE_NP[$d]}); do
#                  path="${name_m[$m]}/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#                  pathout="${OUTDIR[$d]}/${time}/bdy/${name_m[$m]}${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
##                  echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP_S[$d]}+q))]}
#                  echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#                done
#              done
#            fi
#          done
#        elif ((BDYOUT_OPT <= 2)); then
#          for q in $(seq ${SCALE_NP[1]}); do
#            path="mean/bdy_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#            pathout="${OUTDIR[1]}/${time}/bdy/mean${CONNECTOR}boundary$(scale_filename_sfx $((q-1)))"
##            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$q]}
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+q))]}
#          done
#          if ((USE_INIT_FROM_BDY == 1)); then
#            for d in $(seq $DOMNUM); do
#              for q in $(seq ${SCALE_NP[$d]}); do
#                path="mean/init.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#                pathout="${OUTDIR[$d]}/${time}/bdy/mean${CONNECTOR}init_bdy$(scale_filename_sfx $((q-1)))"
##                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP[$d]}+q))]}
#                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((mmean-1)*mem_np+${SCALE_NP[$d]}+q))]}
#              done
#            done
#          fi
#        fi
#      fi
#
#      # fcst
#      #-------------------
#      for mm in $(seq $fmember); do
#        m=$(((c-1) * fmember + mm))
#        for d in $(seq $DOMNUM); do
#          for q in $(seq ${SCALE_NP[$d]}); do
#            path="${name_m[$m]}/hist.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)$(scale_filename_sfx $((q-1)))"
#            pathout="${OUTDIR[$d]}/${time}/fcst/${name_m[$m]}${CONNECTOR}history$(scale_filename_sfx $((q-1)))"
##            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#            if ((OUT_OPT <= 1)); then
#              path="${name_m[$m]}/fcst.d$(printf $DOMAIN_FMT $d)_$(datetime_scale $time)_$(datetime_scale $ftime)$(scale_filename_sfx $((q-1)))"
#              pathout="${OUTDIR[$d]}/${time}/fcst/${name_m[$m]}${CONNECTOR}init_${ftime}$(scale_filename_sfx $((q-1)))"
##              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST_NOLINK}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+q))]}
#            fi
#          done
#        done
#      done
#
#      # log
#      #-------------------
#      if [ "$MPI_TYPE" = 'K' ]; then
#        log_nfmt='.%d'
#      else
#        log_nfmt="-${PROCESS_FMT}"
#      fi
#      if ((LOG_TYPE == 1)); then
#        mlist='1'  ####  the first of the next time??
#        plist='1'
#      else
#        mlist=$(seq $fmember)
#        if ((BDY_ENS == 1)); then
#          mlist_init=$(seq $fmember)
#        elif ((DISK_MODE <= 2)); then # shared run directory: only run one member per cycle
#          mlist_init='1'
#        else # local run directory: run multiple members as needed
#          mlist_init=$(seq $((repeat_mems <= fmember ? $repeat_mems : $fmember)))
#        fi
#        plist=$(seq $totalnp)
#      fi
#
#      if ((BDY_FORMAT != 0 && LOG_OPT <= 2)); then
#        for mm in $mlist_init; do
#          m=$(((c-1) * fmember + mm))
#          for d in $(seq $DOMNUM); do
#            path="log/scale_init.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
#            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale_init/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+1))]}
#          done
#        done
#        if ((c == 1)); then
#          for p in $plist; do
#            if ((nitmax == 1)); then
#              path="log/scale-rm_init_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
#              pathout="${OUTDIR[1]}/${time}/log/fcst_scale_init/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
#            else
#              for it in $(seq $nitmax); do
#                path="log/scale-rm_init_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
#                pathout="${OUTDIR[1]}/${time}/log/fcst_scale_init/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
#                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
#              done
#            fi
#          done
#        fi
#      fi
#      if ((LOG_OPT <= 3)); then
#        for mm in $mlist; do
#          m=$(((c-1) * fmember + mm))
#          for d in $(seq $DOMNUM); do
#            path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).LOG_${time}${SCALE_SFX_NONC_0}"
#            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale/${name_m[$m]}_LOG${SCALE_SFX_NONC_0}"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+1))]}
#            path="log/scale.${name_m[$m]}.d$(printf $DOMAIN_FMT $d).monitor_${time}${SCALE_SFX_NONC_0}"
#            pathout="${OUTDIR[$d]}/${time}/log/fcst_scale/${name_m[$m]}_monitor${SCALE_SFX_NONC_0}"
#            echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${mem2node[$(((m-1)*mem_np+${SCALE_NP[$d]}+1))]}
#          done
#        done
#        if ((c == 1)); then
#          for p in $plist; do
#            if ((nitmax == 1)); then
#              path="log/scale-rm_ens.NOUT_${time}$(printf -- "${log_nfmt}" $((p-1)))"
#              pathout="${OUTDIR[1]}/${time}/log/fcst_scale/NOUT$(printf -- "${log_nfmt}" $((p-1)))"
#              echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
#            else
#              for it in $(seq $nitmax); do
#                path="log/scale-rm_ens.NOUT_${time}_${it}$(printf -- "${log_nfmt}" $((p-1)))"
#                pathout="${OUTDIR[1]}/${time}/log/fcst_scale/NOUT-${it}$(printf -- "${log_nfmt}" $((p-1)))"
#                echo "${pathout}|${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.${proc2node[$p]}
#              done
#            fi
#          done
#        fi
#      fi
#
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
  if (( mem_np_bdy_ < 1 )) && (( BDY_FORMAT < 4 )) ; then
    echo "[Error] $0: Specify DATA_BDY_SCALE_PRC_NUM_X/Y" >&2
    exit 1
  fi
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

if [ "$TOPO_FORMAT" != "prep" ] || [ "$LAND_FORMAT" != "prep" ] ; then
  mkdir -p $OUTDIR/const/topo
  mkdir -p $OUTDIR/const/landuse
  config_file_scale_launcher fcst fcst_scale-rm_pp_ens "f<member>/pp" 1
  OFFLINE_PARENT_BASENAME=

  if ((BDY_FORMAT == 1)); then
    BDYCATALOGUE=${DATA_BDY_SCALE}/const/log/latlon_domain_catalogue.txt
    BDYTOPO=${DATA_BDY_SCALE}/const/topo
  fi

#  if ((BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#    OFFLINE_PARENT_BASENAME="$COPYTOPO"
#  fi

  if [ "$TOPO_FORMAT" != 'prep' ]; then
    CONVERT_TOPO='.true.'
  else
    CONVERT_TOPO='.false.'
  fi
  
  if [ "$LANDUSE_FORMAT" != 'prep' ]; then
    CONVERT_LANDUSE='.true.'
  else
    CONVERT_LANDUSE='.false.'
  fi

  # assume Domain 1

  mkdir -p $OUTDIR/$time/log/fcst_scale_pp

  conf_file_src=$SCRP_DIR/config.nml.scale_pp
  conf="$(cat $conf_file_src | \
           sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$OUTDIR/$time/log/fcst_scale_pp/LOG\"," \
               -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
               -e "/!--TOPOGRAPHY_OUT_BASENAME--/a TOPOGRAPHY_OUT_BASENAME = \"${OUTDIR}/const/topo/topo\"," \
               -e "/!--LANDUSE_OUT_BASENAME--/a LANDUSE_OUT_BASENAME = \"${OUTDIR}/const/landuse/landuse\"," \
               -e "/!--CONVERT_TOPO--/a CONVERT_TOPO = $CONVERT_TOPO," \
               -e "/!--CONVERT_LANDUSE--/a CONVERT_LANDUSE = $CONVERT_LANDUSE," \
               -e "/!--CNVTOPO_name--/a CNVTOPO_name = \"$TOPO_FORMAT\"," \
               -e "/!--GTOPO30_IN_DIR--/a GTOPO30_IN_DIR = \"${SCALE_DB}/topo/GTOPO30/Products\"," \
               -e "/!--DEM50M_IN_DIR--/a DEM50M_IN_DIR = \"${SCALE_DB}/topo/DEM50M/Products\"," \
               -e "/!--CNVLANDUSE_name--/a CNVLANDUSE_name = '$LANDUSE_FORMAT'," \
               -e "/!--GLCCv2_IN_DIR--/a GLCCv2_IN_DIR = \"${SCALE_DB}/landuse/GLCCv2/Products\"," \
               -e "/!--LU100M_IN_DIR--/a LU100M_IN_DIR = \"${SCALE_DB}/landuse/LU100M/Products\"," \
               -e "/!--COPYTOPO_IN_BASENAME--/a COPYTOPO_IN_BASENAME = \"${BDYTOPO}\"," \
               -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${BDYCATALOGUE}\"," \
               -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${OFFLINE_PARENT_BASENAME}\"," \
               -e "/!--OFFLINE_PARENT_PRC_NUM_X--/a OFFLINE_PARENT_PRC_NUM_X = ${DATA_BDY_SCALE_PRC_NUM_X}," \
               -e "/!--OFFLINE_PARENT_PRC_NUM_Y--/a OFFLINE_PARENT_PRC_NUM_Y = ${DATA_BDY_SCALE_PRC_NUM_Y}," \
          )"
   mkdir -p $TMP/f$(printf $MEMBER_FMT 1)
   conf_file="$TMP/f$(printf $MEMBER_FMT 1)/pp.d01_${STIME}.conf"
   echo "$conf" > ${conf_file}

#   for q in $(seq ${SCALE_NP[1]}); do
#      pathin="${DATA_TOPO[$d]}/const/${CONNECTOR_TOPO}topo$(scale_filename_sfx $((q-1)))"
#      path="${TMPROOT}/topo/topo.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#      ln -sf $pathin $path
#
#      pathin2="${DATA_LANDUSE[$d]}/const/${CONNECTOR_LANDUSE}landuse$(scale_filename_sfx $((q-1)))"
#      path2="${TMPROOT}/landuse/landuse.d$(printf $DOMAIN_FMT $d)$(scale_filename_sfx $((q-1)))"
#      ln -sf $pathin2 $path2
#   done
fi


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
                    path="bdy/${name_m[$m]}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
                    #echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                    ln -sf  $pathin $TMP/$path
                  done
                done
              else
                for q in $(seq $mem_np_bdy_); do
                  pathin="${DATA_BDY_SCALE}/${time_bdy}/${BDY_SCALE_DIR}/${BDY_MEAN}${CONNECTOR}history$(scale_filename_bdy_sfx $((q-1)))"
                  path="bdy/mean/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))$(scale_filename_bdy_sfx $((q-1)))"
                  #echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                  ln -sf  $pathin $TMP/$path
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
                    path="bdy/${name_m[$m]}/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
                    #echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                    ln -sf  $pathin $TMP/$path
                  done
                done
              else
                for ifile in $(seq $filenum); do
                  if ((BDY_ROTATING == 1)); then
                    pathin="${data_bdy_i}/${time}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                  else
                    pathin="${data_bdy_i}/${BDY_MEAN}/${filename_prefix[$ifile]}${time_bdy}${filename_suffix[$ifile]}"
                  fi
                  path="bdy/mean/bdyorg_${filenamein_prefix[$ifile]}$(datetime_scale $time_bdy_start_prev)_$(printf %05d $((ibdy-1)))${filenamein_suffix[$ifile]}"
                  #echo "${pathin}|${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
                  ln -sf  $pathin $TMP/$path
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
            LATLON_CATALOGUE_FNAME="${DATA_BDY_SCALE}/const/log/latlon_domain_catalogue.txt"
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
            BASENAME_ORG="${TMPROOT_BDYDATA}/bdy/${mem_bdy}/gradsbdy.conf"
          else
            if ((nbdy <= 1)); then
              bdy_no_suffix="_$(printf %05d 0)"
            else
              bdy_no_suffix=
            fi
            BASENAME_ORG="${TMPROOT_BDYDATA}/bdy/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}"
          fi

          for d in $(seq $DOMNUM); do
            dfmt=$(printf $DOMAIN_FMT $d)

            if ((d == 1)); then
              conf_file_src=$SCRP_DIR/config.nml.scale_init
            else
              conf_file_src=$SCRP_DIR/config.nml.scale_init.d$d
            fi

            mkdir -p $OUTDIR/$time/anal/${mem_bdy}
            mkdir -p $OUTDIR/$time/bdy/${mem_bdy}
            if (( MAKEINIT == 1 )); then
              RESTART_OUT_BASENAME=${OUTDIR[$d]}/$time/anal/${mem_bdy}/init
            else
              RESTART_OUT_BASENAME=${OUTDIR[$d]}/$time/bdy/${mem_bdy}/init_bdy
            fi

            mkdir -p ${OUTDIR[$d]}/$time/log/fcst_scale_init

            conf="$(cat $conf_file_src | \
                sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"${OUTDIR[$d]}/$time/log/fcst_scale_init/${name_m[$m]}_LOG_${time}\"," \
                    -e "/!--FILE_AGGREGATE--/a FILE_AGGREGATE = ${FILE_AGGREGATE}," \
                    -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${time:0:4}, ${time:4:2}, ${time:6:2}, ${time:8:2}, ${time:10:2}, ${time:12:2}," \
                    -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
                    -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${RESTART_OUT_BASENAME}\"," \
                    -e "/!--RESTART_OUT_POSTFIX_TIMELABEL--/a RESTART_OUT_POSTFIX_TIMELABEL = .true.," \
                    -e "/!--TOPOGRAPHY_IN_BASENAME--/a TOPOGRAPHY_IN_BASENAME = \"${INDIR[$d]}/const/topo/topo\"," \
                    -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${INDIR[$d]}/const/landuse/landuse\"," \
                    -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\",")"
            if ((BDY_FORMAT == 1)); then
              conf="$(echo "$conf" | \
                  sed -e "/!--OFFLINE_PARENT_BASENAME--/a OFFLINE_PARENT_BASENAME = \"${TMPROOT_BDYDATA}/bdy/${mem_bdy}/bdyorg_$(datetime_scale $time_bdy_start_prev)_$(printf %05d 0)\"," \
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
                sed -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"bdy/${mem_bdy}/bdy_$(datetime_scale $time)\"," \
                    -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${nbdy}," \
                    -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${ntsteps}," \
                    -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${ntsteps_skip},")"
            else
              #------ before SCALE 5.2
              conf="$(echo "$conf" | \
                sed -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"bdy/${mem_bdy}/bdy.d${dfmt}_$(datetime_scale $time)\"," \
                    -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${nbdy}," \
                    -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${ntsteps}," \
                    -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${ntsteps_skip},")"
              #------ after SCALE 5.3
              #conf="$(echo "$conf" | \
              #  sed -e "/!--MAKE_BOUNDARY--/a MAKE_BOUNDARY = .false.,")"
              #------
            fi
            mkdir -p $TMP/f$(printf $MEMBER_FMT $m)
            conf_file="$TMP/f$(printf $MEMBER_FMT $m)/init.d${dfmt}_${time_s}.conf"
            #echo "  $conf_file"
            echo "$conf" > ${conf_file}

#            if ((stage_config == 1)); then
#              if ((DISK_MODE == 3)); then
#                for q in $(seq $mem_np_); do
#                  echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#                done
#              else
#                echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
#              fi
#            fi
          done # [ d in $(seq $DOMNUM) ]

          #if ((BDY_FORMAT == 4 && (BDY_ENS == 0 || m == 1))); then
          if ((BDY_FORMAT == 4 )) ; then
            mkdir -p $CONFIG_DIR/${mem_bdy}
            conf_file="$TMP/bdy/${mem_bdy}/gradsbdy.conf"
            #echo "  $conf_file"
            if ((nbdy <= 1)); then
              bdy_no_suffix="_$(printf %05d 0)"
            else
              bdy_no_suffix=
            fi
            cat $SCRP_DIR/config.nml.grads_boundary | \
                sed -e "s#--DIR--/bdyatm#${TMPROOT_BDYDATA}/bdy/${mem_bdy}/bdyorg_atm_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                    -e "s#--DIR--/bdysfc#${TMPROOT_BDYDATA}/bdy/${mem_bdy}/bdyorg_sfc_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                    -e "s#--DIR--/bdyland#${TMPROOT_BDYDATA}/bdy/${mem_bdy}/bdyorg_lnd_$(datetime_scale $time_bdy_start_prev)${bdy_no_suffix}#g" \
                > ${conf_file}

#            if ((stage_config == 1)); then
#              echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
#            fi
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
          mkdir -p ${OUTDIR[$d]}/const/log
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

          mkdir -p $OUTDIR/$time/fcst/${name_m[$m]}
          RESTART_IN_BASENAME="${INDIR[$d]}/$time/anal/${name_m[$m]}/init"

          if (( MAKEINIT == 1 )); then
            RESTART_IN_POSTFIX_TIMELABEL=".true."
          else
            RESTART_IN_POSTFIX_TIMELABEL=".false."
          fi

          if ((d == 1)); then
            conf_file_src=$SCRP_DIR/config.nml.scale
          else
            conf_file_src=$SCRP_DIR/config.nml.scale.d$d
          fi

          mkdir -p ${OUTDIR[$d]}/$time/log/fcst_scale

          conf="$(cat $conf_file_src | \
              sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"${OUTDIR[$d]}/$time/log/fcst_scale/${name_m[$m]}_LOG_${time}\"," \
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
                  -e "/!--RESTART_IN_POSTFIX_TIMELABEL--/a RESTART_IN_POSTFIX_TIMELABEL = ${RESTART_IN_POSTFIX_TIMELABEL}," \
                  -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = ${RESTART_OUTPUT}," \
                  -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${OUTDIR[$d]}/$time/fcst/${name_m[$m]}/init\"," \
                  -e "/!--TOPOGRAPHY_IN_BASENAME--/a TOPOGRAPHY_IN_BASENAME = \"${INDIR[$d]}/const/topo/topo\"," \
                  -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${INDIR[$d]}/const/landuse/landuse\"," \
                  -e "/!--FILE_HISTORY_DEFAULT_BASENAME--/a FILE_HISTORY_DEFAULT_BASENAME = \"${OUTDIR[$d]}/$time/fcst/${name_m[$m]}/history\"," \
                  -e "/!--FILE_HISTORY_DEFAULT_TINTERVAL--/a FILE_HISTORY_DEFAULT_TINTERVAL = ${FCSTOUT}.D0," \
                  -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"${OUTDIR[$d]}/$time/log/fcst_scale/${name_m[$m]}_monitor_${time}\"," \
                  -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/land/param.bucket.conf\"," \
                  -e "/!--DOMAIN_CATALOGUE_FNAME--/a DOMAIN_CATALOGUE_FNAME = \"${OUTDIR[$d]}/const/log/latlon_domain_catalogue.txt\"," \
                  -e "/!--DOMAIN_CATALOGUE_OUTPUT--/a DOMAIN_CATALOGUE_OUTPUT = ${DOMAIN_CATALOGUE_OUTPUT}," \
                  -e "/!--URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME--/a  URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/urban/param.kusaka01.dat\"," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAG.29\"," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/PARAPC.29\"," \
                  -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/VARDATA.RM29\"," \
                  -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPROOT_CONSTDB}/dat/rad/cira.nc\"," \
                  -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPROOT_CONSTDB}/dat/rad/MIPAS\",")"
          if ((d == 1)); then
            conf="$(echo "$conf" | \
                sed -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"bdy/${mem_bdy}/bdy_$(datetime_scale $time)\"," \
                    -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = ${bdy_start_time:0:4}, ${bdy_start_time:4:2}, ${bdy_start_time:6:2}, ${bdy_start_time:8:2}, ${bdy_start_time:10:2}, ${bdy_start_time:12:2}," \
                    -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0,")"
          fi
          if [ ! -e "$SCRP_DIR/config.nml.scale_user" ]; then
            if ((OCEAN_INPUT == 1)); then
              if ((OCEAN_FORMAT == 99)); then
                conf="$(echo "$conf" | \
                    sed -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"$OUTDIR/$time/bdy/${mem_bdy}/init_bdy.d${dfmt}_$(datetime_scale $time)\",")"
              fi
            fi
            if ((LAND_INPUT == 1)); then
              if ((LAND_FORMAT == 99)); then
                conf="$(echo "$conf" | \
                    sed -e "/!--LAND_RESTART_IN_BASENAME--/a LAND_RESTART_IN_BASENAME = \"$OUTDIR/$time/bdy/${mem_bdy}/init_bdy.d${dfmt}_$(datetime_scale $time)\",")"
              fi
            fi
          fi
          mkdir -p $CONFIG_DIR/f$(printf $MEMBER_FMT $m)
          mkdir -p $TMP/f$(printf $MEMBER_FMT $m)
          conf_file="$TMP/f$(printf $MEMBER_FMT $m)/run.d${dfmt}_${time_s}.conf"
          #echo "  $conf_file"
          echo "$conf" > ${conf_file}

#          if ((stage_config == 1)); then
#            if ((DISK_MODE == 3)); then
#              for q in $(seq $mem_np_); do
#                echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}.${mem2node[$(((m-1)*mem_np+q))]}
#              done
#            else
#              echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> ${STAGING_DIR}/${STGINLIST}
#            fi
#          fi
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

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=3
stepname[1]='Run SCALE pp'
stepexecdir[1]="$TMPRUN/scale_pp"
stepexecname[1]="scale-rm_pp_ens"
stepname[2]='Run SCALE init'
stepexecdir[2]="$TMPRUN/scale_init"
stepexecname[2]="scale-rm_init_ens"
stepname[3]='Run ensemble forecasts'
stepexecdir[3]="$TMPRUN/scale"
stepexecname[3]="scale-rm_ens"
#stepname[4]='Run verification'
#stepexecdir[4]="$TMPRUN/verify"
#stepexecname[4]="verify"

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run ensemble forecasts and (optional) verifications.

Configuration files:
  config.main
  config.cycle

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP CONF_MODE TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  MEMBERS     List of forecast members ('mean' for ensemble mean)
               all:     Run all members including ensemble mean (default)
               mems:    Run all members but not including ensemble mean
               '2 4 6': Run members 2, 4, 6
  CYCLE       Number of forecast cycles run in parallel
               (default: 1)
  CYCLE_SKIP  Run forecasts every ? cycles
               (default: 1)
  IF_VERF     Run verification? [Not finished!]
               0: No (default)
               1: Yes
              * to run the verification, a shared disk storing observations
                and reference model analyses needs to be used
  IF_EFSO     Use EFSO forecast length and output interval? [Not finished!]
               0: No (default)
               1: Yes
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
MEMBERS=${1:-$MEMBERS}; shift
CYCLE=${1:-$CYCLE}; shift
CYCLE_SKIP=${1:-$CYCLE_SKIP}; shift
IF_VERF=${1:-$IF_VERF}; shift
IF_EFSO=${1:-$IF_EFSO}; shift
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
if [ -z "$MEMBERS" ] || [ "$MEMBERS" = 'all' ]; then
  if (( DET_RUN > 0 )); then
    MEMBERS="mean mdet $(printf "$MEMBER_FMT " $(seq $MEMBER))"
  else
    MEMBERS="mean $(printf "$MEMBER_FMT " $(seq $MEMBER))"
  fi
elif [ "$MEMBERS" = 'mems' ]; then
  MEMBERS=$(printf "$MEMBER_FMT " $(seq $MEMBER))
else
  tmpstr=''
  for m in $MEMBERS; do
    if [ "$m" = 'mean' ] || [ "$m" = 'mdet' ]; then
      tmpstr="$tmpstr$m "
    else
      tmpstr="$tmpstr$(printf $MEMBER_FMT $((10#$m))) "
      (($? != 0)) && exit 1
    fi
  done
  MEMBERS="$tmpstr"
fi
CYCLE=${CYCLE:-0}
CYCLE_SKIP=${CYCLE_SKIP:-1}
IF_VERF=${IF_VERF:-0}
IF_EFSO=${IF_EFSO:-0}
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

if ((BDY_FORMAT >= 1)); then
  if ((BDYCYCLE_INT % BDYINT != 0)); then
    echo "[Error] \$BDYCYCLE_INT needs to be an exact multiple of \$BDYINT" >&2
    exit 1
  fi
  BDY_STARTFRAME_MAX=$((BDYCYCLE_INT / BDYINT))
  if [ -z "$PARENT_REF_TIME" ]; then
    PARENT_REF_TIME=$STIME
    for bdy_startframe in $(seq $BDY_STARTFRAME_MAX); do
      if ((BDY_FORMAT == 1)); then
        BFILE="$DATA_BDY_SCALE/${PARENT_REF_TIME}/hist/${BDY_MEAN}/history${SCALE_SFX_0}"
      elif ((BDY_FORMAT == 2 && BDY_ROTATING == 1)); then
        BFILE="$DATA_BDY_WRF/${PARENT_REF_TIME}/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}" 
      elif ((BDY_FORMAT == 2 && BDY_ROTATING != 1)); then
        BFILE="$DATA_BDY_WRF/${BDY_MEAN}/wrfout_${PARENT_REF_TIME}"
      elif ((BDY_FORMAT == 4 && BDY_ROTATING == 1)); then
        BFILE="$DATA_BDY_GRADS/${PARENT_REF_TIME}/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd"
      elif ((BDY_FORMAT == 4 && BDY_ROTATING != 1)); then 
        BFILE="$DATA_BDY_GRADS/${BDY_MEAN}/atm_${PARENT_REF_TIME}.grd"
      fi

      if [ -s "$BFILE" ] ; then
        break
      fi 

      if ((bdy_startframe == BDY_STARTFRAME_MAX)); then
        echo "[Error] Cannot find boundary files. "$BFILE >&2
        exit 1
      fi

      PARENT_REF_TIME=$(datetime $PARENT_REF_TIME -${BDYINT} s)
    done
  fi
fi

#-------------------------------------------------------------------------------
}

print_setting () {
#-------------------------------------------------------------------------------

for vname in DIR DOMAIN @INDIR @OUTDIR @DATA_TOPO DATA_TOPO_BDY_SCALE @DATA_LANDUSE DATA_BDY_SCALE \
             DATA_BDY_SCALE_PREP DATA_BDY_WRF DATA_BDY_NICAM OBS OBSNCEP DET_RUN TOPO_FORMAT \
             LANDUSE_FORMAT LANDUSE_UPDATE BDY_FORMAT BDY_ENS BDYINT BDYCYCLE_INT PARENT_REF_TIME \
             ENABLE_PARAM_USER OCEAN_INPUT OCEAN_FORMAT LAND_INPUT LAND_FORMAT OBSNUM WINDOW_S WINDOW_E \
             LCYCLE LTIMESLOT MEMBER NNODES NNODES_APPAR PPN PPN_APPAR THREADS @SCALE_NP \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP CONF_MODE \
             FCSTLEN FCSTOUT MAKEINIT OUT_OPT TOPOOUT_OPT LANDUSEOUT_OPT BDYOUT_OPT \
             LOG_OPT LOG_TYPE; do
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

archive_log () {
#-------------------------------------------------------------------------------

if ((LOG_TYPE >= 3)); then
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do
    if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/${time}/log/fcst_scale_pp" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/${time}/log/fcst_scale_pp ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/${time}/log/fcst_scale_pp ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/${time}/log/fcst_scale_pp
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/${time}/log/fcst_scale_pp
        fi
      fi
    fi

    if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/${time}/log/fcst_scale_init" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/${time}/log/fcst_scale_init ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/${time}/log/fcst_scale_init ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/${time}/log/fcst_scale_init
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/${time}/log/fcst_scale_init
        fi
      fi
    fi

    if ((LOG_OPT <= 3)) && [ -d "$OUTDIR/${time}/log/fcst_scale" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/${time}/log/fcst_scale ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/${time}/log/fcst_scale ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/${time}/log -cf $OUTDIR/${time}/log/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/${time}/log/fcst_scale
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/${time}/log -czf $OUTDIR/${time}/log/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/${time}/log/fcst_scale
        fi
      fi
    fi

    time=$(datetime $time $lcycles s)
  done
  if ((TAR_THREAD > 1)); then
    wait
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
