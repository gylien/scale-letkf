#!/bin/bash
#===============================================================================
#
#  Steps of 'cycle.sh'
#
#===============================================================================

staging_list_simple () {
#-------------------------------------------------------------------------------
# TMPDAT

cat >> $STAGING_DIR/stagein.out << EOF
${ENSMODEL_DIR}/scale-rm_pp_ens|scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|scale-rm_ens
${OBSUTIL_DIR}/obsope|obsope
${LETKF_DIR}/letkf|letkf
${COMMON_DIR}/pdbash|pdbash
${COMMON_DIR}/datetime|datetime
${DATADIR}/rad|dat/rad
${DATADIR}/land|dat/land
|sprd/
|log/
EOF

## H08
#  if [ -e "${RTTOV_COEF}" ] && [ -e "${RTTOV_SCCOEF}" ]; then
#    cat >> $STAGING_DIR/stagein.out << EOF
#${RTTOV_COEF}|dat/rttov/rtcoef_himawari_8_ahi.dat
#${RTTOV_SCCOEF}|dat/rttov/sccldcoef_himawari_8_ahi.dat
#EOF
#  fi

if [ "$TOPO_FORMAT" != 'prep' ]; then
  if ((DISK_MODE_TOPO_LANDUSE_DB == 2)); then
    echo "${DATADIR}/topo/${TOPO_FORMAT}/Products|dat/topo/${TOPO_FORMAT}/Products|s" >> $STAGING_DIR/stagein.out
  else
    echo "${DATADIR}/topo/${TOPO_FORMAT}/Products|dat/topo/${TOPO_FORMAT}/Products" >> $STAGING_DIR/stagein.out
  fi
fi
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  if ((DISK_MODE_TOPO_LANDUSE_DB == 2)); then
    echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products|dat/landuse/${LANDUSE_FORMAT}/Products|s" >> $STAGING_DIR/stagein.out
  else
    echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products|dat/landuse/${LANDUSE_FORMAT}/Products" >> $STAGING_DIR/stagein.out
  fi
fi

time=$(datetime $STIME $LCYCLE s)
while ((time <= $(datetime $ETIME $LCYCLE s))); do
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ] && [ -e ${OBS}/${OBSNAME[$iobs]}_${time}.dat ]; then
      echo "${OBS}/${OBSNAME[$iobs]}_${time}.dat|obs.${OBSNAME[$iobs]}_${time}.dat" >> $STAGING_DIR/stagein.out
    fi
  done
  time=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
# TMPOUT

time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0
while ((time <= ETIME)); do
  loop=$((loop+1))
  if ((ONLINE_STGOUT == 1)); then
    stgoutstep="stageout.loop.${loop}"
  else
    stgoutstep='stageout.out'
  fi

  #-------------------
  # stage-in
  #-------------------

  # anal
  #-------------------
  if ((loop == 1 && MAKEINIT != 1)); then
    for m in $(seq $mtot); do
      pathin="${INDIR}/${time}/anal/${name_m[$m]}.init.nc"
      path="${name_m[$m]}/anal.d01_$(datetime_scale $time).nc"
      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+1))]}
    done
  fi

  # topo
  #-------------------
  if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
    pathin="${DATA_TOPO}/const/topo.nc"
    path="topo.d01.nc"
    echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out
  fi

#    # topo (bdy_scale)
#    #-------------------
#    if ((loop == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#      pathin="${DATA_TOPO_BDY_SCALE}.nc"
#      path="bdytopo.nc"
#      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out
#    fi

  # landuse
  #-------------------
  if ((loop == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
    pathin="${DATA_LANDUSE}/const/landuse.nc"
    path="landuse.d01.nc"
    echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out
  fi

  # bdy (prepared)
  #-------------------
  if ((loop == 1 && BDY_FORMAT == 0)); then
    if ((BDY_ENS == 0)); then
      pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${BDY_MEAN}.boundary.nc"
      path="mean/bdy_$(datetime_scale $time).nc"
      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out
    elif ((BDY_ENS == 1)); then
      for m in $(seq $mtot); do
        pathin="${DATA_BDY_SCALE_PREP}/${time}/bdy/${name_m[$m]}.boundary.nc"
        path="${name_m[$m]}/bdy_$(datetime_scale $time).nc"
        echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+1))]}
      done
    fi
  fi

  #-------------------
  # stage-out
  #-------------------

  # anal
  #-------------------
  if ((loop == 1 && MAKEINIT == 1)); then
    for m in $(seq $mtot); do
      path="${name_m[$m]}/anal.d01_$(datetime_scale $time).nc"
      pathout="${OUTDIR}/${time}/anal/${name_m[$m]}.init.nc"
      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
    done
  fi

#    # topo
#    #-------------------
#    if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#      path="topo.d01.nc"
#      pathout="${OUTDIR}/const/topo.nc"
#      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}
#    fi

#    # landuse
#    #-------------------
#    if ((loop == 1 && LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
#      path="landuse.d01.nc"
#      pathout="${OUTDIR}/const/landuse.nc"
#      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}
#    fi

#    # bdy
#    #-------------------
#    if ((loop == 1 && BDY_FORMAT != 0)); then
#      if ((BDY_ENS == 0)); then
#        path="mean/bdy.nc"
#        pathout="${OUTDIR}/${time}/bdy/mean.boundary.nc"
#        echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}
#      elif ((BDY_ENS == 1)); then
#        for m in $(seq $mtot); do
#          path="${name_m[$m]}/bdy.nc"
#          pathout="${OUTDIR}/${time}/bdy/${name_m[$m]}.boundary.nc"
#          echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
#        done
#      fi
#    fi

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
    path="${name_m[$m]}/hist.d01_$(datetime_scale $time).nc"
    pathout="${OUTDIR}/${time}/hist/${name_m[$m]}.history.nc"
    echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
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
    path="${name_m[$m]}/gues.d01_$(datetime_scale $atime).nc"
    pathout="${OUTDIR}/${atime}/gues/${name_m[$m]}.init.nc"
    echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
    if ((m == mmean && SPRD_OUT == 1)); then
      path="sprd/gues.d01_$(datetime_scale $atime).nc"
      pathout="${OUTDIR}/${atime}/gues/sprd.init.nc"
      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
    fi
  done

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
    path="${name_m[$m]}/anal.d01_$(datetime_scale $atime).nc"
    pathout="${OUTDIR}/${atime}/anal/${name_m[$m]}.init.nc"
    echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
    if ((m == mmean && SPRD_OUT == 1)); then
      path="sprd/anal.d01_$(datetime_scale $atime).nc"
      pathout="${OUTDIR}/${atime}/anal/sprd.init.nc"
      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
    fi
  done

#    # diag
#    #-------------------
#    if ((RTPS_INFL_OUT == 1)); then
#      path="rtpsinfl.d01_$(datetime_scale $atime).nc"
#      pathout="${OUTDIR}/${atime}/diag/rtpsinfl.init.nc"
#      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+1))]}
#    fi
#    if ((NOBS_OUT == 1)); then
#      path="nobs.d01_$(datetime_scale $atime).nc"
#      pathout="${OUTDIR}/${atime}/diag/nobs.init.nc"
#      echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mmean-1)*mem_np+1))]}
#    fi

#    # obsgues
#    #-------------------
#    if ((OBSOUT_OPT <= 2)); then
#      for m in $(seq $mtot); do ###### either $mmean or $mmdet ? ######
#        path="${name_m[$m]}/obsgues.d01_${atime}.dat"
#        pathout="${OUTDIR}/${atime}/obsgues/${name_m[$m]}.obsda.dat"
#        echo "${pathout}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((m-1)*mem_np+1))]}
#      done
#    fi

######    # log
######    #-------------------
######    if [ "$MPI_TYPE" = 'K' ]; then
######      log_zeros='0'
######    else
######      log_zeros='000000'
######    fi

######    if ((LOG_OPT <= 2)); then
######      if ((LOG_TYPE == 1)); then
######        path="${time}/log/scale_pp/0001_pp.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale_pp/0001_LOG.pe000000"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale_pp/NOUT.${log_zeros}"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale_init/0001_init.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale_init/0001_gradsbdy.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale_init/0001_LOG.pe000000"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        if ((BDY_ENS == 1)); then
######          path="${time}/log/scale_init/NOUT-1.${log_zeros}"
######        else
######          path="${time}/log/scale_init/NOUT.${log_zeros}"
######        fi
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######      else
######        path="${time}/log/scale_pp"
######        echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
######        path="${time}/log/scale_init"
######        echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
######      fi
######    fi
######    if ((LOG_OPT <= 3)); then
######      if ((LOG_TYPE == 1)); then
######        path="${time}/log/scale/0001_run.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale/0001_LOG.pe000000"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale/NOUT-1.${log_zeros}"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${time}/log/scale/latlon_domain_catalogue.txt"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######      else
######        path="${time}/log/scale"
######        echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
######      fi
######    fi
######    if ((LOG_OPT <= 4)); then
######      if ((LOG_TYPE == 1)); then
######        path="${atime}/log/obsope/obsope.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${atime}/log/obsope/LOG.pe000000"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${atime}/log/obsope/NOUT.${log_zeros}"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${atime}/log/letkf/letkf.conf"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${atime}/log/letkf/LOG.pe000000"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######        path="${atime}/log/letkf/NOUT.${log_zeros}"
######        echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
######      else
######        path="${atime}/log/obsope"
######        echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
######        path="${atime}/log/letkf"
######        echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
######      fi
######    fi

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
}

#===============================================================================

config_file_list () {
#-------------------------------------------------------------------------------

echo "Generate configration files..."

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

  obstime $time

  #-----------------------------------------------------------------------------
  # scale
  #-----------------------------------------------------------------------------

  conf_file="scale-rm_ens_${time}.conf"
  echo "  $conf_file"
  cat $SCRP_DIR/config.nml.ensmodel | \
      sed -e "/!--MEMBER--/a MEMBER = $MEMBER," \
          -e "/!--MEMBER_RUN--/a MEMBER_RUN = $mtot," \
          -e "/!--CONF_FILES--/a CONF_FILES = \"@@@@/run_${time}.conf\"," \
          -e "/!--NNODES--/a NNODES = $NNODES_APPAR," \
          -e "/!--PPN--/a PPN = $PPN_APPAR," \
          -e "/!--MEM_NODES--/a MEM_NODES = $mem_nodes," \
          -e "/!--MEM_NP--/a MEM_NP = $mem_np," \
      > $CONFIG_DIR/${conf_file}
#  echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> $STAGING_DIR/stagein.out

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

    conf_file="${name_m[$m]}/run_${time}.conf"
    echo "  $conf_file"
    mkdir -p $CONFIG_DIR/${name_m[$m]}
    cat $SCRP_DIR/config.nml.scale | \
        sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"log/scale.${name_m[$m]}.LOG_${time}\"," \
            -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = .true.," \
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
            -e "/!--HISTORY_DEFAULT_BASENAME--/a HISTORY_DEFAULT_BASENAME = \"${name_m[$m]}/hist.d01_$(datetime_scale $time)\"," \
            -e "/!--HISTORY_DEFAULT_TINTERVAL--/a HISTORY_DEFAULT_TINTERVAL = ${CYCLEFOUT}.D0," \
            -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"log/scale.${name_m[$m]}.monitor_${time}\"," \
            -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"dat/land/param.bucket.conf\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"dat/rad/PARAG.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"dat/rad/PARAPC.29\"," \
            -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"dat/rad/VARDATA.RM29\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"dat/rad/cira.nc\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"dat/rad/MIPAS\"," \
            -e "/!--TIME_END_RESTART_OUT--/a TIME_END_RESTART_OUT = .false.," \
            -e "/!--RESTART_OUT_ADDITIONAL_COPIES--/a RESTART_OUT_ADDITIONAL_COPIES = ${RESTART_OUT_ADDITIONAL_COPIES}," \
            -e "/!--RESTART_OUT_ADDITIONAL_BASENAME--/a RESTART_OUT_ADDITIONAL_BASENAME = ${RESTART_OUT_ADDITIONAL_BASENAME}" \
        > $CONFIG_DIR/${conf_file}
#    echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> $STAGING_DIR/stagein.out.${mem2node[$(((m-1)*mem_np+1))]}
  done

  #-----------------------------------------------------------------------------
  # letkf
  #-----------------------------------------------------------------------------

  OBS_IN_NAME_LIST=
  for iobs in $(seq $OBSNUM); do
    if [ "${OBSNAME[$iobs]}" != '' ]; then
      OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'obs.${OBSNAME[$iobs]}_${atime}.dat', "
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
          -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = .true.," \
      > $CONFIG_DIR/${conf_file}
  # Most of these parameters are not important for letkf
  cat $SCRP_DIR/config.nml.scale | \
      sed -e "/!--IO_AGGREGATE--/a IO_AGGREGATE = .true.," \
      >> $CONFIG_DIR/${conf_file}
#  echo "$CONFIG_DIR/${conf_file}|${conf_file}" >> $STAGING_DIR/stagein.out




#      NNP=$(cat ${NODEFILE_DIR}/${nodestr} | wc -l)
#      mpiexec -n $NNP -vcoordfile "${NODEFILE_DIR}/${nodestr}" -of-proc log/${stepexecname[$s]}.NOUT_${conf_time} ./${stepexecname[$s]} ${stepexecname[$s]}_${conf_time}.conf || exit $?

#    fi
#  done

  #-------------------
  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
done

#-------------------------------------------------------------------------------
}

#===============================================================================
