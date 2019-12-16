#!/bin/bash
#===============================================================================
#
#  Steps of 'fcst.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=1
stepname[1]='Run SCALE pp'
stepexecdir[1]="$TMPRUN/scale_pp"
stepexecname[1]="scale-rm_pp_ens"

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

CONF_MODE=${1:-$CONF_MODE}; shift
TIME_LIMIT="${1:-$TIME_LIMIT}"

MEMBERS='mdet'

#-------------------------------------------------------------------------------
# assign default values to and standardize the parameters

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

#  if ((MAKEINIT == 1)); then
#    if [ -d "${OUTDIR}/anal" ]; then
#      if [ -n "$(ls ${OUTDIR}/anal 2> /dev/null)" ]; then
#        echo "[Error] $myname: Initial ensemble is to be generated (\$MAKEINIT = 1) at \"${OUTDIR}/anal/\", but existing data are found there;" >&2
#        echo "        Set \$MAKEINIT = 0 or remove \"${OUTDIR}/anal/*\" before running this job." >&2
#        exit 1
#      fi
#    fi
#  fi
fi

#... more detections...

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
${SCRP_DIR}/config.nml.ensmodel|${DAT_SUBDIR}/conf/config.nml.ensmodel
EOF

if [ -e "${SCRP_DIR}/config.nml.grads_boundary" ]; then
  echo "${SCRP_DIR}/config.nml.grads_boundary|${DAT_SUBDIR}/conf/config.nml.grads_boundary" >> ${STAGING_DIR}/${STGINLIST}
fi

#-------------------------------------------------------------------------------
# TMPRUN

cat >> ${STAGING_DIR}/${STGINLIST} << EOF
${ENSMODEL_DIR}/scale-rm_pp_ens|${RUN_SUBDIR}/scale_pp/scale-rm_pp_ens
EOF

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
loop=1
c=1

      # empty directories
      #-------------------

      echo "|${OUT_SUBDIR}/log/fcst_scale_pp/" >> ${STAGING_DIR}/${STGINLIST}

      #-------------------
      # stage-in
      #-------------------

      # topo (bdy_scale)
      #-------------------
      if ((loop == 1 && c == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
        if ((PNETCDF_BDY_SCALE == 1)); then
          pathin="${DATA_TOPO_BDY_SCALE}.nc"
          path="bdytopo/const/topo.nc"
        else
          pathin="${DATA_TOPO_BDY_SCALE}/"
          path="bdytopo/const/"
        fi

        mkdir -p ${TMP}/${DAT_SUBDIR}/${path}
        ln -sf ${pathin}/* ${TMP}/${DAT_SUBDIR}/${path}
#        echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
      fi

      #-------------------
      # stage-out
      #-------------------
#

#      if ((LOG_OPT <= 2)); then
#        if ((LOG_TYPE == 1)); then
#          if ((c == 1)); then
#            path="fcst_scale_pp/${name_m[1]}_pp.conf"
#            echo "${OUTDIR}/log/${path}|${OUT_SUBDIR}/log/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#            path="fcst_scale_pp/${name_m[1]}_LOG${SCALE_SFX_NONC_0}"
#            echo "${OUTDIR}/log/${path}|${OUT_SUBDIR}/log/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#            path="fcst_scale_pp/NOUT.${log_zeros}"
#            echo "${OUTDIR}/log/${path}|${OUT_SUBDIR}/log/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}.1
#          fi
#        else
#          path="fcst_scale_pp/"
#          echo "${OUTDIR}/log/${path}|${OUT_SUBDIR}/log/${path}|${loop}" >> ${STAGING_DIR}/${STGOUTLIST}
#        fi
#      fi

#-------------------
# stage-in
#-------------------

# bdy
#-------------------
if ((BDY_FORMAT >= 1)); then
  if ((BDY_FORMAT == 1)); then
#    if [ -s "$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt" ]; then
#      pathin="$DATA_BDY_SCALE/const/log/latlon_domain_catalogue.txt"
    if [ -s "$DATA_TOPO_BDY_SCALE/../log/latlon_domain_catalogue.txt" ]; then
      pathin="$DATA_TOPO_BDY_SCALE/../log/latlon_domain_catalogue.txt"
      path="bdyorg/latlon_domain_catalogue.txt"
      echo "${pathin}|${DAT_SUBDIR}/${path}" >> ${STAGING_DIR}/${STGINLIST_BDYDATA}
    else
      echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
      echo "        '$DATA_TOPO_BDY_SCALE/../log/latlon_domain_catalogue.txt'" >&2
      exit 1
    fi
  fi
fi # ((BDY_FORMAT >= 1))

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
#  MEMBER_RUN=$rcycle
  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

if (pdrun all $PROC_OPT); then
 echo  " bash $SCRP_DIR/src/pre_scale_pp_node.sh $MYRANK $mem_nodes $mem_np $TMPRUN/scale_pp $MEMBER_RUN $iter fcst"
  bash $SCRP_DIR/src/pre_scale_pp_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale_pp $MEMBER_RUN $iter fcst
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
      if ((DISK_MODE <= 2)); then
        c=$m
      else
        c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
      fi

        bash $SCRP_DIR/src/pre_scale_pp.sh $MYRANK 20191212000000 ${name_m[$m]} \
             $TMPRUN/scale_pp/$(printf $MEMBER_FMT $m) $TMPDAT \
             fcst ${bdytopo} ${bdycatalogue}
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
  MEMBER_RUN=$rcycle
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  if (pdrun $g $PROC_OPT); then
    m=$(((it-1)*parallel_mems+g))
    if ((m >= 1 && m <= MEMBER_RUN)); then
      if ((DISK_MODE <= 2)); then
        c=$m
      else
        c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
      fi

      if [ -n "${stimes[$c]}" ]; then
        bash $SCRP_DIR/src/post_scale_pp.sh $MYRANK ${stimes[$c]} \
             ${name_m[$m]} $TMPRUN/scale_pp/$(printf $MEMBER_FMT $m) $LOG_OPT fcst
      fi
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) end" >&2
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
    if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/log/${time}/fcst_scale_pp" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/log/${time}/fcst_scale_pp ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/log/${time}/fcst_scale_pp ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/log/${time}/fcst_scale_pp
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/log/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/log/${time}/fcst_scale_pp
        fi
      fi
    fi

    if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/log/${time}/fcst_scale_init" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/log/${time}/fcst_scale_init ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/log/${time}/fcst_scale_init ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/log/${time}/fcst_scale_init
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/log/${time}/fcst_scale_init
        fi
      fi
    fi

    if ((LOG_OPT <= 3)) && [ -d "$OUTDIR/log/${time}/fcst_scale" ]; then
      if ((TAR_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= TAR_THREAD)); do
          sleep 1s
        done
        if ((LOG_TYPE == 3)); then
          ( tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/log/${time}/fcst_scale ) &
        elif ((LOG_TYPE == 4)); then
          ( tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/log/${time}/fcst_scale ) &
        fi
      else
        if ((LOG_TYPE == 3)); then
          tar -C $OUTDIR/log/${time} -cf $OUTDIR/log/${time}/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/log/${time}/fcst_scale
        elif ((LOG_TYPE == 4)); then
          tar -C $OUTDIR/log/${time} -czf $OUTDIR/log/${time}/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/log/${time}/fcst_scale
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
