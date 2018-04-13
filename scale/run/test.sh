#!/bin/bash
#===============================================================================

function print_summary () {
#  rm -f config.main
#  rm -f config.${SCPNAME[$i]}
#  rm -f config.nml.letkf
#  rm -f config.nml.obsope
#  rm -f config.nml.scale
#  rm -f config.nml.scale_pp
#  rm -f config.nml.scale_init
#  rm -f config.nml.scale_user
#  rm -f config.nml.grads_boundary

#  rm -f ${SCPNAME[$i]}_job.sh
#  rm -f ${jobname}.?${jobid}

  rm -f test.letkf
  rm -f test.letkf.ref
  rm -f test.log
  rm -f test.time

  if ((NTEST > 0)); then
    echo
    echo "TEST SUMMARY:"
    echo "========================================================================================================"
    echo "Config                               Job    Dom  Preset    MPI_type  Mode     -> Time(s)  Result"
    echo "--------------------------------------------------------------------------------------------------------"
    for j in $(seq $NTEST); do
      printf "%-37s%-7s%-5s%-10s%-10s%-9s-> %-9s%-15s\n" "${CONFIG[$j]}" "${SCPNAME[$j]}" "${DOMNUM[$j]}" "${PRESET[$j]}" "${MPI_TYPE[$j]}" "${CONF_MODE[$j]}" "${Rtime[$j]}" "${Rstatus[$j]}"
    done
    echo "========================================================================================================"

    pass=0
    fail=0
    for i in $(seq $NTEST); do
      if [ "${Rstatus[$i]}" = 'Correct' ]; then
        pass=$((pass+1))
      else
        fail=$((fail+1))
      fi
    done
    echo "PASS: $pass"
    echo "FAIL: $fail"
  fi
}
trap print_summary EXIT

#-------------------------------------------------------------------------------

function print_result () {
  echo
  echo "TEST $i RESULT: ${Rstatus[$i]}"
  echo
}

#-------------------------------------------------------------------------------

NTEST=0

cd "$(dirname "$0")"

#-------------------------------------------------------------------------------

if [ "$1" = '-f' ]; then
  if (($# < 2)); then
    echo "$0: Insufficient arguments" >&2
    exit 1
  fi
  shift
  FILENAME="$1"
  NTEST=0
  while read line; do
    NTEST=$((NTEST+1))
    CONFIG[$NTEST]="$(echo $line | cut -d ' ' -f1)"
    SCPNAME[$NTEST]="$(echo $line | cut -d ' ' -f2)"
    DOMNUM[$NTEST]="$(echo $line | cut -d ' ' -f3)"
    PRESET[$NTEST]="$(echo $line | cut -d ' ' -f4)"
    MPI_TYPE[$NTEST]="$(echo $line | cut -d ' ' -f5)"
    CONF_MODE[$NTEST]="$(echo $line | cut -d ' ' -f6)"
    WTIME_L[$NTEST]="$(echo $line | cut -d ' ' -f7)"
    if [ "${DOMNUM[$NTEST]}" = '-' ]; then
      DOMNUM[$NTEST]=1
    fi
    if [ "${CONF_MODE[$NTEST]}" = '-' ]; then
      CONF_MODE[$NTEST]='dynamic'
    fi
  done < "$FILENAME"
else
  if (($# < 5)); then
    echo "$0: Insufficient arguments" >&2
    exit 1
  fi
  NTEST=1
  CONFIG[1]="$1"; shift
  SCPNAME[1]="$1"; shift
  DOMNUM[1]="$1"; shift
  PRESET[1]="$1"; shift
  MPI_TYPE[1]="$1"; shift
  CONF_MODE[1]="$1"; shift
  WTIME_L[1]="$1"
  if [ "${DOMNUM[1]}" = '-' ]; then
    DOMNUM[1]=1
  fi
  if [ "${CONF_MODE[1]}" = '-' ]; then
    CONF_MODE[1]='dynamic'
  fi
fi

for i in $(seq $NTEST); do
  Rstatus[$i]='Not run'
  Rtime[$i]='0'
done

#-------------------------------------------------------------------------------

for i in $(seq $NTEST); do

  for d in $(seq ${DOMNUM[$i]}); do
    CONFIG_D[$d]="$(echo ${CONFIG[$i]} | cut -d '|' -f${d})"
    if ((d > 1)); then
      CONFIG_D[$d]=${CONFIG_D[$d]:-${CONFIG_D[1]}}
    fi
  done

  echo
  echo "********************************************************************************************"
  echo "* TEST $i: ${CONFIG_D[1]} ${SCPNAME[$i]} ${DOMNUM[$i]} ${PRESET[$i]} ${MPI_TYPE[$i]} ${CONF_MODE[$i]} ${WTIME_L[$i]}"
  echo "********************************************************************************************"
  echo

  if [ "${PRESET[$i]}" = 'K' ] || [ "${PRESET[$i]}" = 'K_rankdir' ]; then
    config_suffix='K'
    script_suffix='_K'
  elif [ "${PRESET[$i]}" = 'K_micro' ]; then
    config_suffix='K'
    script_suffix='_K_micro'
  elif [ "${PRESET[$i]}" = 'OFP' ]; then
    config_suffix='ofp'
    script_suffix='_ofp'
  elif [ "${PRESET[$i]}" = 'Linux_torque' ]; then
    config_suffix='hakushu'
    script_suffix='_torque'
  else
    echo "[Error] Unsupported \$PRESET: ${PRESET[$i]}" >&2
    print_result
    continue
  fi

#-------------------------------------------------------------------------------

  rm -f config.main
  rm -f config.${SCPNAME[$i]}
  rm -f config.nml.letkf
  rm -f config.nml.obsope
  rm -f config.nml.scale*
  rm -f config.nml.grads_boundary

  cat config/${CONFIG_D[1]}/config.main.${config_suffix} | \
      sed -e "s/<DOMNUM>/${DOMNUM[$i]}/g" \
          -e "s/<PRESET>/${PRESET[$i]}/g" \
          -e "s/<MPI_TYPE>/${MPI_TYPE[$i]}/g" \
      > config.main

  cat config/${CONFIG_D[1]}/config.${SCPNAME[$i]} | \
      sed -e "s/<WTIME_L>/${WTIME_L[$i]}/g" \
          -e "s/<CONF_MODE>/${CONF_MODE[$i]}/g" \
      > config.${SCPNAME[$i]}

  ln -fs config/${CONFIG_D[1]}/config.nml.letkf .
  ln -fs config/${CONFIG_D[1]}/config.nml.scale .
  ln -fs config/${CONFIG_D[1]}/config.nml.scale_pp .
  ln -fs config/${CONFIG_D[1]}/config.nml.scale_init .
  if [ -e "config/${CONFIG_D[1]}/config.nml.scale_user" ]; then
    ln -fs config/${CONFIG_D[1]}/config.nml.scale_user .
  fi
  if [ -e "config/${CONFIG_D[1]}/config.nml.obsope" ]; then
    ln -fs config/${CONFIG_D[1]}/config.nml.obsope .
  fi
  if [ -e "config/${CONFIG_D[1]}/config.nml.grads_boundary" ]; then
    ln -fs config/${CONFIG_D[1]}/config.nml.grads_boundary .
  fi

  for d in $(seq 2 ${DOMNUM[$i]}); do
    if [ -s "config/${CONFIG_D[1]}/config.nml.letkf.d${d}" ]; then
      ln -fs config/${CONFIG_D[1]}/config.nml.letkf.d${d} config.nml.letkf.d${d}
    else
      ln -fs config/${CONFIG_D[$d]}/config.nml.letkf config.nml.letkf.d${d}
    fi
    if [ -s "config/${CONFIG_D[1]}/config.nml.scale.d${d}" ]; then
      ln -fs config/${CONFIG_D[1]}/config.nml.scale.d${d} config.nml.scale.d${d}
    else
      ln -fs config/${CONFIG_D[$d]}/config.nml.scale config.nml.scale.d${d}
    fi
    if [ -s "config/${CONFIG_D[1]}/config.nml.scale_init.d${d}" ]; then
      ln -fs config/${CONFIG_D[1]}/config.nml.scale_init.d${d} config.nml.scale_init.d${d}
    else
      ln -fs config/${CONFIG_D[$d]}/config.nml.scale_init config.nml.scale_init.d${d}
    fi
  done

  . config.main || exit $?
  . config.${SCPNAME[$i]} || exit $?
  . src/func_datetime.sh || exit $?

#-------------------------------------------------------------------------------

  { time -p ./${SCPNAME[$i]}${script_suffix}.sh 2>&1 | tee test.log ; } 2> test.time
  rc=${PIPESTATUS[0]}

  Rstatus[$i]='Exit/Error'

#-------------------------------------------------------------------------------

  if [ "${PRESET[$i]}" = 'K' ] || [ "${PRESET[$i]}" = 'K_rankdir' ] || [ "${PRESET[$i]}" = 'K_micro' ]; then
    jobname="${SCPNAME[$i]}_${SYSNAME}"
    jobid=$(grep 'pjsub Job' test.log | cut -d ' ' -f6)
    logdir="$OUTDIR/exp/${jobid}_${SCPNAME[$i]}_${STIME}"
    stdout="$logdir/job.o"
    stderr="$logdir/job.e"
    jobinfo="$logdir/job.i"
    letkflogfmt="NOUT.%d"
  elif [ "${PRESET[$i]}" = 'OFP' ]; then
    jobname="${SCPNAME[$i]}_${SYSNAME}"
    jobid=$(grep 'pjsub Job' test.log | cut -d ' ' -f6)
    logdir="$OUTDIR/exp/${jobid}_${SCPNAME[$i]}_${STIME}"
    stdout="$logdir/job.o"
    stderr="$logdir/job.e"
    letkflogfmt="NOUT-${PROCESS_FMT}"
  elif [ "${PRESET[$i]}" = 'Linux_torque' ]; then
    jobname="${SCPNAME[$i]}_job.sh"
    jobid=$(grep 'qsub Job' test.log | cut -d ' ' -f3)
    logdir="$OUTDIR/exp/${jobid}_${SCPNAME[$i]}_${STIME}"
    stdout="$logdir/job.o"
    stderr="$logdir/job.e"
    letkflogfmt="NOUT-${PROCESS_FMT}"
#  elif [ "${PRESET[$i]}" = 'Linux' ]; then
#    ...
#    Rtime[$i]=$(grep 'real' test.time | cut -d ' ' -f2)
  fi

  if [ -e "$jobinfo" ]; then
    Rtime[$i]=$(grep 'ELAPSE TIME (USE)' $jobinfo | cut -d '(' -f3 | cut -d ')' -f1)
  else
    Rtime[$i]=$(grep 'real' test.time | cut -d ' ' -f2)
  fi

  if ((rc != 0)); then
    print_result
    continue
  fi

  if [ ! -e "$stdout" ] || [ ! -e "$stderr" ]; then
    print_result
    continue
  fi

  echo
  echo "STANDARD ERROR MESSAGE:"
  echo "========================================"
  cat $stderr
  echo "========================================"

  if [ -z "$(tail -n 1 $stderr | grep "Finish")" ]; then
    print_result
    continue
  fi

  Rstatus[$i]='Done/Incorrect'

  rm -f ${SCPNAME[$i]}_job.sh
  rm -f ${jobname}.?${jobid}

#-------------------------------------------------------------------------------

  if [ "${SCPNAME[$i]}" = 'cycle' ]; then
    eatime=$(datetime $ETIME $LCYCLE s)
    str_search="OBSERVATIONAL DEPARTURE STATISTICS.*GLOBAL\|GLOBAL.*OBSERVATIONAL DEPARTURE STATISTICS"
    str_grep="OBSERVATIONAL DEPARTURE STATISTICS"

    cdiff=0
    for d in $(seq ${DOMNUM[$i]}); do
      sed -n "/${str_search}/,+7p" ${OUTDIR[1]}/${eatime}/log/letkf/$(printf $letkflogfmt ${SCALE_NP_S[$d]}) | grep -v "$str_grep" > test.letkf
      sed -n "/${str_search}/,+7p" ${OUTDIR[1]}/results/${eatime}/log/letkf/$(printf $letkflogfmt ${SCALE_NP_S[$d]}) | grep -v "$str_grep" > test.letkf.ref

      if [ -n "$(diff test.letkf test.letkf.ref)" ]; then
        echo
        echo "REFERENCE RESULT [Domain $d]:"
        echo
        cat test.letkf.ref
        echo
        echo "THIS RESULT [Domain $d]:"
        echo
        cat test.letkf
        echo
        cdiff=1
      fi
    done
    if ((cdiff == 1)); then
      print_result
      continue
    fi
  fi

  Rstatus[$i]='Correct'
  print_result

done # [ i in $(seq $NTEST) ]

#===============================================================================
