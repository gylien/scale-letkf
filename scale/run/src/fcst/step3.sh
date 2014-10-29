#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#  August 2014, Modified from GFS-LETKF, Guo-Yuan Lien
#
#===============================================================================

cd "$(dirname "$0")"
if [ ! -f configure.sh ]; then
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi

. configure.sh
. distribute.sh
#. stageinout.sh
. datetime.sh
. util.sh

nsteps=4
declare -a stepname
stepname[1]='Initialization (stage in)'
stepname[2]='Run ensemble forecasts'
stepname[3]='Run verification'
stepname[4]='Finalization (stage out)'

#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[fcst.sh] Run ensemble forecasts.
          *use settings in 'configure.sh'

Usage: $0 STIME [MEMBERS] [CYCLE] [CYCLE_SKIP] [IF_VERF] [IF_EFSO] [STEP]

  STIME       Start time of the cycle
  MEMBERS     List of forecast members ('mean' for ensemble mean)
              'all':   run all members including ensemble mean
              'mems':  run all members but not including ensemble mean
              '2 4 6': run members 2, 4, 6
              (default: 'all')
  CYCLE       Number of forecast cycles run in parallel
              (default: 1)
  CYCLE_SKIP  Run forecasts every ? cycles
              (default: 1)

  ... belows are not finished yet ...

  IF_VERF     Run verification?
              0: Do not run verification
              1: Run verification
              (default: 0)
              * to run the verification, a shared disk storing observations
                and reference model analyses needs to be used
  IF_EFSO     Use EFSO forecast length and output interval?
              0: No
              1: Yes
              (default: 0)

  ......

  STEP        Steps of the script to be executed
EOF
for s in $(seq $nsteps); do
  echo "              ${s}:     ${stepname[$s]}" 1>&2
done
  cat 1>&2 << EOF
              'all': Run all steps
              '2-3': Run steps 2 to 3
              '2-' : Run steps after 2
              (default: all)

EOF
  exit 1
fi

#-------------------------------------------------------------------------------

MEMBERS="${2:-all}"
CYCLE=${3:-1}
CYCLE_SKIP=${4:-1}
IF_VERF=${5:-0}
IF_EFSO=${6:-0}
STEP=${7:-all}

for c in $(seq $CYCLE); do
  if [ "$c" -eq 1 ]; then
    STIME[$c]=$(datetime $1)
  else
    STIME[$c]=$(datetime ${STIME[$((c-1))]} $((LCYCLE*CYCLE_SKIP)) s)
  fi
#  Syyyymmddhh[$c]=${STIME[$c]:0:10}
  stimef[$c]="$(datetime_fmt ${STIME[$c]})"
done

#if [ "$IF_EFSO" = '1' ]; then
#  FCSTLEN=$EFSOFLEN
#  FCSTOUT=$EFSOFOUT
#fi

#-------------------------------------------------------------------------------

mkdir -p $TMP1
tmpnode="$TMP1/node"
tmpstagein="$TMP1/stagein"
tmpstageout="$TMP1/stageout"

# Files that only a copy is needed
tmpprog="$STMP1/prog"
tmpdata="$STMP1/data"

# Files that are different in each node
if ((LTMPUSE == 0)); then
  tmpscale="$STMP1/scale"
  tmpverify="$STMP1/verify"
  tmpout="$STMP1/out"
  if ((STMPHOST == 0)); then
    ln -s $OUTDIR $tmpout
  fi
elif ((LTMPUSE == 1)); then
  tmpscale="$LTMP2/scale"
  tmpverify="$LTMP2/verify"
  tmpout="$LTMP2/out"
#  if [ "$IF_VERF" = '1' ]; then
#    echo "[Warning] $0: \$OBS, \$ANLGRDP, \$ANLGRDP2 need to be in a shared disk that can be seen from all computing nodes." 1>&2
#  fi
else
  echo "[Error] $0: Unsupported \$LTMPUSE setting." 1>&2
  exit 1
fi

echo "[$(datetime_now)] Start fcst.sh $@" 1>&2

#-------------------------------------------------------------------------------

if [ "$MEMBERS" = 'all' ] || [ "$MEMBERS" = 'mems' ]; then
  if [ "$MEMBERS" = 'all' ]; then
    MEMBERS="mean "
  else
    MEMBERS=''
  fi
  for m in $(seq $MEMBER); do
    MEMBERS="${MEMBERS}$(printf $MEMBER_FMT $m) "
  done
else
  memberstmp=''
  for m in $MEMBERS; do
    if [ "$m" = 'mean' ] || [ "$m" = 'sprd' ]; then
      memberstmp="$memberstmp$m "
    else
      memberstmp="$memberstmp$(printf $MEMBER_FMT $((10#$m))) "
      (($? != 0)) && exit 1
    fi
  done
  MEMBERS="$memberstmp"
fi

#-------------------------------------------------------------------------------

declare -a run_step
parse_steps $nsteps "$STEP"

#for i in $(seq $nsteps); do
#  echo "$i ${run_step[$i]}"
#done

#===============================================================================
# Print job information

echo
echo " +----------------------------------------------------------------+"
echo " |                        SCALE-Forecasts                         |"
echo " +----------------------------------------------------------------+"
for s in $(seq $nsteps); do
  if (("${run_step[$s]}" == 1)); then
    printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
  fi
done
echo " +----------------------------------------------------------------+"
echo
echo "  Number of cycles:         $CYCLE"
echo "  Forecast start time:"
for c in $(seq $CYCLE); do
  printf "    Cycle %-5s %s\n" "$c:" "${stimef[$c]}"
done
echo
echo "  Forecast length:          $FCSTLEN s"
echo "  Output interval:          $FCSTOUT s"
echo

declare -a procs
declare -a mem2proc
declare -a node
declare -a name_m
declare -a node_m
distribute_fcst "$MEMBERS" $CYCLE $MTYPE 1 machinefile

#echo
#echo "  # forecast members:   $fmember"

echo
echo "===================================================================="




exit





#===============================================================================
# 1. Copy necessary files
p=1
step_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstagein
rm -f $tmpstagein/*

# Fix files
cat > $tmpstagein/run.2 << EOF
ck|${FIXGLOBAL}|fix
EOF

# Executable files
cat > $tmpstagein/run.1 << EOF
ck|${EXECGLOBAL}/global_fcst|prog/global_fcst
ck|${EXECGLOBAL}/global_sfchdr|prog/global_sfchdr
ck|${EXECGLOBAL}/global_sighdr|prog/global_sighdr
ck|${DIR}/ssio/ss2grd|prog/ss2grd
ck|${DIR}/ssio/ss2grdp|prog/ss2grdp
ck|${DIR}/verify/verify|prog/verify
EOF

if [ "$SHAREDISK" = '0' ]; then
# Forecast initial analyses
  for c in `seq $CYCLES`; do
    for m in `seq $fmember`; do
      mt=$(((c-1) * fmember + m))
      echo "ne|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig" >> $tmpstagein/out.${node_m[$mt]}
      echo "ne|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc" >> $tmpstagein/out.${node_m[$mt]}
    done
  done
fi

stagein $ltmprun1 $ltmprun2 $ltmpout 1  # clean stagein

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 2. Ensemble forecasts
p=2
step_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/run_gfs.sh .
for c in `seq $CYCLES`; do
  cf=`printf '%04d' $c`
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    $MPIBIN/mpiexec -host ${node_m[$mt]} bash run_gfs.sh \
                    $ltmpout/anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig \
                    $ltmpout/anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc \
                    $FCSTLEN $FCSTOUT $ltmpgfs/${cf}_${name_m[$mt]} \
                    $ltmpprog $ltmpfix &
    sleep 0.05s
  done
done
wait

#-------------------------------------------------------------------------------

echo
ppnl=$((ppn*mem_nodes_gfs))
pcount=0
for c in `seq $CYCLES`; do
  if [ "$LOG_OPT" -le 2 ]; then
    mkdir -p $OUTDIR/log/gfsfcst/${Syyyymmddhh[$c]}
  fi
  cf=`printf '%04d' $c`
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    np=`cat $tmpnode/machinefile.gfs.${cf}.${name_m[$mt]} | wc -l`
    if [ "${node_m[$mt]}" = "${node[1]}" ]; then
      pcount=$((pcount+np))
      if [ "$pcount" -gt "$ppnl" ]; then
        echo "    wait..."
        wait
        pcount=$np
      fi
    fi
    echo "  run GFS forecast: ${stimef[$c]}, member ${name_m[$mt]} on node ${nodes_gfs_m[$mt]}"
    if [ "$LOG_OPT" -le 2 ]; then
      $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gfs.${cf}.${name_m[$mt]} -n $np \
                      -wdir $ltmpgfs/${cf}_${name_m[$mt]} \
                      ./global_fcst > $OUTDIR/log/gfsfcst/${Syyyymmddhh[$c]}/gfs_${name_m[$mt]}.log 2>&1 &
    else
      $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.gfs.${cf}.${name_m[$mt]} -n $np \
                      -wdir $ltmpgfs/${cf}_${name_m[$mt]} \
                      ./global_fcst > /dev/null 2>&1 &
    fi
    sleep 0.05s
  done
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 3. Convert GFS sig/sfc format into grid format
p=3
step_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/datetime.sh .

cat > fcst_31.sh << EOF
. datetime.sh
mem="\$1"
cyc="\$2"
stime="\$3"
Syyyymmddhh=\${stime:0:10}
mkdir -p $ltmpout/fcst/\${Syyyymmddhh}/\${mem}
mkdir -p $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}
mkdir -p $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}

mkdir -p $ltmpssio/\${cyc}_\${mem}
cd $ltmpssio/\${cyc}_\${mem}
ln -fs $ltmpprog/ss2grd .
ln -fs $ltmpprog/ss2grdp .
fh=0
while [ "\$fh" -le "$FCSTLEN" ]; do
  fhh="\$(printf '%02d' \$fh)"
  fhhh="\$(printf '%03d' \$fh)"
  Fyyyymmddhh=\$(datetime \$stime \$fh h | cut -c 1-10)
  cd $ltmpssio/\${cyc}_\${mem}
  rm -f fort.*
  if [ "\$fh" -eq 0 ]; then
    ln -fs $ltmpgfs/\${cyc}_\${mem}/sig_ini fort.11
    ln -fs $ltmpgfs/\${cyc}_\${mem}/sfc_ini fort.12
  else
    ln -fs $ltmpgfs/\${cyc}_\${mem}/SIG.F\${fhh} fort.11
    ln -fs $ltmpgfs/\${cyc}_\${mem}/SFC.F\${fhh} fort.12
  fi
  ./ss2grd
  mv -f fort.31 $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd
EOF
if [ "$FOUT_OPT" -le 2 ]; then
  cat >> fcst_31.sh << EOF
  ./ss2grdp
  mv -f fort.31 $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd
EOF
fi
if [ "$FOUT_OPT" -le 1 ]; then
  cat >> fcst_31.sh << EOF
  if [ "\$fh" -eq 0 ]; then
    cp -fL $ltmpgfs/\${cyc}_\${mem}/sig_ini $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sig
    cp -fL $ltmpgfs/\${cyc}_\${mem}/sfc_ini $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sfc
  else
    mv -f $ltmpgfs/\${cyc}_\${mem}/SIG.F\${fhh} $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sig
    mv -f $ltmpgfs/\${cyc}_\${mem}/SFC.F\${fhh} $ltmpout/fcst/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.sfc
  fi
EOF
fi
cat >> fcst_31.sh << EOF
fh=\$((fh+$FCSTOUT))
done
EOF

#-------------------------------------------------------------------------------

echo
ppnl=$((ppn*2))
np=1
pcount=0
for c in `seq $CYCLES`; do
  cf=`printf '%04d' $c`
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    if [ "${node_m[$mt]}" = "${node[1]}" ]; then
      pcount=$((pcount+np))
      if [ "$pcount" -gt "$ppnl" ]; then
        echo "    wait..."
        wait
        pcount=$np
      fi
    fi
    echo "  ${stimef[$c]}, member ${name_m[$mt]} on node '${node_m[$mt]}'"
    $MPIBIN/mpiexec -host ${node_m[$mt]} bash fcst_31.sh ${name_m[$mt]} $cf "${STIME[$c]}" &
    sleep 0.05s
  done
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================
# 4. Verification
if [ "$IF_VERF" = '1' ]; then
p=4
step_timing $p
#-------------------------------------------------------------------------------

cd $TMPMPI
cp -f $RUNDIR/datetime.sh .

cat > fcst_41.sh << EOF
. datetime.sh
mem="\$1"
cyc="\$2"
stime="\$3"
Syyyymmddhh=\${stime:0:10}

mkdir -p $ltmpverify/\${cyc}_\${mem}
cd $ltmpverify/\${cyc}_\${mem}
ln -fs $ltmpprog/verify .
fh=0
while [ "\$fh" -le "$FCSTLEN" ]; do
  fhhh=\$(printf '%03d' \$fh)
  Fyyyymmddhh=\$(datetime \$stime \$fh h | cut -c 1-10)

  if [ -s "$ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd" ] &&
     [ -s "$ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd" ]; then
    cd $ltmpverify/\${cyc}_\${mem}
    rm -f fcst.grd fcstp.grd obs??.dat ana??.grd
    ln -s $ltmpout/fcstg/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd fcst.grd
    ln -s $ltmpout/fcstgp/\${Syyyymmddhh}/\${mem}/\${Fyyyymmddhh}.grd fcstp.grd
###### only support shared disk
    cat $OBS/obs\${Fyyyymmddhh}/t.dat > obs01.dat
    cat $OBS/obs\${Fyyyymmddhh}/t-1.dat >> obs01.dat
    cat $OBS/obs\${Fyyyymmddhh}/t+1.dat >> obs01.dat
    ln -s $ANLGRDP/\${Fyyyymmddhh}.grd ana01.grd
#    ln -s $ANLGRDP2/\${Fyyyymmddhh}.grd ana02.grd
######
    ./verify > /dev/null 2>&1

    mkdir -p $ltmpout/verfo1/\${fhhh}/\${mem}
    mkdir -p $ltmpout/verfa1/\${fhhh}/\${mem}
    mkdir -p $ltmpout/verfa2/\${fhhh}/\${mem}
    mv -f vrfobs01.dat $ltmpout/verfo1/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
    mv -f vrfana01.dat $ltmpout/verfa1/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
#    mv -f vrfana02.dat $ltmpout/verfa2/\${fhhh}/\${mem}/\${Fyyyymmddhh}.dat
  fi
fh=\$((fh+$FCSTOUT))
done
EOF

#-------------------------------------------------------------------------------

echo
ppnl=$ppn
#ppnl=$((ppn*2))
np=1
pcount=0
for c in `seq $CYCLES`; do
  cf=`printf '%04d' $c`
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    if [ "${node_m[$mt]}" = "${node[1]}" ]; then
      pcount=$((pcount+np))
      if [ "$pcount" -gt "$ppnl" ]; then
        echo "    wait..."
        wait
        pcount=$np
      fi
    fi
    echo "  ${stimef[$c]}, member ${name_m[$mt]} on node '${node_m[$mt]}'"
    $MPIBIN/mpiexec -host ${node_m[$mt]} bash fcst_41.sh ${name_m[$mt]} $cf "${STIME[$c]}" &
    sleep 0.05s
  done
done
echo "    wait..."
wait

#-------------------------------------------------------------------------------
echo
echo "===================================================================="
fi
#===============================================================================
# 5. Collect outputs (stageout)
p=5
step_timing $p
#-------------------------------------------------------------------------------

for c in `seq $CYCLES`; do
  STIMEgrads=$(datetimegrads ${STIME[$c]})
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    if [ "$FOUT_OPT" -le 1 ]; then
      mkdir -p $OUTDIR/fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}
    fi
    mkdir -p $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}
    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 x > \
                     $OUTDIR/fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhx.ctl
    if [ "$FOUT_OPT" -le 2 ]; then
      mkdir -p $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}
      $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${FCSTOUT}hr 10000 p > \
                       $OUTDIR/fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/yyyymmddhhp.ctl
    fi

    fh=0
    while [ "$fh" -le "$FCSTLEN" ]; do
      fhhh=`printf '%03d' $fh`
      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
      mkdir -p $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
      cd $OUTDIR/fcstv/${fhhh}/${name_m[$mt]}
      ln -fs ../../../fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
      if [ ! -s 'yyyymmddhhx.ctl' ]; then
        $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
                         yyyymmddhhx.ctl
      fi
      if [ "$FOUT_OPT" -le 2 ]; then
        mkdir -p $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
        cd $OUTDIR/fcstvp/${fhhh}/${name_m[$mt]}
        ln -fs ../../../fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd .
        if [ ! -s 'yyyymmddhhp.ctl' ]; then
          $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
                           yyyymmddhhp.ctl
        fi
      fi
    fh=$((fh + FCSTOUT))
    done
  done
done

#-------------------------------------------------------------------------------
if [ "$SHAREDISK" = '0' ]; then
#-------------------------------------------------------------------------------

cd $TMPMPI
mkdir -p $tmpstageout
rm -f $tmpstageout/*

for c in `seq $CYCLES`; do
  for m in `seq $fmember`; do
    mt=$(((c-1) * fmember + m))
    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sig" >> $tmpstageout/out.${node_m[$mt]}
    echo "rm|anal/${name_m[$mt]}/${Syyyymmddhh[$c]}.sfc" >> $tmpstageout/out.${node_m[$mt]}
    fh=0
    while [ "$fh" -le "$FCSTLEN" ]; do
      fhhh=`printf '%03d' $fh`
      Fyyyymmddhh=$(datetime ${STIME[$c]} $fh h | cut -c 1-10)
      if [ "$FOUT_OPT" -le 1 ]; then
        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sig" >> $tmpstageout/out.${node_m[$mt]}
        echo "mv|fcst/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.sfc" >> $tmpstageout/out.${node_m[$mt]}
      fi
      echo "mv|fcstg/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
      if [ "$FOUT_OPT" -le 2 ]; then
        echo "mv|fcstgp/${Syyyymmddhh[$c]}/${name_m[$mt]}/${Fyyyymmddhh}.grd" >> $tmpstageout/out.${node_m[$mt]}
      fi
      echo "mv|verfo1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
      echo "mv|verfa1/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
      echo "mv|verfa2/${fhhh}/${name_m[$mt]}/${Fyyyymmddhh}.dat" >> $tmpstageout/out.${node_m[$mt]}
    fh=$((fh + FCSTOUT))
    done
  done
done

#-------------------------------------------------------------------------------

stageout $ltmpout 0  # clean stageout
$MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes \
                rm -fr $LTMP1/${tmpsubdir} $LTMP2/${tmpsubdir} &
wait

#-------------------------------------------------------------------------------
elif [ "$SHAREDISK" = '1' ]; then
#-------------------------------------------------------------------------------

rm -fr $ltmprun1 $ltmprun2

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------
echo
echo "===================================================================="
#===============================================================================

#rm -fr $TMPMPI

echo
echo " +----------------------------------------------------------------+"
echo " |              GFS-Forecasts successfully completed              |"
echo " +----------------------------------------------------------------+"
echo

echo "[$(now_format)] Finish fcst.sh $@" 1>&2
exit 0

#===============================================================================
