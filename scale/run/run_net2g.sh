#!/bin/bash

MODE=debug-flat
#MODE=regular-flat
ELAPSE="00:30:00"

EXP=2000m_InSnd_LT_SN14_Mac_0605_PAWR
EXP=2000m_InSnd_LT_SN14_Mac_0605_NODA
EXP=2000m_LT_SN14_0710_DA_TH02_RTPS0.95_0802
EXP=2000m_LT_SN14_0710_DA_TH02_RTPS0.95_0802

EXP=500m_WK1982_LT_SN14_NATURE

#EXP=2000m_WK1982_LT_SN14_NATURE
#EXP=2000m_WK1982_LT_SN14_DA_AKSOY_UVP_0704_MEM01

#EXP=2000m_WK1982_NOLT_TOMITA_RHO
#EXP=2000m_WK1982_NOLT_TOMITA_NORHO

#EXP=2000m_WK1982_NOLT_SN14
#SMEM=2
#EMEM=5 # 12min

# 4mem ~ 12min
# 15mem > 24min
# 5mem ~ 15min


SMEM=0
EMEM=0

NMEM=$((EMEM - SMEM + 1))

#PRC_MEM=64 # PRC for each member
PRC_MEM=1024 # PRC for each member
TPRC=$((PRC_MEM * NMEM))
TPRC=$PRC_MEM # TEST
PPN=64

if (( TPRC < PPN )) ; then
  PPN=$TPRC
fi
NNODE=$((TPRC / PPN))


NET2G=/work/hp150019/f22013/SCALE-LETKF/scale-LT/scale_lt_devel_20190522/scale-rm/util/netcdf2grads_h/net2g
RUNDIR=`pwd`

OUTPUT="/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/${EXP}"


TS=1
TE=61
DT=2 # 1min
STIME=20010101012000

STIME=20010101011000
TS=1
TE=21
DT=2 # 1min

STIME=20010101010000
TS=1
TE=16
DT=2 # 1min

T_MERGE_OUT=T # CHECK

#TS=41
#TE=41
#T_MERGE_OUT=F # CHECK

# Variable
#VNAME='"U", "V", "W", "PRES", "T", "QV", "QC", "QR", "QI", "QS", "QG", "QCRG_C"'
VNAME='"QHYD"'
VNAME='"QI", "QS","QG", "QHYD"'
VNAME='"QI", "U", "V"'

VNAME='"QHYD", "QCRG_TOT", "NegFLASH", "PosFLASH"'
VNAME='"QHYD"'
VNAME='"FlashPoint"'
#VNAME='"QCRG_TOT"'
VNAME='"NegFLASH", "PosFLASH", "QCRG_G"'
#VNAME='"NegFLASH"'
VNAME='"LTpath"'

VNAME='"QCRG_TOT", "QHYD", "LTpath"'
VNAME='"QCRG_TOT", "QHYD"'
VNAME='"LTpath"'
VNAME='"PosFLASH"'

#VNAME='"Ex", "Ey", "Ez", "FlashPoint"'

# Level
#ZMAX=40
#ZMIN=1
#DZ=5
#TARGET_ZLEV=`seq -s , ${ZMIN} ${DZ} ${ZMAX}`
#NUM_ZLEV=$((ZMAX / DZ))
TARGET_ZLEV="1, 5, 10, 15, 20, 25, 30, 35, 40"
NUM_ZLEV=9

TARGET_ZLEV=`seq 1 40`
NUM_ZLEV=40

# Make a temporal directory
TMPRUN=${RUNDIR}/../tmp_net2g/${EXP}_${SMEM}_${EMEM}
rm -rf ${TMPRUN}
mkdir -p ${TMPRUN}

cp ${NET2G} ${TMPRUN}



# Make a job script for net2g

JOBSH_MEM=${TMPRUN}/net2g.sh
cat << EOF >> ${JOBSH_MEM}
#!/bin/sh

#PJM -N NET2G
#PJM -L rscgrp=${MODE}
#PJM -L node=${NNODE}
#PJM -L elapse=${ELAPSE}
#PJM --mpi proc=${TPRC}
#PJM --omp thread=1
#PJM -g hp150019

ulimit -s unlimited

export FORT_FMT_RECL=400

#rm -f machinefile
#for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
#  for ippn in \$(seq $PPN); do
#    echo "\$inode" >> machinefile
#  done
#done

module load hdf5/1.8.17
module load netcdf/4.4.1
module load netcdf-fortran/4.4.3
ulimit -s unlimited
export OMP_STACKSIZE=128m

CTIME=\$(date +"%Y/%m/%d %H:%M:%S")
echo "Start "\$CTIME

EOF

# Member loop
for MEM in `seq $SMEM $EMEM`
do
  mmmm=$(printf %04d $MEM)
  if (( MEM == 0 )) ; then
    mmmm="mean"
  fi
  echo $MEM" "$mmmm

  IDIR=$OUTPUT/$STIME/fcst/${mmmm}
  ODIR=$OUTPUT/$STIME/fcst/${mmmm}

  # run.conf file
  RMEM=${mmmm}
  CONF_RMEM=${RMEM}_net2g.conf

  RUNCONF_RMEM=${RMEM}_run.conf
  RUNCONF_RMEM_ORG=${OUTPUT}/${STIME}/log/fcst_scale/${RUNCONF_RMEM}

  cp ${RUNCONF_RMEM_ORG} ${TMPRUN}/${RUNCONF_RMEM}

  CONF_MEM=${mmmm}_net2g.conf
cat << EOF >> ${TMPRUN}/${CONF_MEM}
!#--------------------------------------------------------------------------------------------------------------
!#   Namelist for netcdf2grads-H (under the nickname of POPSCA)
!#   author: Team SCALE
!#--------------------------------------------------------------------------------------------------------------

&LOGOUT
 LOG_BASENAME   = "LOG_${mmmm}",                   !# [C] name of log file; when "STDOUT" is specified,
                                               !#     log will be outputed to standard out.
 LOG_ALL_OUTPUT = .false.,                     !# [L] log file output in all mpi processes
 LOG_LEVEL      = 1,                           !# [I] log level; 0=debug, 1=normal ( default=1 )
/

&INFO
 START_TSTEP = ${TS},                              !# [I] No. of timestep for start step
 END_TSTEP   = ${TE},                             !# [I] No. of timestep for end step
 INC_TSTEP   = ${DT},                              !# [I] increment size of timestep
 DOMAIN_NUM  = 1,                              !# [I] number of domain
 ZCOUNT      = ${NUM_ZLEV},                              !# [I] the number of target z-levels
 ZSTART      = 1,                              !# [I] the number of start grid for z-level (default=1).
 CONFFILE    = "${RUNCONF_RMEM}",          !# [C] path to the config file for the run of scale-rm
 IDIR        = "${IDIR}",                       !# [C] path to the directory of history.pe*.nc
 ODIR        = "${ODIR}",                    !# [C] path to the directory of grads files (output)
 Z_LEV_TYPE  = "original",                         !# [C] output type ("plev","zlev","original","anal")
 Z_MERGE_OUT = .true.,                         !# [L] data output as an array vertically merged
 T_MERGE_OUT = "${T_MERGE_OUT}",
/

!&EXTRA
! EXTRA_TINTERVAL = 600.0,                      !# [F] specify irregular time interval (i.e. time1, 2...)
! EXTRA_TUNIT     = "SEC",                      !# [C] time unit for extra_tinterval
!/

!&ANAL
! ANALYSIS    = "ave"                           !# [C] analysis options = max, min, sum, ave (if Z_LEV_TYPE="anal")
!/

&VARI
 VNAME       = ${VNAME},        !# [C] the name of variable to convert
                                               !#     in this version, 3D variables and 2D variables cannot
                                               !#     be executed together. plz, execute separately.
 TARGET_ZLEV = ${TARGET_ZLEV},                 !# [I] array of target z-levels:
                                               !#     set vertical grid number or height(m) (if Z_LEV=TYPE="original")
                                               !#     set height(m)                        (if Z_LEV=TYPE="zlev")
                                               !#     set pressure(hPa)                    (if Z_LEV=TYPE="plev")
/

&GRADS
 DELT        = "1mn"                           !# [C] for grads ctl file. Defalt is "1mn".
 STIME       = "00:00Z01JAN2000"               !# [C] for grads ctl file. Defalt is "00:00Z01JAN2000".
/

EOF


cat << EOF >> ${JOBSH_MEM}
mpiexec -n ${PRC_MEM} ./net2g ${CONF_MEM}&
EOF

done # member loop


cat << EOF >> ${JOBSH_MEM}
wait
CTIME=\$(date +"%Y/%m/%d %H:%M:%S")
echo "End "\$CTIME
EOF


cd ${TMPRUN}
pjsub net2g.sh
cd -


exit











