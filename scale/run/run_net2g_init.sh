#!/bin/bash

EXP=2000m_InSnd_LT_SN14_Mac_0605_PAWR
EXP=2000m_InSnd_LT_SN14_Mac_0605_NODA

NNODE=1
PPN=1
PRC=1

NET2G_INIT=/work/hp150019/f22013/SCALE-LETKF/scale-LT/scale_lt_devel_20190522/scale-rm/util/netcdf2grads/net2g_init
RUNDIR=`pwd`

OUTPUT="/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/${EXP}"

TS=1
TE=61
DT=4
STIME=20000101011000
MEM=mean


#VNAME='"U", "V", "W", "PRES", "T", "QV", "QC", "QR", "QI", "QS", "QG", "QCRG_C"'
#VNAME='"QHYD", "QCRG_TOT"'
VNAME='"QV"'

ZMAX=40
DZ=4
TARGET_ZLEV=`seq -s , 1 ${DZ} ${ZMAX}`
NUM_ZLEV=$((ZMAX / DZ))

DIR=anal

IDIR=$OUTPUT/$STIME/${DIR}/${MEM}
ODIR=$OUTPUT/$STIME/${DIR}/${MEM}

TMPRUN=${RUNDIR}/../tmp_net2g
rm -rf ${TMPRUN}
mkdir -p ${TMPRUN}

TMPRUN_MEM=${TMPRUN}/${MEM}
mkdir -p $TMPRUN_MEM

cp ${NET2G_INIT} ${TMPRUN_MEM}

CONF_MEM=namelist.in

INITCONF_MEM=${MEM}_init.conf
INITCONF_MEM_ORG=/work/hp150019/f22013/SCALE-LETKF/scale-LT/OUTPUT/2000m_InSnd_LT_SN14_Mac_0605_NODA/init.conf

cp ${INITCONF_MEM_ORG} ${TMPRUN_MEM}/${INITCONF_MEM}

cat << EOF >> ${TMPRUN_MEM}/${CONF_MEM}
&info
timestep=1,
conffile="${INITCONF_MEM}"
odir="${output_grads_dir}"
vcount=1
&end
&vari
vname="${ielem}",
&end
&grads
delt="10mn"
stime="${stime}"
&end

EOF
exit



!#--------------------------------------------------------------------------------------------------------------
!#   Namelist for netcdf2grads-H (under the nickname of POPSCA)
!#   author: Team SCALE
!#--------------------------------------------------------------------------------------------------------------

&LOGOUT
 LOG_BASENAME   = "LOG_d01",                   !# [C] name of log file; when "STDOUT" is specified,
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
 CONFFILE    = "${RUNCONF_MEM}",          !# [C] path to the config file for the run of scale-rm
 IDIR        = "${IDIR}",                       !# [C] path to the directory of history.pe*.nc
 ODIR        = "${ODIR}",                    !# [C] path to the directory of grads files (output)
 Z_LEV_TYPE  = "original",                         !# [C] output type ("plev","zlev","original","anal")
 Z_MERGE_OUT = .true.,                         !# [L] data output as an array vertically merged
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

# Make a job script for net2g


JOBSH_MEM=${TMPRUN_MEM}/net2g.sh
cat << EOF >> ${JOBSH_MEM}
#!/bin/sh

#PJM -N NET2G
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L node=${NNODE}
#PJM -L elapse=00:30:00
#PJM --mpi proc=${PRC}
#PJM --omp thread=1
#PJM -g hp150019

ulimit -s unlimited

export FORT_FMT_RECL=400

rm -f machinefile
for inode in \$(cat \$I_MPI_HYDRA_HOST_FILE); do
  for ippn in \$(seq $PPN); do
    echo "\$inode" >> machinefile
  done
done

module load hdf5/1.8.17
module load netcdf/4.4.1
module load netcdf-fortran/4.4.3
ulimit -s unlimited
export OMP_STACKSIZE=128m

mpiexec -n ${PRC} ./net2g ${CONF_MEM}

EOF

exit
