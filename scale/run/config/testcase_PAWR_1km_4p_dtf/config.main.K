#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================

#PRESET='<PRESET>'                           # 'K' or 'K_rankdir' or 'K_micro'
PRESET='K'

#===============================================================================

DIR="$(cd "$(pwd)/.." && pwd)"              # Root directory of the SCALE-LETKF

USERNAME=gylien
USERNAME=honda

INDIR="/volume63/data/hp150019/gylien/exp/scale-letkf/benchmark/d4_1km_initbdy_v180531"
OUTDIR="/volume63/data/hp150019/${USERNAME}/exp/scale-letkf/benchmark/testcase_PAWR_1km_4p_dtf" # Directory for SCALE-LETKF output

LIBDTF_PATH=

DTF_MPMD=1 # 1: MPMD launch mode (e.g., mpiexec -n ./scale-rm_ens: -n XX ./letkf)
           # 0: Launch two background jobs (e.g., mpiexec -n XX ./scale-rm_ens& mpiexec -n XX ./letkf&)
SPLIT_WRAP="/volume63/data/hp150019/honda/SCALE-LETKF/scale-productrun-CREST0222/scale-productrun-crest-7f2c3dbd_20181023_stat/split_world_wrap/libsplitworld.so"

#===============================================================================
# Location of model/data files

SCALEDIR="$DIR/../.."                                                # Directory of the SCALE model
DATADIR="/volume63/data/hp150019/gylien/scale-letkf/database"        # Directory of the SCALE database

DATA_TOPO=                                                           # Directory of the prepared topo files
#DATA_TOPO_BDY_SCALE="/volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/d3_1km_312/const/topo"
DATA_LANDUSE=                                                        # Directory of the prepared landuse files
#DATA_BDY_SCALE="/volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/d3_1km_312"
DATA_BDY_SCALE_PRC_NUM_X=3
DATA_BDY_SCALE_PRC_NUM_Y=3
DATA_BDY_SCALE_PREP=                                                 # Directory of the prepared SCALE boundary files
DATA_BDY_WRF=                                                        # Directory of the boundary data in WRF format
DATA_BDY_NICAM=                                                      # Directory of the boundary data in NICAM format (not finished)

OBS="/volume63/data/hp150019/gylien/obs/20130713/raw_v2a"             # Observation data in LETKF observation format
OBSNCEP=
OBS_USE_JITDT=0

#===============================================================================
# model/data file options

PNETCDF=1               # 0: Do not use PnetCDF
                        # 1: Use PnetCDF

PNETCDF_BDY_SCALE=0     # Use PNETCDF I/O for the parent domain SCALE data? (default: $PNETCDF)
                        # 0: Do not use PnetCDF
                        # 1: Use PnetCDF

DET_RUN=1               # 0: Disable the deterministic run
                        # 1: Enable  the deterministic run
DET_RUN_CYCLED=1

TOPO_FORMAT='prep'      # 'prep': Use prepared topo files in $DATA_TOPO
                        # 'GTOPO30' (requires compatible 'config.nml.scale_pp')
                        # 'DEM50M'  (requires compatible 'config.nml.scale_pp')

LANDUSE_FORMAT='prep'   # 'prep': Use prepared landuse files in $DATA_LANDUSE
                        # 'GLCCv2' (requires compatible 'config.nml.scale_pp')
                        # 'LU100M' (requires compatible 'config.nml.scale_pp')
LANDUSE_UPDATE=0        # 0: Time-invariant landuse files
                        # 1: Time-variant landuse files

BDY_FORMAT=0            # 0: SCALE boundary files (with exactly same domain settings; do not need additional preprocessing)
                        # 1: SCALE history (requires compatible 'config.nml.scale_init')
                        # 2: WRF           (requires compatible 'config.nml.scale_init')
                        # 3: NICAM         (requires compatible 'config.nml.scale_init')
                        # 4: GrADS         (requires compatible 'config.nml.scale_init')
BDY_SINGLE_FILE=0       # 0: Length of a boundary file = $BDYCYCLE_INT (e.g., files made by data assimilation cycles)
                        # 1: Length of a boundary file is long enough so that only a single boundary file is used for each forecast
BDY_SCALE_DIR='hist'    # Directory name of the SCALE history files when $BDY_FORMAT = 1

BDY_ENS=1               # 0: Fixed boundary files for all memebers
                        # 1: Ensemble boundary files
BDY_ROTATING=0          # 0: Use a same series of boundary files for all initial time
                        # 1: Use different series of boundary files for different initial time

BDYINT=300
BDYCYCLE_INT=300

PARENT_REF_TIME=20130713060000

ENABLE_PARAM_USER=1     # 0: Do not enable the 'PARAM_USER' section of the SCALE namelist
                        # 1: Enable the 'PARAM_USER' section of the SCALE namelist (require 'config.nml.scale_user' and customized version of SCALE)

OCEAN_INPUT=0           # 0: No ocean input (use cycling ocean variables)
                        # 1: Update the ocean variables every cycle
OCEAN_FORMAT=0          # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)
LAND_INPUT=0            # 0: No land input (use cycling land variables)
                        # 1: Update the land variables every cycle
LAND_FORMAT=0           # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)

OBSNUM=1
OBSNAME[1]=pawr
OBS_FORMAT[1]=PAWR_TOSHIBA

#===============================================================================
# Cycling settings

WINDOW_S=30        # SCALE forecast time when the assimilation window starts (second)
WINDOW_E=30        # SCALE forecast time when the assimilation window ends (second)
LCYCLE=30          # Length of a DA cycle (second)
LTIMESLOT=30       # Timeslot interval for 4D-LETKF (second)

#===============================================================================
# Parallelization settings

#MEMBER=100         # Ensemble size
#
#NNODES=408         # Number of nodes
MEMBER=50         # Ensemble size

NNODES=208         # Number of nodes
PPN=1              # Number of processes per node

THREADS=8          # Number of threads per process

SCALE_NP=4         # Number of processes to run SCALE

BGJOB_INT='0.1s'   # Interval of multiple background job submissions

#===============================================================================
# Temporary directories to store runtime files

DISK_MODE_OBS=2

ONLINE_STGOUT=0             # Stage out right after each cycle (do not wait until the end of the job)?
                            #  0: No
                            #  1: Yes

SYSNAME="$(basename $OUTDIR)"                # A unique name in the machine
TMPSUBDIR="scale-letkf_${SYSNAME}"           # (used to identify multiple runs in the same time)

TMP="/scratch/$(id -ng)/${USERNAME}/$TMPSUBDIR" # Temporary directory shared among all nodes
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary directory only on the server node
#TMPL=
#TMPSL="/tmp/${USERNAME}/$TMPSUBDIR"  # Temporary directory only on the server node and on a local file system
TMP_JITDATA="$TMP/jit-data"

CLEAR_TMP=0                 # Clear temporary directories after the completion of job?
                            #  0: No
                            #  1: Yes

#===============================================================================
# Environmental settings

MPIRUN="mpiexec"
if (which $MPIRUN > /dev/null 2>&1); then
  MPIRUN=$(which $MPIRUN)
fi

SCP='cp -L'
SCP_HOSTPREFIX=''
#SCP="scp -q"
#SCP_HOSTPREFIX="XXXX:"

STAGE_THREAD=8
TAR_THREAD=8

PYTHON="python"

#BUFRBIN=

#===============================================================================
# Machine-independent source file

. config.rc

#===============================================================================
