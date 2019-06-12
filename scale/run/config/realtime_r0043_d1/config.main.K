#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================

PRESET='<PRESET>'                           # 'K' or 'K_rankdir' or 'K_micro'

EXPNAME=r0043
EXPDOMAIN=1

#===============================================================================

DIR="$(cd "$(pwd)/.." && pwd)"              # Root directory of the SCALE-LETKF

OUTDIR="/volume63/data/hp150019/gylien/realtime/${EXPNAME}/exp_d${EXPDOMAIN}" # Directory for SCALE-LETKF output

#===============================================================================
# Location of model/data files

SCALEDIR="$DIR/../.."                                                # Directory of the SCALE model
DATADIR="/volume63/data/hp150019/gylien/scale-letkf/database"        # Directory of the SCALE database

DATA_TOPO=                                                           # Directory of the prepared topo files
DATA_TOPO_BDY_SCALE=
DATA_LANDUSE=                                                        # Directory of the prepared landuse files
DATA_BDY_SCALE=                                                      # Directory of the boundary data in SCALE history format (parent domain)
DATA_BDY_SCALE_PRC_NUM_X=
DATA_BDY_SCALE_PRC_NUM_Y=
DATA_BDY_SCALE_PREP=                                                 # Directory of the prepared SCALE boundary files
DATA_BDY_WRF=                                                        # Directory of the boundary data in WRF format
DATA_BDY_NICAM=                                                      # Directory of the boundary data in NICAM format (not finished)
DATA_BDY_GRADS="/volume63/data/hp150019/gylien/realtime/${EXPNAME}/<DATA_BDY_GRADS>" # Directory of the boundary data in GrADS format

OBS="/volume63/data/hp150019/gylien/realtime/${EXPNAME}/ncepobs"     # Directory of observation data in LETKF obs format
OBSNCEP=                                                             # Directory of observation data in NCEP BUFR format

#===============================================================================
# model/data file options

DET_RUN=1               # 0: Disable the deterministic run
                        # 1: Enable  the deterministic run

TOPO_FORMAT='prep'      # 'prep': Use prepared topo files in $DATA_TOPO
                        # 'GTOPO30' (requires compatible 'config.nml.scale_pp')
                        # 'DEM50M'  (requires compatible 'config.nml.scale_pp')

LANDUSE_FORMAT='prep'   # 'prep': Use prepared landuse files in $DATA_LANDUSE
                        # 'GLCCv2' (requires compatible 'config.nml.scale_pp')
                        # 'LU100M' (requires compatible 'config.nml.scale_pp')
LANDUSE_UPDATE=0        # 0: Time-invariant landuse files
                        # 1: Time-variant landuse files

BDY_FORMAT=4            # 0: SCALE boundary files (with exactly same domain settings; do not need additional preprocessing)
                        # 1: SCALE history (requires compatible 'config.nml.scale_init')
                        # 2: WRF           (requires compatible 'config.nml.scale_init')
                        # 3: NICAM         (requires compatible 'config.nml.scale_init')
                        # 4: GrADS         (requires compatible 'config.nml.scale_init')
BDY_SINGLE_FILE=0       # 0: Length of a boundary file = $BDYCYCLE_INT (e.g., files made by data assimilation cycles)
                        # 1: Length of a boundary file is long enough so that only a single boundary file is used for each forecast
BDY_SCALE_DIR='hist'    # Directory name of the SCALE history files when $BDY_FORMAT = 1

BDY_ENS=0               # 0: Fixed boundary files for all memebers
                        # 1: Ensemble boundary files
BDY_ROTATING=1          # 0: Use a same series of boundary files for all initial time
                        # 1: Use different series of boundary files for different initial time

BDYINT=21600
BDYCYCLE_INT=21600

PARENT_REF_TIME=

ENABLE_PARAM_USER=1     # 0: Do not enable the 'PARAM_USER' section of the SCALE namelist
                        # 1: Enable the 'PARAM_USER' section of the SCALE namelist (require 'config.nml.scale_user' and customized version of SCALE)

OCEAN_INPUT=1           # 0: No ocean input (use cycling ocean variables)
                        # 1: Update the ocean variables every cycle
OCEAN_FORMAT=99         # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)
LAND_INPUT=1            # 0: No land input (use cycling land variables)
                        # 1: Update the land variables every cycle
LAND_FORMAT=99          # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)

OBSNUM=1
OBSNAME[1]=obs

#===============================================================================
# Cycling settings

WINDOW_S=10800     # SCALE forecast time when the assimilation window starts (second)
WINDOW_E=32400     # SCALE forecast time when the assimilation window ends (second)
LCYCLE=21600       # Length of a DA cycle (second)
LTIMESLOT=3600     # Timeslot interval for 4D-LETKF (second)

#===============================================================================
# Parallelization settings

MEMBER=50          # Ensemble size

NNODES=<NNODES>    # Number of nodes
PPN=1              # Number of processes per node

THREADS=8          # Number of threads per process

SCALE_NP=48        # Number of processes to run SCALE

BGJOB_INT='0.1s'   # Interval of multiple background job submissions

#===============================================================================
# Temporary directories to store runtime files

ONLINE_STGOUT=0             # Stage out right after each cycle (do not wait until the end of the job)?
                            #  0: No
                            #  1: Yes

SYSNAME="$(basename $OUTDIR)"                # A unique name in the machine
TMPSUBDIR="scale-letkf_${SYSNAME}"           # (used to identify multiple runs in the same time)

TMP="/scratch/$(id -ng)/gylien/$TMPSUBDIR" # Temporary directory shared among all nodes
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary directory only on the server node
#TMPL=

CLEAR_TMP=1                 # Clear temporary directories after the completion of job?
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
