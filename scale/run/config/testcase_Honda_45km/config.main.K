#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================

#PRESET='<PRESET>'       # 'K' or 'K_rankdir' or 'K_micro'
PRESET='K_micro'       # 'K' or 'K_rankdir' or 'K_micro'

#===============================================================================

# -- H08 --
# -- RTTOV_DIR --
DIR_RTTOV=/data/share025/honda/RTTOV
RTTOV_COEF=${DIR_RTTOV}/rtcoef_rttov11/rttov7pred54L/rtcoef_himawari_8_ahi.dat
RTTOV_SCCOEF=${DIR_RTTOV}/rtcoef_rttov11/cldaer/sccldcoef_himawari_8_ahi.dat

#===============================================================================


EXP="testcase_Him8_45km"

DIR="$(cd "$(pwd)/.." && pwd)"              # Root directory of the SCALE-LETKF

DDIR="$(cd "$(pwd)/../../../.." && pwd)"    # Root directory of the test suite

OUTDIR="${DDIR}/OUTPUT/${EXP}" # Directory for SCALE-LETKF output

#===============================================================================
# Location of model/data files

SCALEDIR="$DIR/../.."                                                # Directory of the SCALE model
DATADIR="/volume63/data/hp150019/gylien/scale-letkf/database"        # Directory of the SCALE database

PEST_PMAX=2
PEST_NAME[1]=PARAM_ATMOS_PHY_MP_TOMITA08_Cr
PEST_NAME[2]=PARAM_ATMOS_PHY_MP_TOMITA08_Cs


SCALEDIR="$DIR/../.."                                                # Directory of the SCALE model
DATADIR="/volume63/data/hp150019/gylien/database"                    # Directory of the SCALE database

DATA_TOPO="/volume63/data/hp150019/gylien/scale-letkf-test-suite/exp/testcase_45km_4p_l36_data"       # Directory of the prepared topo files
DATA_TOPO_BDY_SCALE=
DATA_LANDUSE="/volume63/data/hp150019/gylien/scale-letkf-test-suite/exp/testcase_45km_4p_l36_data"     # Directory of the prepared landuse files
DATA_BDY_SCALE="/volume63/data/hp150019/honda/SCALE-LETKF/scale-5.2.2/OUTPUT/KYUSHU2017_D1"                                                      # Directory of the boundary data in SCALE history format (parent domain)
DATA_BDY_SCALE_PREP="/volume63/data/hp150019/gylien/scale-letkf-test-suite/exp/testcase_45km_4p_l36_data"  # Directory of the prepared SCALE boundary files
#DATA_BDY_WRF="/volume63/data/hp150019/gylien/scale-letkf-test-suite/ncepfnl/wrfout"              # Directory of the boundary data in WRF format
DATA_BDY_WRF="/volume63/data/hp150019/honda/SCALE-LETKF/scale-5.2.2/OUTPUT/testcase_Him8_45km/wrfbdy"
DATA_BDY_NICAM=                                                      # Directory of the boundary data in NICAM format (not finished)
DATA_BDY_GRADS="/volume63/data/hp150019/gylien/scale-letkf-test-suite/ncepfnl_grads"             # Directory of the boundary data in GrADS format


OBS="/volume63/data/hp150019/honda/SCALE-LETKF/scale-5.2.2/OUTPUT/testcase_Him8_45km/obs"
OBSNCEP=                               # Directory of observation data in NCEP BUFR format

#===============================================================================
# model/data file options

#TOPO_FORMAT='GTOPO30'   # 'prep': Use prepared topo files in $DATA_TOPO
TOPO_FORMAT='prep'   # 'prep': Use prepared topo files in $DATA_TOPO
                        # 'GTOPO30' (requires compatible 'config.nml.scale_pp')
                        # 'DEM50M'  (requires compatible 'config.nml.scale_pp')

#LANDUSE_FORMAT='GLCCv2' # 'prep': Use prepared landuse files in $DATA_LANDUSE
LANDUSE_FORMAT='prep' # 'prep': Use prepared landuse files in $DATA_LANDUSE
                        # 'GLCCv2' (requires compatible 'config.nml.scale_pp')
                        # 'LU100M' (requires compatible 'config.nml.scale_pp')
LANDUSE_UPDATE=0        # 0: Time-invariant landuse files
                        # 1: Time-variant landuse files

BDY_FORMAT=2            # 0: SCALE boundary files (with exactly same domain settings; do not need additional preprocessing)

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

BDYINT=21600
BDYCYCLE_INT=21600


PARENT_REF_TIME=


OCEAN_INPUT=0           # 0: No ocean input (use cycling ocean variables)
                        # 1: Update the ocean variables every cycle
OCEAN_FORMAT=99         # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)
LAND_INPUT=0            # 0: No land input (use cycling land variables)
                        # 1: Update the land variables every cycle
LAND_FORMAT=99          # 0: SCALE init files (with exactly same domain settings; do not need additional preprocessing)
                        # 99: From the same file as used in generating the boundary conditions ($BDY_FORMAT)

OBSNUM=1
#OBSNAME[1]=H08                                    
OBSNAME[1]=obs                                 

OBSOPE_SEPARATE[1]=0

#===============================================================================
# Cycling settings

WINDOW_S=3600     # SCALE forecast time when the assimilation window starts (second)
WINDOW_E=3600     # SCALE forecast time when the assimilation window ends (second)
LCYCLE=3600       # Length of a DA cycle (second)
LTIMESLOT=3600     # Timeslot interval for 4D-LETKF (second)

#===============================================================================
# Parallelization settings

MEMBER=3          # Ensemble size

NNODES=16        # Number of nodes

PPN=1              # Number of processes per node
THREADS=8          # Number of threads per process

SCALE_NP=4        # Number of processes to run SCALE

BGJOB_INT='0.1s'   # Interval of multiple background job submissions

#===============================================================================
# Temporary directories to store runtime files

ONLINE_STGOUT=0             # Stage out right after each cycle (do not wait until the end of the job)?
                            #  0: No
                            #  1: Yes

SYSNAME="$(basename $OUTDIR)"                # A unique name in the machine
TMPSUBDIR="scale-letkf_${SYSNAME}"           # (used to identify multiple runs in the same time)

TMP="/scratch/$(id -ng)/honda/$TMPSUBDIR" # Temporary directory shared among all nodes
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary directory only on the server node
#TMPL=
#TMPL="/tmp/honda/$TMPSUBDIR"

CLEAR_TMP=0                 # Clear temporary directories after the completion of job?
                            #  0: No
                            #  1: Yes

#===============================================================================
# Environmental settings

MPIBIN=$(dirname $(which mpiexec))
MPIRUN="$MPIBIN/mpiexec"

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
