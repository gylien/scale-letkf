#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================

PRESET='K_rankdir'                          # 'K' or 'K_rankdir' or 'K_micro'

USERNAME=honda
#===============================================================================

DIR="$(cd "$(pwd)/.." && pwd)"              # Root directory of the SCALE-LETKF

DDIR="$(cd "$(pwd)/../../../.." && pwd)"    # Root directory of the test suite

INDIR=                                      # Directory for SCALE-LETKF input (default: $OUTDIR)
OUTDIR="$DDIR/exp/testcase_45km_4p_l36"     # Directory for SCALE-LETKF output

#===============================================================================
# Location of model/data files

MODELDIR="$DIR/../../bin/${SCALE_SYS}"                               # Directory of the SCALE model executables
DATADIR="$DDIR/database"                                             # Directory of the SCALE database

DATA_TOPO=                                                           # Directory of the prepared topo files
DATA_TOPO_BDY_SCALE=
DATA_LANDUSE=                                                        # Directory of the prepared landuse files
DATA_BDY_SCALE=                                                      # Directory of the boundary data in SCALE history format (parent domain)
DATA_BDY_SCALE_PREP=                                                 # Directory of the prepared SCALE boundary files
DATA_BDY_WRF="$DDIR/ncepfnl/wrfout"                                  # Directory of the boundary data in WRF format
DATA_BDY_NICAM=                                                      # Directory of the boundary data in NICAM format (not finished)
DATA_BDY_GRADS="$DDIR/ncepfnl_grads"                                 # Directory of the boundary data in GrADS format

OBS="$DDIR/obs/prepbufr_obs_eastasia"                                # Observation data in LETKF observation format
OBSNCEP=

#===============================================================================
# model/data file options


PARENT_REF_TIME=

ENABLE_PARAM_USER=0     # 0: Do not enable the 'PARAM_USER' section of the SCALE namelist
                        # 1: Enable the 'PARAM_USER' section of the SCALE namelist (require 'config.nml.scale_user' and customized version of SCALE)


OBSNUM=1
OBSNAME[1]=obs                                    
OBSNAME[2]=radar
OBSOPE_SEPARATE[1]=0
OBSOPE_SEPARATE[2]=0

#===============================================================================
# Cycling settings

WINDOW_S=600     # SCALE forecast time when assimilation window starts (second)
WINDOW_E=600     # SCALE forecast time when assimilation window ends (second)
LCYCLE=600       # Length of a GFS-LETKF cycle (second)
LTIMESLOT=600     # Timeslot interval for 4D-LETKF (second)

#===============================================================================
# Parallelization settings

MEMBER=10          # Ensemble size

NNODES=44          # Number of nodes
PPN=1              # Number of processes per node

THREADS=8          # Number of threads per process

SCALE_NP=4         # Number of processes to run SCALE

BGJOB_INT='0.1s'   # Interval of multiple background job submissions

#ENABLE_SET=1       ######

#===============================================================================
# Temporary directories to store runtime files

DATA_BDY_TMPLOC=1           # Location of the temporary directory for DATA_BDY
                            #  1: in $TMPDAT
                            #  2: in $TMPOUT

DISK_MODE_TOPO_LANDUSE_DB=2

DISK_MODE_DATA_BDY=2

DISK_MODE_DATA_TOPO=2
DISK_MODE_DATA_LANDUSE=2
DISK_MODE_DATA_BDYPREP=2

ONLINE_STGOUT=0             # Stage out right after each cycle (do not wait until the end of the job)?
                            #  0: No
                            #  1: Yes

#SIMPLE_STGOUT=1             # Stage out the entire 'out' temporary directory instead of specifying each file?
#                            #  0: No
#                            #  1: Yes

SYSNAME="$(basename $OUTDIR)"                # A unique name in the machine
TMPSUBDIR="scale-letkf_${SYSNAME}" # (used to identify multiple runs in the same time)

TMP="/scratch/$(id -ng)/${USERNAME}/$TMPSUBDIR" # Temporary directory shared among all nodes
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary directory only on the server node
#TMPL=

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

SCP_THREAD=8
TAR_THREAD=8

PYTHON="python"

#BUFRBIN=

#===============================================================================
# Machine-independent source file

. config.rc

#===============================================================================
