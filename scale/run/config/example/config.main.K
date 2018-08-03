#!/bin/bash
#===============================================================================
#
#  Main settings for SCALE-LETKF scripts
#
#===============================================================================
# Top-level settings

PRESET='K_rankdir'      # Pre-determined combination of settings for particular system environments
                        # - 'K':  K computer regular (staged) job
                        # - 'K_rankdir':  K computer regular (staged) job with the rank directory enabled
                        # - 'K_micro':  K computer micro job
                        # - 'OFP':  Oakforest-PACS
                        # - 'Linux':  Linux (single machine)
                        # - 'Linux_torque':  Linux cluster with a Torque job scheduling system

RUN_LEVEL=              # [SHOULD NOT SET BY USER] Level to control how many parts of the code will be run 
                        # - 0:  Run everything
                        # - 1:  Run everything but skipping detecting errors
                        # - 2:  [NOT IMPLEMENTED] Staging list files are ready; skip generating them
                        # - 3:  Staging-in has been done; skip staging
                        # - 4:  Staging-in/out is done outside this script; skiping staging
                        #  (default: 0)

#===============================================================================
# Location of the root and input/output directories

DIR="$(cd "$(pwd)/.." && pwd)"              # Root directory of the SCALE-LETKF source code

DDIR="$(cd "$(pwd)/../../../.." && pwd)"    ## Root directory of the test suite

INDIR=                                      # SCALE-LETKF experiment directory for input data
                                            #  (The initial condition data need to exist in this directory unless $MAKEINIT = 1 in config.cycle or config.fcst)
                                            #  (default: $OUTDIR)
OUTDIR="$DDIR/exp/testcase_45km_4p_l36"     # SCALE-LETKF experiment directory for output data

#===============================================================================
# Location of model/data files

SCALEDIR="$DIR/../.."                       # Directory of the SCALE model source code
DATADIR="$DDIR/database"                    # Directory of the SCALE database

DATA_TOPO=                                  # SCALE-LETKF experiment directory under which prepared topo files exist
                                            #  (effective only when $TOPO_FORMAT = 'prep') (default: $INDIR)
DATA_LANDUSE=                               # SCALE-LETKF experiment directory under which prepared landuse files exist
                                            #  (effective only when $LANDUSE_FORMAT = 'prep') (default: $INDIR)
DATA_BDY_SCALE_PREP=                        # SCALE-LETKF experiment directory under which prepared SCALE boundary files exist
                                            #  (effective only when $BDY_FORMAT = 0) (default: $INDIR)
DATA_BDY_SCALE=                             # Parent domain's SCALE-LETKF experiment directory for boundary data
                                            #  (effective only when $BDY_FORMAT = 1: offline-nesting run)
DATA_TOPO_BDY_SCALE=                        # Directory of the parent domain's topo files
                                            #  (effective only when $BDY_FORMAT = 1: offline-nesting run)
DATA_BDY_WRF="$DDIR/ncepfnl/wrfout"         # Directory of boundary data in WRF format; see [[File naming convention of data directories]]
                                            #  (effective only when $BDY_FORMAT = 2)
DATA_BDY_NICAM=                             ## [NOT IMPLEMENTED] Directory of boundary data in NICAM format; see [[File naming convention of data directories]]
                                            ##  (effective only when $BDY_FORMAT = 3)
DATA_BDY_GRADS="$DDIR/ncepfnl_grads"        # Directory of boundary data in GrADS format; see [[File naming convention of data directories]]
                                            #  (effective only when $BDY_FORMAT = 4)
DATA_ADDINFL=                               # Directory of additive inflation files
                                            #  (effective only when $ADDINFL = 1 in config.cycle) (default: $INDIR)

OBS="$DDIR/obs/prepbufr_obs_eastasia"       # Directory of observation data files; see [[File naming convention of data directories]]
OBSNCEP=                                    ## (not implemented)

#===============================================================================
# Model/data file options

PNETCDF=1               # Use PnetCDF for single-file I/O?
                        # - 0:  No
                        # - 1:  Yes
                        #  (default: 0)

PNETCDF_BDY_SCALE=      # Used PnetCDF I/O in the parent domain's SCALE-LETKF experiment?
                        #  (effective only when $BDY_FORMAT = 1: offline-nesting run)
                        # - 0:  No
                        # - 1:  Yes
                        #  (default: $PNETCDF)

DET_RUN=0               # Enable the deterministic run (Schraff et al. 2016 QJRMS)?
                        # - 0:  No
                        # - 1:  Yes
                        #  (default: 0)

TOPO_FORMAT='GTOPO30'   # Topography data file format
                        # - 'prep':  Use prepared topo files in $DATA_TOPO
                        # - 'GTOPO30':  Use 'GTOPO30' dataset (requires compatible 'config.nml.scale_pp')
                        # - 'DEM50M':  Use 'DEM50M' dataset (requires compatible 'config.nml.scale_pp')

LANDUSE_FORMAT='GLCCv2' # Land-use data file format
                        # - 'prep':  Use prepared landuse files in $DATA_LANDUSE
                        # - 'GLCCv2':  Use 'GLCCv2' dataset (requires compatible 'config.nml.scale_pp')
                        # - 'LU100M':  Use 'LU100M' dataset (requires compatible 'config.nml.scale_pp')
LANDUSE_UPDATE=0        # Use time-variant landuse files
                        # - 0:  No (time-invariant landuse files)
                        # - 1:  Yes
                        #  (default: 0)

BDY_FORMAT=2            # Boundary data format
                        # - 0:  Use prepared SCALE boundary files with exactly the same domain settings; no additional pre-processing is required ('scale_init' step will be skipped)
                        # - 1:  SCALE history file: offline-nesting run (requires compatible 'config.nml.scale_init')
                        # - 2:  WRF (requires compatible 'config.nml.scale_init')
                        ## - 3:  [NOT IMPLEMENTED] NICAM (requires compatible 'config.nml.scale_init')
                        # - 4:  GrADS (requires compatible 'config.nml.scale_init')
BDY_SINGLE_FILE=0       # Boundary data for each cycle are contained in a single file of a sufficient length?
                        # - 0:  No - Assume the length of a single boundary file = $BDYCYCLE_INT which may be shorter than the required forecast length; more than one boundary files may be used for each cycle (e.g., files made by data assimilation cycles)
                        # - 1:  Yes - Assume the length of a single boundary file >= the required forecast length (i.e., $WINDOW_E for 'cycle' job; $FCSTLEN for 'fcst' job); always only one single boundary file is used for each cycle
                        #  (default: 0)
BDY_SCALE_DIR='hist'    # Sub-directory name of the SCALE history files in the parent domain's SCALE-LETKF experiment directory
                        #  (effective only when $BDY_FORMAT = 1: offline-nesting run)
                        #  (default: 'hist')

BDY_MEAN='mean'         # Sub-directory name representing the ensemble mean in the boundary data directory
                        #  (default: 'mean')
BDY_ENS=0               # Use ensemble boundary conditions?
                        # - 0:  No - Use a fixed boundary condition (under $BDY_MEAN) for all memebers
                        # - 1:  Yes - Use ensemble boundary conditions
                        #  (default: 0)
BDY_ROTATING=0          # Use different series of boundary data for different cycles?
                        # - 0:  No - Use the same series of boundary data files for all cycles
                        # - 1:  Yes - Use different series of boundary data files (which can overlap in time) for different cycles
                        # Note that this setting affects the file naming convention in the $DATA_BDY_* boundary data directory; see [[File naming convention of data directories]]
                        #  (default: 0)

BDYINT=21600            # Time interval (second) of each time frame in the boundary data with multiple time frames in a file.
                        # If there is only a single time frame in the boundary data files, set this value to be equal to $BDYCYCLE_INT, the time interval of each boundary data file.
                        # In the offline-nesting run, this should be equal to the time interval of the history files in the parent domain's experiment.
                        #  (default: $LCYCLE)
BDYCYCLE_INT=           # Time interval (second) of each boundary data file (regardless of the multiple time frames in a file).
                        # In the offline-nesting run, this should be equal to $LCYCLE in the parent domain's experiment.
                        #  (default: $BDYINT)

PARENT_REF_TIME=        # [OPTIONAL] Any time (format: YYYYMMDDhhmmss) at which a boundary data file exists (for boundary data with multiple time frames in a file, use the start time in a file).
                        # This is used as a reference time to search the boundary data files; the program will search available boundary data files by plus/minus increments of $BDYCYCLE_INT based on this reference time.
                        # This setting is only important when $LCYCLE < $BDYCYCLE_INT, likely in a offline-nesing run. In the case that $LCYCLE >= $BDYCYCLE_INT, this setting can be left blank.

ENABLE_PARAM_USER=0     # Use '&PARAM_USER' namelist section in the 'config.nml.scale_user' configuration file?
                        # - 0:  No
                        # - 1:  Yes - This setting can only be used with a customized version of SCALE for SCALE-LETKF and it requires 'config.nml.scale_user' configuration file
                        #  (default: 0)

OCEAN_INPUT=1           # In every cycle, update ocean variables from the ocean input files?
                        # - 0:  No - Cycle the ocean variables same as for atmospheric variables
                        # - 1:  Yes - The list of ocean variables to be updated can be controlled in 'config.nml.scale_user' if using a customized version of SCALE for SCALE-LETKF; otherwise, all ocean variables will be updated
                        #  (default: 0)
OCEAN_FORMAT=99         # The source of the ocean input files
                        #  (effective only when $OCEAN_INPUT = 1)
                        # - 0:  Use prepared SCALE restart (init) files with exactly the same domain settings; no additional pre-processing is required
                        # - 99:  Use the same data source of the boundary data (based on $BDY_FORMAT setting)
                        #  (default: 99)
LAND_INPUT=1            # In every cycle, update land variables from the land input files?
                        # - 0:  No - Cycle the land variables same as for atmospheric variables
                        # - 1:  Yes - The list of land variables to be updated can be controlled in 'config.nml.scale_user' if using a customized version of SCALE for SCALE-LETKF; otherwise, all land variables will be updated
                        #  (default: 0)
LAND_FORMAT=99          # The source of the land input files
                        #  (effective only when $OCEAN_INPUT = 1)
                        # - 0:  Use prepared SCALE restart (init) files with exactly the same domain settings; no additional pre-processing is required
                        # - 99:  Use the same data source of the boundary data (based on $BDY_FORMAT setting)
                        #  (default: 99)

OBSNUM=1                # Number of observation data files (in $OBS directory) to be assimilated at the same time
OBSNAME[1]=obs          # Array of the file name prefix for each observation data file [1..$OBSNUM]; see [[File naming convention of data directories]]
OBSNAME[2]=radar
OBSOPE_SEPARATE[1]=0    # Array of values indicating whether each observation data file should be processed using a separate observation operator program [1..$OBSNUM]
                        # - 0:  No - Use built-in observation operator in the LETKF program
                        # - 1:  Yes - Use a separate observation operator program
                        #  (default: 0)
OBSOPE_SEPARATE[2]=0

#===============================================================================
# Cycling settings

WINDOW_S=10800     # SCALE forecast time when the assimilation window starts (second)
                   #  (Set this value to be equal to $LCYCLE for 3D-LETKF)
WINDOW_E=32400     # SCALE forecast time when the assimilation window ends (second)
                   #  (Set this value to be equal to $LCYCLE for 3D-LETKF)
LCYCLE=21600       # Length of a data assimilation cycle (second)
LTIMESLOT=3600     # Time slot interval for 4D-LETKF (second)
                   #  (Set this value to be equal to $LCYCLE for 3D-LETKF)

#===============================================================================
# Parallelization settings

MEMBER=10          # Ensemble size

NNODES=44          # Number of computing nodes used to run the job
PPN=1              # Number of processes per node (parallelized by MPI with $NNODES*$PPN processes)

THREADS=8          # Number of threads per MPI process (parallelized by OpenMP or automatic parallelization with this number of threads)

SCALE_NP=4         # Number of MPI processes to run a single member of the SCALE forecast

BGJOB_INT='0.1s'   ## [DEPRECATED] Interval of multiple background job submissions

#===============================================================================
# Use of the temporary runtime directories

DISK_ALL_RANK_LOCAL=0       # [USUALLY SET BY PRESET] Local disks for each MPI process in the same node are independent?
                            # - 0:  No
                            # - 1:  Yes - Choose this when using the rank-directory on the K computer
                            #  (default: 0)

DISK_MODE=3                 # Type of disk and the way of using the disk for general input/output files (except for the database, boundary, and observation files) in runtime
                            # - 1:  Use a shared disk; only create symbolic links of the input/output files to $TMP directory for computation
                            # - 2:  Use a shared disk; copy the input files to $TMP directory for computation and copy the output files back (staging)
                            # - 3:  Use local disks; copy the input files to $TMPL directory for computation and copy the output files back (staging)
                            #  (default: 2)
DISK_MODE_CONSTDB=2         # Type of disk and the way of using the disk for database (topo/landuse) files in runtime
                            # - 1-3:  Same meaning as for $DISK_MODE
                            #  (default: $DISK_MODE)
DISK_MODE_BDYDATA=2         # Type of disk and the way of using the disk for boundary data files in runtime
                            # - 1-3:  Same meaning as for $DISK_MODE
                            #  (default: $DISK_MODE)
DISK_MODE_OBS=2             # Type of disk and the way of using the disk for observation data files in runtime
                            # - 1-3:  Same meaning as for $DISK_MODE
                            #  (default: $DISK_MODE)

ONLINE_STGOUT=0             ## [NOT TESTED] Stage out right after each cycle (do not wait until the end of the job)?
                            ## - 0:  No
                            ## - 1:  Yes

IO_ARB=0                    ## [NOT COMPLETE] Use I/O arbitrator?
                            ## - 0:  No
                            ## - 1:  Yes

SYSNAME="$(basename $OUTDIR)"      # A string to be used as a part of the temporary directory name
TMPSUBDIR="scale-letkf_${SYSNAME}" # A part of the temporary directory name (should be unique in the machine to avoid any potential conflict when running multiple jobs at the same time)

TMP="/scratch/$(id -ng)/gylien/$TMPSUBDIR" # Temporary runtime directory path on a shared disk (available on both the head node and all computing nodes)
                                           #  (On the K computer, effective only for the 'micro' job)
TMPS="$DIR/tmp/$TMPSUBDIR"  # Temporary runtime directory path available only on the head node
TMPSL=                      # [OPTIONAL] Temporary runtime directory path available only on the head node and on a local file system
                            #  (default: $TMPS)
TMPL=                       # Temporary runtime directory path on local disks (available on all computing nodes; do not need to be available on the head node)
                            #  (On the K computer, do not need to set this)

CLEAR_TMP=0                 # Clear the temporary directories after the completion of the job?
                            # - 0:  No
                            # - 1:  Yes
                            #  (default: 1)

#===============================================================================
# Environmental settings

MPI_TYPE=                   # [USUALLY SET BY PRESET] Type of MPI software
                            # - 'K':  FUJITSU MPI on the K computer
                            # - 'openmpi'  Open MPI
                            # - 'impi'  Intel MPI
                            # - 'sgimpt':  SGI MPT

MPIRUN="mpiexec"            # Path of the 'mpiexec' or 'mpirun' program
if (which $MPIRUN > /dev/null 2>&1); then
  MPIRUN=$(which $MPIRUN)
fi

SCP='cp -L'                 # Command (on computing nodes) used for built-in file staging (to copy runtime files from input/to output directories)
SCP_HOSTPREFIX=''           # [OPTIONAL] When using scp-like commands for built-in file staging, the prefix for the hostname of the head node

#SCP="scp -q"
#SCP_HOSTPREFIX="XXXX:"

STAGE_THREAD=8              # Number of threads for built-in parallel file staging
TAR_THREAD=8                # Number of threads for parallel archiving log files after the completion of the job
                            #  (effective only when $LOG_TYPE > 3)

PYTHON="python"             # Command for the 'python' program

#BUFRBIN=

#===============================================================================
# Machine-independent source file

. config.rc

#===============================================================================
