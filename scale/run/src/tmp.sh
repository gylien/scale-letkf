#!/bin/bash
#===============================================================================
#
#  Script to prepare the SCALE model run directory.
#  August  2014  created,  Guo-Yuan Lien
#  October 2014  modified, Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage: $0 INIT BDY TOPO LANDUSE STIME FCSTLEN FCSTINT TMPDIR MODELDIR DATADIR
#
#    INIT      Basename of SCALE initial files
#    BDY       Basename of SCALE boundary files
#    TOPO      Basename of SCALE topography files
#    LANDUSE   Basename of SCALE land use files
#    STIME     Start time (format: YYYYMMDDHHMMSS)
#    FCSTLEN   Forecast length (second)
#    FCSTINT   Output interval (second)
#    TMPDIR    Temporary directory to run the model
#    MODELDIR  Directory of SCALE model executable files
#    DATADIR   Directory of SCALE data files
#
#===============================================================================

cd "$(dirname "$0")/.."
if [ ! -f configure.sh ]; then
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi
. configure.sh

#-------------------------------------------------------------------------------

INIT="$1"
BDY="$2"
TOPO="$3"
LANDUSE="$4"
STIME="$5"
FCSTLEN="$6"
FCSTINT="$7"
TMPDIR="$8"
MODELDIR="${9:-$MODELDIR}"
DATADIR="${10:-$DATADIR}"

#-------------------------------------------------------------------------------

scale_exec="${MODELDIR}/scale-les"

s_yyyy=${STIME:0:4}
s_mm=${STIME:4:2}
s_dd=${STIME:6:2}
s_hh=${STIME:8:2}
s_ii=${STIME:10:2}
s_ss=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*
cd $TMPDIR

ln -fs $scale_exec .

ln -fs $DATADIR/rad/PARAG.29 .
ln -fs $DATADIR/rad/PARAPC.29 .
ln -fs $DATADIR/rad/VARDATA.RM29 .
ln -fs $DATADIR/rad/cira.nc .
ln -fs $DATADIR/rad/MIPAS/day.atm .
ln -fs $DATADIR/rad/MIPAS/equ.atm .
ln -fs $DATADIR/rad/MIPAS/sum.atm .
ln -fs $DATADIR/rad/MIPAS/win.atm .

ipe=0

$mem_np

while [ -s "$INIT$(printf $SCALE_SFX $ipe)" ]; do
  ln -fs $INIT$(printf $SCALE_SFX $ipe) init$(printf $SCALE_SFX $ipe)
ipe=$((ipe+1))
done
  
ipe=0
while [ -s "$BDY$(printf $SCALE_SFX $ipe)" ]; do
  ln -fs $BDY$(printf $SCALE_SFX $ipe) boundary$(printf $SCALE_SFX $ipe)
ipe=$((ipe+1))
done

ipe=0
ipef=`printf '%06d' $ipe`
while [ -s "${TOPO}.pe${ipef}.nc" ]; do
  ln -fs ${TOPO}.pe${ipef}.nc topo.pe${ipef}.nc
ipe=$((ipe+1))
ipef=`printf '%06d' $ipe`
done

ipe=0
ipef=`printf '%06d' $ipe`
while [ -s "${LANDUSE}.pe${ipef}.nc" ]; do
  ln -fs ${LANDUSE}.pe${ipef}.nc landuse.pe${ipef}.nc
ipe=$((ipe+1))
ipef=`printf '%06d' $ipe`
done

#===============================================================================

cat << EOF > run.conf
&PARAM_IO
 IO_LOG_ALLNODE  = .false.,
/

&PARAM_PRC
 PRC_NUM_X       = 6,
 PRC_NUM_Y       = 6,
 PRC_PERIODIC_X = .false.,
 PRC_PERIODIC_Y = .false.,
/

&PARAM_INDEX
 KMAX = 30,
 IMAX = 16,
 JMAX = 16,
/

&PARAM_GRID
 DX        = 2000.D0,
 DY        = 2000.D0,
 FZ(:) =     120.00D0,   240.00D0,   360.00D0,   480.00D0,   600.00D0,
             720.00D0,   840.00D0,   960.00D0,  1074.90D0,  1198.70D0,
            1345.40D0,  1529.00D0,  1763.50D0,  2062.90D0,  2441.20D0,
            2912.40D0,  3490.50D0,  4189.50D0,  5023.40D0,  6006.20D0,
            7000.00D0,  8000.00D0,  9000.00D0, 10000.00D0, 11000.00D0,
           12000.00D0, 13000.00D0, 14000.00D0, 15000.00D0, 16000.00D0,
 BUFFER_DZ =  5000.D0,
 BUFFER_DX =  12000.D0,
 BUFFER_DY =  12000.D0,
 BUFFFACT  =  1.00D0,
/

&PARAM_LAND_INDEX
 LKMAX = 5,
/

&PARAM_URBAN_INDEX
 UKMAX = 4,
/

&PARAM_LAND_GRID
 LDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0, 1.00D0,
/

&PARAM_URBAN_GRID
 UDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0,
/

&PARAM_TIME
 TIME_STARTDATE             = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = ${FCSTLEN}.D0,
 TIME_DURATION_UNIT         = "SEC",
 TIME_DT                    = 1.0D0,
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = 0.5D0,
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
 TIME_DT_ATMOS_PHY_TB       = 0.5D0,
 TIME_DT_ATMOS_PHY_TB_UNIT  = "SEC",
 TIME_DT_ATMOS_PHY_SF       = 2.0D0,
 TIME_DT_ATMOS_PHY_SF_UNIT  = "SEC",
 TIME_DT_ATMOS_PHY_MP       = 2.0D0,
 TIME_DT_ATMOS_PHY_MP_UNIT  = "SEC",
 TIME_DT_ATMOS_PHY_RD       = 300.0D0,
 TIME_DT_ATMOS_PHY_RD_UNIT  = "SEC",
 !TIME_DT_OCEAN              = 5.D0,
 !TIME_DT_OCEAN_UNIT         = "SEC",
 !TIME_DT_LAND               = 5.D0,
 !TIME_DT_LAND_UNIT          = "SEC",
 !TIME_DT_URBAN              = 5.D0,
 !TIME_DT_URBAN_UNIT         = "SEC",
/

&PARAM_RESTART
 RESTART_OUTPUT       = .true.,
 RESTART_IN_BASENAME  = "init",
 RESTART_OUT_BASENAME = "restart",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_TOPO
 TOPO_IN_BASENAME = "topo",
/

&PARAM_LANDUSE
 LANDUSE_IN_BASENAME = "landuse",
/

&PARAM_MAPPROJ
 MPRJ_basepoint_lon = 134.6D0,
 MPRJ_basepoint_lat =  34.0D0,
/

&PARAM_TRACER
 TRACER_TYPE = 'TOMITA08',
/

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "HEVI",
 ATMOS_PHY_TB_TYPE = "SMAGORINSKY",  !SMAGORINSKY
 ATMOS_PHY_SF_TYPE = "COUPLE",       !COUPLE
 ATMOS_PHY_MP_TYPE = "TOMITA08",     !TOMITA08 or SN14
 ATMOS_PHY_RD_TYPE = "MSTRNX",       !MSTRNX
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_IN_BASENAME  = "init",
 ATMOS_RESTART_OUT_BASENAME = "restart",
 ATMOS_RESTART_OUTPUT       = .true.,
 ATMOS_VARS_CHECKRANGE      = .true.,
/

&PARAM_ATMOS_PHY_MP_VARS
 ATMOS_PHY_MP_RESTART_IN_BASENAME  = "init",
 ATMOS_PHY_MP_RESTART_OUT_BASENAME = "restart",
 ATMOS_PHY_MP_RESTART_OUTPUT       = .true.,
/

&PARAM_ATMOS_PHY_RD_VARS
 ATMOS_PHY_RD_RESTART_IN_BASENAME  = "init",
 ATMOS_PHY_RD_RESTART_OUT_BASENAME = "restart",
 ATMOS_PHY_RD_RESTART_OUTPUT       = .true.,
/

&PARAM_ATMOS_PHY_SF_VARS
 ATMOS_PHY_SF_RESTART_IN_BASENAME  = "init",
 ATMOS_PHY_SF_RESTART_OUT_BASENAME = "restart",
 ATMOS_PHY_SF_RESTART_OUTPUT       = .true.,
/

!&PARAM_ATMOS_PHY_TB_VARS
! ATMOS_PHY_TB_RESTART_IN_BASENAME  = "init",
! ATMOS_PHY_TB_RESTART_OUT_BASENAME = "restart",
! ATMOS_PHY_TB_RESTART_OUTPUT       = .false.,
!/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE       = "INIT",
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE         = "REAL",
 ATMOS_BOUNDARY_IN_BASENAME  = 'boundary',
! ATMOS_BOUNDARY_OUT_BASENAME = 'boundary_check',
 ATMOS_BOUNDARY_USE_VELZ     = .true.,
 ATMOS_BOUNDARY_USE_VELX     = .true.,
 ATMOS_BOUNDARY_USE_VELY     = .true.,
! ATMOS_BOUNDARY_USE_POTT     = .true.,
 ATMOS_BOUNDARY_USE_DENS     = .true.,
 ATMOS_BOUNDARY_USE_QV       = .true.,
 ATMOS_BOUNDARY_FRACZ        = 1.0D0,
 ATMOS_BOUNDARY_FRACX        = 1.0D0,
 ATMOS_BOUNDARY_FRACY        = 1.0D0,
 ATMOS_BOUNDARY_TAUZ         = 10.D0,
 ATMOS_BOUNDARY_TAUX         = 10.D0,
 ATMOS_BOUNDARY_TAUY         = 10.D0,
 ATMOS_BOUNDARY_UPDATE_DT    = 600.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF_COEF = 1.D-4,
/

&PARAM_ATMOS_PHY_RD_MSTRN
 ATMOS_PHY_RD_MSTRN_KADD                  = 30,
 ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = "PARAG.29",
 ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = "PARAPC.29",
 ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = "VARDATA.RM29",
 ATMOS_PHY_RD_MSTRN_NBAND                 = 29,
/

&PARAM_ATMOS_PHY_RD_PROFILE
 ATMOS_PHY_RD_PROFILE_TOA                   = 100.D0,
 ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME    = "cira.nc",
 ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = ".",
 DEBUG                                      = .true.,
/

!&PARAM_ATMOS_PHY_SF_VARS
! ATMOS_PHY_SF_DEFAULT_SFC_TEMP   = 300.D0,
! ATMOS_PHY_SF_DEFAULT_SFC_albedo = 0.04D0,
!/

&PARAM_OCEAN
 OCEAN_TYPE = "SLAB",
/

&PARAM_OCEAN_VARS
 OCEAN_RESTART_IN_BASENAME  = "init",
 OCEAN_RESTART_OUT_BASENAME = "restart",
 OCEAN_RESTART_OUTPUT       = .true.,
 OCEAN_VARS_CHECKRANGE      = .true.,
/

&PARAM_OCEAN_SLAB
 OCEAN_PHY_SLAB_DEPTH = 10.D0,
/

&PARAM_OCEAN_ROUGHNESS
 OCEAN_ROUGHNESS_TYPE = 'MOON07'
/

&PARAM_LAND
 LAND_TYPE = "BUCKET",
/

&PARAM_LAND_VARS
 LAND_RESTART_IN_BASENAME  = "init",
 LAND_RESTART_OUT_BASENAME = "restart",
 LAND_RESTART_OUTPUT       = .true.,
 LAND_VARS_CHECKRANGE      = .true.,
/

&PARAM_LAND_BUCKET
 LAND_PHY_UPDATE_BOTTOM_WATER = .true.,
 LAND_PHY_UPDATE_BOTTOM_TEMP  = .false.,
/

&PARAM_URBAN
 URBAN_TYPE = "UCM",
/

&PARAM_URBAN_VARS
 URBAN_RESTART_IN_BASENAME  = "init",
 URBAN_RESTART_OUT_BASENAME = "restart",
 URBAN_RESTART_OUTPUT       = .true.,
 URBAN_VARS_CHECKRANGE      = .true.,
/

&PARAM_CPL
 CPL_TYPE_AtmOcn = "BULK",
 CPL_TYPE_AtmLnd = "BULK",
 CPL_TYPE_AtmUrb = "BULK",
/

&PARAM_CPL_BULKCOEF
 CPL_BULKCOEF_TYPE = 'BH91',
/

&PARAM_CPL_ATMURB_BULK
 URBAN_UCM_STRGR = 0.24D0,
! URBAN_UCM_STRGB = 0.009D0,
! URBAN_UCM_STRGG = 0.24D0,
/

&PARAM_HISTORY
 HISTORY_DEFAULT_BASENAME  = "history",
 HISTORY_DEFAULT_TINTERVAL = ${FCSTINT}.D0,
 HISTORY_DEFAULT_TUNIT     = "SEC",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL4",
 HISTORY_DEFAULT_ZINTERP   = .true.,
/

&HISTITEM item='T'    /
&HISTITEM item='PRES' /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='PT'   /
&HISTITEM item='RH'   /
&HISTITEM item='DENS' /

&HISTITEM item='TKE'  /

&HISTITEM item='QV'   /
&HISTITEM item='QC'   /
&HISTITEM item='QI'   /
&HISTITEM item='QR'   /
&HISTITEM item='QS'   /
&HISTITEM item='QG'   /
&HISTITEM item='QHYD' /

&HISTITEM item='PREC' /
&HISTITEM item='RAIN' /
&HISTITEM item='SNOW' /

&HISTITEM item='OLR'  /
&HISTITEM item='SLR'  /
&HISTITEM item='OSR'  /
&HISTITEM item='SSR'  /

&HISTITEM item='U10'  /
&HISTITEM item='V10'  /
&HISTITEM item='T2'   /
&HISTITEM item='Q2'   /

&HISTITEM item='SFC_TEMP'   /
&HISTITEM item='OCEAN_TEMP' /
&HISTITEM item='LAND_TEMP'  /
&HISTITEM item='LAND_WATER' /

&PARAM_MONITOR
 MONITOR_STEP_INTERVAL = 90,
/

&MONITITEM item='QDRY' /
&MONITITEM item='QTOT' /
&MONITITEM item='ENGT' /
&MONITITEM item='ENGP' /
&MONITITEM item='ENGK' /
&MONITITEM item='ENGI' /

&PARAM_LAND_DATA
 index       = 1,
 description = "bare ground",
 STRGMAX     =  0.20D0,
 STRGCRT     =  0.15D0,
 TCS         =  0.25D0,
 HCS         =  1.30D+6,
 DFW         =  3.38D-6,
 Z0M         =  0.01D0,
/

&PARAM_LAND_DATA
 index       = 2,
 description = "grassland",
 STRGMAX     =  0.20D0,
 STRGCRT     =  0.10D0,
 TCS         =  0.25D0,
 HCS         =  1.30D+6,
 DFW         =  3.38D-6,
 Z0M         =  0.10D0,
/

&PARAM_LAND_DATA
 index       = 3,
 description = "deciduous forest",
 STRGMAX     =  0.20D0,
 STRGCRT     =  0.05D0,
 TCS         =  0.25D0,
 HCS         =  1.30D+6,
 DFW         =  3.38D-6,
 Z0M         =  0.30D0,
/

&PARAM_LAND_DATA
 index       = 4,
 description = "paddy",
 STRGMAX     =  0.30D0,
 STRGCRT     =  0.03D0,
 TCS         =  0.25D0,
 HCS         =  2.00D+6,
 DFW         =  3.38D-6,
 Z0M         =  0.01D0,
/
EOF

#===============================================================================

if [ "$MPIEXEC" != 'no' ]; then
  if [ -s "$MACHINE" ]; then
    $MPIEXEC -n $NP -machinefile $MACHINE $SCALEEXEC run.conf > scale.log 2>&1
  else
    $MPIEXEC -n $NP $SCALEEXEC run.conf > scale.log 2>&1
  fi
fi

#===============================================================================

exit 0
