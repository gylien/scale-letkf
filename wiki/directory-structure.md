# Directory structure
## Overview

```
$TOPDIR/                                 ### top directory for SCALE-LETKF realtime system (this branch)
├── admin/                                # control scripts
│   └── monitor/                          # monitoring scripts
│       ├── letkf/                        #   -- for LETKF
│       └── realtime/                     #   -- for realtime operation
├── external/                            ### external data
│   ├── ncepgfs/                          # NCEP GFS grib2 raw data and download/convert/plot scripts 
│   ├── ncepgfs_grads/                    #       -- grads format data
│   ├── ncepobs_gdas/                     # NCEP PREPBUFR raw data and download/decode scripts
│   └── ncepobs_gdas_letkf/               #            -- decoded data
│
├── scale_ope/                           ### SCALE model
│   ├── bin/                              # binary  (double precision)
│   ├── bin_single/                       #         (single)
│   ├── include/                          # include (double precision)
│   ├── include_single/                   #         (single)
│   ├── lib/                              # library (double precision)
│   ├── lib_single/                       #         (single)
│   ├── scale-rm/                         # SCALE-RM main source files  
|   ├── scalelib/                         # SCALE common library source files
|   ├── sysdep/                           # system-dependent configuration
│   ├── utils/                             
|   |
│   ├── scale-letkf_ope/                 ### SCALE-LETKF for domain 1-3 (branch 'realtime_fcst_D1-3_ope')  
│   |   ├── common/                        
│   |   └── scale/                         
│   |       ├── arch/                      
│   |       ├── common/                    
│   |       ├── ensmodel/                 # main programs 
│   |       ├── letkf/                     
│   |       ├── obs/                       
│   |       ├── run/                      # run scripts for D1 cycle
│   |       ├── run_d1_ext/               #                 D1 extended forecast
│   |       ├── run_d1-2/                 #                 D1 and D2 forecast (online nesting)
│   |       ├── run_d3/                   #                 D3 forecast
│   |       ├── run_d4_init/              #                 D4 initial and boundary conditions
│   |       └── run_pp/                   # create topo and landuse files
|   |
│   └── scale-letkf_ope_d4/              ### SCALE-LETKF for PAWR DA in domain 4 (branch 'dacycle_OFP_ope')  
│       ├── common/                        
│       └── scale/                         
│           ├── arch/                      
│           ├── common/                    
│           ├── dacycle/                  # main program of PAWR-SCALE-LETKF using direct transfer
│           ├── ensmodel/                 # not used
│           ├── letkf/                    
│           ├── obs/                       
│           │   └── radar/                # including realtime processing of PAWR observation  
│           └── run/                      # run scripts for D4 realtime DA cycle and extended forecasts
│ 
└── result/                              ### Results
    ├── ope/                              # Double precision
    │   ├── d1/                           # Domain 1 analysis
    │   ├── ncepgfs_grads_da/             #   temporal paths used in D1 cycle
    │   └── ncepobs/                      #   temporal paths used in D1 cycle
    └── ope_single/                       # Single precision
        ├── d1/                           # Domain 1 forecast (online)
        ├── d1_ext/                       #          extended forecast 
        ├── d2/                           # Domain 2 forecast (online)
        ├── d3/                           # Domain 3 forecast Kansai area
        ├── d3_Kanto/                     #                   Kanto area
        ├── d4_1km/                       # Domain 4 analysis/forecast (test)
        ├── d4_500m/                      #                            (ope)
        ├── d4_250m/                      #                            (champaign)
        ├── ncepgfs_grads/                #   temporal paths used in D1-2 forecast 
        └── ncepgfs_grads_ext/            #   temporal paths used in D1 extended forecast

```


## Output directories
(edited from [https://github.com/gylien/scale-letkf/wiki/SCALE-LETKF-experiment-directory](https://github.com/gylien/scale-letkf/wiki/SCALE-LETKF-experiment-directory))

For D1 and D2 : such as `$OUTDIR = $TOPDIR/result/ope/d1`

```
$OUTDIR/                                  # SCALE-LETKF output directory
├── YYYYMMDDhhmmss/                       # date and time;
│   │                                     # for 'fcst' jobs: forecast initial time
│   ├── gues/                             # first-guess files (for 'cycle' jobs only)
│   │   ├── <member>/                     # 4-digit ensemble member number
│   │   │   └── init.pe<process>.nc
│   │   ├── mean/                         # ensemble mean
│   │   │   └── init.pe<process>.nc
│   │   ├── mdet/                         # deterministic run
│   │   │   └── init.pe<process>.nc
│   │   └── sprd/                         # ensemble spread
│   │       └── init.pe<process>.nc
│   ├── anal/                             # analysis files
│   │   │                                 # (used for initial conditions)
│   │   ├── <member>/                     # 4-digit ensemble member number
│   │   │   ├── init.pe<process>.nc       # analysis files
│   │   │   ├── init_ocean.pe<process>.nc # initial conditions for ocean data
│   │   │   │                             # when $OCEAN_INPUT=1 and $OCEAN_FORMAT=0
│   │   │   └── init_land.pe<process>.nc  # initial conditions for land data
│   │   │                                 # when $LAND_INPUT=1 and $LAND_FORMAT=0
│   │   ├── mean/                         # ensemble mean
│   │   │   ├── init.pe<process>.nc       #
│   │   │   ├── init_ocean.pe<process>.nc # (same as for ensemble members)
│   │   │   └── init_land.pe<process>.nc  #
│   │   ├── mdet/                         # deterministic run
│   │   │   ├── init.pe<process>.nc       #
│   │   │   ├── init_ocean.pe<process>.nc # (same as for ensemble members)
│   │   │   └── init_land.pe<process>.nc  #
│   │   └── sprd/                         # ensemble spread
│   │       └── init.pe<process>.nc
│   ├── hist/                             # history files in cycling DA runs saved
│   │   │                                 # for future offline-nesting run
│   │   │                                 # (for 'cycle' jobs only)
│   │   ├── <member>/                     # 4-digit ensemble member number
│   │   │   └── history.pe<process>.nc
│   │   ├── mean/                         # ensemble mean
│   │   │   └── history.pe<process>.nc
│   │   └── mdet/                         # deterministic run
│   │       └── history.pe<process>.nc
│   ├── fcst/                             # extended forecasts (for 'fcst' jobs only)
│   │   ├── <member>/                              # 4-digit ensemble member number
│   │   │   ├── history.pe<process>.nc             # forecast history files
│   │   │   └── init_YYYYMMDDHHmmss.pe<process>.nc # restart files in the end of the
│   │   │                                          # forecasts;
│   │   │                                          # YYYYMMDDHHmmss: forecast end time
│   │   ├── mean/                                  # ensemble mean
│   │   │   ├── history.pe<process>.nc             # (same as for ensemble members)
│   │   │   └── init_YYYYMMDDHHmmss.pe<process>.nc #
│   │   └── mdet/                                  # deterministic run
│   │       ├── history.pe<process>.nc             # (same as for ensemble members)
│   │       └── init_YYYYMMDDHHmmss.pe<process>.nc #

│   ├── fcstgp/                             # forecast grads files
│   ├── fcstgpi/                           # plot files 

│   ├── bdy/                              # boundary files saved for repeated use
│   │   ├── <member>/                     # 4-digit ensemble member number
│   │   │   └── boundary.pe<process>.nc
│   │   ├── mean/                         # ensemble mean
│   │   │   └── boundary.pe<process>.nc
│   │   └── mdet/                         # deterministic run
│   │       └── boundary.pe<process>.nc
│   ├── log/                              # log files (may be archived based on 
│   │   │                                 #            $LOG_TYPE setting)
│   │   ├── scale_pp/                     # for 'scale_pp' steps in 'cycle' jobs
│   │   ├── scale_init/                   # for 'scale_init' steps in 'cycle' jobs
│   │   ├── scale/                        # for 'scale' steps in 'cycle' jobs
│   │   ├── obsope/                       # for 'obsope' steps in 'cycle' jobs
│   │   ├── letkf/                        # for 'letkf' steps in 'cycle' jobs
│   │   ├── fcst_scale_pp/                # for 'scale_pp' steps in 'fcst' jobs
│   │   ├── fcst_scale_init/              # for 'scale_init' steps in 'fcst' jobs
│   │   └── fcst_scale/                   # for 'scale' steps in 'fcst' jobs
│   └── obs/                              # observation-space diagnostic files
│       └── obsdep.dat                    # observation departures
└── const/                                # time-invariant files
    ├── topo/                             # topography files
    │   └── topo.pe<process>.nc
    ├── landuse/                          # land-use files
    │   └── landuse.pe<process>.nc
    └── log/
        └── latlon_domain_catalogue.txt   # SCALE "domain catalogue" files
                                          # (for offline nesting)
```


For D3 : `$OUTDIR = $TOPDIR/result/ope_single/d3`
- identified by the initial times of its own and of a forecast in the parent domain (D2) for init/boundary data
- downscaled init/boundary data for D4 is stored in subdirectories

```
$OUTDIR/                                               ### SCALE-LETKF output directory
├── ref_YYYYMMDDhhmmss/                                 # PARENT DOMAIN (D2) forecast init time 
│   └── YYYYMMDDhhmmss/                                 # forecast init time
│       ├── anal/                                       # analysis files
│       │   └── <member>/                               # 4-digit ensemble member number and/or mean/mdet
│       │      └── init.pe<process>.nc                  # analysis files
│       ├── fcst/                                       # extended forecasts (for 'fcst' jobs only)
│       │   └── <member>/                               # 4-digit ensemble member number and/or mean/mdet
│       │       ├── history.pe<process>.nc              # forecast history files
│       │       └── init_YYYYMMDDHHmmss.pe<process>.nc  # restart files
│       ├── log/                                        # log files
│       │   ├── fcst_scale_init/                        # for 'scale_init' steps in 'fcst' jobs
│       │   └── fcst_scale/                             # for 'scale' steps in 'fcst' jobs
│       ├── plot/                                       # quicklook png files 
│       │   └── <member>/                               #   4-digit ensemble member number and/or mean/mdet
│       │       └── xxxxxx.png                           
│       └── d4_XXX/                                    ### Init/bdy data for Domain 4 (XXX = 1km/500m/250m)
│           ├── anal/                                   # Initial files with timestamp
│           │   └── <member>/                           # 4-digit ensemble member number and/or mean/mdet
│           │       └── init_YYYYMMDD-HHmmss.000.pe<process>.nc   
│           ├── bdy/                                    # Boundary files with timestamp
│           │   └── <member>/                              # 4-digit ensemble member number and/or mean/mdet
│           │       └── boundary_YYYYMMDD-HHmmss.000.pe<process>.nc   
│           ├── exp/                                    # runtime log files archive 
│           └── log/                                    # log files of creating D4 init/bdy archive
└── const/                                              # time-invariant files
    ├── topo/                                           # topography files
    │   └── topo.pe<process>.nc
    ├── landuse/                                        # land-use files
    │   └── landuse.pe<process>.nc
    └── log/
        └── latlon_domain_catalogue.txt                 # SCALE "domain catalogue" files
                                                        # (for offline nesting)
```

