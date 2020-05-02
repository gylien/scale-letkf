## Create topography and landuse files

To create topograpy and landuse files, use scripts in the directory `run_pp`.

### Configuration

To generate new topo and landuse files, set `config.main` as follows (D3 case for example). 
```
DOMAIN=d3
NX=192 ### Number of grid
NY=192 ### Number of grid
PRC_NUM_X=16 ### Number of process (subdomain) in X direction
PRC_NUM_Y=16 ### Number of process (subdomain) in Y direction 
DATA_TOPO_BDY_SCALE="${TOPDIR}/result/ope/d2/const/topo"   ### Parent topo files (for domain nesting) 
DATA_BDY_SCALE_PRC_NUM_X=16 ### Parent number of process  
DATA_BDY_SCALE_PRC_NUM_Y=16 ### 
```

Source data files for topography and landuse are specified in `config.main` as follows. 
DEM50M and LU100M may be better for D4 as they have finer resolution. 

```
TOPO_FORMAT='GTOPO30'     # 'prep': Use prepared topo files in $DATA_TOPO
                        # 'GTOPO30' (requires compatible 'config.nml.scale_pp')
                        # 'DEM50M'  (requires compatible 'config.nml.scale_pp')
LANDUSE_FORMAT='GLCCv2'   # 'prep': Use prepared landuse files in $DATA_LANDUSE
                        # 'GLCCv2' (requires compatible 'config.nml.scale_pp')
                        # 'LU100M' (requires compatible 'config.nml.scale_pp')
```

Execute `pp_ofp.sh` to generate topo and landuse files (and `const/log/latlon_domain_catalogue.txt`). 
(The name of a log file should not be changed, as it is used for job end detection. )

The job takes only a few nodes and a few minutes. If not specified, `debug-flat` resource group is used. 
```
./pp_ofp.sh &> pp_ofp.log &
```
By default, the destination of output files are `**your SCALE-LETKF top directory**/result/test_const/const`.
Check the validity of files and move the directory `const` to your target directory. 
