# Getting started

## Get the source

```
cd (directory-to-place-the-system) 
git clone https://github.com/gylien/scale-letkf.git -b realtime-interface SCALE-LETKF-rt
cd SCALE-LETKF-rt
git clone https://gitlab.com/scale-met/scale.git -b LETKF-5.3.4 scale_ope
cd scale_ope
git clone https://github.com/gylien/scale-letkf.git -b realtime_fcst_D1-3_ope scale-letkf_ope 
git clone https://github.com/gylien/scale-letkf.git -b dacycle_OFP_ope  scale-letkf_ope_d4
```
See [directory structure](directory-structure.md) for the entire structure of the system.

## Build SCALE-LETKF

(See also [usage-on-OFP](usage-on-OFP.md))

### Build SCALE-RM

exec files for both single-precision and double-precision need to be created. 

```
cd scale_ope
cp scale-letkf_ope/Mkinclude_scale ./Mkinclude   ### Special config to create multiple bin/lib/include sets
export SCALE_SYS=OFP
export SCALE_ENABLE_OPENMP=T 
export SCALE_USE_SINGLEFP=T ### single precision 
cd scale-rm/src
make 
make allclean
export SCALE_USE_SINGLEFP=F ### double precision
make
cd -
```
This makes pairs of bin/lib/include directories respectively for single and double precision execusion.

### Build LETKF

Two separate LETKF source codes respectively for D1-3 and D4(realtime) are needed to be compiled.

If you need to perform simulation with both single and double precision, you need to compile twice with `$SCALE_USE_SINGLEFP=T` and `F`.

By default, D1 DA cycle is performed with double precision and D1-3 ensemble forecasts are performed with single precision. 

```
cd scale-letkf_ope/scale
ln -s arch/configure.user.OFP ./configure.user
export SCALE_USE_SINGLEFP=T ### single precision 
make 
export SCALE_USE_SINGLEFP=F ### double precision 
make 
cd -
cd scale-letkf_ope_d4/scale
ln -s arch/configure.user.OFP ./configure.user
export SCALE_USE_SINGLEFP=T ### single precision 
make 
```
The one for D1-3 creates exec files as follows. If `$SCALE_USE_SINGLEFP=T`, the prefix `_single` is added to each file name.  
```
ensmodel/scale-rm_ens
ensmodel/scale-rm_init_ens
ensmodel/scale-rm_pp_ens
```
The one for D4 creats a combined exec file as follows.
```
dacycle/dacycle
```

## Prepare initial data

Initial background data and ensemble perturbations are required to launch DA cycle for Domain 1.
See [initialization](initialize.md).

## Prepare boundary and observation data

Boundary and observation data need to be downloaded from NCEP websites. The scripts in `$TOPDIR/external/ncepgfs` and `$TOPDIR/external/ncepobs_gdas` do that job. To start the realtime operation, the scripts need to be automatically executed every 6 hours. 

See [Prepare NCEP GFS and PREPBUFR data](prepare-ncep-gfs-and-prepbufr-data.md) for the detail. 

## Prepare topography and landuse files 

Before you start a DA cycle or a forecast, you also need topography and landuse files using the same grid and MPI domain-decomposition with the forecast domain. A file `latlon_domain_catalogue.txt` is also needed for an offline domain nesting.

See [create topography and landuse files](create-topography-and-landuse-files.md).

## Prepare configure files for physical processes 
Physical processes in SCALE-RM such as radiation and ground soil model need parameter table files. 
You don't need to prepare anything unless you need to change the default settings of physical processes.
The necessary files are copied from the default database in `scale-rm/test`. 

## Here we go
Go to `admin` and [start DA cycle for domain 1](control.md).  

