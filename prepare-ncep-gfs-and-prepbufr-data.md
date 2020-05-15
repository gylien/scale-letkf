# Prepare NCEP GFS and PREPBUFR data

## NCEP GFS and PREPBUFR

NCEP GFS and PREPBUFR data sets are currently used for SCALE-LETKF domain 1 as boundary conditions and observation. 
6-hourly automatic download is necessary for realtime operation (see crontab_smp.txt).

Original data is available at NCEP NOMAD ([https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod)) in wgrib2 format. 
No user registration is required.

For realtime continuous operation, be sure to keep up with the information:[NWS - Service change notice](https://www.weather.gov/notification/#scn)


### NCEP GFS model data

*(copyed from [Guo-Yuan Lien's description](https://github.com/aamemiya/scale-letkf-rt-interface/blob/previous/scale-letkf-rt/README.md))*

In `ncepgfs` directory, create a `mtime` file indicating the date and time (in UTC) of the completely downloaded data 
(can be fake: the data do not need to be actually there); e.g.,
~~~
echo "2018-06-30 18" > mtime
~~~
The next run of `get_ncep_gfs.sh` script will try to download the next 6-h data. If the download process is successful, 
the script will modify the `mtime` file to be the date and time of the newly downloaded data.

There are several data processing scripts called in `get_ncep_gfs.sh` after downloading the data:

<!--
- **`run/wps/convert.sh` and `run/wrf/convert.sh`: Conversion from the NCEP GFS grid2-format files to the NetCDF-format 
  WRF files**

    To conduct this step, copy a complied WRF Preprocessing System (WPS) code into `ncepgfs/run/wps/`, and copy the run 
    directory (`WRFV3/run`) of a compiled WRF model into `ncepgfs/run/wrf/`. Otherwise, comment out to disable these 
    `run/wps/convert.sh` and `run/wrf/convert.sh` execution sections. *This step is not necessary if using GrADS-format 
    input for the SCALE-LETKF.*
    
    The WRF-format data are saved into `<DATA_DIR_FOR_SCALE-LETKF>/ncepgfs_wrf/`.
-->

- **`run/grads/convert.sh`: Generate GrADS ctl file for the NCEP grib2-format files**

    This step just runs `g2ctl` and `gribmap` to generate the GrADS ctl files for plotting, without generating new data 
    files.

- **`run/plot/plot.sh`: Plot the GFS analysis data**

    The graphic files are saved into `YYYYMMDDHH/plot/`.    
<!--
; in addition, they are also sent to a web server. Edit 
    `web_host` and `web_remote_dir` in `run/plot/plot.sh` for the host name and destination directory of the web server.
-->
- **`run/grads_scale/convert.sh`: Conversion from the NCEP GFS grib2-format files to the GrADS-format files for 
  SCALE-LETKF input**

    This step extracts a certain required variables from the NCEP GFS data and write them into GrADS-format binary 
    files to be used for SCALE-LETKF input. <!--*This step is not necessary if using WRF-format input for the SCALE-LETKF.*-->

    *See the resent change in `convert_prmsl.sh` due to [the GFS format change](https://www.weather.gov/media/notification/scn19-40gfs_v15_1.pdf) from 2019.6.12. See also [https://gitlab.com/scale-met/scale/issues/1297](https://gitlab.com/scale-met/scale/issues/1297). 

    The GrADS-format data are saved into `<DATA_DIR_FOR_SCALE-LETKF>/data/ncepgfs_grads/`.

### PREPBUFR decoder

In `ncepobs_gdas` directory, the `get_ncep_obs.sh` script is used to download the NCEP GDAS PREPBUFR data. These data 
are very similar to the NCEP GFS PREPBUFR data, but they are available later than the GFS observation data, and have 
more complete data within the data assimilation window. We can choose to use either the GFS or GDAS observation data to
run the SCALE-LETKF near-real-time system. *GDAS is used in the current SCALE-LETKF system.*


A data processing script is called in `get_ncep_gfs.sh` after downloading the data:

- **`run/dec_prepbufr/convert.sh`: Decode BUFR data to generate a binary file for SCALE-LETKF**

    The binary data is saved into `<DATA_DIR_FOR_SCALE-LETKF>/data/ncepobs_gdas/`.

    This script is dependent on following two decoder tools.
1. [NCEP BUFRLIB](https://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/)
2. fortran program dec_prepbufr (source available in [scale-letkf/scale/obs](https://github.com/gylien/scale-letkf/tree/master/scale/obs))

*The crutial bug fix was performed in 2019.3. The error in decoding virtual temperature was fixed. See the detail [here](https://github.com/gylien/scale-letkf/commit/3ad0aba7169037bcb4c9cf0e4c7fe17fe3778f94#diff-972592ed1fe45ba69df5ed75317fc650R173). 


- **`run/plot/map`: generate a map of PREPBUFR observation distribution** (optional)

## Automatic download

As NOMAD provides data of only latest two weeks, it is necessary to keep updating NCEP data necessary for boundary conditions and data assimilation. 

Normally NCEP GFS data is updated around 4,10,16,22 UTC (13,19,1,7 JST) (about 4 hours delay from the reality). PREPBUFR data comes later, around 6, 12, 18, 0 UTC (15,21,3,9 JST) (about 6 hours delay). 
The cron job for these scripts can be like:

~~~
0 2,8,14,20 * * * <PATH_TO_SCRIPTS>/data/ncepgfs/get_ncep_gfs.sh &> get_ncep_gfs.log &
0 4,10,16,22 * * * <PATH_TO_SCRIPTS>/data/ncepobs_gdas/get_ncep_obs.sh &> get_ncep_obs.log &
~~~
