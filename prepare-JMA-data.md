
## JMA data (optional)

JMA model forecast and precipitation observation and nowcasting data are avaiable and can be used for realtime comparison. 

### MSM forecast data 

MSM is the operational regional mesoscale model developed by JMA.
The model configuration is summarized [here(pdf file)](http://www.jma.go.jp/jma/kishou/books/nwptext/51/appendix_A.pdf).
The state of the art of the JMA models is descibed [here](http://www.jma.go.jp/jma/kishou/books/nwptext/51/No51.pdf) in detail.

MSM forecast data can be downloaded, converted and plotted automatically by `get_msm.sh` in the same way as NCEP GFS. 
The original data is downloaded from [Kyoto Univ. RISH website](http://database.rish.kyoto-u.ac.jp/arch/jmadata/) in NetCDF format and converted to GrADS. 
Surface data such as precipitaion, 10m-level wind, surface temprature is avaiable hourly. 

The original NetCDF and converted GrADS-format data is saved into `<DATA_DIR_FOR_SCALE-LETKF>/data/msm/`.

### Precipitation data

**Since the beginning of FY2020, the account of Japan Meteorological Businesss Support Center has been expired. JMA precipitation data is available through [Meteorology Research Consotium](https://www.mri-jma.go.jp/Project/cons/index.html) with a few hours delay.**

JMA provides a couple of precipitation observation/forecast data sets through [Japan Meteorological Business Support Center](http://www.jmbsc.or.jp/jp/).  
For verifications of precipitation forecast in different timescales, it is benefitable to use two different JMA datasets below. 

 - **analysis / forecast**  
The JMA precipitation analysis/forecast data is suitable for a reference for precipitation **in hourly timescale.** 
The analysis is obtained from multiple observation including radar and ground-based observation network. The forecast is the combination of nowcasting based on motion vector and a 1-km resolution NWP model forecast (LFM). 
This data corresponds to [this page] (https://www.jma.go.jp/jp/kaikotan/index.html) in the JMA website.

| analysis/forecast |  |  
| --- | --- |  
| Spatial resolution | 1 km |  
| Temporal resolution | *30 min |   
| Forecast length | 6 hours |  
| Refresh | *30 min |  

The data is automatically downloaded in `/jmbsc/jmbsc_anal_1km/`.   
The scripts `data/JMA_precip/draw/grads/plot_{anal,fcst}.sh` make plots of analysis/forecast precipitation of input time in png format. 

*currently the data and plot are updated hourly. 

The current plot is avaibable at [realtime fcst_d2 quicklook](http://daweb.r-ccs27.riken.jp/~amemiya/scale/fcst_d2.php).

- **radar / nowcast**  
The JMA precipitation radar/nowcast data is suitable for a reference for localized precipitation in **a timescale shorter than an hour.** 
The radar product is made from multiple JMA radar data and QC-ed. The data has been converted to precipitation rate (mm/h) unit.

This data corresponds to [this page] (https://www.jma.go.jp/jp/highresorad/index.html) in the JMA website.

| radar/nowcast |  |  
| --- | --- |  
| Spatial resolution | 1 km |
| Temporal resolution | 5 min | 
| Forecast length | 1 hour |
| Refresh | 5 min |

The data is automatically downloaded in `/jmbsc/jmbsc_radar_5min/` and `/jmbsc/jmbsc_nowcast_1km_5min/`.   
The scripts `data/JMA_precip/draw/fortran/plot_{radar,nowcast}.sh` make plots of radar/nowcast precipitation of input time in png format. 


The cron job for automatic download and plot for these data is:
~~~
30 * * * * <PATH_TO_SCRIPTS>/data/JMA_precip/draw/auto.sh &> /data9/amemiya/scale-letkf-rt-interface/data/JMA_precip/draw/auto.log 
0,5,10,15,20,25,30,35,40,45,50,55 * * * * <PATH_TO_SCRIPTS>/data/JMA_precip/draw/auto_nowcast.sh &> /data9/amemiya/scale-letkf-rt-interface/data/JMA_precip/draw/auto_nowcast.log
~~~
