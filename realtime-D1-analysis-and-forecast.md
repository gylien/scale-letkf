## SCALE-LETKF realtime system domain 1
### Overview

* Domain setting  

| | |
| --- | --- |
| Target area | Japan |
| Map projection | Lambert conformal conic projection|
| dx | 18 km |
| Model top | 28.8 km |
| levels | 36 |
| dt | 40 sec |
| Number of processes | 192 |
| Number of subdomain in X | 12 |
| Number of subdomain in Y | 16 |


* Boundary data

NCEP GFS analysis and forecast data is used for the boundary condition.   
The data is obtained from [NOMAD](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/) and converted every 6 hours (See [here](prepare-ncep-gfs-and-prepbufr-data.md)).

### Analysis cycle

* Start from the previous step analysis data 

### Forecast 

Extended forecasts for 5 days ahead are performed every 6 hours, initialized by `mdet` of D1 analysis. 
