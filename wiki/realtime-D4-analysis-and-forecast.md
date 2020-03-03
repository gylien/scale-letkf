## SCALE-LETKF realtime system domain 4
### Overview

* Domain setting  

| | | 
| --- | --- |
| Target area | Kansai / Kanto |
| Map projection | Mercator |
| dx | 250m / 500 m / 1 km | 
| Model top | 28.8 km |
| levels | 60 |
| dt | 4 sec |
| Number of processes | 256 |
| Number of subdomain in X | 16 |
| Number of subdomain in Y | 16 |
<img src="https://github.com/aamemiya/shared_image/blob/master/D4_Kobe.png" height="250px">
<img src="https://github.com/aamemiya/shared_image/blob/master/D4_Tokyo.png" height="250px">

### Analysis cycle

* Start from the previous step analysis data 

### Forecast 

Extended forecasts for 5 days ahead are performed every 6 hours, initialized by `mdet` of D1 analysis. 
