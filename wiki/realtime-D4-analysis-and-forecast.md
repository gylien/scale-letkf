## SCALE-LETKF realtime system domain 4

See [wiki - PAWR DA with the merged version](https://github.com/gylien/scale-letkf/wiki/PAWR-DA-with-the-merged-version) for a more detailed instruction. 

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
| Number of processes | 256/1024/4096 |
| Number of subdomain in X | 16/32/64 |
| Number of subdomain in Y | 16/32/64 |

<img src="https://github.com/aamemiya/shared_image/blob/master/D4_Kobe.png" height="250px"><img src="https://github.com/aamemiya/shared_image/blob/master/D4_Tokyo.png" height="250px">

### Analysis cycle

* Start from the previous step analysis data 

### Forecast 

Extended forecasts for 30 minutes ahead are performed, initiated with every 30 seconds analysis fields. 
