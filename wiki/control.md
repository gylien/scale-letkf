## Control SCALE-LETKF analyses and forecasts 

The scripts in `admin` control realtime operation of the SCALE-LETKF system. 

The following scripts run a corresponding single analysis or forecast. 
```
admin_cycle.sh
admin_fcst_d1_ext.sh
admin_fcst_d1-2.sh
admin_fcst_d3.sh \*
```
\* Including D4 initial and boundary conditions downscaled from resultant D3 forecast.

(TODO: scripts for D4)

The following scripts launch automatic realtime execusion of corresponding analysis or forecast. 
```
auto_cycle.sh
auto_fcst_d1_ext.sh
auto_fcst_d1-2.sh
auto_fcst_d3.sh
```

### Configuration

`admin.rc` sets the common configuration. 

```
--admin.rc--

scale_ver="ope"
letkf_ver="ope"
expname="ope"
#expname="test"

nmem=50

nmem_d1_ext='mdet' ### deterministic run
nmem_d2=2 ### disk saving mode
nmem_d3=2 ### disk saving mode
#nmem_d2=50
#nmem_d3=50

LCYCLE=21600
FCSTLEN=432000
FCSTLEN_d2=86400

realtimebase="/work/${groupname}/share/SCALE-LETKF-rt"
###realtimebase="/work/${groupname}/${username}/HPCC_SCALE-LETKF-rt"
```

The version of source code of the SCALE model and SCALE-LETKF can be altered by changing `scale_ver` and `letkf_ver`. 

The number of members to run LETKF or forecasts in each domain can be modified by `nmem` and `nmem_{domain}`.
The forecast length for domain 1 (extended) and domain 1-2 (used for D3 init/boundary) can be modified by `FCSTLEN` and `FCSTLEN_d2`.

**(LCYCLE should not be changed)**

If more detailed configuration needs to be modified (such as changing output variables), go to `scale_ope/scale-letkf_ope/scale/run{domain}` and edit scripts or namelist files. 


### D1 analysis cycle

The latest time of analysis data is recorded in `admin_cycle.time` in a "YYYY-MM-DD HH:MM:SS" format.
If a target time needs to be specified, for example,
```
echo '2020-01-01 00:00:00' > admin_cycle.time
```
and the next D1 analysis creates an analysis at "2020-01-01 06:00:00". 

To run a single analysis cycle, 
```
nohup ./admin_cycle.sh &
```
The progress will be shown in `admin_cycle.log`. While the program is running a temporary file `admin_cycle.lock` appears to prevent accidental multiple execution. When it ends successfully, the time stamp in `admin_cycle.time` is updated. 

To run a realtime analysis cycle automatically, 
```
nohup ./auto_cycle.sh &> auto_cycle.log &
```
this repeatedly call `admin_cycle.sh` or wait until NCEP data for next time step is avaiable.

### D1 extended forecast

D1 5-day extended forecast from an arbitrary initial time is run by 
```
nohup ./admin_fcst_d1_ext.sh "2020-01-01 06:00:00" &
```
as long as initial files (D1 analysis) and NCEP GFS files are ready.  
The progress will be shown in `admin_fcst.log`.

To run it automatically, specify initial time to start. 
```
nohup ./auto_fcst_d1_ext.sh "2020-01-01 06:00:00" & > auto_fcst_d1_ext.log &
```
this automatically run forecast from "2020-01-01 06:00:00" first and from "2020-01-01 12:00:00" next, and so on.

Multiple forecasts from different initial times can be run simulateneously (see also [here](misc.md)).


### D1 and D2 forecast

D1 and D2 online-nesting forecast from an arbitrary initial time is run in the similar way. 
```
nohup ./admin_fcst_d1-2.sh "2020-01-01 06:00:00" &
```
as long as initial files (D1 analysis) and NCEP GFS files are ready.  
The progress will be shown in `admin_fcst_d1-2.log`.

To run it automatically, specify initial time to start. 
```
nohup ./auto_fcst_d1-2.sh "2020-01-01 06:00:00" & > auto_fcst_d1-2.log &
```
this automatically run forecast from "2020-01-01 06:00:00" first and from "2020-01-01 12:00:00" next, and so on.

Multiple forecasts from different initial times can be run simulateneously (see also [here](misc.md)).

### D3 forecast and creating D4 init/boundary files

To run a D3 forecast, an initial time of *D2 forecast* to be used as initial/boundary condition, an initial time of D3 forecast, and a forecast length need to be specified.   
```
nohup ./admin_fcst_d1-2.sh "2020-01-01 06:00:00" "2020-01-01 10:00:00" 21600 &
```
as long as initial and boundary files (D2 forecast initialized at 2020-01-01 06Z and includes forecast range 2020-01-01 10Z-16Z) are ready.  
The progress will be shown in `admin_fcst_d3.log`.

The automatic run of D3 is a bit more complicated than D1-2 because it is designed to automatically adjust the init time to current real time.
```
nohup ./auto_fcst_d3.sh & > auto_fcst_d3.log &
```
(description to be added)

Multiple forecasts from different initial times can be run simulateneously (see also [here](misc.md)).



### Contenious run of automatic scripts
While automatic run scripts like `auto_cycle.sh` is running, a temporary file such as `running_cycle` appears. It stores **the hostname (ofpXX or obcxXX) and PID** of the corresponding automatic program currently running. 

To manually stop automatic scripts, **log in to the corresponding host** and kill the process with the corresponding PID.


