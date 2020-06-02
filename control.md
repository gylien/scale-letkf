## Control SCALE-LETKF analyses and forecasts 

The scripts in `admin` control realtime operation of the SCALE-LETKF system. 

```
admin_cycle.sh           ### run a single D1 DA cycle
admin_fcst_d1_ext.sh     ### run a single D1 extended forecast
admin_fcst_d1-2.sh       ### run a single D1-2 forecast
admin_fcst_d3.sh         ### run a single D3 forecast and generate D4 initial/boundary conditions from it 

auto_cycle.sh            ### automatically run D1 DA cycle for multiple times  
auto_fcst_d1_ext.sh      ### automatically run D1 extended forecasts 
auto_fcst_d1-2.sh        ### automatically run D1-2 forecasts
auto_fcst_d3.sh          ### automatically run realtime D3 forecasts 

auto_cycle_d4.sh         ### automatically run realtime D4 DA cycle and extended forecasts 

time_offset.txt          ### time offset in second for past mode (optional)
```
### Quick start

**Before running the system, you should make sure no jobs are running.**  
If running* files exist at the admin directory, re-login to the login node where pre-existing jobs are submitted and kill all processes.  
Information on nodes and process IDs is included in running* files.  
  
After that, remove running* files and clean up run directories.
```
rm -f ../scale_ope/scale-letkf_ope/scale/run_d1-2/waiting_list
rm -f ../scale_ope/scale-letkf_ope/scale/run_d1-2/*stat*

rm -f ../scale_ope/scale-letkf_ope/scale/run_d3/waiting_list
rm -f ../scale_ope/scale-letkf_ope/scale/run_d3/*stat*
```
  
To launch the realtime system from an initial date, say, 00UTC 1/1/2020, execute the following commands from a login node on OFP.  
**To monitor the bash processes, run 'ps auxf |grep {account}' on the "same" login node**  
```
echo '2020-01-01 00' > admin_cycle.time
nohup ./auto_cycle.sh &> auto_cycle.log &
nohup ./auto_fcst_d1-2.sh "2020-01-01 00:00:00" &> auto_fcst_d1-2.log &
nohup ./auto_fcst_d3.sh &> auto_fcst_d3.log &
```
To run D4 cycle and forecast, re-login to ofp02 and run the following command.
```
nohup ./auto_cycle_d4.sh &> auto_cycle_d4.log &
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

The latest time of analysis data is recorded in `admin_cycle.time` in a "YYYY-MM-DD HH" format.
If a target time needs to be specified, for example,
```
echo '2020-01-01 00' > admin_cycle.time
```
and the next D1 analysis creates an analysis at "2020-01-01 06". 

To run a single analysis cycle, 
```
nohup ./admin_cycle.sh &
```
The progress will be shown in `admin_cycle.log`. While the program is running a temporary file `admin_cycle.lock` appears to prevent accidental multiple execution. When `admin_cycle.sh` ends successfully, the time stamp in `admin_cycle.time` is updated. 

To run a realtime analysis cycle automatically, 
```
nohup ./auto_cycle.sh &> auto_cycle.log &
```
this repeatedly call `admin_cycle.sh` or wait until when NCEP GFS and PREPBUFR data for next time step is avaiable.

### D1 extended forecast

D1 5-day extended forecast from an arbitrary initial time is run by 
```
nohup ./admin_fcst_d1_ext.sh "2020-01-01 06:00:00" &
```
as long as initial files (D1 analysis) and NCEP GFS files are ready.  
The progress will be shown in `admin_fcst.log`.

To run it automatically, specify initial time to start. 
```
nohup ./auto_fcst_d1_ext.sh "2020-01-01 06:00:00" &> auto_fcst_d1_ext.log &
```
this automatically run forecast from "2020-01-01 06:00:00" first and from "2020-01-01 12:00:00" next, and so on.

Simutaneous run of multiple forecasts from different initial times is supported (see also [here](misc.md)). This also applies for D1-2 and D3. 

### D1 and D2 forecast

D1 and D2 online-nesting forecast from an arbitrary initial time is run in the similar way. 
```
nohup ./admin_fcst_d1-2.sh "2020-01-01 06:00:00" &
```
as long as initial files (D1 analysis) and NCEP GFS files are ready.  
The progress will be shown in `admin_fcst_d1-2.log`.

To run it automatically, specify initial time to start. 
```
nohup ./auto_fcst_d1-2.sh "2020-01-01 06:00:00" &> auto_fcst_d1-2.log &
```
this automatically run forecast from "2020-01-01 06:00:00" first and from "2020-01-01 12:00:00" next, and so on.

### D3 forecast and creating D4 init/boundary files

To run a D3 forecast, an initial time of *D2 forecast* to be used as initial/boundary condition, an initial time of D3 forecast, and a forecast length need to be specified.   
```
nohup ./admin_fcst_d1-2.sh "2020-01-01 06:00:00" "2020-01-01 10:00:00" 21600 &
```
as long as initial and boundary files (D2 forecast initialized at 2020-01-01 06Z and includes forecast range 2020-01-01 10Z-16Z) are ready.  
The progress will be shown in `admin_fcst_d3.log`.

The automatic run of D3 is a bit more complicated than D1-2 because it is designed to automatically adjust the init time to current real time.
```
nohup ./auto_fcst_d3.sh &> auto_fcst_d3.log &
```

### Continuous run of automatic scripts
While automatic run scripts like `auto_cycle.sh` is running, a temporary file such as `running_cycle` appears. It stores **the hostname (ofpXX or obcxXX) and PID** of the corresponding automatic program currently running. 

To manually stop automatic scripts, **log in to the corresponding host** and kill the process with the corresponding PID.

### Past mode 
The scripts supports the "past mode", in which the test case in a specific time in the past can be treated as a virtual realtime operation. 
To activate past mode, set time offset in `time_offset.txt` in a unit of second.
For example,  if you want to treat the past "2019 August 24th 15 UTC" as if it is happening now, set the offset as follows.  
```
nows= `date -u +%s`
pasts=`date -ud "2019-08-24 15:00:00" +%s`
echo `expr $pasts - $nows` > time_offset.txt
```
Then the time offset applies for the scripts `auto_fcst_d3.sh` and `auto_cycle_d4.sh` as well as for monitoring tools.  

Set the latest analysis time in `admin_cycle.time`
```
echo "2019-08-24 00" > admin_cycle.time
```

To imitate 6-hourly download of realtime NCEP data for past events, `external/past/get_past.sh` can be used instead of `get_ncep_gfs.sh` and `get_ncep_obs.sh`.
Set crontab to call it periodically. 
```
### crontab
0 * * * * $EXTDIR/past/get_past.sh > $EXTDIR/past/get_past.log
```
