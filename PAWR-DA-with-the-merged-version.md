(copied from [wiki/PAWR-DA-with-the-merged-version](https://github.com/gylien/scale-letkf/wiki/PAWR-DA-with-the-merged-version))

# Getting source
## SCALE 
https://gitlab.com/scale-met/scale/-/tree/LETKF-5.3.6/scale-rm/src  
or copy from   
/work/hp150019/share/honda/SCALE-LETKF/src/scale-LETKF-5.3.6.tar.gz (on OFP)
  
## LETKF
`git clone https://github.com/gylien/scale-letkf.git`  
`cd scale-letkf`  
`git checkout dacycle_OFP`  
If you want to update the above repository to the newest remote, `git pull`  
  
or simply   
`cp /work/hp150019/share/honda/SCALE-LETKF/src/<newest LETKF zip file below> YOUR-DIRECTORY` and unzip  
1bce6a0dacycle_OFP (as of 3/16/2020)  
2349224dacycle_OFP.zip (as of 4/20/2020)  
ef56421dacycle_OFP.zip (as of 4/28/2020)   
bc6908ddacycle_OFP.zip (as of 5/11/2020)   
b2ce4a9dacycle_OFP.zip (as of 5/12/2020)   

   
#  Compile the codes

## Intel compiler
`module unload intel`  
`module load intel/2019.5.281`  

## SCALE
`export SCALE_SYS=OFP`  
`export SCALE_ENABLE_PNETCDF=F`  
`export SCALE_ENABLE_OPENMP=T`  
`export SCALE_USE_SINGLEFP=T`  
   
`cd scale-rm/src`   
`make -j8`  
  
## LETKF
`cd scale-letkf/scale`  
`make -j2`  
  
When uses run a test case near Tokyo (e.g., D4_500m_testcases or D4_500m_20190910), please make sure ENABLE_SAITAMA_MPW=T in arch/configure.user.OFP  
  
# Run a job  
`cd run`  
`ln -sf config/$EXPNAME/config.* .`  (e.g., EXPNAME= D4_500m_testcases or D4_500m_20190910)
`ln -sf config.main.ofp_honda config.main`  
`ln -sf config.cycle_honda config.cycle`      
Edit config.main (e.g., set OUTDIR, OBS, etc.)  
Edit config.cycle (e.g., STIME, ETIME, etc.)  
Edit config.nml.letkf (e.g., MIN_RADAR_REF_MEMBER_OBSRAIN, RADAR_SO_SIZE_HORI, etc.)  
`nohup ./cycle_ofp.sh >log&`  
  
LETKF logs can be monitored by  
`tail scale-letkf/scale/tmp/scale-letkf_$EXPNAME/cycle_job.o*`  
    
# Check results
## Analysis & forecast imagery
Two-dimensional plots drawing radar reflectivity are stored at $OUTDIR/$STIME/dafcst  
     
## GrADS data
Analysis and forecast ensemble means (state variables and radar reflectivity) are stored at $OUTDIR/$STIME/mean_grads
Examples of GrADS control files are stored in config/D4_500m_testcases   
  
fcst_ref.ctl: Radar reflectivity and radial velocity simulated from a SCALE forecast  
mean_pawr_grads.ctl:  Radar reflectivity analysis/guess simulated from the ensemble mean  
mean_grads.ctl:  Analyzed/first-guess state variables  
Users need to modify "dset" and "tdef" as your experiment  

# Tunable parameters
## PARAM_LETKF
* RELAX_ALPHA_SPREAD: Coefficient for RTPS  
* POSITIVE_DEFINITE_Q: Negative qv in analysis is replaced with zero
* POSITIVE_DEFINITE_QHYD: Negative qh (hydrometeors) in analysis are replaced with zero

## PARAM_LETKF_OBS
* HORI_LOCAL & VERT_LOCAL: Localization scales 
* MAX_NOBS_PER_GRID: Maximum number of observations assimilated in each grid (so-called "observation number limit")  

## PARAM_LETKF_RADAR
* RADAR_SO_SIZE_HORI & RADAR_SO_SIZE_VERT: Super-ob scale for PAWR  
* MIN_RADAR_REF_MEMBER_OBSRAIN: Required rainy members if PAWR observation is rainy
* MIN_RADAR_REF_MEMBER_OBSNORAIN: Required rainy members if PAWR observation is clear  
* RADAR_THIN_HORI: Thinning interval in horizontal (# of grids, RADAR_THIN_HORI=1 corresponds no horizontal thinning)  
* RADAR_THIN_VERT: Thinning interval in vertical (# of grids, RADAR_THIN_VERT=1 corresponds no vertical thinning)  

## PARAM_OBS_ERROR
* OBSERR_RADAR_REF & OBSERR_RADAR_VR: Observation errors  
  
# Important setting
## config.cycle
* STIME: Start time of analysis (In D4_500m_testcases, 4 cases are available as of 4/28/2020)  
* ETIME: End time of analysis
* TIME_LIMIT: Elapse time limit  
* DACYCLE_RUN_FCST: Run extended (dafcst) forecasts from the analysis (if DACYCLE_RUN_FCST=1)    
* MAX_DACYCLE_RUN_FCST: Maximum number of dacycle forecasts  
* NUM_DACYCLE_FCST_MEM: Number of dafcst members  
* DACYCLE_RUN_FCST_TIME: Forecast time for dacycle forecasts    
* ICYC_DACYCLE_RUN_FCST: Cycle # of 1st dacycle-forecast start (if this is 10, dafcst forecasts are not initiated before the 10th analysis)  
* ICYC_DACYCLE_RUN_ANALYSIS: Cycle # of 1st analysis (if this is 10, no obs are assimilated before the 10th cycle)  
* OUT_GRADS_DAFCST: GrADS output from dafcst  
* OUT_MEAN_GRADS: GrADS output (ensemble-mean analysis/guess)  
* USE_MDET_FCST: if this is 1(0), dafcst foreacsts are initiated from the deterministic member (ensemble mean)  

