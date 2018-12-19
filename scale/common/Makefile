include ../configure.user

OBJS    = SFMT.o common.o common_mpi.o common_mtx.o netlib.o common_sort.o common_rand.o common_letkf.o \
	  common_ncio.o common_scale.o common_obs_scale.o common_mpi_scale.o \
          common_nml.o
          #common_precip.o

OBS_OBJS   = common.o common_nml.o common_scale.o 

ifeq ($(H08),T)
OBJS    += scale_H08_fwd.o
OBS_OBJS   += scale_H08_fwd.o
endif
ifeq ($(GPR),T)
OBJS	+= common_geosatpr.o module_simulator.o module_opt_micro.o \
	   simulator_radar.o module_fov.o
OBS_OBJS	+= common_geosatpr.o module_simulator.o module_opt_micro.o \
	   	   simulator_radar.o module_fov.o
endif

all:	$(OBJS) datetime pdbash 

datetime:	datetime.o SFMT.o common.o
	$(SFC) $(FOPTS) -o datetime datetime.o SFMT.o common.o
pdbash:	pdbash.o SFMT.o
	$(FC) $(FOPTS) -o pdbash pdbash.o

common.o:	../../common/common.f90
	$(SFC) $(FOPTS) -c $<
common_mpi.o:	../../common/common_mpi.f90 common.o
	$(FC) $(FOPTS) -c $<
common_mtx.o:	../../common/common_mtx.f90 common.o
	$(SFC) $(FOPTS) -c $<
SFMT.o: ../../common/SFMT.f90
	$(SFC) $(FOPTS) -c $<
netlib.o:	../../common/netlib.f
	$(SFC) $(FOPTS) -c $<
netlibblas.o:	../../common/netlibblas.f
	$(SFC) $(FOPTS) -c $<
common_sort.o:	../../common/common_sort.f90 common.o
	$(SFC) $(FOPTS) -c $<
common_rand.o:	../../common/common_rand.f90
	$(SFC) $(FOPTS) -c $<
common_letkf.o:	../../common/common_letkf.f90 common.o common_mtx.o
	$(SFC) $(FOPTS) -c $<
common_ncio.o:	common_ncio.f90 common.o
	$(SFC) $(FOPTS) -c $< $(NETCDF_INC)
common_scale.o:	common_scale.f90 common.o common_nml.o common_mpi.o common_ncio.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
ifeq ($(H08),T)
scale_H08_fwd.o:        scale_H08_fwd.F90 common.o 
	$(FC) $(FOPTS) -c $< $(RT_INCLUDE) $(RT_LDFLAGS) $(SCALE_INC)
endif
ifeq ($(GPR),T)
common_geosatpr.o: ../jsim/common_geosatpr.f90 common.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
module_simulator.o: ../jsim/module_simulator.f90 common.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
module_opt_micro.o: ../jsim/module_opt_micro.f90 common.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
simulator_radar.o: ../jsim/simulator_radar.f90 common.o module_opt_micro.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
module_fov.o: ../jsim/module_fov.f90 common.o
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
endif
common_obs_scale.o:	common_obs_scale.f90 $(OBS_OBJS)
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
common_mpi_scale.o:	common_mpi_scale.f90 common.o common_nml.o common_mpi.o \
	common_scale.o common_obs_scale.o 
	$(FC) $(FOPTS) -c $< $(NETCDF_INC) $(SCALE_INC)
#common_precip.o:	common_precip.f90 common.o common_scale.o common_letkf.o
#	$(SFC) $(FOPTS) -c $<
common_nml.o:	common_nml.f90 common.o
	$(FC) $(FOPTS) -c $< $(INCS) $(SCALE_INC)

datetime.o:	datetime.f90 common.o
	$(SFC) $(FOPTS) -c $<
pdbash.o:	pdbash.f90
	$(FC) $(FOPTS) -c $<

clean:
	rm -f *.o *.mod *.lst datetime pdbash
