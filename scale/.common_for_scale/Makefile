include ../configure.user

OBJS    = SFMT.o common.o common_mpi.o common_mtx.o netlib.o common_letkf.o \
	  common_ncio.o common_scale.o common_obs_scale.o common_mpi_scale.o \
          common_nml.o \
          mod_admin_time.o

all:	$(OBJS) 

common.o:	../../common/common.f90
	$(FC) $(FOPTS_SCALE) -c ../../common/common.f90
common_mpi.o:	../../common/common_mpi.f90 common.o
	$(FC) $(FOPTS_SCALE) -c ../../common/common_mpi.f90
common_mtx.o:	../../common/common_mtx.f90 common.o
	$(FC) $(FOPTS_SCALE) -c ../../common/common_mtx.f90
SFMT.o: ../../common/SFMT.f90
	$(FC) $(FOPTS_SCALE) -c ../../common/SFMT.f90
netlib.o:	../../common/netlib.f
	$(FC) $(FOPTS_SCALE) -c ../../common/netlib.f
netlibblas.o:	../../common/netlibblas.f
	$(FC) $(FOPTS_SCALE) -c ../../common/netlibblas.f
common_letkf.o:	../../common/common_letkf.f90 common.o common_mtx.o
	$(FC) $(FOPTS_SCALE) -c ../../common/common_letkf.f90
common_ncio.o:	../common/common_ncio.f90 common.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_ncio.f90 $(NETCDF_INC)
common_scale.o:	../common/common_scale.f90 common.o common_nml.o common_mpi.o common_ncio.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_scale.f90 $(NETCDF_INC) $(SCALE_INC)
common_obs_scale.o:	../common/common_obs_scale.f90 common.o common_scale.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_obs_scale.f90 $(NETCDF_INC) $(SCALE_INC)
common_mpi_scale.o:	../common/common_mpi_scale.f90 common.o common_nml.o common_mpi.o \
	common_scale.o common_obs_scale.o mod_admin_time.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_mpi_scale.f90 $(NETCDF_INC) $(SCALE_INC)
common_nml.o:	../common/common_nml.f90 common.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_nml.f90 $(INCS) $(SCALE_INC)

mod_admin_time.o:
	$(FC) $(FOPTS_SCALE) -c $(SCALE)/scale-les/src/admin/mod_admin_time.f90 $(SCALE_INC)

clean:
	rm -f *.o *.mod *.lst
