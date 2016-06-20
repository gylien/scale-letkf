include ../configure.user

OBJS    = SFMT.o common.o common_mpi.o common_mtx.o netlib.o common_letkf.o \
	  common_ncio.o common_scale.o common_obs_scale.o common_mpi_scale.o \
          common_nml.o \
          farb_connect.o farb_util.o farb_conf_dict.o farb_conf_parse.o

ifeq ($(H08),T)
OBJS    += scale_H08_fwd.o
endif

all:	$(OBJS) 

farb_connect.o:	../common/farb_connect.c ../common/farb_connect.h
	$(CC) $(COPTS) -fPIC -c ../common/farb_connect.c
farb_util.o:	../common/farb_util.c ../common/farb_util.h
	$(CC) $(COPTS) -fPIC -c ../common/farb_util.c
farb_conf_dict.o:	../common/farb_conf_dict.c ../common/farb_conf_dict.h
	$(CC) $(COPTS) -fPIC -c ../common/farb_conf_dict.c
farb_conf_parse.o:	../common/farb_conf_parse.c ../common/farb_conf_parse.h
	$(CC) $(COPTS) -fPIC -c ../common/farb_conf_parse.c

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
ifeq ($(H08),T)
scale_H08_fwd.o: ../common/scale_H08_fwd.F90 common.o
	$(FC) $(FOPTS) -c ../common/scale_H08_fwd.F90 $(RT_INCLUDE) $(RT_LDFLAGS) $(SCALE_INC)
common_obs_scale.o:	../common/common_obs_scale.f90 common.o common_scale.o scale_H08_fwd.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_obs_scale.f90 $(NETCDF_INC) $(SCALE_INC)
else
common_obs_scale.o:	../common/common_obs_scale.f90 common.o common_scale.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_obs_scale.f90 $(NETCDF_INC) $(SCALE_INC)
endif
common_mpi_scale.o:	../common/common_mpi_scale.f90 common.o common_nml.o common_mpi.o \
	common_scale.o common_obs_scale.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_mpi_scale.f90 $(NETCDF_INC) $(SCALE_INC)
common_nml.o:	../common/common_nml.f90 common.o
	$(FC) $(FOPTS_SCALE) -c ../common/common_nml.f90 $(INCS) $(SCALE_INC)

clean:
	rm -f *.o *.mod *.lst
