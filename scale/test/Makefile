include ../configure.user

PROGS   = group obsope

OBJS    =

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o
#              $(COMMON_DIR)/common_mtx.o     \
#              $(COMMON_DIR)/netlib.o         \
#              $(COMMON_DIR)/common_letkf.o   \
#              $(COMMON_DIR)/common_ncio.o    \
#              $(COMMON_DIR)/common_scale.o     \
#              $(COMMON_DIR)/common_obs_scale.o \
#              $(COMMON_DIR)/common_mpi_scale.o
#              $(COMMON_DIR)/common_precip.o

INCS    = -I/data1/gylien/scale_ens/dc_utils \
          -I/data1/gylien/scale_ens/gtool \
          -I/data1/gylien/scale_ens/include
#          -I/data1/gylien/scale/scalelib/include

LIBS    = -L/data1/gylien/scale_ens/lib -lscale -lgtool -ldcutils \
          $(NETCDF_LIB) # $(LIBBLAS)

all:	$(PROGS)

group:	group.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o group group.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)

group.o:	group.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c group.f90 -I$(COMMON_DIR) $(INCS)

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod $(PROGS)

mod_ioapi.o:	common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -c mod_ioapi.f90 -I$(COMMON_DIR) $(INCS)

obsope.o: mod_ioapi.o common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -c obsope.f90 -I$(COMMON_DIR) $(INCS)

obsope: obsope.o mod_ioapi.o common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o obsope obsope.o mod_ioapi.o $(COMMON_OBJS) $(LIBS)
#-mkl
