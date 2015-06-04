include ../configure.user

PROGS   = scale-les_ens

OBJS    = $(shell ls obj/*.o)

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o     \
              $(COMMON_DIR)/mod_admin_time.o
#              $(COMMON_DIR)/common_scalelib.o
#              $(COMMON_DIR)/common_precip.o

INCS    = $(SCALE_INC)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LIBBLAS)

all:	$(PROGS)

scale-les_ens:	scale-les_ens.o $(OBJS) common_modules
	$(FC) $(FOPTS) $(FMCMODEL) -o scale-les_ens scale-les_ens.o \
	$(COMMON_OBJS) $(OBJS) $(LIBS)

scale-les_ens.o:	scale-les_ens.f90 $(OBJS) common_modules
	$(FC) $(FOPTS) -c scale-les_ens.f90 -I$(COMMON_DIR) $(INCS) \
	-Iobj

common_modules:
	( cd $(COMMON_DIR) ; make )

clean:
	rm -f *.o *.mod *.lst $(PROGS) letkf mean efso
