include ../configure

PROGS   = obsope obssim #obsmake dec_prepbufr

OBJS    = obsope_tools.o

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
              $(COMMON_DIR)/common_nml.o

INCS    = $(SCALE_INC) $(NETCDF_INC) -I$(SCALE_RM_OBJDIR)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

.PHONY:	all clean FORCE

all:	$(PROGS)

obsope:	obsope.o $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obsope obsope.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS) -I$(SCALE_RM_OBJDIR)
obsmake:	obsmake.o $(OBS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obsmake obsmake.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
obssim:	obssim.o $(OBS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obssim obssim.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
dec_prepbufr:	dec_prepbufr.o $(COMMON_DIR)
	$(FC) $(FOPTS) -o dec_prepbufr dec_prepbufr.o $(COMMON_OBJS) \
	$(LIBS) $(BUFR_LIB)

obsope.o:	obsope.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
obsmake.o:	obsmake.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
obssim.o:	obssim.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
dec_prepbufr.o:	dec_prepbufr.f90 $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

obsope_tools.o:	obsope_tools.f90 $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

$(COMMON_DIR):	FORCE
	$(MAKE) -C $@

clean:
	rm -f *.o *.mod *.lst $(PROGS) obsope

FORCE:
