include ../configure

PROGS   = dacycle

OBJS    = 

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_sort.o    \
              $(COMMON_DIR)/common_rand.o    \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o   \
              $(COMMON_DIR)/common_scalerm.o \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o

INCS    = $(SCALE_INC) $(NETCDF_INC) -I$(SCALE_RM_OBJDIR)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

OBS_DIR = ../obs

OBS_OBJS = $(OBS_DIR)/obs_tools.o \
           $(OBS_DIR)/obsope_tools.o \
           $(OBS_DIR)/radar_obs.o \
           $(OBS_DIR)/radar_tools.o \
           $(OBS_DIR)/read_toshiba.o \
           $(OBS_DIR)/read_toshiba_f.o

LETKF_DIR = ../letkf

LETKF_OBJS = $(LETKF_DIR)/letkf_obs.o $(LETKF_DIR)/letkf_tools.o

SCALE_RM_OBJS = $(shell ls $(SCALE_RM_OBJDIR)/*.o | grep -v "scale-rm.o\|scale-rm_pp.o\|scale-rm_init.o")

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

.PHONY:	all clean FORCE

all:	$(PROGS)

dacycle:	dacycle.o $(OBJS) $(COMMON_DIR) $(OBS_DIR) $(LETKF_DIR)
	$(FC) $(FOPTS_EXE) $(FMCMODEL) -o dacycle dacycle.o \
	$(COMMON_OBJS) $(OBJS) $(OBS_OBJS) $(LETKF_OBJS) $(SCALE_RM_OBJS) $(LIBS)

dacycle.o:	dacycle.f90 $(OBJS) $(COMMON_DIR) $(OBS_DIR) $(LETKF_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) -I$(OBS_DIR) -I$(LETKF_DIR) $(INCS)

$(COMMON_DIR):	FORCE
	$(MAKE) -C $@

$(OBS_DIR):	FORCE
	$(MAKE) -C $@ obj

$(LETKF_DIR):	FORCE
	$(MAKE) -C $@ letkf_obs.o letkf_tools.o

clean:
	rm -f *.o *.mod *.lst $(PROGS) dacycle

FORCE:
