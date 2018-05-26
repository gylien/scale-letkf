include ../configure

PROGS   = letkf

OBJS    = letkf_obs.o letkf_tools.o

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
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_scalerm.o   \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o

INCS    = $(SCALE_INC) $(NETCDF_INC)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

OBS_DIR = ../obs

OBS_OBJS = $(OBS_DIR)/obs_tools.o \
           $(OBS_DIR)/obsope_tools.o \
           $(OBS_DIR)/radar_obs.o \
           $(OBS_DIR)/radar_tools.o \
           $(OBS_DIR)/read_toshiba.o \
           $(OBS_DIR)/read_toshiba_f.o

SCALE_RM_OBJS = $(shell ls $(SCALE_RM_OBJDIR)/*.o | grep -v "scale-rm.o\|scale-rm_pp.o\|scale-rm_init.o")

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

ifeq ($(ENABLE_JITDT),T)
OBS_OBJS += $(OBS_DIR)/jitdt_read_toshiba.o \
            $(OBS_DIR)/jitdt_read_toshiba_f.o
INCS    += $(JITDT_INC)
LIBS    += $(JITDT_LIB)
endif

.PHONY:	all clean FORCE

all:	$(PROGS)

letkf:	letkf.o $(OBJS) $(COMMON_DIR) $(OBS_DIR)
	$(FC) $(FOPTS_EXE) $(FMCMODEL) -o letkf letkf.o \
	$(COMMON_OBJS) $(OBJS) $(OBS_OBJS) $(SCALE_RM_OBJS) $(LIBS)

letkf.o:	letkf.f90 $(OBJS) $(COMMON_DIR) $(OBS_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) -I$(OBS_DIR) $(INCS)

letkf_obs.o:	letkf_obs.f90 $(COMMON_DIR) $(OBS_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) -I$(OBS_DIR) $(INCS)
letkf_tools.o:	letkf_tools.f90 letkf_obs.o $(COMMON_DIR) $(OBS_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) -I$(OBS_DIR) $(INCS)

$(COMMON_DIR):	FORCE
	$(MAKE) -C $@

$(OBS_DIR):	FORCE
	$(MAKE) -C $@ obj

clean:
	rm -f *.o *.mod *.lst $(PROGS) letkf

FORCE:
