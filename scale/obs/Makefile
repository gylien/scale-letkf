include ../configure

PROGS   = obsope #obsmake obssim #dec_prepbufr

OBJS    = obs_tools.o obsope_tools.o \
          radar_obs.o radar_tools.o \
	  read_toshiba$(TOSHIBA_OPT).o read_toshiba$(TOSHIBA_OPT)_f.o
 
COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_mtx.o     \
              $(COMMON_DIR)/netlib.o         \
              $(COMMON_DIR)/common_letkf.o   \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_scalerm.o   \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o 

INCS    = $(SCALE_INC) $(NETCDF_INC)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

SCALE_RM_OBJS = $(shell ls $(SCALE_RM_OBJDIR)/*.o | grep -v "scale-rm.o\|scale-rm_pp.o\|scale-rm_init.o")

ifeq ($(H08),T)
COMMON_OBJS += $(COMMON_DIR)/scale_H08_fwd.o
LIBS    += $(RT_INCLUDE) $(RT_LDFLAGS)
endif

ifeq ($(ENABLE_JITDT),T)
OBJS    += jitdt_read_toshiba$(TOSHIBA_OPT).o jitdt_read_toshiba_f.o
INCS    += $(JITDT_INC)
LIBS    += $(JITDT_LIB)
endif

ifeq ($(ENABLE_PLOT_DCL),T)
COMMON_OBJS += $(COMMON_DIR)/plot_dbz_DCL_obs.o
COMMON_OBJS += $(COMMON_DIR)/plot_dbz_DCL.o
LIBS += $(LIBDCL)
endif

.PHONY:	all clean FORCE

all:	$(PROGS)

obsope:	obsope.o $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS_EXE) $(FMCMODEL) -o obsope obsope.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
obsmake:	obsmake.o $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS_EXE) $(FMCMODEL) -o obsmake obsmake.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
obssim:	obssim.o $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS_EXE) $(FMCMODEL) -o obssim obssim.o \
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

obj:	$(OBJS)

obs_tools.o:	obs_tools.f90 radar_obs.o $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS) -I$(SCALE_RM_OBJDIR)
obsope_tools.o:	obsope_tools.f90 obs_tools.o $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

ifeq ($(ENABLE_JITDT),T)
ifeq ($(ENABLE_PLOT_DCL),T)
radar_obs.o:	radar/radar_obs.f90 radar_tools.o read_toshiba$(TOSHIBA_OPT)_f.o jitdt_read_toshiba_f.o ../common/plot_dbz_DCL_obs.o $(COMMON_DIR)
else
radar_obs.o:	radar/radar_obs.f90 radar_tools.o read_toshiba$(TOSHIBA_OPT)_f.o jitdt_read_toshiba_f.o  $(COMMON_DIR)
endif
else
ifeq ($(ENABLE_PLOT_DCL),T)
radar_obs.o:	radar/radar_obs.f90 radar_tools.o read_toshiba$(TOSHIBA_OPT)_f.o ../common/plot_dbz_DCL_obs.o $(COMMON_DIR)
else
radar_obs.o:	radar/radar_obs.f90 radar_tools.o read_toshiba$(TOSHIBA_OPT)_f.o $(COMMON_DIR)
endif
endif
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
radar_tools.o:	radar/radar_tools.f90 $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

read_toshiba.o:	radar/read_toshiba.c radar/read_toshiba.h
	$(CC) $(COPTS) -c $<
read_toshiba_f.o:	radar/read_toshiba_f.f90 read_toshiba.o
	$(FC) $(FOPTS) -c $<
read_toshiba$(TOSHIBA_OPT).o:	radar/read_toshiba$(TOSHIBA_OPT).c radar/read_toshiba$(TOSHIBA_OPT).h
	$(CC) $(COPTS) -c $<
read_toshiba$(TOSHIBA_OPT)_f.o:	radar/read_toshiba$(TOSHIBA_OPT)_f.f90 read_toshiba$(TOSHIBA_OPT).o
	$(FC) $(FOPTS) -c $<
jitdt_read_toshiba$(TOSHIBA_OPT).o:	radar/jitdt_read_toshiba$(TOSHIBA_OPT).c radar/jitdt_read_toshiba$(TOSHIBA_OPT).h
	$(CC) $(COPTS) -c $< $(JITDT_INC)
jitdt_read_toshiba_mpr.o:	radar/jitdt_read_toshiba_mpr.c radar/jitdt_read_toshiba_mpr.h
	$(CC) $(COPTS) -c $< $(JITDT_INC)

jitdt_read_toshiba_f.o:	radar/jitdt_read_toshiba_f.f90 jitdt_read_toshiba$(TOSHIBA_OPT).o read_toshiba$(TOSHIBA_OPT)_f.o
	$(FC) $(FOPTS) -c $<

$(COMMON_DIR):	FORCE
	$(MAKE) -C $@

clean:
	rm -f *.o *.mod *.lst $(PROGS)

FORCE:
