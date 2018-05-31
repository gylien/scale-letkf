# On-off settings
#-----------------

SCALE_USER_VER = T

# -- for Himawari-8 assimilatin on the K computer (T.Honda)
# -- If you want to compile without RTTOV (Himawari assimilation),
# -- you should set "H08 = F".
H08        = F

ENABLE_JITDT = T

ENABLE_DTF = T

# Libraries
#-----------------

#NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
#NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
NETCDF_INC = $(SCALE_NETCDF_INCLUDE)
NETCDF_LIB = $(SCALE_NETCDF_LIBS)

LAPACK_LIB = -SSL2BLAMP

BUFR_LIB   =

# Himawari-8
DIR_RTTOV  = /data/share005/honda/RTTOV
RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov11.2.0_test -lrttov11.2.0_other \
             -lrttov11.2.0_coef_io -lrttov11.2.0_parallel -lrttov11.2.0_main
RT_INCLUDE = -I$(DIR_RTTOV)/mod -I$(DIR_RTTOV)/include

# JIT-DT
JITDT_INC  = -I/home/hp150019/k01977/opt/nowcast/include
JITDT_LIB  = -L/home/hp150019/k01977/opt/nowcast/lib -lk-jitdt-pk

# Compiler flags
#-----------------

#SFC        = frtpx
#FC         = mpifrtpx
#CC         = mpifccpx

FOPTS_SCALE = $(FFLAGS)
COPTS_SCALE = $(CFLAGS)

FOPTS      = $(FOPTS_SCALE) \
             -Knoparallel,noeval,nopreex \
             -Kopenmp
COPTS      = $(COPTS_SCALE) \
             -Knoparallel,noeval,nopreex

FOPTS_EXE  = $(FOPTS_SCALE) \
             -Kopenmp

FMCMODEL   =
FBYTESWAP  = 
FFREE      = -Free
FFIXED     = -Fixed
