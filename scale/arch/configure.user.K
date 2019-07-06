TOPDIR    = ../../..
SCALE_INC = -I$(TOPDIR)/dc_utils \
            -I$(TOPDIR)/gtool    \
            -I$(TOPDIR)/include
SCALE_LIB = -L$(TOPDIR)/lib -lscale -lgtool -ldcutils
SCALE_RM_OBJDIR = $(TOPDIR)/scale-rm/src/.libs

include $(TOPDIR)/sysdep/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

#NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
#NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
NETCDF_INC = $(SCALE_NETCDF_INCLUDE)
NETCDF_LIB = $(SCALE_NETCDF_LIBS)

LAPACK_LIB = -SSL2BLAMP

BUFR_LIB  =

#######

#SFC       = frtpx
#FC        = mpifrtpx

FOPTS_SCALE = $(FFLAGS) 

#FOPTS     = -Kfast,ocl,openmp,noeval -V -Qa,d,i,p,t,x -Koptmsg=2 \
#            -x-    \
#            -Ksimd \
#            -Kauto,threadsafe
FOPTS     = $(FOPTS_SCALE) \
            -Knoparallel,noeval,nopreex \
            -Kopenmp

FMCMODEL  =
FBYTESWAP = 
FFREE     = -Free
FFIXED    = -Fixed

# TC vital DA
TCV       = F

# -- for Himawari-8 assimilatin on the K computer (T.Honda)
# -- If you want to compile without RTTOV (Himawari assimilation),
# -- you should set "H08 = F".
H08       = T

#DIR_RTTOV = /data/share005/honda/RTTOV
#DIR_RTTOV = /data/share025/honda/RTTOV
DIR_RTTOV = /volume63/data/hp120282_2/honda/SCALE-LETKF/RTTOV
RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov11.2.0_test -lrttov11.2.0_other \
             -lrttov11.2.0_coef_io -lrttov11.2.0_parallel -lrttov11.2.0_main
RT_INCLUDE = -I$(DIR_RTTOV)/mod -I$(DIR_RTTOV)/include

ifeq ($(TCV),T)
FOPTS     += -DTCV
endif


ifeq ($(H08),T)
FOPTS     += -DH08
endif
