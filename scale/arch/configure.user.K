SCALE     = ../../..
SCALE_INC = -I$(SCALE)/dc_utils \
            -I$(SCALE)/gtool    \
            -I$(SCALE)/include
SCALE_LIB = -L$(SCALE)/lib -lscale -lgtool -ldcutils
SCALE_RM_OBJDIR = $(SCALE)/scale-rm/src/.libs

include $(SCALE)/sysdep/Makedef.$(SCALE_SYS)
include $(SCALE)/Mkinclude

#NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
#NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
NETCDF_INC = $(NETCDF_INCLUDE)
NETCDF_LIB = $(NETCDF_LIBS)

NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm
#NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm

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

# -- for Himawari-8 assimilatin on the K computer (T.Honda)
# -- If you want to compile without RTTOV (Himawari assimilation),
# -- you should set "H08 = F".
H08       = T

HDF5_PREFIX  = /opt/klocal/hdf5-1.10.1-fujitsu

DIR_RTTOV = /volume64/data/ra001011/honda/RTTOV12.2_netcdf
#RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov12_wrapper -lrttov12_mw_scatt -lrttov12_brdf_atlas -lrttov12_emis_atlas -lrttov12_other -lrttov12_parallel -lrttov12_coef_io -lrttov12_hdf -lrttov12_main -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5

RT_LDFLAGS = -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
             -L$(DIR_RTTOV)/lib -lrttov12_brdf_atlas -lrttov12_emis_atlas -lrttov12_mw_scatt -lrttov12_other -lrttov12_coef_io -lrttov12_hdf -lrttov12_parallel -lrttov12_main 

##-lrttov12_other \
#             -lrttov12_coef_io -lrttov12_parallel -lrttov12_main -lrttov12_hdf \
#             -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5


#RT_LDFLAGS = -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
#             -L$(DIR_RTTOV)/lib -lrttov12_other -lrttov12_coef_io -lrttov12_hdf -lrttov12_parallel -lrttov12_main 

#RT_INCLUDE = -D_RTTOV_HDF -I$(HDF5_PREFIX)/include \
#             -I$(DIR_RTTOV)/mod -I$(DIR_RTTOV)/include


#RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov12_emis_atlas -lrttov12_brdf_atlas -lrttov12_other -lrttov12_mw_scatt -lrttov12_coef_io -lrttov12_hdf -lrttov12_parallel -lrttov12_main \
#             -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 

RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov12_wrapper -lrttov12_mw_scatt -lrttov12_brdf_atlas -lrttov12_emis_atlas -lrttov12_other -lrttov12_parallel -lrttov12_coef_io -lrttov12_hdf -lrttov12_main -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
#-L$(DIR_RTTOV)/lib -lrttov12_emis_atlas -lrttov12_brdf_atlas -lrttov12_other -lrttov12_mw_scatt -lrttov12_coef_io -lrttov12_hdf -lrttov12_parallel -lrttov12_main \
#             -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 


RT_INCLUDE = -I$(DIR_RTTOV)/include \
             -D_RTTOV_HDF -I/opt/aics/netcdf/k-serial-noszip/include -I$(DIR_RTTOV)/mod
RT_LDFLAGS = -L$(DIR_RTTOV)/lib -lrttov12_wrapper -lrttov12_mw_scatt -lrttov12_brdf_atlas -lrttov12_emis_atlas -lrttov12_other -lrttov12_parallel -lrttov12_coef_io -lrttov12_hdf -lrttov12_main -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm

ifeq ($(H08),T)
FOPTS     += -DH08
FOPTS_SCALE += -DH08
endif
