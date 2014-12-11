SFC       = frtpx
FC        = mpifrtpx
FOPTS     = -Kfast,parallel,ocl,openmp \
# preex,array_private,noalias=s,mfunc=2 \
            -Koptmsg=1
#            -Qi -Qt -X03 -Ncompdisp -Cpp   \
#            -Kdynamic_iteration -Ksimd     \
#            -Kauto,threadsafe
FOPTS2    = $(FOPTS)
FMCMODEL  =
FFREE     = -Free
FFIXED    = -Fixed

LIBBLAS   = -SSL2BLAMP

NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm

SCALE_INCS = -I/volume63/data/ra000015/gylien/scale_ens/dc_utils \
             -I/volume63/data/ra000015/gylien/scale_ens/gtool    \
             -I/volume63/data/ra000015/gylien/scale_ens/include
SCALE_LIBS = -L/volume63/data/ra000015/gylien/scale_ens/lib -lscale -lgtool -ldcutils

BUFR_LIB  = 
