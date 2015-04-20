SFC       = frtpx
FC        = mpifrtpx
#FOPTS     = -Kfast,parallel,ocl,preex,array_private,noalias=s,mfunc=2
#FOPTS     = -Kfast,parallel,ocl,preex,array_private \
#            -Qi -Qt -X03 -Ncompdisp -Koptmsg=1 -Cpp                   \
#            -x-                                                       \
#            -Ksimd                                                    \
#            -Kauto,threadsafe
FOPTS     = -Kfast,ocl,openmp,noeval -V -Qa,d,i,p,t,x -Koptmsg=2
#FOPTS     = -Kfast,ocl,noeval -V -Qa,d,i,p,t,x -Koptmsg=2
#FOPTS     = -Kfast,parallel,ocl,noeval -V -Qa,d,i,p,t,x -Koptmsg=2
#FOPTS     = -Kfast,parallel,ocl,openmp,noeval -V -Qa,d,i,p,t,x -Koptmsg=2
#FOPTS     = -O0                                                       \
#            -Qi -Qt -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp \
#            -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global        \
#            -DDEBUG

FOPTS2    = $(FOPTS)
FMCMODEL  =
FFREE     = -Free
FFIXED    = -Fixed

LIBBLAS   = -SSL2BLAMP

NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
#NETCDF_INC = -I/opt/aics/netcdf/k-noszip/include
#NETCDF_LIB = -L/opt/aics/netcdf/k-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm

SCALE_INC = -I/volume63/data/hp150019/gylien/scale-letkf/scale_ens/dc_utils \
            -I/volume63/data/hp150019/gylien/scale-letkf/scale_ens/gtool    \
            -I/volume63/data/hp150019/gylien/scale-letkf/scale_ens/include
SCALE_LIB = -L/volume63/data/hp150019/gylien/scale-letkf/scale_ens/lib -lscale -lgtool -ldcutils

BUFR_LIB  = 
