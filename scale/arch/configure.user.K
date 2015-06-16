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
FOPTS3    = $(FOPTS) -Kparallel
FMCMODEL  =
FFREE     = -Free
FFIXED    = -Fixed

LIBBLAS   = -SSL2BLAMP

NETCDF_INC = -I/opt/aics/netcdf/k-serial-noszip/include
NETCDF_LIB = -L/opt/aics/netcdf/k-serial-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
#NETCDF_INC = -I/opt/aics/netcdf/k-noszip/include
#NETCDF_LIB = -L/opt/aics/netcdf/k-noszip/lib-static -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm

SCALE     = /volume63/data/hp150019/gylien/scale-letkf/update/scale_ens

SCALE_INC = -I$(SCALE)/dc_utils \
            -I$(SCALE)/gtool    \
            -I$(SCALE)/include
SCALE_LIB = -L$(SCALE)/lib -lscale -lgtool -ldcutils
SCALE_LES_OBJDIR = /volume63/data/hp150019/gylien/scale-letkf/update/scale/scale-les/test/case_real/check_japan/run/.libs

BUFR_LIB  = 
