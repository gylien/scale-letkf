#!/bin/sh

module load hdf5
module load netcdf
module load netcdf-fortran

NCDFPATH=/work/opt/local/apps/intel/2018.1.163/netcdf-fotran/4.4.3

/work/hp150019/c24140/Lib/dcl/bin/dclfrt draw_uvw.f90 -o draw_uvw -L${NCDFPATH}/lib -lnetcdff  -I${NCDFPATH}/include ### -CU -CB -traceback
