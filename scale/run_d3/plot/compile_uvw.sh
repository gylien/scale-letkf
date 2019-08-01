#!/bin/sh

module load hdf5
module load netcdf-fortran

NCDFPATH=/work/opt/local/apps/intel/2018.1.163/netcdf-fotran/4.4.3

dclfrt draw_uvw.f90 -o draw_uvw -L${NCDFPATH}/lib -lnetcdff  -I${NCDFPATH}/include -CU -CB -traceback