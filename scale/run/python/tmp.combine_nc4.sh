#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=debug-flat
#PJM -L node=1
#PJM -L elapse=0:10:00
#PJM -g hp150019
#PJM -j
#------- Program execution -------#

MPIE=/work/hp150019/share/sw/anaconda3/5.2.0/bin/mpiexec
SPATH=/work/hp150019/f22013/SCALE-LETKF/scale-LT/scale_lt_devel_20190522/scale-letkf/scale/run/python
${MPIE} -n 1 python3  ${SPATH}/combine_nc4.py
