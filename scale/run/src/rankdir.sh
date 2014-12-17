#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE topo creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[pre_scale_pp_topo.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK STIME TMPDIR EXECDIR DATADIR

  MYRANK   My rank number (not used)

EOF
  exit 1
fi

MYRANK="$1"

#===============================================================================

pwd > ../

#===============================================================================

exit 0
