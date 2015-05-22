#!/bin/bash

wkdir="$(cd "$( dirname "$0" )" && pwd)"
cd ${wkdir}

#------

if (($# < 3)); then
  echo "$0: Insufficient arguments" >&2
  exit 1
fi

STIME="$1"
TIME_DT="$2"
TIME_DT_DYN="$3"

#------

cat config/EastAsia_18km_48p/config.main.K | \
    sed -e "s/<STIME>/${STIME}/g" \
    > config.main

cat config/EastAsia_18km_48p/config.cycle | \
    sed -e "s/<STIME>/${STIME}/g" \
    > config.cycle

cat config/EastAsia_18km_48p/config.nml.scale | \
    sed -e "s/<TIME_DT>/${TIME_DT}/g" | \
    sed -e "s/<TIME_DT_DYN>/${TIME_DT_DYN}/g" \
    > config.nml.scale

./cycle_K.sh
