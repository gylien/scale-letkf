#!/bin/bash

#===============================================================================

. config.main

MEM="$1"; shift

echo ""
echo "${INIT_LOC_ENS}"

IFS=','

#-----
ifile=${FULL_NAME}_${STIME}.txt
while read line
do
  # divide by comma
  set -- $line 
  echo $0

done < ${INIT_LOC_ENS}

#===============================================================================

exit 0
