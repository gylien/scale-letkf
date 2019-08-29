#!/bin/bash
#
# JIT-DT offline test launcher
#

# $1: Number of cycles

# load config.main
. config.main

rm -f ${TMP_JITDATA}/job.running

CNT=$1

watch -g ls ${TMP_JITDATA}/job.running > /dev/null
echo "Launch offline JIT-DT TEST script!"
#cd ${TMP_JITDATA}/../../bin/
#echo $pwd
#./src/jitdt-offline-gen.sh ${LCYCLE} ${CNT} ${JIT_TOP} ${JIT_TARFILE} ${JIT_OFFLOG}
./src/jitdt-offline-gen.sh ${LCYCLE} ${CNT} ${TMP_JITDATA}/../../ ${JIT_TARFILE} ${JIT_OFFLOG}
#cd -

exit 0

