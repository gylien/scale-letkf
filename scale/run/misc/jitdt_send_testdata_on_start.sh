#!/bin/bash
#===============================================================================
#
#  (Re)Start jit-kwatch daemon, and send test data when the job starts
#
#===============================================================================

. config.main

TEST_DIR='/volume63/data/hp150019/gylien/obs/20130713'
TEST_TARFILE='osaka_20130713060030-070000_v2a.tar'
TEST_LOG='osaka_20130713060030-070000.log'
TEST_INTERVAL=30
TEST_ITERATION=120

echo
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Restart daemon..."
echo

mkdir -p ${TMP_JITDATA}
rm -f ${TMP_JITDATA}/job.running
mkdir -p /tmp/$(id -nu)

jit-kwatch stop /tmp/$(id -nu)
jit-kwatch start /tmp/$(id -nu) ${TMP_JITDATA}

echo
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Wait until job is running..."
echo

while [ ! -e "${TMP_JITDATA}/job.running" ]; do
  sleep 1s
done

echo
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Start sending data..."
echo

ktestgen.pl ${TEST_DIR}/${TEST_LOG} \
            ${TEST_DIR}/${TEST_TARFILE} \
            /tmp/$(id -nu) \
            ${TMP_JITDATA} \
            3 ${TEST_INTERVAL} ${TEST_ITERATION}

sleep $((TEST_INTERVAL + 60))s

echo
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Stop daemon..."
echo

jit-kwatch stop /tmp/$(id -nu)

echo
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done"
echo
