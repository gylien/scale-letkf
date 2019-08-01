#!/bin/bash
#
# copied from /work/hp150019/share/nowcast/bin/offline-gen.sh on 07/31/2019
#
# $1: Transfer time interval (second)
# $2: Transfer cycle number
# $3: Top directory
# $4: tar file name
# $5: Log file name
#


INT=$1 # sec
CNT=$2
TOP=$3
TARFILE=$4
LOGFILE=$5

# default
if [ -z "$1" ]; then
  INT=30 # sec
fi
if [ -z "$2" ]; then
  CNT=2
fi
if [ -z "$3" ]; then
  echo "[Warning] Top directory for JIT-DT was not specified"
  exit 1
fi
if [ -z "$4" ]; then
  echo "[Warning] tar file for JIT-DT was not specified"
  exit 1
fi
if [ -z "$5" ]; then
  echo "[Warning] Log file for JIT-DT was not specified"
  exit 1
fi
#####

OFFLINE_TOP=$TOP/offline
DATADIR=$TOP/testdata
#BIN=$TOP/bin
#LOGFILE=testdata.txt


#/work/hp150019/share/nowcast/bin/ltestgen.pl $DATADIR/$LOGFILE $DATADIR/$TARFILE $OFFLINE_TOP/data ${INT} ${CNT}
#/work/hp150019/f22013/SCALE-LETKF/scale-5.3.2/scale-5.3.x_LETKF_sngl/scale-letkf-dacycle/scale/run/src/ltestgen.pl $DATADIR/$LOGFILE $DATADIR/$TARFILE $OFFLINE_TOP/data ${INT} ${CNT}
./src/ltestgen.pl $LOGFILE $TARFILE $OFFLINE_TOP/data ${INT} ${CNT}
