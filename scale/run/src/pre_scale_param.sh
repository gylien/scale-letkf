#!/bin/bash
#===============================================================================
#
#  Script to modify run.conf for parameter estimation
#
#===============================================================================

. config.main

if (($# < 6)); then
  cat >&2 << EOF

[pre_scale_param.sh] Modify a run.conf for SCALE model run.

Usage: $0 MYRANK MEM STIME FCSTLEN FCSTINT HISTINT TMPDIR OUT_OPT [SCPCALL]

  MYRANK   My rank number (not used)
  MEM      Name of the ensemble member
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval of restart files (second)
  HISTINT  Output interval of history files (second)
  TMPDIR   Temporary directory to run the model
  OUT_OPT
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
HISTINT="$1"; shift
TMPDIR="$1"; shift
OUT_OPT="$1"; shift
SCPCALL="${1:-cycle}"; shift

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

# file list

#files=`ls ${TMPOUT}/param/*_${STIME}.txt`
#echo "files param: "${files}

echo ""
echo  ${TMPOUT}/param/PARAM_LIST.txt

IFS=','


#-----
# Parameter loop
while read line
do
  set -- $line 
  FULL_NAME=$1

  ifile=${FULL_NAME}_${STIME}.txt
  while read line
  do
    # divide by comma
    set -- $line 
    PNAME=$1
    PMEM=$2
    PVAL=$3

    PMEM=$(printf '%04d' $PMEM)

    if [ $MEM == $PMEM ] || [ "${PMEM}" = 'mean' ] || [ "${PMEM}" = 'mdet' ]; then
      sed -i "/!--${FULL_NAME}--/a ${PNAME} = ${PVAL}," $TMPDIR/run.conf
      break
    fi

  done < ${TMPOUT}/param/${ifile}

done < ${TMPOUT}/param/PARAM_LIST.txt # Parameter loop

#===============================================================================

exit 0
