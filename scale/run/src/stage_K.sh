#!/bin/bash
#===============================================================================
#
#  Create stage-in and stage-out scripts for the K-computer
#
#===============================================================================

. config.all
#. src/func_util.sh

STAGING_DIR="$1"

#-------------------------------------------------------------------------------
# stage-in: Files in TMPDAT directory

if [ -s "$STAGING_DIR/stagein.dat" ]; then
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if [ -z "$TMPL" ]; then
          echo "#PJM --stgin-dir \"${source} ${TMPDAT}/${destin}\""
        else
          if ((TMPDAT_MODE == 2)); then
            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPDAT}/${destin}\""
          else
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPDAT}/${destin}\""
          fi
        fi
      else
        if [ -z "$TMPL" ]; then
          echo "#PJM --stgin \"${source} ${TMPDAT}/${destin}\""
        else
          if ((TMPDAT_MODE == 2)); then
            echo "#PJM --stgin \"rank=0 ${source} 0:${TMPDAT}/${destin}\""
          else
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPDAT}/${destin}\""
          fi
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.dat"
fi

#-------------------------------------------------------------------------------
# stage-in: Files in TMPOUT directory

i=0
while [ -s "$STAGING_DIR/stagein.out.$((i+1))" ]; do
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -z "$TMPL" ]; then
        echo "#PJM --stgin \"${source} ${TMPDAT}/${destin}\""
      else
        echo "#PJM --stgin \"rank=${i} ${source} ${i}:${TMPDAT}/${destin}\""
      fi
    fi
  done < "$STAGING_DIR/stagein.out.$((i+1))"
  i=$((i+1))
done

#-------------------------------------------------------------------------------
# stage-in: nodefiles

if [ -z "$TMPL" ]; then
  echo "#PJM --stgin-dir \"$TMPS/node ./node\""
else
  echo "#PJM --stgin-dir \"rank=* $TMPS/node %r:./node\""
fi

#-------------------------------------------------------------------------------
# stage-in: scripts

if [ -z "$TMPL" ]; then
  echo "#PJM --stgin \"$SCRP_DIR/fcst.sh ./runscp/fcst.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/src/* ./runscp/src/\""
  echo "#PJM --stgin \"$TMPS/config.all ./runscp/config.all\""
  echo "#PJM --stgin \"$TMPS/config.fcst ./runscp/config.fcst\""
  echo "#PJM --stgin \"$SCRP_DIR/scale_pp_topo.conf ./runscp/scale_pp_topo.conf\""
  echo "#PJM --stgin \"$SCRP_DIR/scale_pp_landuse.conf ./runscp/scale_pp_landuse.conf\""
  echo "#PJM --stgin \"$SCRP_DIR/scale_init.conf ./runscp/scale_init.conf\""
  echo "#PJM --stgin \"$SCRP_DIR/scale.conf ./runscp/scale.conf\""
else
  echo "#PJM --stgin \"rank=* $SCRP_DIR/fcst.sh %r:./runscp/fcst.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/src/* %r:./runscp/src/\""
  echo "#PJM --stgin \"rank=* $TMPS/config.all %r:./runscp/config.all\""
  echo "#PJM --stgin \"rank=* $TMPS/config.fcst %r:./runscp/config.fcst\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/scale_pp_topo.conf %r:./runscp/scale_pp_topo.conf\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/scale_pp_landuse.conf %r:./runscp/scale_pp_landuse.conf\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/scale_init.conf %r:./runscp/scale_init.conf\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/scale.conf %r:./runscp/scale.conf\""
fi

#===============================================================================

exit 0
