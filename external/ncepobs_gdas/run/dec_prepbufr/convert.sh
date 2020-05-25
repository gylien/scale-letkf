#!/bin/bash -l

source ~/.bashrc

tt="$1"
prepbufrdir="$2"
outdir="$3"

#tt='2015-05-07 06'
#prepbufrdir='/data7/gylien/realtime/ncepobs_gdas/2015050706'
#outdir='/data7/gylien/realtime/ncepobs_gdas_letkf/2015050706'

BUFRBIN="/work/hp150019/c24140/scale-letkf-rt/external/lib/bufrlib/10.1.0_intel/bin"

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

timef=$(date -ud "$tt" '+%Y-%m-%d %H')
timef2=$(date -ud "$tt" '+%Y%m%d%H%M%S')
timef3=$(date -ud "$tt" '+%Y%m%d%H')
echo "[$timef]"

rm -f prepbufr.in fort.90


#file="$prepbufrdir/prepbufr.${timef3}"
file="../../${timef3}/prepbufr.${timef3}"

#echo $file
#exit

wc -c "$file" | $BUFRBIN/grabbufr "$file" prepbufr.in

[ -f fort.90 ] && rm fort.90


###rm exec_dec_prepbufr.sh.e*
###rm exec_dec_prepbufr.sh.o*

###cp template.sh  exec_dec_prepbufr.sh
###echo "mkdir -p $outdir" >> exec_dec_prepbufr.sh
###echo "mv fort.90 $outdir/obs_${timef2}.dat" >> exec_dec_prepbufr.sh
##### must be executed on calc node ; otherwise AVX512 is not supported
###pjsub ./exec_dec_prepbufr.sh

module load hdf5/1.10.5
module load netcdf/4.7.0
module load netcdf-fortran/4.4.5

DECBIN="/work/hp150019/share/SCALE-LETKF-rt/scale_noAVX/scale-letkf/scale/obs/dec_prepbufr"

[ ! -f dec_prepbufr_fixed ] && ln -s $DECBIN dec_prepbufr_fixed

./dec_prepbufr_fixed

touch fort.90
mkdir -p $outdir
mv fort.90 $outdir/obs_${timef2}.dat
