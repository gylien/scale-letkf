#!/bin/sh

grads_version=grads


gradsfile_dir="$1"
gradsfile_base_s=msm_s
gradsfile_base_p=msm_p

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}
YYYYMMDDHH=`basename $gradsfile_dir`
YYYY=${YYYYMMDDHH:0:4}
MM=${YYYYMMDDHH:4:2}
DD=${YYYYMMDDHH:6:2}
HH=${YYYYMMDDHH:8:2}
MMSS=0000
MF=`date -R -d "$YYYY$MM$DD" | awk '{print $3}'`


TIMEF=${HH}Z${DD}${MF}${YYYY}
ncfile_s=MSM${YYYYMMDDHH}S.nc
ncfile_p=MSM${YYYYMMDDHH}P.nc
gradsfile_s=$gradsfile_base_s.${YYYYMMDDHH}${MMSS}
gradsfile_p=$gradsfile_base_p.${YYYYMMDDHH}${MMSS}


sed  -e "s/<--GRADSFILE_S-->/${gradsfile_s}/g" -e "s/<--NCFILE_S-->/${ncfile_s}/g" template_make_s.gs  > $gradsfile_dir/make_s.gs
sed  -e "s/<--GRADSFILE_P-->/${gradsfile_p}/g" -e "s/<--NCFILE_P-->/${ncfile_p}/g" template_make_p.gs  > $gradsfile_dir/make_p.gs

sed -e "s/<--GRADSFILE_S-->/${gradsfile_s}/g" -e "s/<--TIMEF-->/${TIMEF}/g" template_ctl_s.ctl  > $gradsfile_dir/$gradsfile_s.ctl
sed -e "s/<--GRADSFILE_P-->/${gradsfile_p}/g" -e "s/<--TIMEF-->/${TIMEF}/g" template_ctl_p.ctl  > $gradsfile_dir/$gradsfile_p.ctl
#sed -e "" template_ctl_p.ctl  > $gradsfile_dir/$gradsfile_p.ctl

cd $gradsfile_dir
$grads_version -bcl make_s.gs 
$grads_version -bcl make_p.gs 


