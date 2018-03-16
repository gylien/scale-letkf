#/bin/bash

#list="
#2015summer_18km
#2015summer_18km_init
#BDA_d1_15km_48p_bf10
#BDA_d2_5km_64p_bf16
#BDA_d3_1km_9p_bf20
#BDA_d4_100m_64p_bf40
#BDA_d4_1km_4p_bf20
#BDA_d4_250m_36p_bf40
#BDA_d4_250m_720p_bf40
#BDA_d4_500m_9p_bf30
#realtime_r0023_d1
#realtime_r0043_d1
#realtime_r0050_d1
#testcase_45km_4p_l36_pnetcdf
#"
# DET_RUN_CYCLED = .true.

#list="
#testcase_PAWR_5km_4p_pnetcdf
#"
# DET_RUN_CYCLED = .false.

list="
BDA_d3_100m_256p_bf40
BDA_d3_100m_64p_bf40
BDA_d3_100m_720p_bf40
realtime_r0020_d1
realtime_v160405_d1
testcase_15km_48p_l36
testcase_45km_4p_l36_bdyprep
testcase_45km_4p_l36_gfsinput
testcase_45km_4p_l36_pnetcdf_bdyprep
testcase_45km_4p_l36_pnetcdf_ioarb
testcase_45km_4p_l72
testcase_PAWR_1km_4p
testcase_PAWR_1km_4p_pnetcdf
"
# No DET_RUN_CYCLED

for j in $list; do
  if [ "$j" != "./replace.sh" ]; then

    for i in $(ls ../config/${j}/config.nml.letkf ../config/${j}/config.nml.obsope); do
      echo $i

      sed -i '/&PARAM_ENSEMBLE/,/\//d' $i
      sed -i '/&PARAM_LETKF_PRC/,/\//d' $i
    done

    cat > ../config/${j}/config.nml.ensmodel <<EOF
&PARAM_ENSEMBLE
!--MEMBER--
!--MEMBER_RUN--
!--MEMBER_ITER--
!--CONF_FILES--
!--CONF_FILES_SEQNUM--
!--DET_RUN--
! DET_RUN_CYCLED = .true.,
/

&PARAM_PROCESS
!--PPN--
!--MEM_NODES--
!--NUM_DOMAIN--
!--PRC_DOMAINS--
! COLOR_REORDER = .true.,
/

EOF

  fi
done
