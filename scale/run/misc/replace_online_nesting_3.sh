#/bin/bash

list="
2015summer_18km
2015summer_18km_init
BDA_d1_15km_48p_bf10
BDA_d2_5km_64p_bf16
BDA_d3_100m_256p_bf40
BDA_d3_100m_64p_bf40
BDA_d3_100m_720p_bf40
BDA_d3_1km_9p_bf20
BDA_d4_100m_64p_bf40
BDA_d4_1km_4p_bf20
BDA_d4_250m_36p_bf40
BDA_d4_250m_720p_bf40
BDA_d4_500m_9p_bf30
realtime_r0020_d1
realtime_r0023_d1
realtime_r0043_d1
realtime_r0050_d1
realtime_v160405_d1
testcase_15km_48p_l36
testcase_45km_4p_l36_bdyprep
testcase_45km_4p_l36_gfsinput
testcase_45km_4p_l36_pnetcdf_bdyprep
testcase_45km_4p_l36_pnetcdf
testcase_45km_4p_l36_pnetcdf_ioarb
testcase_45km_4p_l72
testcase_PAWR_1km_4p
testcase_PAWR_1km_4p_pnetcdf
testcase_PAWR_5km_4p_pnetcdf
"

for j in $list; do
  if [ "$j" != "./replace.sh" ]; then

#    for i in $(ls ../config/${j}/config.nml.scale); do
#      echo $i

#      sed -i '/TIME_END_RESTART_OUT/a \/\n\n&PARAM_COMM_CARTESC_NEST\n!--ONLINE_DOMAIN_NUM--\n!--ONLINE_IAM_PARENT--\n!--ONLINE_IAM_DAUGHTER--\n ONLINE_BOUNDARY_USE_QHYD = .true.,\n ONLINE_AGGRESSIVE_COMM   = .true.,\n ONLINE_SPECIFIED_MAXRQ   = 10000,' $i

#      sed -i '/ RESTART_OUTPUT/c !--RESTART_OUTPUT--' $i
#    done

    for i in $(ls ../config/${j}/config.nml.scale_init); do
      echo $i

      sed -i '/BOUNDARY_UPDATE_DT/a !--MAKE_BOUNDARY--' $i
    done

  fi
done
