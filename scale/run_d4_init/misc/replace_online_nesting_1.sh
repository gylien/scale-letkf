#/bin/bash

#list="
#2015summer_18km
#2015summer_18km_init
#BDA_d1_15km_48p_bf10
#realtime_r0020_d1
#realtime_r0023_d1
#realtime_r0043_d1
#realtime_r0050_d1
#realtime_v160405_d1
#"

list="
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
testcase_PAWR_1km_4p
testcase_PAWR_1km_4p_pnetcdf
testcase_PAWR_5km_4p_pnetcdf
"

for j in $list; do
  if [ "$j" != "./replace.sh" ]; then

    for i in $(ls ../config/${j}/config.nml.scale_pp ../config/${j}/config.nml.scale_init); do
      echo $i

      if [ "$i" = "../config/${j}/config.nml.scale_init" ]; then
        for k in $(ls ../config/${j}/config.main*); do
#          sed -i "/DATA_BDY_SCALE=/a DATA_BDY_SCALE_PRC_NUM_X=\nDATA_BDY_SCALE_PRC_NUM_Y=" $k

          OFFLINE_PARENT_PRC_NUM_X=$(grep 'OFFLINE_PARENT_PRC_NUM_X' $i | cut -d '=' -f2 | cut -d ',' -f1 | tr -d '[:space:]')
          OFFLINE_PARENT_PRC_NUM_Y=$(grep 'OFFLINE_PARENT_PRC_NUM_Y' $i | cut -d '=' -f2 | cut -d ',' -f1 | tr -d '[:space:]')
          sed -i "/DATA_BDY_SCALE=/a DATA_BDY_SCALE_PRC_NUM_X=${OFFLINE_PARENT_PRC_NUM_X}\nDATA_BDY_SCALE_PRC_NUM_Y=${OFFLINE_PARENT_PRC_NUM_Y}" $k
        done
      fi
      sed -i "/OFFLINE_PARENT_PRC_NUM_X/c !--OFFLINE_PARENT_PRC_NUM_X--" $i
      sed -i "/OFFLINE_PARENT_PRC_NUM_Y/c !--OFFLINE_PARENT_PRC_NUM_Y--" $i

    done

  fi
done
