#/bin/bash

#list=$(ls */config.nml.scale*)

#list="
#testcase_PAWR_1km_4p_pnetcdf_bdyprep/config.nml.scale
#testcase_PAWR_1km_4p_pnetcdf_ioarb/config.nml.scale
#testcase_PAWR_500m_9p/config.nml.scale
#testcase_PAWR_500m_9p_pnetcdf/config.nml.scale
#testcase_PAWR_500m_9p_pnetcdf_bdyprep/config.nml.scale
#testcase_PAWR_500m_9p_pnetcdf_ioarb/config.nml.scale
#testcase_PAWR_100m_720p_pnetcdf_ioarb/config.nml.scale
#testcase_PAWR_250m_720p_pnetcdf_ioarb/config.nml.scale
#"

list="
testcase_PAWR_1km_4p_pnetcdf_bdyprep/config.nml.scale*
testcase_PAWR_1km_4p_pnetcdf_ioarb/config.nml.scale*
testcase_PAWR_500m_9p/config.nml.scale*
testcase_PAWR_500m_9p_pnetcdf/config.nml.scale*
testcase_PAWR_500m_9p_pnetcdf_bdyprep/config.nml.scale*
testcase_PAWR_500m_9p_pnetcdf_ioarb/config.nml.scale*
testcase_PAWR_100m_720p_pnetcdf_ioarb/config.nml.scale*
testcase_PAWR_250m_720p_pnetcdf_ioarb/config.nml.scale*
"

for i in $list; do
  if [ "$i" != "./replace.sh" ]; then
    echo $i

#    sed -i '/OBSNUM=/,/OBSNAME\[2/c OBSNUM=1\nOBSNAME[1]=obs' $i
#    sed -i '/OBSNUM=/,/OBSNAME\[2/c OBSNUM=1\nOBSNAME[1]=radar' $i
#    sed -i "/OBS_IN_FORMAT =/c \ OBS_IN_FORMAT = 'PREPBUFR'," $i/config.nml.obsope
#    sed -i "/OBS_IN_FORMAT =/c \ OBS_IN_FORMAT = 'RADAR'," $i/config.nml.obsope
#    sed -i "/!--IO_AGGREGATE--/c /\n\n&PARAM_FILE\n!--FILE_AGGREGATE--" $i
#    DT=$(grep 'TIME_DT_ATMOS_PHY_TB ' $i | cut -d '=' -f2)
#    sed -i "/TIME_DT_ATMOS_PHY_TB_UNIT/a \ TIME_DT_ATMOS_PHY_BL       =${DT}\n\ TIME_DT_ATMOS_PHY_BL_UNIT  = \"SEC\"," $i
#    sed -i "/ATMOS_PHY_TB_TYPE/c \ ATMOS_PHY_TB_TYPE = \"NONE\",\n\ ATMOS_PHY_BL_TYPE = \"MYNN\"," $i
#    sed -i "/ATMOS_PHY_TB_TYPE/a \ ATMOS_PHY_BL_TYPE = \"NONE\"," $i
#    sed -i "/ATMOS_PHY_TB_TYPE/c \ ATMOS_PHY_TB_TYPE = \"SMAGORINSKY\",\n\ ATMOS_PHY_BL_TYPE = \"MYNN\"," $i

    sed -i "/ATMOS_DYN_enable_coriolis/c \ ATMOS_DYN_coriolis_type              = \"SPHERE\"," $i
    sed -i "s/HISTORY_OUTPUT_STEP0/FILE_HISTORY_OUTPUT_STEP0/" $i
    sed -i "/&PARAM_HISTORY/c &PARAM_FILE_HISTORY" $i
    sed -i "s/HISTORY_DEFAULT_BASENAME/FILE_HISTORY_DEFAULT_BASENAME/" $i
    sed -i "s/HISTORY_DEFAULT_TINTERVAL/FILE_HISTORY_DEFAULT_TINTERVAL/" $i
    sed -i "s/HISTORY_DEFAULT_TUNIT/FILE_HISTORY_DEFAULT_TUNIT/" $i
    sed -i "s/HISTORY_DEFAULT_TAVERAGE/FILE_HISTORY_DEFAULT_TAVERAGE/" $i
    sed -i "s/HISTORY_DEFAULT_DATATYPE/FILE_HISTORY_DEFAULT_DATATYPE/" $i
    sed -i "s/History_DEFAULT_ZCOORD/FILE_HISTORY_DEFAULT_ZCOORD/" $i
    sed -i "s/HISTORY_OUTPUT_STEP0/FILE_HISTORY_OUTPUT_STEP0/" $i
    sed -i "/&PARAM_HIST/,+3d" $i
    sed -i "s/HISTITEM item/HISTORY_ITEM name/g" $i
    sed -i "/&HISTORY_ITEM name='RHL'/d" $i
    sed -i "/&HISTORY_ITEM name='RHI'/d" $i
    sed -i "/&PARAM_NEST/c &PARAM_COMM_CARTESC_NEST" $i
    sed -i "s/NEST_INTERP_LEVEL/COMM_CARTESC_NEST_INTERP_LEVEL/g" $i

    sed -i "/&PARAM_GRID/c &PARAM_ATMOS_GRID_CARTESC" $i
    sed -i "/&PARAM_INDEX/c &PARAM_ATMOS_GRID_CARTESC_INDEX" $i
    sed -i "/&PARAM_LAND_GRID/c &PARAM_LAND_GRID_CARTESC" $i
    sed -i "/&PARAM_LAND_INDEX/c &PARAM_LAND_GRID_CARTESC_INDEX" $i
    sed -i "/&PARAM_URBAN_GRID/c &PARAM_URBAN_GRID_CARTESC" $i
    sed -i "/&PARAM_URBAN_INDEX/c &PARAM_URBAN_GRID_CARTESC_INDEX" $i
    sed -i "/&PARAM_MAPPROJ/c &PARAM_MAPPROJECTION" $i
    sed -i "s/MPRJ_/MAPPROJECTION_/g" $i

  fi
done
