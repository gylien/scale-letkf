#/bin/bash

list="
2015summer_18km
2015summer_18km_init
realtime_r0023_d1
realtime_r0043_d1
realtime_r0050_d1
testcase_45km_4p_l36_gfsinput
"

for j in $list; do
  if [ "$j" != "./replace.sh" ]; then

    for i in $(ls ../config/${j}/config.nml.grads_boundary); do
      echo $i

      sed -i 's#--DIR--/bdysfc#--FNAME_SFC--#g' $i
      sed -i 's#--DIR--/bdyatm#--FNAME_ATMOS--#g' $i
      sed -i 's#--DIR--/bdyland#--FNAME_LAND--#g' $i
    done

  fi
done
