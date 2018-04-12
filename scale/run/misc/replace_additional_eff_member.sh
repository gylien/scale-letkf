#/bin/bash

list="$(grep -rn 'RESTART_OUT_ADDITIONAL_BASENAME' ../config | cut -d ':' -f1)"

for i in $list; do
  echo $i

  sed -i '/RESTART_OUT_ADDITIONAL_BASENAME/a !--RESTART_OUT_ADDITIONAL_EFF_MEMBER--' $i
done
