#!/bin/sh

ctime=$1

if [ -z $ctime ];then
 echo "specify YYYYMMDDHHMMSS" 
 exit
fi

for item in u v tv q ;do
  ./draw_profile $ctime $item
done
mv location_*.txt map
mogrify -trim map_*.png 
mogrify -trim profile_*.png 
mv map_*.png map

rm figs/*.png
mv profile_*.png figs
