#!/bin/bash -l

SRCDIR="/work/jh150019/share/SCALE-LETKF-rt/result/ope/d4_500m/dafcst_img"

hostname="weather"
www_DATA="/srv/www/html/nowcast/scale-letkf_saitama/"
www_ARCHIVE="/srv/www_rw/data/scale-letkf_saitama/archive"
json_list="dataformat.json dataformat.json_en"

name="rainrate-gmap"
nmax_anal=21  ### 10 min
nmax_fcst=60  ### 30 min

fcstlen=1800

limit_waitsec=300

startTimef=$1
stopTimef=$2

function array2json(){
  arrayin=("$@")
  jout="[`printf '"%q",' "${arrayin[@]}"`"
  jout="${jout%,}]"
  jout=`echo $jout | sed -e 's#"##g' `
  echo "$jout"
}
function array2jsons(){
  arrayin=("$@")
  len=${#arrayin[@]}
  len=`expr $len - 1`
  jout='["'${arrayin[0]}'"'
  for i in `seq $len`;do
    jout="$jout"',"'${arrayin[$i]}'"'
  done
  echo "$jout"]
}


targetTimef=$startTimef

nows=`date -u +%s`

rm -rf work
mkdir -p work
cd work

nowb=$nows
nows=`date -u +%s`
echo 'prep workdir : '`expr $nows - $nowb`

while [ `date -ud "$targetTimef" +%s` -le `date -ud "$stopTimef" +%s` ] ;do

rm ${name}_*.png 

iwaitsec=0
while [ $iwaitsec -lt $limit_waitsec ] ;do

 echo "try" ${targetTimef} ...
 
  ffile_prefix="fcst_dbz_"
  ffile_suffix="_FT`printf %04d $fcstlen`s_z01957m.png"
  ffiles=`ls -x $SRCDIR/${ffile_prefix}*${ffile_suffix}`

  testTimef=`date -ud "- 30 second $targetTimef" +"%Y-%m-%d %H:%M:%S"` 

  iflag=0
  for fitem in $ffiles;do
    ffile=`basename $fitem`
    len=${#ffile_prefix}
    testtime=${ffile:${len}:15}
    testTimef="${testtime:0:4}-${testtime:4:2}-${testtime:6:2} ${testtime:9:2}:${testtime:11:2}:${testtime:13:2}"
    tests=`date -ud "${testTimef}" +%s`
  if [ -f $SRCDIR/$ffile ] && [ $tests -ge `date -ud "$targetTimef" +%s` ] ; then
    iflag=1
    timestamp=`date -ud "$testTimef" +"%Y%m%d-%H%M%S"`
    indx=0
    list_flagforecast=()
    list_image_timestamp=()
    list_image_timestamp_en=()
    list_image_unixtime=()

  for i in `seq $nmax_anal`;do
    diff=`expr 30 \* \( $nmax_anal - $i \)`
    timestampf=`date -ud "- $diff second $testTimef" +"%Y-%m-%d %H:%M:%S"`
    timestamp=`date -ud "$timestampf" +"%Y%m%d-%H%M%S"`
    timestamps=`date -ud "$timestampf" +"%s"`
    timestampl=`date -ud "9 hour $timestampf" +"%Y/%m/%d %H:%M:%S"`
    afile=$SRCDIR/anal_dbz_${timestamp}_z01957m.png 
    if [ -f $afile ];then
      indx=`expr $indx + 1`
      ln -sf $afile ${name}_`printf %04d $indx`.png 
      list_flagforecast=("${list_flagforecast[@]}" false)
      list_image_timestamp=("${list_image_timestamp[@]}" "解析 $timestampl")
      list_image_timestamp_en=("${list_image_timestamp[@]}" "analysis $timestampl")
      list_image_unixtime=("${list_image_unixtime[@]}" $timestamps)
    fi
  done

    timestamp=`date -ud "$testTimef" +"%Y%m%d-%H%M%S"`
  for i in `seq $nmax_fcst`;do
    diff=`expr 30 \* $i`
    timestampf=`date -ud "$diff second $testTimef" +"%Y-%m-%d %H:%M:%S"`
    timestamps=`date -ud "$timestampf" +"%s"`
    timestampl=`date -ud "9 hour $timestampf" +"%Y/%m/%d %H:%M:%S"`
    ffile=$SRCDIR/fcst_dbz_${timestamp}_FT`printf %04d $diff`s_z01957m.png 

    if [ -f $ffile ];then
    indx=`expr $indx + 1`
      ln -sf $ffile ${name}_`printf %04d $indx`.png 
      list_flagforecast=("${list_flagforecast[@]}" true)
      list_image_timestamp=("${list_image_timestamp[@]}" "予測 $timestampl")
      list_image_timestamp_en=("${list_image_timestamp[@]}" "forecast $timestampl")
      list_image_unixtime=("${list_image_unixtime[@]}" $timestamps)
    fi
  done

#nowb=$nows
#nows=`date -u +%s`
#echo 'make link : '`expr $nows - $nowb`

    timestampf=`date -ud "$testTimef" +"%Y/%m/%d %H:%M:%S"`
    timestamps=`date -ud "$testTimef" +"%s"`
    timestampl=`date -ud "9 hour $timestampf" +"%Y/%m/%d %H:%M:%S"`

json=$(cat <<EOS
  {
"forecast_flag_list": `array2json "${list_flagforecast[@]}"`,
"image_timestamp_list": `array2jsons "${list_image_timestamp_en[@]}"`,
"image_unixtime": `array2json "${list_image_unixtime[@]}"`,
"config": 
{"NUM_IMAGES":${#list_flagforecast[@]},"GRAYOUT_COLOR_OBS":"#000000","GRAYOUT_COLOR_FCST":"#0000FF","GRAYOUT_COLOR_INIT":"#000000","GRAYOUT_OPACITY":0.2,"RAINRATE_IMAGE_PREFIX":"/nowcast/scale-letkf_saitama/raw/rainrate-gmap_","COORDINATES_SW":[35.500955985933494,139.16510747108475],"COORDINATES_NE":[36.22041101372449,140.052849071884],"MAP_CENTER":[35.84349792480469,139.60897827148438],"TURNOFF_AUTO_UPDATE":1800000,"AUTO_UPDATE_FREQENCY":30000,"INITTIME":"${timestampl}","INITUNIXTIME":${timestamps},"OBSERVATORY":"SAITAMA"
  }
  }
EOS
  )
  echo "$json" > dataformat.json_en

json=$(cat <<EOS
  {
"forecast_flag_list": `array2json "${list_flagforecast[@]}"`,
"image_timestamp_list": `array2jsons "${list_image_timestamp[@]}"`,
"image_unixtime": `array2json "${list_image_unixtime[@]}"`,
"config": 
{"NUM_IMAGES":${#list_flagforecast[@]},"GRAYOUT_COLOR_OBS":"#000000","GRAYOUT_COLOR_FCST":"#0000FF","GRAYOUT_COLOR_INIT":"#000000","GRAYOUT_OPACITY":0.2,"RAINRATE_IMAGE_PREFIX":"/nowcast/scale-letkf_saitama/raw/rainrate-gmap_","COORDINATES_SW":[35.500955985933494,139.16510747108475],"COORDINATES_NE":[36.22041101372449,140.052849071884],"MAP_CENTER":[35.84349792480469,139.60897827148438],"TURNOFF_AUTO_UPDATE":1800000,"AUTO_UPDATE_FREQENCY":30000,"INITTIME":"${timestampl}","INITUNIXTIME":${timestamps},"OBSERVATORY":"SAITAMA"
  }
  }
EOS
  )
  echo "$json" > dataformat.json
  
#nowb=$nows
#nows=`date -u +%s`
#echo 'make json : '`expr $nows - $nowb`

  chmod 644 *.png 
  chmod 664 *json* 

  img_list=`ls -x ${name}*`
  NowTime=$timestamps

  echo "send" $timestampl "..." 

  scp -p ${img_list} ${hostname}:${www_DATA}/raw/
  scp -p ${json_list}  ${hostname}:${www_DATA}/json/

  ssh ${hostname} "mkdir -p ${www_ARCHIVE}/${NowTime}"

  scp -p ${img_list} ${hostname}:${www_ARCHIVE}/${NowTime}/
  scp -p ${json_list} ${hostname}:${www_ARCHIVE}/${NowTime}/

  ssh ${hostname} "echo ${NowTime} > ${www_ARCHIVE}/../latestInitUnixTime"
#  ssh ${hostname} "echo ${NowTime} > ${www_ARCHIVE}/../lastAcceptedUnixTime"
 iwaitsec=$limit_waitsec

#nowb=$nows
#nows=`date -u +%s`
#echo 'send : '`expr $nows - $nowb`

  fi

done ### ffile

if [ $iflag == 0 ];then
  sleep 2
  iwaitsec=`expr $iwaitsec + 2`
fi
done ### iwaitsec < limit_waitsec

targetTimef=`date -ud "30 second $testTimef" +"%Y-%m-%d %H:%M:%S"`
done ### while targetTimef < stopTimef

echo "== finished =="
cd -

