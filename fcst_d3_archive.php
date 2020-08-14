<?php

date_default_timezone_set('UTC');

$n_set = 4;
$sets[0] = 'd3';
$sets[1] = 'd2_d3';
$sets[2] = 'msm_d3';
$sets[3] = 'JMA_precip/nowcast_d3_grads/archive';
#$sets[4] = 'JMA_precip/nowcast_d3/realtime';
$set_name[0] = 'SCALE-LETKF Domain 3';
$set_name[1] = 'SCALE-LETKF Domain 2';
$set_name[2] = 'JMA MSM';
$set_name[3] = 'JMA Radar';

#$set_name[1] = 'JMA precip anal/fcst';
#$set_name[2] = 'JMA precip radar/nowcast';

$n_prod = 8;
$n_prod_msm = 2;
$prods[0] = 'sfc_prcp';
$prods[1] = 'sfc_uvt';
$prods[2] = '1000m_dbz';
$prods[3] = '1000m_uvw';
$prods[4] = '3000m_dbz';
$prods[5] = '3000m_uvw';
$prods[6] = '5000m_dbz';
$prods[7] = '5000m_uvw';
$prods_f[0] = 'Surface [Prcp]';
$prods_f[1] = 'Surface [Wind+Temp]';
$prods_f[2] = '1000m [Radar ref]';
$prods_f[3] = '1000m [Wind]';
$prods_f[4] = '3000m [Radar ref]';
$prods_f[5] = '3000m [Wind]';
$prods_f[6] = '5000m [Radar ref]';
$prods_f[7] = '5000m [Wind]';


$datadir = 'data';

$n_mem=2;
$cmems[0]='mdet' ;
$cmems[1]='mean' ;

$path_nodata_msm = $datadir."/nodata_Tokyo.png";
###

 $xmargin = 5 ;
 $xsize = 10 ;

###
#$set = 0;
#if (isset($_GET['set']) && ($_GET['set'] <= $n_set))
#  $set = $_GET['set'];

#--

#$prod_index = 0;
#$prod = $prods[$prod_index];
#if (isset($_GET['prod'])) {
#  for ($n = 0; $n < $n_prod; $n++) {
#    if ($_GET['prod'] == $prods[$n]) {
#      $prod = $_GET['prod'];
#      $prod_index = $n;
#   }
#  }
#}

#--

###


$btime_oldest = "2020-05-01";
$tstep = 43200 ;

$time_sec = strtotime($btime_oldest);

$n_btime=0 ;
while ($time_sec < date("U")) {
  $btimetest =date("YmdHis", $time_sec);
  if (file_exists('monitor/archive/monitor_plot_d3_'.$btimetest.'.png')) {
    $btimes[$n_btime] = $btimetest ;
    $btimes_f[$n_btime] = date("Y-m-d H", $time_sec)."Z" ;
    $n_btime++;
  }
    $time_sec=$time_sec+$tstep;
}
$btimes = array_reverse($btimes);
$btimes_f = array_reverse($btimes_f);


$b_index=0;
$basetime=$btimes[$b_index]; 
if (isset($_GET['basetime'])) {
  $basetime = $_GET['basetime'];
 for ($n=0;$n < $n_btime;$n++){
if( $basetime == $btimes[$n] ){
$b_index = $n ;
break;
}
}
}


 $json_xloc = file_get_contents('./json/archive/xlocs_d3_'.$basetime.'.json' );
 $json_yloc = file_get_contents('./json/archive/ylocs_d3_'.$basetime.'.json' );
 $json_itime = file_get_contents('./json/archive/itimes_d3_'.$basetime.'.json' );
 $json_ptime = file_get_contents('./json/archive/ptimes_d3_'.$basetime.'.json' );
 $XlocMarker = json_decode($json_xloc,true);
 $YlocMarker = json_decode($json_yloc,true);
 $InitTimes = json_decode($json_itime,true);
 $ParentTimes = json_decode($json_ptime,true);


$nowU = strtotime($basetime);
$nowdate = date("Y n/j",$nowU+9*3600);
$nowY = date("Y",$nowU);

$Uright = $nowU + 27 * 3600;

$s = 0;
$n_rtimes[0]=count($InitTimes[0]);

$n_imgs[$s]=0;
for ($r_index = 0; $r_index < $n_rtimes[$s]; $r_index++) {
if (isset($InitTimes[$s][$r_index])) {
for ($p = 0; $p < $n_prod; $p++) {
$path = "$datadir/" . $sets[$s] . "/ref_".$nowY.$ParentTimes[$s][$r_index]."0000/".$nowY.$InitTimes[$s][$r_index]. "0000/" .  $cmems[0]."/".$prods[$p] ; ###

$imgs = scandir($path);

$n_img = 0;

$img_index = 0;
foreach ($imgs as $img_file) {
if (ereg ("$prods[$p]_f([0-9]{6}).png", $img_file, $regs)) {
    $files[$n_img] = $img_file;
    $n_img++ ;
  }
}

for ($n = 0; $n < $n_img; $n++) {
  for ($c = 0;$c < $n_mem; $c++){
  $imgsrc[$c][$s][$r_index][$n][$p] =  "$datadir/" . $sets[$s] . "/ref_".$nowY.$ParentTimes[$s][$r_index]."0000/".$nowY.$InitTimes[$s][$r_index]. "0000/" .  $cmems[$c]."/".$prods[$p]. "/" . $files[$n];
}
}

if ($n_img > $n_imgs[$s]){
 $n_imgs[$s]=$n_img ;
}

}
}
}

  ### D2 and MSM

$s = 1 ;
$n_rtimes[$s]=count($InitTimes[$s]);

$n_imgs[$s]=0;
for ($r_index = 0; $r_index < $n_rtimes[$s]; $r_index++) {
if (isset($InitTimes[$s][$r_index])) {
for ($p = 0; $p < $n_prod_msm; $p++) {
$path = "$datadir/" . $sets[$s] . "/".$nowY.$InitTimes[$s][$r_index]. "0000/" .$cmems[0]."/".  $prods[$p] ; ###
$imgs = scandir($path);

$n_img = 0;

$img_index = 0;
foreach ($imgs as $img_file) {
if (ereg ("$prods[$p]_f([0-9]{6}).png", $img_file, $regs)) {
    $files[$n_img] = $img_file;
    $n_img++ ;
  }
}

for ($n = 0; $n < $n_img; $n++) {
  for ($c = 0;$c < $n_mem; $c++){
  $imgsrc[$c][$s][$r_index][$n][$p] = "$datadir/" . $sets[$s] . "/".$nowY.$InitTimes[$s][$r_index]. "0000/" .$cmems[$c]."/".  $prods[$p]."/" . $files[$n];
  }
}

if ($n_img > $n_imgs[$s]){
 $n_imgs[$s]=$n_img ;
}

}
}
}


$s = 2 ;
$n_rtimes[$s]=count($InitTimes[$s]);

$n_imgs[$s]=0;
for ($r_index = 0; $r_index < $n_rtimes[$s]; $r_index++) {
if (isset($InitTimes[$s][$r_index])) {
for ($p = 0; $p < $n_prod_msm; $p++) {
$path = "$datadir/" . $sets[$s] . "/".$nowY.$InitTimes[$s][$r_index]. "0000/" . $prods[$p] ; ###
$imgs = scandir($path);

$n_img = 0;

$img_index = 0;
foreach ($imgs as $img_file) {
if (ereg ("$prods[$p]_f([0-9]{6}).png", $img_file, $regs)) {
    $files[$n_img] = $img_file;
    $n_img++ ;
  }
}

for ($n = 0; $n < $n_img; $n++) {
  $imgsrc[0][$s][$r_index][$n][$p] = "$datadir/" . $sets[$s] . "/".$nowY.$InitTimes[$s][$r_index]. "0000/".  $prods[$p]."/" . $files[$n];
  $imgsrc[1][$s][$r_index][$n][$p] = "$datadir/" . $sets[$s] . "/".$nowY.$InitTimes[$s][$r_index]. "0000/".  $prods[$p]."/" . $files[$n];
}

if ($n_img > $n_imgs[$s]){
 $n_imgs[$s]=$n_img ;
}

}
}
}


for ($r_index = 0; $r_index < $n_rtimes[$s]; $r_index++) {
for ($n = 0; $n < $n_imgs[$s]; $n++) { 
for ($p = $n_prod_msm; $p < $n_prod; $p++) {
  $imgsrc[0][$s][$r_index][$n][$p] = "$path_nodata_msm";
  $imgsrc[1][$s][$r_index][$n][$p] = "$path_nodata_msm";
}
}
}


### radar/nowcast 
$s = 3 ;
$n_rtimes[$s]=1;
$p = 0 ;
$dpath = "$datadir/" . "$sets[$s]";
$dates = array_slice(scandir($dpath),2);
$n = 0;
foreach ($dates as $date){
#echo $date;
#echo substr($date,4,4);
#echo $InitTimes[$s][0];
#exit ;

if ( intval(substr($date,4,4)) >= intval(substr($InitTimes[$s][0],0,4)) ){
$path = "$datadir/" . "$sets[$s]". "/$date" ;
$imgs = scandir($path);


foreach ($imgs as $img_file) {
  if (ereg ("radar_([0-9]{11})000.png", $img_file, $regs) ) {  ### reduced to 10-min -- tori aezu --
    if (intval(substr($regs[1],4,6)) >= intval($InitTimes[$s][0]) ) {
    $imgsrc[0][$s][0][$n][$p] = "$path/" .$img_file;
    $imgsrc[1][$s][0][$n][$p] = "$path/" .$img_file;
    $n++ ;
#    echo $n.' '.count($XlocMarker[$s][0]).' '.$img_file."\n";
    if ( $n == count($XlocMarker[$s][0]) ) break 2 ;
  }
  }
}
}
}
 $n_imgs[$s]=$n ;

/*

#--

### anal/fcst 
$s = 1 ;
$p = 0;
$path = "$datadir/" . "$sets[$s]";
$imgs = scandir($path);
$n = 0;
foreach ($imgs as $img_file) {
  if (ereg ("anal_([0-9]{14}).png", $img_file, $regs) || ereg  ("fcst_([0-9]{14}).png", $img_file, $regs)) {
   if (intval(substr($regs[1],4,6)) >= intval($InitTimes[$s][0]) ) {
    $imgsrc[$s][0][$n][0] = "$path/" .$img_file;
    $n++ ;
  }
  }
}
$n_rtimes[$s] = 1 ;
$n_imgs[$s] = $n ;

### radar/nowcast 
$s = 2 ;
$p = 0 ;
$path = "$datadir/" . "$sets[$s]";
$imgs = scandir($path);
$n = 0;
foreach ($imgs as $img_file) {
  if (ereg ("radar_([0-9]{14}).png", $img_file, $regs)) {
   if (intval(substr($regs[1],4,6)) >= intval($InitTimes[$s][0]) ) {
    $imgsrc[$s][0][$n][0] = "$path/" .$img_file;
    $n++ ;
  }
  }
}
foreach ($imgs as $img_file) {
  if (ereg  ("nowcast_([0-9]{14}).png", $img_file, $regs)) {
    if (intval(substr($regs[1],4,6)) >= intval($InitTimes[$s][0]) ) {
    $imgsrc[$s][0][$n][0] = "$path/" .$img_file;
    $n++ ;
  }
  }
}
$n_rtimes[$s] = 1 ;
$n_imgs[$s] = $n ;

*/

?>


<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Cache-Control" content="no-cache">
<meta http-equiv="Refresh" content="300">
<title>Model forecasts quicklook for domain 3</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
<!--
var interval = 400 ;
var now_index = 0 ;
var now_indexb = 0;
var now_set = 1;
var now_setb = 1;
var now_mem = 0;
var now_memb = 0;
<?php 
if ( $n_imgs[0] > 0 ){
echo "var now_set  = 0;\n";
echo "var now_setb  = 0;\n";
}
?>
var now_nr = 0;
var now_nrb = 0;
var now_prod  = 0;
var now_prodb = 0;
var isload = false ;
var is_animation = false ;

var set_name = new Array();

<?php
  for ($s = 0; $s < $n_set; $s++) {
       echo "set_name[$s] = '" . $set_name[$s] . "';\n";  
}
?>


var imgsrc = new Array();

<?php
for ($c = 0; $c < $n_mem; $c++) {
 echo "imgsrc[$c] = new Array();\n";
  for ($s = 0; $s < $n_set; $s++) {
   if ($n_imgs[$s] > 0) {
    echo "imgsrc[$c][$s] = new Array();\n";
    for ($r_index=0; $r_index < $n_rtimes[$s]; $r_index++) {
      echo "imgsrc[$c][$s][$r_index] = new Array();\n";
     for ($n_index=0; $n_index < $n_imgs[$s]; $n_index++) {
       echo "imgsrc[$c][$s][$r_index][$n_index] = new Array();\n";
     foreach ($imgsrc[$c][$s][$r_index][$n_index] as $n => $isrc) {
       echo "imgsrc[$c][$s][$r_index][$n_index][$n] = '" . $isrc . "';\n";
    }

   }
  }
}
}
}
?>


var XlocMarker = new Array();
var YlocMarker = new Array();
var InitTimes = new Array();
var ParentTimes = new Array();

<?php

  for ($s = 0; $s < $n_set; $s++) {
    echo "XlocMarker[$s] = new Array();\n";
    echo "YlocMarker[$s] = new Array();\n";
    echo "InitTimes[$s] = new Array();\n";
    echo "ParentTimes[$s] = new Array();\n";
       foreach ($InitTimes[$s] as $n => $isrc) {
      echo "InitTimes[$s][$n] = '" . $isrc . "';\n";
}
       foreach ($ParentTimes[$s] as $n => $isrc) {
      echo "ParentTimes[$s][$n] = '" . $isrc . "';\n";
}
    for ($r_index=0; $r_index < $n_rtimes[$s]; $r_index++) {
      echo "XlocMarker[$s][$r_index] = new Array();\n";
      echo "YlocMarker[$s][$r_index] = new Array();\n";
     foreach ($XlocMarker[$s][$r_index] as $n => $isrc) {
       echo "XlocMarker[$s][$r_index][$n] = Number(" . $isrc . ");\n";
    }
     foreach ($YlocMarker[$s][$r_index] as $n => $isrc) {
       echo "YlocMarker[$s][$r_index][$n] = '" . $isrc . "';\n";
    }
  }
}


?>




function onload_script() {
  isload = true ;
  isover = false ;

  document.getElementById('title_set1').textContent = set_name[now_set] ;
  document.getElementById('title_set2').textContent = set_name[now_setb] ;
  document.getElementById('img_area').src = imgsrc[now_mem][now_set][now_nr][now_index][now_prod] ;
  document.getElementById('img_area2').src = imgsrc[now_memb][now_setb][now_nrb][now_indexb][now_prodb] ;
  marker_red();
  marker_blue();


<?php
if ($n_img > 0)
  echo "  change_image() ;";
?>



<?php
echo "  set = 0;\n";
if (isset($_GET['anim']) && ($_GET['anim'] == 'true') && ($n_img > 0)) {
  echo "  animation_image() ;\n";
}
?>
}

function change_prod() {
  if (isload) {
    var prod_obj = document.getElementById('prod') ;
    var prod2_obj = document.getElementById('prod2') ;
    now_prod = prod_obj.selectedIndex ;
    now_prodb = prod2_obj.selectedIndex ;
    change_image();
  }
}


function change_image() {
  if (isload) {
console.log(now_mem + ' ' + now_set + ' ' + now_nr + ' ' + now_index + ' ' + now_prod );
    document.getElementById('img_area').src = imgsrc[now_mem][now_set][now_nr][now_index][now_prod] ;
    document.getElementById('img_area2').src = imgsrc[now_memb][now_setb][now_nrb][now_indexb][now_prodb] ;
    marker_red();
    marker_blue();
  }
}

function animation_image() {
  if (isload) {
    var ani_obj = document.getElementById('animation') ;
    if (is_animation == false) {
      is_animation = true ;
      ani_obj.value = 'Stop' ;
      play() ;
    }
    else {
      is_animation = false ;
      ani_obj.value = 'Animation' ;
    }
  }
}

function animation_image() {
  if (isload) {
    var ani_obj = document.getElementById('animation') ;
    if (is_animation == false) {
      is_animation = true ;
      ani_obj.value = 'Stop' ;
      play() ;
    }
    else {
      is_animation = false ;
      ani_obj.value = 'Animation' ;
    }
  }
}


function play() {
  if (is_animation == true) {
    if (XlocMarker[now_set][now_nr][now_index+1] === undefined || XlocMarker[now_setb][now_nrb][now_indexb+1] === undefined ){
    now_index = 0 ;
    now_indexb = 0 ;
    while (XlocMarker[now_set][now_nr][now_index] <  XlocMarker[now_setb][now_nrb][now_indexb] ) {
      now_index++ ;   
}
    while (XlocMarker[now_set][now_nr][now_index] >  XlocMarker[now_setb][now_nrb][now_indexb] ) {
      now_indexb++ ;   
}
}
console.log('indx' + XlocMarker[now_set][now_nr][now_index] + ' indx2 ' + XlocMarker[now_setb][now_nrb][now_indexb] );
    if (XlocMarker[now_set][now_nr][now_index] < XlocMarker[now_setb][now_nrb][now_indexb] ) {
        now_index++ ;
} else if (XlocMarker[now_set][now_nr][now_index] > XlocMarker[now_setb][now_nrb][now_indexb] ) {
        now_indexb++ ;
} else {
    now_index++ ;
    now_indexb++ ;
}
    
    change_image() ;
    setTimeout('play();', interval) ;
  }
}


function change_btime() {
  if (isload) {
    var btime_obj = document.getElementById('btime') ;
    location.href = '?basetime=' + btime_obj.value  ;
  }
}

function change_mem(ele) {
  var myid = ele.id ;
  nums = myid.split("_");
  now_indx = nums[2];
  if (now_indx == 1){
  now_mem = nums[3];
  }
  if (now_indx == 2){
  now_memb = nums[3];
  }
  change_image();
 <?php 
  for ($i=0;$i<$n_mem;$i++){
   echo "document.getElementById('button_mem_' + now_indx + '_' + ".$i.").style.backgroundColor = '#FCFCFC';";
  }
?>
  ele.style.backgroundColor = "#A4A4A4";
}

//function play() {
//  if (is_animation == true) {

//   if (XlocMarker[now_set][now_nr][now_index] < XlocMarker[now_setb][now_nrb][now_indexb] ) {
//        now_index++ ;
//} else if (XlocMarker[now_set][now_nr][now_index] > XlocMarker[now_setb][now_nrb][now_indexb] ) {
//        now_indexb++ ;
//} else {
//    now_index++ ;
//    now_indexb++ ;
//}
//    
//    if (XlocMarker[now_set][now_nr][now_index+1] === undefined || XlocMarker[now_setb][now_nrb][now_indexb+1] === undefined ){
//    now_index = 0 ;
//    now_indexb = 0 ;
//    while (XlocMarker[now_set][now_nr][now_index] <  XlocMarker[now_setb][now_nrb][now_indexb] ) {
//      now_index++ ;   
//}
//    while (XlocMarker[now_set][now_nr][now_index] >  XlocMarker[now_setb][now_nrb][now_indexb] ) {
//      now_indexb++ ;   
//}
//}
//    change_image() ;
//    setTimeout('play();', interval) ;
//  }
//}

function marker_red() {
      xmargin = 5 ;
      document.getElementById('marker_red').style.left = Number(XlocMarker[now_set][now_nr][now_index])-xmargin + 'px';
      document.getElementById('marker_red').style.top =  YlocMarker[now_set][now_nr][now_index] + 'px';
if ( now_set == 1 ){
      document.getElementById('marker_red').style.height = '60px';
}else{
      document.getElementById('marker_red').style.height = '30px';
}
}

function marker_blue() {
      xmargin = 5 ;
      if ( now_set == now_setb && now_nr && now_nrb && now_index == now_indexb      ) { xmargin = 4 }
      document.getElementById('marker_blue').style.left =
      Number(XlocMarker[now_setb][now_nrb][now_indexb]) -xmargin + 'px';
      document.getElementById('marker_blue').style.top =  YlocMarker[now_setb][now_nrb][now_indexb] + 'px';
if ( now_setb == 1 ){
      document.getElementById('marker_blue').style.height = '60px';
}else{
      document.getElementById('marker_blue').style.height = '30px';
}
}

function changeMapImage(ele,imgPath) {
  var myid = ele.id ;
  document.getElementById(myid).src = imgPath;
}

function change_Image_byMarker(ele) { 
  myid = ele.id ;
  nums = myid.split("_") ;
  now_set = nums[1] ;    
  now_nr = nums[2] ;    
  now_index = nums[3] ;    
  document.getElementById('title_set1').textContent = set_name[now_set] ;
  document.getElementById('img_area').src = imgsrc[now_mem][now_set][now_nr][now_index][now_prod] ;
  marker_red(); 
}
function change_Image2_byMarker(ele) { 
  myid = ele.id ;
  nums = myid.split("_") ;
  now_setb = nums[1] ;    
  now_nrb = nums[2] ;    
  now_indexb = nums[3] ;    
  document.getElementById('title_set2').textContent = set_name[now_setb] ;
  document.getElementById('img_area2').src = imgsrc[now_mem][now_setb][now_nrb][now_indexb][now_prodb] ;
  marker_blue(); 
}

// -->
</script>
</head>

<body onLoad="onload_script();">

<div align="left">

<?php echo $res." ".$out[0];?>


  <p class="p_title">Model forecasts quicklook for domain 3</p>
<?php
 if (intval($time_offset) < 0){
 echo  "<p class=\"p_title\"><span style=\"color:red\">past mode : ".$nowdate." </p>\n";
 }
 if (intval($time_offset) > 0){
 echo  "<p class=\"p_title\"><span style=\"color:red\">future mode : offset +".$time_offset." s </p>\n";
 } 
?>
      <p class="p_menu_title">Base Time (JST):</p>
      <select name="btime" id='btime' onChange="change_btime();" class="menu" style="background-color: #FFE0C6">
<?php
for ($n = 0; $n < $n_btime; $n++) {
  if ($n == $b_index)
    echo "        <option value=\"" . $btimes[$n] . "\" selected=\"selected\">" . $btimes_f[$n] . "</option>\n";
  else
    echo "        <option value=\"" . $btimes[$n] . "\">" . $btimes_f[$n] . "</option>\n";
}
?>
      </select>
</div>

<div align="left">
<table border="0" cellpadding="4" cellspacing="0"><tr>
  <td valign="top" style="text-align: left" width="300px">
    <form name="view" action="#">
      <p class="p_menu_title">Products:</p>
      <select name="prod" id='prod' onChange="change_prod();" class="menu_d3" style="background-color: #C0FFF0">
<?php
for ($n = 0; $n < $n_prod; $n++) {
  if ($n == $prod_index)
    echo "        <option value=\"" . $prods[$n] . "\" selected=\"selected\">" . $prods_f[$n] . "</option>\n";
  else
    echo "        <option value=\"" . $prods[$n] . "\">" . $prods_f[$n] . "</option>\n";
}
?>
      </select>
</td>
<td valign="top" style="text-align: right" width="300px">
<p class="p_menu_title" align="left">
<br><?php
for ($i=0;$i<$n_mem;$i++){
if ($i==$now_mem){
echo "<input id=\"button_mem_1_$i\" type=\"button\" value=\"$cmems[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #A4A4A4;\" onclick=\"change_mem(this)\"><br>";
}else{
echo "<input id=\"button_mem_1_$i\" type=\"button\" value=\"$cmems[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #FCFCFC;\" onclick=\"change_mem(this)\"><br>";
}
}
?>
</p>
</td>

  <td valign="bottom" style="text-align: center" width="200px">
<?php
if (array_sum($n_imgs) == 0)
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" disabled>\n";
else
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" onClick=\"animation_image();\">\n";
?>
    </form>
  </td>
<td valign="top" style="text-align: right" width="300px">
<p class="p_menu_title" align="right">
<br><?php
for ($i=0;$i<$n_mem;$i++){
if ($i==$now_memb){
echo "<input id=\"button_mem_2_$i\" type=\"button\" value=\"$cmems[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #A4A4A4;\" onclick=\"change_mem(this)\"><br>";
}else{
echo "<input id=\"button_mem_2_$i\" type=\"button\" value=\"$cmems[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #FCFCFC;\" onclick=\"change_mem(this)\"><br>";
}
}
?>
</p>
</td>

<td valign="top" style="text-align: right" width="300px">
    <form name="view" action="#">
      <p class="p_menu_title">Products:</p>
      <select name="prod" id='prod2' onChange="change_prod();" class="menu" style="background-color: #C0FFF0">
<?php
for ($n = 0; $n < $n_prod; $n++) {
  if ($n == $prod2_index)
    echo "        <option value=\"" . $prods[$n] . "\" selected=\"selected\">" . $prods_f[$n] . "</option>\n";
  else
    echo "        <option value=\"" . $prods[$n] . "\">" . $prods_f[$n] . "</option>\n";
}
?>
      </select>
</td>
</tr></table>



<table border="0" cellpadding="4" cellspacing="0"><tr>
  <td valign="top" style="text-align: center; padding-left: 0px"
  width="600px"><p id="title_set1" class="p_model" style="border-style: ridge;
  border-color:red; border-width:2px ;">SCALE-LETKF d3 mdet</p></td>
  <td valign="top" style="text-align: center; padding-left: 0px"
  width="600px"><p id="title_set2" class="p_model" style="border-style: ridge;
  border-color:blue; border-width:2px ;">SCALE-LETKF d3 mdet</p></td>
</tr>
<tr>

  <td valign="middle" style="text-align: left; padding-left: 0px" width="600px" >
<?php
$filepath = $imgsrc[0][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area\" src=\"$filepath\" border=\"0\" width=\"700px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area\" border=\"0\" width=\"700px\">\n";
?>
  </td>

  <td valign="middle" style="text-align: left; padding-left: 0px">

<?php
$filepath = $imgsrc[1][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area2\" src=\"$filepath\" border=\"0\" width=\"700px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area2\" border=\"0\" width=\"700px\">\n";
?>

  </td>
</tr></table>
</div>

<br>
<br>


<div align="left" style="position:relative; top:0px; left:0px;">
  <img id="monitor_campus" src="monitor/archive/monitor_plot_d3_<?php echo $basetime;?>.png" style="margin:0px;"> 



 <?php
  echo "<img id=\"marker_red\" src=\"./img/marker-d2.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:30px; width:".$xsize."px\";>\n";
  echo "<img id=\"marker_blue\" src=\"./img/marker-d2_blue.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:30px; width:".$xsize."px\";>\n";

 for ($s = $n_set-1 ; $s >= 0; $s--) {
 $hgtmrk="30px";
 if ($s == 1) $hgtmrk="60px" ;
 for ($nr = 0; $nr < $n_rtimes[$s]; $nr++) {
 for ($n = 0; $n < $n_imgs[$s]; $n++) {
  if (isset($XlocMarker[$s][$nr][$n]) ) {  
  echo  " <img id=\"markerposition_".$s."_".$nr."_".$n.   "\" src=\"./img/marker_trans_new.png\" ; style=\"position:absolute;  top: ".$YlocMarker[$s][$nr][$n]."px; left:".($XlocMarker[$s][$nr][$n]-$xmargin)."px; margin:0px; height:".$hgtmrk."; width:".$xsize."px;\" onclick=\"change_Image_byMarker(this);\" oncontextmenu=\"change_Image2_byMarker(this); return false\" onmouseover=\"changeMapImage(this,'./img/marker-d2_grey.png');\" onmouseout=\"changeMapImage(this,'./img/marker_trans_new.png');\" >  \n"; 
  };  
  };
  };
  };
 ?>

</div>

</body>
</html>
