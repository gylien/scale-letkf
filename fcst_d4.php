<?php

date_default_timezone_set('UTC');
$max_rtime = 30;

$n_set = 2;
$sets[0] = 'd4/realtime';
$sets[1] = 'd4/realtime';
#$sets[1] = 'PAWR_nowcast/realtime';
$sets[2] = 'JMA_precip/nowcast_d4/realtime';

$set_name[0] = 'PAWR SCALE-LETKF';
$set_name[1] = 'PAWR observation';
#$set_name[1] = 'PAWR nowcast';
$set_name[2] = 'JMA radar / nowcast';

$n_prod = 0;
$prods[0] = 'radar_ref';
$prods[1] = 'sfc_wind';
$prods[2] = 'sfc_2mtemp';
$prods[3] = 'sfc_temp';
$prods_f[0] = 'Radar refrectivity';
$prods_f[1] = 'Surface [Wind]';
$prods_f[2] = 'Surface [T2m]';
$prods_f[3] = 'Surface [Tskin]';

$datadir = 'data';

###

 $json_xloc = file_get_contents('./json/xlocs_d4.json' );
 $json_yloc = file_get_contents('./json/ylocs_d4.json' );
 $json_itime = file_get_contents('./json/itimes_d4.json' );
 $XlocMarker = json_decode($json_xloc,true);
 $YlocMarker = json_decode($json_yloc,true);
 $InitTimes = json_decode($json_itime,true);

 $xmargin = 5 ;
 $xsize = 10 ;

 $n_step = count($XlocMarker[0][0]);


 $InitTimeU = strtotime(substr($InitTimes[0][0],0,4)."-".substr($InitTimes[0][0],4,2)."-".substr($InitTimes[0][0],6,2)."
 ".substr($InitTimes[0][0],8,2).":".substr($InitTimes[0][0],10,2).":".substr($InitTimes[0][0],12,2));

 for ( $i = 0 ; $i <= $n_step-1 ; $i++ ) {
 $TimeMarker[$i] = date('YmdHis',$InitTimeU + $i * 30) ;
}
#--
  $time_offset=file_get_contents('time_offset.txt');
  $nowU = date("U") + intval($time_offset);
  $nowdate = date("Y n/j",$nowU+9*3600);
  $nowY = date("Y",$nowU);
 

$s = 0 ;
$path = "$datadir/" . "$sets[$s]";
$levels = array_diff(scandir($path),array('.','..','.gitignore'));

$nlev=count($levels);

for ($lv = 0 ; $lv < $nlev ; $lv++ ) {
$n = 0;

$levels[$lv]=$levels[$lv+3]; ######
$levels_f[$lv]=substr($levels[$lv],2,5) ;

$imgs = scandir($path."/".$levels[$lv]);

foreach ($imgs as $img_file) {
  if (ereg ("anal_([0-9]{14}).png", $img_file, $regs) || ereg  ("fcst_([0-9]{14}).png", $img_file, $regs)) {
    if (intval(substr($regs[1],0,14)) >= intval($InitTimes[0][0]) ) {
    $vtime[$s][$n] = substr($regs[1],0,14);
    $imgsrc[$s][$lv][$n] = $path."/".$levels[$lv]."/".$img_file;
    $n++ ;
  }
  }
}
}
$n_imgs[$s] = count($imgsrc[$s][0]) ;



for ($s = 1; $s < $n_set ; $s++) {
for ($lv = 0 ; $lv < $nlev ; $lv++ ) {
$path = "$datadir/" . "$sets[$s]";
$imgs = scandir($path."/".$levels[$lv]);
$n = 0;
foreach ($imgs as $img_file) {
#  if (ereg ("radar_([0-9]{14}).png", $img_file, $regs)) {
  if (ereg ("obs_([0-9]{14}).png", $img_file, $regs)) {
    if (intval(substr($regs[1],0,14)) >= intval($InitTimes[0][0]) ) {
    $vtime[$s][$n] = substr($regs[1],0,14);
    $imgsrc[$s][$lv][$n] = $path."/".$levels[$lv]."/".$img_file;
    $n++ ;
  }
  }
}



#foreach ($imgs as $img_file) {
#  if (ereg  ("nowcast_([0-9]{14}).png", $img_file, $regs)) {
#    if (intval(substr($regs[1],0,14)) >= intval($InitTimes[0][0]) ) {
#    $vtime[$s][$n] = substr($regs[1],0,14);
#    $imgsrc[$s][$n] = "$path/" .$img_file;
#    $n++ ;
#  }
#  }
#}

} # lv

$n_imgs[$s] = $n ;
}

?>

<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Pragma" content="no-cache">
<meta http-equiv="Cache-Control" content="no-cache">
<meta http-equiv="Refresh" content="60">
<title>Model forecasts quicklook for domain 4 realtime </title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
<!--
var interval = 400 ;
var isload = false ;
var is_animation = false ;

var now_set = 0 ;

var n_step = <?php echo $n_step; ?> ;
var now_step = n_step - 1 ;
var n_level = <?php echo $nlev; ?> ;
var now_level = 0 ;
var levels = new Array() ;
<?php
  for ($lv = 0; $lv < $nlev; $lv++) {
       echo "levels[$lv] = '" . $levels[$lv] . "';\n";  
}
?>

var n_imgs = new Array() ;
<?php
  for ($s = 0; $s < $n_set; $s++) {
       echo "n_imgs[$s] = '" . $n_imgs[$s] . "';\n";  
}
?>

var now_index = new Array() ;
now_index[0] = 0 ;
now_index[1] = 0 ;
now_index[2] = 0 ;

var set_name = new Array();

<?php
  for ($s = 0; $s < $n_set; $s++) {
       echo "set_name[$s] = '" . $set_name[$s] . "';\n";  
}
?>


var imgsrc = new Array();

var imgnodata = new Array();
imgnodata[0]="data/d4/nodata.png";
imgnodata[1]="data/d4/nodata.png";
imgnodata[2]="data/d4/nodata.png";

var vtime = new Array();

var XlocMarker = new Array();
var YlocMarker = new Array();
var InitTimes = new Array();


<?php
  for ($s = 0; $s < $n_set; $s++) {
   if ($n_imgs[$s] > 0) {
    echo "imgsrc[$s] = new Array();\n";
    echo "vtime[$s] = new Array();\n";
    for ($lv = 0; $lv < $nlev; $lv++) {
    echo "imgsrc[$s][$lv] = new Array();\n";
      foreach ($imgsrc[$s][$lv] as $n => $isrc) {
          echo "imgsrc[$s][$lv][$n] = '" . $isrc . "';\n";
      }
    }
     foreach ($vtime[$s] as $n => $isrc) {
       echo "vtime[$s][$n] = '" . $isrc . "';\n";
    }
}
}

    echo "InitTimes[0] = new Array();\n";
    echo "InitTimes[0][0] = '" . $InitTimes[0][0] . "';\n";

    echo "XlocMarker[0] = new Array();\n";
    echo "YlocMarker[0] = new Array();\n";
    echo "XlocMarker[0][0] = new Array();\n";
    echo "YlocMarker[0][0] = new Array();\n";
     foreach ($XlocMarker[0][0] as $n => $isrc) {
       echo "XlocMarker[0][0][$n] = '" . $isrc . "';\n";
    }
     foreach ($YlocMarker[0][0] as $n => $isrc) {
       echo "YlocMarker[0][0][$n] = '" . $isrc . "';\n";
    }
    echo "TimeMarker = new Array();\n";
  foreach ($TimeMarker as $n => $time) {
   echo "TimeMarker[$n] =".$time.";\n";
}

?>



function onload_script() {
  isload = true ;
  isover = false ;
  document.getElementById('img_area0').src = imgsrc[0][now_level][now_index[0]] ;
//  document.getElementById('img_area1').src = imgsrc[1][now_index[1]] ;
  document.getElementById('img_area2').src = imgsrc[1][now_level][now_index[2]] ;
  marker_blue();


<?php
if ($n_imgs[$set] > 0)
  echo "  change_image() ;";
?>


<?php
if (isset($_GET['anim']) && ($_GET['anim'] == 'true')) {
  echo "  animation_image() ;\n";
}
?>
}


function change_image() {
  if (isload) {
  var srcs = new Array();
    for ( s = 0 ; s < 3 ; s++) {   
  i_img=9999;  
  intv=100000;
  for (n=0;n<n_imgs[s];n++){

//   if ( TimeMarker[now_step] == vtime[s][n] ) {
   if (( Number(TimeMarker[now_step])-Number(vtime[s][n]) <= intv ) && (Number(TimeMarker[now_step]) >=Number(vtime[s][n]))) {
   i_img = n ;
   intv = Number(TimeMarker[now_step]) - Number(vtime[s][n]) ;
  }

  }
  if (i_img < n_imgs[s] ) {
  srcs[s] = imgsrc[s][now_level][i_img] ;
 }else{
  srcs[s] = imgnodata[s];
}
}

  document.getElementById('img_area0').src = srcs[0] ;
//  document.getElementById('img_area1').src = srcs[1] ;
  document.getElementById('img_area2').src = srcs[1] ;
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

function play() {
  if (is_animation == true) {
    now_step++ ;
    if (now_step === n_step ) now_step = 0 ;
    change_image() ;
    setTimeout('play();', interval) ;
  }
}

function marker_blue() {
      xmargin = 5 ;
      document.getElementById('marker_blue').style.left = Number(XlocMarker[0][0][now_step]) -xmargin + 'px';
      document.getElementById('marker_blue').style.top =  YlocMarker[0][0][now_step] + 'px';
}


function changeMapImage(ele,imgPath) {
  var myid = ele.id ;
  document.getElementById(myid).src = imgPath;
}

function change_Image_byMarker(ele) { 
  var srcs = new Array();
  myid = ele.id ;
  nums = myid.split("_") ;
  now_step = Number(nums[1]) ;    

  change_image();
//  
//    for ( s = 0 ; s < 3 ; s++) {   
//
//  i_img=9999;  
//  intv=100000;
//  for (n=0;n<n_imgs[s];n++){

//   if ( TimeMarker[now_step] == vtime[s][n] ) {
//   if (( Number(TimeMarker[now_step])-Number(vtime[s][n]) <= intv ) && (Number(TimeMarker[now_step]) >=Number(vtime[s][n]))) {
//   i_img = n ;
//   intv = Number(TimeMarker[now_step]) - Number(vtime[s][n]) ;
//  }
//
//  }
//  if (i_img < n_imgs[s] ) {
//  srcs[s] = imgsrc[s][now_level][i_img] ;
// }else{
//  srcs[s] = imgnodata[s];
//}
//}

//  document.getElementById('img_area0').src = srcs[0] ;
////  document.getElementById('img_area1').src = srcs[1] ;
//  document.getElementById('img_area2').src = srcs[2] ;
//  marker_blue(); 
//
}

function change_level(ele) {
  var myid = ele.id ;
  nums = myid.split("_");
  now_level = nums[2];
  change_image();
 <?php 
  for ($i=0;$i<$nlev;$i++){
   echo "document.getElementById('button_level_$i').style.backgroundColor = '#FCFCFC';";
  }
?>
  ele.style.backgroundColor = "#A4A4A4";
}

// -->
</script>
</head>

<body onLoad="onload_script();">

<div align="left">
  <p class="p_title">Model forecasts quicklook for domain 4 : realtime </p>
<?php
 if (intval($time_offset) < 0){
 echo  "<p class=\"p_title\"><span style=\"color:red\">past mode : ".$nowdate." </p>\n";
 }
 if (intval($time_offset) > 0){
 echo  "<p class=\"p_title\"><span style=\"color:red\">future mode : offset +".$time_offset." s </p>\n";
 } 
?>
</div>


<div align="left">
<table border="0" cellpadding="4" cellspacing="0"><tr>

  <td valign="top" style="text-align: left" width="100px">
<?php
if ($n_step == 0)
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" disabled>\n";
else
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" onClick=\"animation_image();\">\n";
?>

  </td>
<td valign="top" style="text-align: right" width="1300px">
<p class="p_menu_title" align="left">
Level
<?php
for ($i=0;$i < $nlev;$i++){
if ($i==0){
echo "<input id=\"button_level_$i\" type=\"button\" value=\"$levels_f[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #A4A4A4;\" onclick=\"change_level(this)\">";
}else{
echo "<input id=\"button_level_$i\" type=\"button\" value=\"$levels_f[$i]\" style=\"font-size:16px; width:100px; margin:10px; padding:0px; background-color: #FCFCFC;\" onclick=\"change_level(this)\">";
}
}
?>
</p>
</td>
</tr></table>

<table border="0" cellpadding="4" cellspacing="0"><tr>
  <td valign="top" style="text-align: center; padding-left: 0px"
  width="600px"><p id="title_set1" class="p_model" style="border-style: ridge;
  border-color:red; border-width:2px ;">PAWR-SCALE-LETKF</p></td>
  <td valign="top" style="text-align: center; padding-left: 0px"
  width="600px"><p id="title_set2" class="p_model" style="border-style: ridge;
  border-color:blue; border-width:2px ;">PAWR OBS</p></td>
</tr>
<tr>
  <td valign="top" style="text-align: right; padding-left: 0px" width="600px">
<?php
$filepath = $imgsrc[0][0][$n_imgs[0]-1];
if (($n_imgs[0] > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area0\" src=\"$filepath\" border=\"0\" width=\"700px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area0\" border=\"0\"  width=\"700px\">\n";
?>
  </td>

  <td valign="top" style="text-align: left; padding-left: 0px">
<?php
$filepath = $imgsrc[2][0][$n_imgs[2]-1];
if (($n_imgs[2] > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area2\" src=\"$filepath\"
  border=\"0\" width=\"700px\">\n";
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

<p style="width:1200px; overflow:auto; position:relative; top:0px; left:0px;" >
  <img id="monitor_campus" src="./monitor/monitor_plot_d4.png" style="margin:0px;"> 
 <?php
  echo "<img id=\"marker_blue\" src=\"./img/marker-d2_blue.png\" style=\"position:absolute; top:".$YlocMarker[0][0][$n_step-1]."px; left:".($XlocMarker[0][0][$n_step-1]-$xmargin)."px; margin:0px; height:200px; width:".$xsize."px\";>\n";

 for ($n = 0; $n < $n_step; $n++) {
  echo  " <img id=\"markerposition_$n\" src=\"./img/marker_trans_new.png\" ;
 style=\"position:absolute;  top: ".$YlocMarker[0][0][$n]."px;
 left:".($XlocMarker[0][0][$n]-$xmargin)."px; margin:0px; height:200px; width:".$xsize."px;\"
 onclick=\"change_Image_byMarker(this);\" onmouseover=\"changeMapImage(this,'./img/marker-d2_grey.png');\" onmouseout=\"changeMapImage(this,'./img/marker_trans_new.png');\" >  \n"; 
  };
 ?>
</p>

</div>


</body>
</html>
