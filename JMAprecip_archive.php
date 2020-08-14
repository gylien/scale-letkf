<?php

date_default_timezone_set('UTC');
$max_rtime = 500;

$n_area = 4;
$areas[0] = 'nowcast_japan/archive';
$areas[1] = 'nowcast_d2/archive';
$areas[2] = 'nowcast_d3/archive';
$areas[3] = 'nowcast_d4/archive';
$area_name[0] = 'Japan';
$area_name[1] = 'Domain 2';
$area_name[2] = 'Domain 3';
$area_name[3] = 'Domain 4';

$datadir = 'data/JMA_precip/';
$subdir_scale = 'scale_verify/'; 
$subdir_msm = 'msm_verify/'; 
$subdir_JMAprecip = 'JMA_precip/verify/'; 
$subdir_sounding = 'scale_verify_sounding/'; 

$times = array_slice(scandir($datadir.$areas[0]),2);
$n_time = count($times);

for ($itime=0; $itime < $n_time; $itime++){
  $ryr = substr($times[$itime], 0, 4);
  $rmn = substr($times[$itime], 4, 2);
  $rdy = substr($times[$itime], 6, 2);
  $times_f[$itime] = date("Y-m-d", mktime($rhr, 0, 0, $rmn, $rdy, $ryr)) ;

$vfiles[$itime] = array_slice(scandir($datadir.$areas[0].'/'.$times[$itime]),2);
$n_vtime[$itime] = count($vfiles[$itime]);
for($ivtime=0;$ivtime<$n_vtime[$itime];$ivtime++){
$vtime_png= explode('_',$vfiles[$itime][$ivtime]);
$vtimes_png[$itime][$ivtime]=substr($vtime_png[1],6,14);
$vtimes_f[$itime][$ivtime] = date("H:i:s", mktime($rhr, 0, 0, $rmn, $rdy, $ryr)) ;
}
}

$now_itime = $n_time - 1 ;
$now_ivtime = 0 ;
$now_area0 = 0;
$now_area1 = 1;

$now_stime = '2019060100';
?>

<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<title>SCALE-LETKF forecast domain 2 comparison: 2019 Summer</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">
<script language="javascript" type="text/javascript">
<!--
var ImgPathSCALE = "<?php echo $datadir.$subdir_scale; ?>";
var ImgPathMSM = "<?php echo $datadir.$subdir_msm; ?>";
var ImgPathJMAprecip = "<?php echo $datadir.$subdir_JMAprecip; ?>";
var ImgPathsounding = "<?php echo $datadir.$subdir_sounding; ?>";

var inteval = 400 ;
var now_index = 0 ;
var isload = false ;

var now_itime= <?php echo $now_itime; ?>;
var now_ivtime= <?php echo $now_ivtime; ?>;
var now_prod= <?php echo $now_prod; ?>;
var now_area0= <?php echo $now_area0; ?>;
var now_area1= <?php echo $now_area1; ?>;

var now_location= 0;
var now_stime= '2019060100';

var times = new Array();
<?php
 foreach ($times as $n => $time){
  echo "times[$n] = '".$time."';\n";
}
?>

var locations = new Array();
<?php
 foreach ($location as $n => $time){
  echo "locations[$n] = '".$time."';\n";
}
?>


var vtimes_png    = new Array();
var vtimes_f  = new Array();
var vtimes_f2 = new Array();
var stimes = new Array();
<?php
 for ($itime=0;$itime<$n_time; $itime++){
 echo "vtimes_png[$itime] = new Array();\n";
 echo "vtimes_f[$itime] = new Array();\n";
 echo "vtimes_f2[$itime] = new Array();\n";
 echo "stimes[$itime] = new Array();\n";
 foreach ($vtimes_png[$itime] as $n => $val){
  echo "vtimes_png[$itime][$n] = '".$val."';\n";
 }
 foreach ($vtimes_f[$itime] as $n => $val){
  echo "vtimes_f[$itime][$n] = '".$val."';\n";
 }
 foreach ($vtimes_f2[$itime] as $n => $val){
  echo "vtimes_f2[$itime][$n] = '".$val."';\n";
 }
 foreach ($stimes[$itime] as $n => $val){
  echo "stimes[$itime][$n] = '".$val."';\n";
 }
}
?>
var prods = new Array();
<?php
 foreach ($prods as $n => $prod){
  echo "prods[$n] = '".$prod."';\n";
}
?>

var areas = new Array();
<?php
 foreach ($areas as $n => $area){
  echo "areas[$n] = '".$area."';\n";
}
?>

var now_time =  times[now_itime];

function onload_script(){
 change_rtime();
 changeMap();
}
 
function changeMap(){
  now_vtime_f2 = vtimes_f2[now_itime][now_ivtime]
  now_stime = stimes[now_itime][now_ivtime]
  imgPath0 = ImgPathSCALE  + now_time  + '/' + areas[now_area0] + prods[now_prod] + '/' + areas[now_area0] + prods[now_prod] + '_' + vtimes_png[now_itime][now_ivtime];
  imgPath1 = ImgPathMSM  + areas[now_area0] + prods[now_prod] + '/'+ areas[now_area0] + prods[now_prod] + '_' + now_vtime_f2 + '.png';
  imgPath2 = ImgPathJMAprecip + 'd2/anal_'+ now_vtime_f2 + '0000.png';
  imgPath3 = ImgPathSCALE  + now_time  + '/' + areas[now_area1] + prods[now_prod] + '/' + areas[now_area1] + prods[now_prod] + '_' + vtimes_png[now_itime][now_ivtime];
  imgPath4 = ImgPathMSM  + areas[now_area1] + prods[now_prod] + '/'+ areas[now_area1] + prods[now_prod] + '_' + now_vtime_f2 + '.png';
  imgPath5 = ImgPathJMAprecip + areas[now_area1].slice(0,-1) + '/anal_'+ now_vtime_f2 + '0000.png';
  imgPath6 = ImgPathsounding + now_time.slice(0,-4)+ '/' + locations[now_location] + '_'+ now_stime + '.png';
  document.getElementById('img_area0').src = imgPath0;
  document.getElementById('img_area1').src = imgPath1;
  document.getElementById('img_area2').src = imgPath2;
  document.getElementById('img_area3').src = imgPath3;
  document.getElementById('img_area4').src = imgPath4;
  document.getElementById('img_area5').src = imgPath5;
  document.getElementById('img_area6').src = imgPath6;
}

function change_area1(ele) {
  var myid = ele.id ;
  nums = myid.split("_");
  now_area1 = nums[2];
  changeMap();
 <?php 
  for ($i=1;$i<$n_area;$i++){
   echo "document.getElementById('button_area_$i').style.backgroundColor = '#FCFCFC';";
  }
?>
  ele.style.backgroundColor = "#A4A4A4";
}

function change_location(ele) {
  var myid = ele.id ;
  nums = myid.split("_");
  now_location = nums[2];
  changeMap();
 <?php 
  for ($i=0;$i<$n_location;$i++){
   echo "document.getElementById('button_loc_$i').style.backgroundColor = '#FCFCFC';";
  }
?>
  ele.style.backgroundColor = "#A4A4A4";
}


  
function change_prod(ele) {
  var myid = ele.id ;
  nums = myid.split("_");
  now_prod = nums[2];
  changeMap();
 <?php 
  for ($i=0;$i<$n_prod;$i++){
   echo "document.getElementById('button_prod_$i').style.backgroundColor = '#FCFCFC';";
  }
?>
  ele.style.backgroundColor = "#A4A4A4";
}

function change_rtime(){
    var rtime_obj = document.getElementById('rtime') ;
    now_itime = rtime_obj.selectedIndex; 
    now_time = rtime_obj.value;
    var vtime_obj = document.getElementById('vtime') ;
    now_ivtime = vtime_obj.selectedIndex;

    time_options = document.getElementById('vtime').options;
 for (n=0; n < time_options.length ; n++){
    document.getElementById('vtime').options[n].text  = vtimes_f[now_itime][n];
    document.getElementById('vtime').options[n].value = vtimes_png[now_itime][n];
 }
 changeMap();
}

function change_vtime(){
    var vtime_obj = document.getElementById('vtime') ;
    now_ivtime = vtime_obj.selectedIndex;
    changeMap();
}
// -->
</script>
</head>

<body onLoad="onload_script();">

<div align="left">
  <p class="p_title">SCALE-LETKF forecast domain 2 comparison : 2019 Summer</p>
</div>

<div align="left">
<table border="0" cellpadding="4" cellspacing="0"><tr>
  <td valign="top" style="text-align: left">
    <form name="view" action="#">

      <p class="p_menu_title">Base Time:</p>
      <select name="rtime" id='rtime' size="16" onChange="change_rtime();" class="menu" style="background-color: #FFE0C6">
<?php
for ($itime = 0; $itime < $n_time; $itime++) {
  if ($itime == $now_itime)
    echo "        <option value=\"" . $times[$itime] . "\" selected=\"selected\">" . $times_f[$itime] . "</option>\n";
  else
    echo "        <option value=\"" . $times[$itime] . "\">" . $times_f[$itime] . "</option>\n";
}
?>
      </select>

      <p class="p_menu_title">Valid Time:</p>
      <select name="vtime" id='vtime' size="16" onChange="change_vtime();" class="menu" style="background-color: #C0FFF0">>
<?php
  for ($ivtime = 0; $ivtime < $n_vtime[$now_itime]; $ivtime++) {
   if ($ivtime == $now_ivtime)
    echo "        <option value=\"" . $vfiles[$now_itime][$ivtime] . "\" selected=\"selected\">" . $vtimes_f[$now_itime][$ivtime] . "</option>\n";
   else
    echo "        <option value=\"" . $vfiles[$now_itime][$ivtime] . "\">" . $vtimes_f[$now_itime][$ivtime] . "</option>\n";

  }
?>
      </select>
   </form>
<br>
<p class="p_menu_title" align="left"> JMA daily charts</p>
<a href="http://www.data.jma.go.jp/fcd/yoho/data/hibiten/2019/201906.pdf" target='_blank'>201906</a><br>
<a href="http://www.data.jma.go.jp/fcd/yoho/data/hibiten/2019/201907.pdf" target='_blank'>201907</a><br>
<a href="http://www.data.jma.go.jp/fcd/yoho/data/hibiten/2019/201908.pdf" target='_blank'>201908</a><br>
<a href="http://www.data.jma.go.jp/fcd/yoho/data/hibiten/2019/201909.pdf" target='_blank'>201909</a>
  </td>

  <td valign="middle" style="text-align: left; padding-left: 0px">
<table border="0" cellpadding="4" cellspacing="0">
<tr><td>
<p class="p_model" align="left">Product
<?php
for ($i=0;$i<count($prods);$i++){
if ($i==$now_prod){
echo "<input id=\"button_prod_$i\" type=\"button\" value=\"$prods_f[$i]\" style=\"font-size:16px; width:160px; margin:20px; padding:0px; background-color: #A4A4A4;\" onclick=\"change_prod(this)\">";
}else{
echo "<input id=\"button_prod_$i\" type=\"button\" value=\"$prods_f[$i]\" style=\"font-size:16px; width:160px; margin:20px; padding:0px; background-color: #FCFCFC;\" onclick=\"change_prod(this)\">";
}
}
?>
</p>
</td></tr>
<tr><td>
<p class="p_model" align="left">Domain 2
</p>
</td></tr>

<tr><td>
<table border="0" cellpadding="4" cellspacing="0" style="table-layout:fixed">
 <tr><td width="620px"><p id="title_set1" class="p_model" style="border-style: ridge; border-color:red; border-width:2px ;">SCALE-LETKF d2 mdet</p></td>
     <td><p id="title_set2" class="p_model" style="border-style: ridge; border-color:blue; border-width:2px ;">MSM Analysis</p></td>
     <td><p id="title_set3" class="p_model" style="border-style: ridge; border-color:green; border-width:2px ;">JMA precip analysis</p></td>
</tr>
 <tr> <td valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = $datadir.$subdir_scale.$times[$now_itime]."/".$areas[$now_area0].$prods[$now_prod]."/".$vfiles[$now_itime][$now_ivtime];
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area0\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area0\" border=\"0\">\n";
?>
  </td>

  <td valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = $datadir.$subdir_msm.$times[$now_itime]."/".$areas[$now_area0].$prods[$now_prod]."/".$prods[0]."_f000000.png";
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area1\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area1\" border=\"0\">\n";
?>
  </td>
<td width="600px">
<?php
$filepath = $datadir.$subdir_JMAprecip."anal_".$vtimes_f2[$now_itime][$now_ivtime]."0000.png";
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area2\" src=\"$filepath\" border=\"0\" width=\"600px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area2\" border=\"0\" width=\"600px\">\n";
?>
 
</td>
</tr>
</table>
</td></tr>
<tr><td>
<p class="p_model" align="left">Domain 3
<input id="button_area_1" type="button" value="Kansai" style="font-size:16px; width:140px; margin:10px; padding:0px; background-color: #A4A4A4;" onclick="change_area1(this)">
<input id="button_area_2" type="button" value="Kanto" style="font-size:16px; width:140px; margin:10px; padding:0px; background-color: #FCFCFC;" onclick="change_area1(this)">
</p>
</td></tr>

<tr><td>
<table border="0" cellpadding="4" cellspacing="0" style="table-layout:fixed">
 <tr> <td width="620px" valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = $datadir.$subdir_scale.$times[$now_itime]."/".$areas[$now_area1].$prods[$now_prod]."/".$vfiles[$now_itime][$now_ivtime];
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area3\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area3\" border=\"0\">\n";
?>
  </td>

  <td width="620px" valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = $datadir.$subdir_msm.$times[$now_itime]."/".$areas[$now_area1].$prods[$now_prod]."/".$prods[0]."_f000000.png";
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area4\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area4\" border=\"0\">\n";
?>
  </td>
<td width="620px">
<?php
$filepath = $datadir.$subdir_JMAprecip."anal_".$vtimes_f2[$now_itime][$now_ivtime]."0000.png";
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area5\" src=\"$filepath\" border=\"0\" width=\"600px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area5\" border=\"0\" width=\"600px\">\n";
?>
 
</td>
</tr>
</table>
</td></tr>
</table>
</td>

<td valign='top' width="400px">
<table border="0" cellpadding="4" cellspacing="0"  style="table-layout:fixed" width='400px'>
<tr><td valign='bottom' height="160px" ><p class="p_model" align="left">Sounding</p></td></tr>
<tr><td width='400px'>
<?php
for ($i=0;$i<count($location);$i++){
if ($i==$now_location){
echo "<input id=\"button_loc_$i\" type=\"button\" value=\"$location[$i]\" style=\"font-size:16px; width:140px; margin:20px; padding:0px; background-color: #A4A4A4;\" onclick=\"change_location(this)\">";
}else{
echo "<input id=\"button_loc_$i\" type=\"button\" value=\"$location[$i]\" style=\"font-size:16px; width:140px; margin:20px; padding:0px; background-color: #FCFCFC;\" onclick=\"change_location(this)\">";
}
}
?>
</td></tr>
<tr><td width="400px">
<?php
$filepath = $datadir.$subdir_sounding.$now_location."_".$now_stime.".png";
if (file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area6\" src=\"$filepath\" border=\"0\" width=\"400px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area6\" border=\"0\" width=\"400px\">\n";
?>
</td></tr>

</table>

</td>

</tr></table>
</div>



</body>
</html>
