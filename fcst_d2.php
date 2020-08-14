<?php

/* JMA precip turned off since 2020/04/01 */

date_default_timezone_set('UTC');
$max_rtime = 7;

$n_set = 3;
### $n_set = 2;  

$sets[0] = 'd2';
$sets[1] = 'msm';
#$sets[2] = 'JMA_precip/anal_d2/realtime';
$sets[2] = 'JMA_precip/nowcast_d2/realtime';

$set_name[0] = 'SCALE-LETKF 6km mdet';
$set_name[1] = 'JMA MSM';
$set_name[2] = 'JMA precip radar';

$n_prod = 9;
$prods[0] = 'sfc_prcp';
$prods[1] = 'sfc_wind';
$prods[2] = 'sfc_2mtemp';
$prods[3] = '925_temp';
$prods[4] = '925_theq';
$prods[5] = '850_temp';
$prods[6] = '850_theq';
$prods[7] = '700_temp';
$prods[8] = '500_temp';
#$prods[10] = '300_wspd';
#$prods[11] = 'olr_ir';
#*$prods[12] = 'max_ref';
#$prods[3] = 'sfc_temp';
#$prods[4] = 'sfc_aprcp';
#$prods[5] = 'sfc_asnow';
#$prods[12] = 'stn_63518_temp';
#$prods[6] = 'prcp_cwb';
#$prods[7] = 'prcp_cwb_1h';
#$prods[8] = 'radar_cwb';
#$prods[9] = 'sfc_temp_cwb';
#$prods[10] = 'sfc_wind_tw';
$prods_f[0] = 'Surface [SLP+Prcp]';
$prods_f[1] = 'Surface [Wind]';
$prods_f[2] = 'Surface [T2m]';
$prods_f[3] = '925hPa [T+RH+Wind]';
$prods_f[4] = '925hPa [Theta_eq]';
$prods_f[5] = '850hPa [T+RH+Wind]';
$prods_f[6] = '850hPa [Theta_eq]';
$prods_f[7] = '700hPa [T+RH+Wind]';
$prods_f[8] = '500hPa [T+RH+Wind]';
#$prods_f[10] = '300hPa [Hgt+Wspd]';
#$prods_f[11] = 'Cloud-IR [OLR]';
#$prods_f[12] = 'Radar [MaxRef]';
#$prods_f[12] = 'Station-Kobe [T]';
#$prods_f[6] = '(TW-CWB) Rainfall';
#$prods_f[7] = '(TW-CWB) Rainfall-1H';
#$prods_f[8] = '(TW-CWB) Radar';
#$prods_f[9] = '(TW-CWB) Surface-T';
#$prods_f[10] = '(TW) Surface-Wind';
#$prods_f[3] = 'Surface [Tskin]';
#$prods_f[4] = 'Surface [PrcpTot]';
#$prods_f[5] = 'Surface [SnowTot]';

$datadir = 'data';
$cmem='mdet' ;

###

 $json_xloc = file_get_contents('json/xlocs_d2.json' );
 $json_yloc = file_get_contents('json/ylocs_d2.json' );
 $json_itime = file_get_contents('json/itimes_d2.json' );
 $XlocMarker = json_decode($json_xloc,true);
 $YlocMarker = json_decode($json_yloc,true);
 $InitTimes = json_decode($json_itime,true);

 $xmargin = 5 ;
 $xsize = 10 ;


$set = 0;
if (isset($_GET['set']) && ($_GET['set'] <= $n_set))
  $set = $_GET['set'];
$set2 = 1;
if (isset($_GET['set2']) && ($_GET['set2'] <= $n_set))
  $set2 = $_GET['set2'];

#--
 $time_offset=file_get_contents('time_offset.txt');
 $nowtime = date("U") + intval($time_offset)+ 9*3600;
 $nowdate = date("Y n/j",$nowtime);
 $nowY = date("Y",$nowtime);

$prod_index = 0;
$prod = $prods[$prod_index];
if (isset($_GET['prod'])) {
  for ($n = 0; $n < $n_prod; $n++) {
    if ($_GET['prod'] == $prods[$n]) {
      $prod = $_GET['prod'];
      $prod_index = $n;
    }
  }
}

$prod2_index = 0;
$prod2 = $prods[$prod2_index];
if (isset($_GET['prod2'])) {
  for ($n = 0; $n < $n_prod; $n++) {
    if ($_GET['prod2'] == $prods[$n]) {
      $prod2 = $_GET['prod2'];
      $prod2_index = $n;
    }
  }
}
#--

for ($s = 0; $s < $n_set ; $s++) {
 for ($r_index = 0; $r_index < $max_rtime ; $r_index++) {
 if (isset($InitTimes[$s][$r_index])) {
  $rtimes[$s][$r_index] = $nowY.strval($InitTimes[$s][$r_index])."0000";
  $n_rtimes[$s] = $r_index + 1;
}
}
}
###

for ($s = 0; $s < $n_set-1; $s++) {
for ($r_index = 0; $r_index < $n_rtimes[$s]; $r_index++) {

if (isset($InitTimes[$s][$r_index])) {


for ($p = 0; $p < $n_prod; $p++) {
$prod_l = $prods[$p] ;

$path = "$datadir/" . $sets[$s] . "/".$rtimes[$s][$r_index]."/$prod_l";
if  ($s == 0){
$path = "$datadir/" . $sets[$s] . "/".$rtimes[$s][$r_index]."/".$cmem."/$prod_l";
}
#echo $path;
#exit;

$imgs = scandir($path);
$n_img = 0;
foreach ($imgs as $img_file) {
  if (ereg ("${prod_l}.gif", $img_file, $regs) || ereg ("${prod_l}.png", $img_file, $regs)) {
    $files[$n_img] = $img_file;
    $n_img++ ;
  }
  else if (ereg ("${prod_l}_f([0-9]{6}).gif", $img_file, $regs) || ereg ("${prod_l}_f([0-9]{6}).png", $img_file, $regs)) {
    $files[$n_img] = $img_file;
    $n_img++ ;
  }
}

for ($n = 0; $n < $n_img; $n++){
  $imgsrc[$s][$r_index][$n][$p] = "$path/" . $files[$n];
}

}

}
}
}



$s = $n_set-1 ;
$p = 0;
$path = "$datadir/" . "$sets[$s]";
$imgs = scandir($path);
$n = 0;
foreach ($imgs as $img_file) {
#  if (ereg ("anal_([0-9]{14}).png", $img_file, $regs) || ereg  ("fcst_([0-9]{14}).png", $img_file, $regs)) {
  if (ereg ("radar_([0-9]{10})0000.png", $img_file, $regs) || ereg  ("nowcast_([0-9]{10})0000.png", $img_file, $regs)) {
    if (intval(substr($regs[1],4,6)) >= intval($InitTimes[$s][0]) ) {
    $imgsrc[$s][0][$n][$p] = "$path/" .$img_file;
    $n++ ;
  }
  }
}




$n_imgs[$s] = $n ;

$n_imgs[0] = 25 ; 
$n_imgs[1] = 34 ;

?>

<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Refresh" content="60">
<title>Model forecasts quicklook for domain 2</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
<!--
var interval = 400 ;
var isload = false ;
var is_animation = false ;

<?php
 if ( $n_rtimes[$set] > 0 ) { 
echo "var now_set = 0 ;\n";
echo "var now_nr = ". ($n_rtimes[$set]-1) ." ;\n" ;
 }else{
echo "var now_set = 1 ;\n";
echo "var now_nr = ". ($n_rtimes[$set2]-2) . " ;\n";
}
?>
var now_index = 0 ;
var now_prod = 0 ;

var now_setb = 1 ;
var now_nrb = <?php echo $n_rtimes[$set2]-1; ?> ;
var now_indexb = 0 ;
var now_prodb = 0 ;

var set_name = new Array();

<?php
  for ($s = 0; $s < $n_set; $s++) {
       echo "set_name[$s] = '" . $set_name[$s] . "';\n";  
}
?>

var imgsrc = new Array();

var XlocMarker = new Array();
var YlocMarker = new Array();
var InitTimes = new Array();


<?php
  for ($s = 0; $s < $n_set; $s++) {
   if ($n_imgs[$s] > 0) {
    echo "imgsrc[$s] = new Array();\n";
    for ($r_index=0; $r_index < $n_rtimes[$s]; $r_index++) {
      echo "imgsrc[$s][$r_index] = new Array();\n";
     for ($n_index=0; $n_index < $n_imgs[$s]; $n_index++) {
       echo "imgsrc[$s][$r_index][$n_index] = new Array();\n";
     foreach ($imgsrc[$s][$r_index][$n_index] as $n => $isrc) {
       echo "imgsrc[$s][$r_index][$n_index][$n] = '" . $isrc . "';\n";
    }

   }
  }
}
}

  for ($s = 0; $s < $n_set; $s++) {
    echo "XlocMarker[$s] = new Array();\n";
    echo "YlocMarker[$s] = new Array();\n";
    echo "InitTimes[$s] = new Array();\n";
       foreach ($InitTimes[$s] as $n => $isrc) {
      echo "InitTimes[$s][$n] = '" . $isrc . "';\n";
}
    for ($r_index=0; $r_index < $n_rtimes[$s]; $r_index++) {
      echo "XlocMarker[$s][$r_index] = new Array();\n";
      echo "YlocMarker[$s][$r_index] = new Array();\n";
     foreach ($XlocMarker[$s][$r_index] as $n => $isrc) {
       echo "XlocMarker[$s][$r_index][$n] = '" . $isrc . "';\n";
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
  document.getElementById('img_area').src = imgsrc[now_set][now_nr][now_index][now_prod] ;
  document.getElementById('img_area2').src = imgsrc[now_setb][now_nrb][now_indexb][now_prodb] ;
  marker_red();
  marker_blue();


<?php
if ($n_imgs[$set] > 0)
  echo "  change_image() ;";
?>


<?php
echo "  set = ${set};\n";
if (isset($_GET['anim']) && ($_GET['anim'] == 'true') && ($n_imgs[$set] > 0)) {
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
    document.getElementById('img_area').src = imgsrc[now_set][now_nr][now_index][now_prod] ;
    document.getElementById('img_area2').src = imgsrc[now_setb][now_nrb][now_indexb][now_prodb] ;
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

function marker_red() {
      xmargin = 5 ;
      document.getElementById('marker_red').style.left = Number(XlocMarker[now_set][now_nr][now_index])-xmargin + 'px';
      document.getElementById('marker_red').style.top =  YlocMarker[now_set][now_nr][now_index] + 'px';
}

function marker_blue() {
      xmargin = 5 ;
      if ( now_set == now_setb && now_nr && now_nrb && now_index == now_indexb      ) { xmargin = 4 }
      document.getElementById('marker_blue').style.left =
      Number(XlocMarker[now_setb][now_nrb][now_indexb]) -xmargin + 'px';
      document.getElementById('marker_blue').style.top =  YlocMarker[now_setb][now_nrb][now_indexb] + 'px';
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
  document.getElementById('img_area').src = imgsrc[now_set][now_nr][now_index][now_prod] ;
  marker_red(); 
}
function change_Image2_byMarker(ele) { 
  myid = ele.id ;
  nums = myid.split("_") ;
  now_setb = nums[1] ;    
  now_nrb = nums[2] ;    
  now_indexb = nums[3] ;    
  document.getElementById('title_set2').textContent = set_name[now_setb] ;
  document.getElementById('img_area2').src = imgsrc[now_setb][now_nrb][now_indexb][now_prodb] ;
  marker_blue(); 
}
// -->
</script>
</head>

<body onLoad="onload_script();">

<div align="left">
  <p class="p_title">Model forecasts quicklook for domain 2</p>
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

  <td valign="top" style="text-align: left" width="600px">
    <form name="view" action="#">
      <p class="p_menu_title">Products:</p>
      <select name="prod" id='prod' onChange="change_prod();" class="menu" style="background-color: #F8E0EC">
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
  <td valign="bottom" style="text-align: center" width="200px">
<?php
if ($n_imgs[$set] == 0)
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" disabled>\n";
else
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" onClick=\"animation_image();\">\n";
?>
    </form>
  </td>
<td valign="top" style="text-align: right" width="600px">
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
  border-color:red; border-width:2px ;">SCALE-LETKF d2 mdet</p></td>
  <td valign="top" style="text-align: center; padding-left: 0px"
  width="600px"><p id="title_set2" class="p_model" style="border-style: ridge;
  border-color:blue; border-width:2px ;">SCALE-LETKF d2 mdet</p></td>
</tr>
<tr>
  <td valign="top" style="text-align: right; padding-left: 0px" width="600px">
<?php
$filepath = $imgsrc[0][$n_rtime[0]][0][0];
if (($n_imgs[0] > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "    <img class=\"img_show\" id=\"img_area\" src=\"$filepath\" border=\"0\" width=\"700px\">\n";
}
else
  echo "    <img class=\"img_show\" id=\"img_area\" border=\"0\" width=\"700px\">\n";
?>
  </td>

  <td valign="top" style="text-align: left; padding-left: 0px">
<?php
$filepath = $imgsrc2[1][$n_rtime[1]][0][0];
if (($n_imgs[1] > 0) && file_exists($filepath)) {
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
  <img id="monitor_campus" src="./monitor/monitor_plot_d2.png" style="margin:0px;"> 

 <?php
  echo "<img id=\"marker_red\" src=\"./img/marker-d2.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:15px; width:".$xsize."px\";>\n";
  echo "<img id=\"marker_blue\" src=\"./img/marker-d2_blue.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:15px; width:".$xsize."px\";>\n";

 for ($s = 0; $s < $n_set; $s++) {
 for ($nr = 0; $nr < $n_rtimes[$s]; $nr++) {
 for ($n = 0; $n < $n_imgs[$s]; $n++) {
  if (isset($XlocMarker[$s][$nr][$n]) ) {  
  echo  " <img id=\"markerposition_".$s."_".$nr."_".$n.   "\" src=\"./img/marker_trans_new.png\" ;
 style=\"position:absolute;  top: ".$YlocMarker[$s][$nr][$n]."px;
 left:".($XlocMarker[$s][$nr][$n]-$xmargin)."px; margin:0px; height:15px; width:".$xsize."px;\"
 onclick=\"change_Image_byMarker(this);\" oncontextmenu=\"change_Image2_byMarker(this); return false\" onmouseover=\"changeMapImage(this,'./img/marker-d2_grey.png');\" onmouseout=\"changeMapImage(this,'./img/marker_trans_new.png');\" >  \n"; 
  };  
  };
  };
  };
 ?>

</div>


</body>
</html>
