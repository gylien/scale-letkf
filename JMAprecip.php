<?php

date_default_timezone_set('UTC');
$max_rtime = 24;

$n_set = 2;
$sets[0] = 'JMA_precip/anal_d3/realtime';
$sets[1] = 'JMA_precip/nowcast_d3/realtime';
$set_name[0] = 'JMA precip anal/fcst';
$set_name[1] = 'JMA precip radar/nowcast';

$n_prod = 1;
$prods[0] = 'rain';
$prods_f[0] = 'Surface [Prcp]';

$datadir = 'data';
$cmem='mean' ;
###

 $json_xloc = file_get_contents('json/xlocs_JMA.json' );
 $json_yloc = file_get_contents('json/ylocs_JMA.json' );
 $json_itime = file_get_contents('json/itimes_JMA.json' );
 $XlocMarker = json_decode($json_xloc,true);
 $YlocMarker = json_decode($json_yloc,true);
 $InitTimes = json_decode($json_itime,true);

 $xmargin = 5 ;
 $xsize = 10 ;

#--

$set = 0;
$n_rtimes[0]=count($InitTimes[0]);

###


### anal/fcst 
$s = 0 ;
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
$s = 1 ;
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
?>


<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Refresh" content="300">
<title>JMA precipitation products quicklook</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
<!--
var interval = 400 ;
var now_index = 0 ;
var now_indexb = 0;
var now_set  = 0;
var now_setb = 1;
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
?>


var XlocMarker = new Array();
var YlocMarker = new Array();
var InitTimes = new Array();

<?php

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
if ($n_img > 0)
  echo "  change_image() ;";
?>



<?php
echo "  set = ${set};\n";
if (isset($_GET['anim']) && ($_GET['anim'] == 'true') && ($n_img > 0)) {
  echo "  animation_image() ;\n";
}
?>
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

    if (XlocMarker[now_set][now_nr][now_index] < XlocMarker[now_setb][now_nrb][now_indexb] ) {
        now_index++ ;
} else if (XlocMarker[now_set][now_nr][now_index] > XlocMarker[now_setb][now_nrb][now_indexb] ) {
        now_indexb++ ;
} else {
    now_index++ ;
    now_indexb++ ;
}
    
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
  <p class="p_title">JMA precipitation products</p>
</div>

<div align="left">
<?php
if ($n_imgs[1] == 0)
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" disabled>\n";
else
  echo "        <input type=\"button\" name=\"animation\" id=\"animation\" value=\"Animation\" onClick=\"animation_image();\">\n";
?>
<br>

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
  <img id="monitor_campus" src="./monitor/monitor_plot_JMAprecip.png" style="margin:0px;"> 



 <?php
  echo "<img id=\"marker_red\" src=\"./img/marker-d2.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:50px; width:".$xsize."px\";>\n";
  echo "<img id=\"marker_blue\" src=\"./img/marker-d2_blue.png\" style=\"position:absolute; top:0px; left:0px; margin:0px; height:50px; width:".$xsize."px\";>\n";

 for ($s = 0; $s < $n_set; $s++) {
 for ($nr = 0; $nr < $n_rtimes[$s]; $nr++) {
 for ($n = 0; $n < $n_imgs[$s]; $n++) {
  if (isset($XlocMarker[$s][$nr][$n]) ) {  
  echo  " <img id=\"markerposition_".$s."_".$nr."_".$n.   "\" src=\"./img/marker_trans_new.png\" ; style=\"position:absolute;  top: ".$YlocMarker[$s][$nr][$n]."px; left:".($XlocMarker[$s][$nr][$n]-$xmargin)."px; margin:0px; height:50px; width:".$xsize."px;\" onclick=\"change_Image_byMarker(this);\" oncontextmenu=\"change_Image2_byMarker(this); return false\" onmouseover=\"changeMapImage(this,'./img/marker-d2_grey.png');\" onmouseout=\"changeMapImage(this,'./img/marker_trans_new.png');\" >  \n"; 
  };  
  };
  };
  };
 ?>

</div>

</body>
</html>
