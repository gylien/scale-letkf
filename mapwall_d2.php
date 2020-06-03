<?php

date_default_timezone_set('UTC');
$max_rtime = 20;
$refresh_int = 10;

$n_prod = 3;
$prods[0] = 'sfc_prcp';
$prods[1] = 'sfc_wind';
$prods[2] = 'sfc_2mtemp';
#$prods[3] = '850_theq';
#$prods[3] = 'sfc_temp';
#$prods[1] = 'max_ref';
#$prods[2] = 'sfc_2mtemp';
#$prods[3] = 'sfc_temp';
#$prods[4] = '850_temp';
#$prods[5] = '700_vvel';
#$prods[7] = '300_wspd';
#$prods[8] = 'olr_ir';


$prods_f[0] = 'Surface [SLP+Rain]';
$prods_f[1] = 'Surface [Wind]';
$prods_f[2] = 'Surface [T2m]';
#$prods_f[3] = '850hPa [theq+Wind] 700hPa [RH]';
#$prods_f[2] = '500hPa [Hgt+Vort]';
#$prods_f[3] = '850hPa [T+RH+Wind]';
#$prods_f[1] = 'Radar [MaxRef]';
#$prods_f[2] = 'Surface [T2m]';
#$prods_f[3] = 'Surface [Tskin]';
#$prods_f[4] = '850hPa [T+RH+Wind]';
#$prods_f[5] = '700hPa [Hgt+Vvel]';
#$prods_f[7] = '300hPa [Hgt+Wspd]';
#$prods_f[8] = 'Cloud-IR [OLR]';


$datadir = 'data';

$model = 'r0051_nest_d2';
$modelv = 'gfs';

#$prod_index = 0;
#$prod = $prods[$prod_index];
#if (isset($_GET['prod'])) {
#  for ($n = 0; $n < $n_prod; $n++) {
#    if ($_GET['prod'] == $prods[$n]) {
#      $prod = $_GET['prod'];
#      $prod_index = $n;
#    }
#  }
#}

$dts = scandir("$datadir/$model", 1);
$n_rtime = 0;
foreach ($dts as $dt_dir) {
  if (ereg ("([0-9]{10})", $dt_dir)) {
    $rtimes[$n_rtime] = $dt_dir;
    $ryr = substr($dt_dir, 0, 4);
    $rmn = substr($dt_dir, 4, 2);
    $rdy = substr($dt_dir, 6, 2);
    $rhr = substr($dt_dir, 8, 2);
    $rtimes_f[$n_rtime] = date("Y-m-d_H", mktime($rhr, 0, 0, $rmn, $rdy, $ryr)) . 'Z (' .
                          date("d_H", mktime($rhr+9, 0, 0, $rmn, $rdy, $ryr)) . 'L)' ;
    $n_rtime++;
    if ($n_rtime >= $max_rtime) break;
  }
}
$rtimes = array_reverse($rtimes);
$rtimes_f = array_reverse($rtimes_f);

$r_index = -1;
if (isset($_GET['rtime'])) {
  for ($n = 0; $n < $n_rtime; $n++) {
    if ($_GET['rtime'] == $rtimes[$n]) {
      $rtime = $_GET['rtime'];
      $r_index = $n;
    }
  }
}
if ($r_index == -1) {
  $r_index = $n_rtime - 1;
  $rtime = $rtimes[$r_index];
}

$ryr = substr($rtime, 0, 4);
$rmn = substr($rtime, 4, 2);
$rdy = substr($rtime, 6, 2);
$rhr = substr($rtime, 8, 2);
#if ($rhr <= 16)
#  $rtime_verify_start = date("YmdH", mktime(17, 0, 0, $rmn, $rdy, $ryr));
#else
#  $rtime_verify_start = date("YmdH", mktime(17, 0, 0, $rmn, $rdy+1, $ryr));

$path = "$datadir/$model/$rtime";
$imgs = scandir("$path/" . $prods[0]);
$n_img = 0;
$img_index = 0;
foreach ($imgs as $img_file) {
  if (ereg ($prods[0] . "_f([0-9]{6}).gif", $img_file, $regs) || ereg ($prods[0] . "_f([0-9]{6}).png", $img_file, $regs)) {
    $files[0][$n_img] = $prods[0] . "/" . $img_file;

    ######
    $files[1][$n_img] = $prods[1] . "/" . $prods[1] . "_f" . $regs[1] . ".png";
    $files[2][$n_img] = $prods[2] . "/" . $prods[2] . "_f" . $regs[1] . ".png";
    $files[3][$n_img] = $prods[3] . "/" . $prods[3] . "_f" . $regs[1] . ".png";
    ######

    $vtimes[$n_img] = date("YmdHis", mktime($rhr, 0, 0+$regs[1], $rmn, $rdy, $ryr));
    $vtimes_f[$n_img] = date("Y-m-d_H", mktime($rhr, 0, 0+$regs[1], $rmn, $rdy, $ryr)) . 'Z (' .
                        date("d_H", mktime($rhr+9, 0, 0+$regs[1], $rmn, $rdy, $ryr)) . 'L)' ;
    if (isset($_GET['vtime']) && ($_GET['vtime'] == $vtimes[$n_img]))
      $img_index = $n_img;
    $n_img++ ;
  }
}

######

#$verf_num = 0;
#for ($n = 0; $n < $n_img; $n++)
#  $verf[0][$n] = "$path/" . $files[$n];
#$verf_name[0] = 'SCALE';

#if ($prod_index <= 7) {
#  $verf_num = 2;
#  for ($n = 0; $n < $n_img; $n++) {
#    $verf[1][$n] = "$datadir/$modelv/$rtime/" . $files[$n];
#    $verf[2][$n] = "$datadir/$modelv/" . $vtimes[$n] . "/${prod}_f000000.png";
#  }
#  $verf_name[1] = 'GFS (F)';
#  $verf_name[2] = 'GFS (A)';
#}

#$verify = 0;
#if (isset($_GET['verify']) && ($_GET['verify'] <= $verf_num))
#  $verify = $_GET['verify'];


######

?>
<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Refresh" content="600">
<title>SCALE-LETKF experimental real-time forecasts domain 2</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
<!--
var inteval = 400 ;
var now_index = 0 ;
var replay_num = 0 ;
var isload = false ;
var is_animation = true ;
//var verify = 0 ;

var imgsrc = new Array();
<?php
if ($n_img > 0) {
  for ($p = 0; $p < $n_prod; $p++) {
    echo "imgsrc[$p] = new Array();\n";
    foreach ($files[$p] as $n => $ifile) {
      echo "imgsrc[$p][$n] = '" . $ifile . "';\n";
    }
  }
}
?>

function onload_script() {
  isload = true ;
  now_index = <?php echo $img_index; ?> ;
<?php
if ($n_img > 0)
  echo "  change_image() ;";
?>

  if (document.images) {
    var imgo = new Array();
<?php
if ($n_img > 0) {
  for ($p = 0; $p < $n_prod; $p++) {
echo "    imgo[$p] = new Array();\n" ;
    foreach ($files[$p] as $n => $ifile) {
      echo "    imgo[$p][$n] = new Image();\n";
      echo "    imgo[$p][$n].src = imgsrc[$p][$n];\n";
    }
  }
}
?>
  }

<?php
#echo "  verify = ${verify};\n";
#if (isset($_GET['anim']) && ($_GET['anim'] == 'true') && ($n_img > 0)) {
  echo "  animation_image() ;\n";
#}
?>
}

//function change_prod() {
//  if (isload) {
//    var rtime_obj = document.getElementById('rtime') ;
//    var vtime_obj = document.getElementById('vtime') ;
//    var prod_obj = document.getElementById('prod') ;
//    location.href = 'index.php?prod=' + prod_obj.value + '&rtime=' + rtime_obj.value + '&vtime=' + vtime_obj.value + '&verify=' + verify + '&anim=' + is_animation ;
//  }
//}

function change_rtime() {
  if (isload) {
    var rtime_obj = document.getElementById('rtime') ;
    var vtime_obj = document.getElementById('vtime') ;
    var prod_obj = document.getElementById('prod') ;
//    location.href = 'index.php?prod=' + prod_obj.value + '&rtime=' + rtime_obj.value + '&vtime=' + vtime_obj.value + '&verify=' + verify + '&anim=' + is_animation ;
    location.href = 'mapwall_d2.php?rtime=' + rtime_obj.value + '&vtime=' + vtime_obj.value ;
  }
}

//function change_verify(iv) {
//  if (isload) {
//    var rtime_obj = document.getElementById('rtime') ;
//    var vtime_obj = document.getElementById('vtime') ;
//    var prod_obj = document.getElementById('prod') ;
//    verify = iv
//    location.href = 'index.php?prod=' + prod_obj.value + '&rtime=' + rtime_obj.value + '&vtime=' + vtime_obj.value + '&verify=' + verify + '&anim=' + is_animation ;
//  }
//}

//function change_verify_over(iv) {
//  if (isload) {
//    document.getElementById('img_area').src = imgsrc[iv][now_index] ;
//    if (iv != <?php echo $verify; ?>) {
//      document.getElementById('verify' + String(iv)).style.backgroundColor = '#FFFF33';
//      document.getElementById('verify' + String(iv)).style.borderStyle = 'outset';
//      document.getElementById('verify<?php echo $verify; ?>').style.backgroundColor = '#FFFFFF';
//      document.getElementById('verify<?php echo $verify; ?>').style.borderStyle = 'inset';
//    }
//  }
//}

//function change_verify_out(iv) {
//  if (isload) {
//    document.getElementById('img_area').src = imgsrc[<?php echo $verify; ?>][now_index] ;
//    if (iv != <?php echo $verify; ?>) {
//      document.getElementById('verify' + String(iv)).style.backgroundColor = '#FFFFFF';
//      document.getElementById('verify' + String(iv)).style.borderStyle = 'inset';
//      document.getElementById('verify<?php echo $verify; ?>').style.backgroundColor = '#FFFF33';
//      document.getElementById('verify<?php echo $verify; ?>').style.borderStyle = 'outset';
//    }
//  }
//}

//function change_verify_up(iv) {
//  change_verify(iv);
//}

//function change_verify_down(iv) {
//  if (isload) {
//    document.getElementById('verify' + String(iv)).style.borderStyle = 'inset';
//  }
//}

//function update_link() {
//  if (isload) {
//    var rtime_obj = document.getElementById('rtime') ;
//    var vtime_obj = document.getElementById('vtime') ;
//    var prod_obj = document.getElementById('prod') ;
//    location.href = 'index.php?prod=' + prod_obj.value + '&rtime=' + rtime_obj.value + '&vtime=' + vtime_obj.value + '&verify=' + verify + '&anim=' + is_animation ;
//  }
//}

function change_image() {
  if (isload) {
    var vtime_obj = document.getElementById('vtime') ;
    vtime_obj.options[now_index].selected = true ;
<?php
for ($p = 0; $p < $n_prod; $p++) {
  echo "    document.getElementById('img_area_${p}').src = '${path}/' + imgsrc[${p}][now_index] ;\n" ;
}
?>
  }
}

//function change_image_click() {
//  if (isload) {
//    if (is_animation == false) {
//      var vtime_obj = document.getElementById('vtime') ;
//      now_index = vtime_obj.selectedIndex ;
//      document.getElementById('img_area').src = imgsrc[<?php echo $verify; ?>][now_index] ;
//    }
//  }
//}

function animation_image() {
  if (isload) {
//    var ani_obj = document.getElementById('animation') ;
//    var vtime_obj = document.getElementById('vtime') ;
//    if (is_animation == false) {
//      is_animation = true ;
//      ani_obj.value = 'Stop' ;
//      vtime_obj.disabled = true ;
      play() ;
//    }
//    else {
//      is_animation = false ;
//      ani_obj.value = 'Animation' ;
//      vtime_obj.disabled = false ;
//    }
  }
}

function play() {
  if (is_animation == true) {
    now_index++ ;
    if (now_index >= <?php echo $n_img; ?>) {
      now_index = 0 ;
      replay_num++ ;
      if (replay_num >= <?php echo $refresh_int; ?>) {
        location.href = 'mapwall_d2.php' ;
      }
    }
    change_image() ;
    setTimeout('play();', inteval) ;
  }
}

// -->
</script>
</head>

<body onLoad="onload_script();">

<div align="left">
  <p class="p_title">SCALE-LETKF experimental real-time forecasts domain 2</p>
</div>

<div align="left">
<table border="0" cellpadding="4" cellspacing="0">
  <tr>

  <td valign="top" style="text-align: left" rowspan="2">
    <form name="view" action="#">
      <p class="p_menu_title">Base Time:</p>
      <select name="rtime" id='rtime' onChange="change_rtime();" class="menu" style="background-color: #FFE0C6">
<?php
for ($n = 0; $n < $n_rtime; $n++) {
  if ($n == $r_index)
    echo "        <option value=\"" . $rtimes[$n] . "\" selected=\"selected\">" . $rtimes_f[$n] . "</option>\n";
  else
    echo "        <option value=\"" . $rtimes[$n] . "\">" . $rtimes_f[$n] . "</option>\n";
}
?>
      </select>

      <p class="p_menu_title">Valid Time:</p>
      <select name="vtime" id='vtime' size="22" class="menu" disabled>
<?php
if ($n_img == 0)
  echo "        <option value=\"\">Not available</option>\n";
else {
  for ($n = 0; $n < $n_img; $n++) {
    echo "        <option value=\"" . $vtimes[$n] . "\">" . $vtimes_f[$n] . "</option>\n";
  }
}
?>
      </select>
    </form>
  </td>

    <td valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = "$path/" . $files[0][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "      <img class=\"img_show\" id=\"img_area_0\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "      <img class=\"img_show\" id=\"img_area_0\" border=\"0\">\n";
?>
    </td>
    <td valign="middle" style="text-align: left; padding-left: 0px">
<?php
$filepath = "$path/" . $files[1][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "      <img class=\"img_show\" id=\"img_area_1\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "      <img class=\"img_show\" id=\"img_area_1\" border=\"0\">\n";
?>
    </td>
  </tr>
  <tr>

    <td valign="middle" style="text-align: left; padding-left: 0px">

<?php
$filepath = "$path/" . $files[2][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "      <img class=\"img_show\" id=\"img_area_2\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "      <img class=\"img_show\" id=\"img_area_2\" border=\"0\">\n";
?>

    </td>
    <td valign="middle" style="text-align: left; padding-left: 0px">
<!--
<?php
$filepath = "$path/" . $files[3][$img_index];
if (($n_img > 0) && file_exists($filepath)) {
  $size = getimagesize($filepath);
  echo "      <img class=\"img_show\" id=\"img_area_3\" src=\"$filepath\" border=\"0\">\n";
}
else
  echo "      <img class=\"img_show\" id=\"img_area_3\" border=\"0\">\n";
?>
-->
    </td>

  </tr>
</table>
</div>

</body>
</html>
