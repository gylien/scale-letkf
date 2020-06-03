<?php
$cvar='u';
if (isset($_GET['var']) )  $cvar = $_GET['var'];
if ($cvar=='tv'){
$cvar_map='t' ;
}else{
$cvar_map=$cvar ; 
}
?>

<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Cache-Control" content="no-cache">
<title>LETKF monitor test (sounding)</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
//<!--

var ImgPathBase = "./monitor/letkf_monitor/sounding/" ;

var now_var = "<?php echo $cvar ;?>" ; 
var now_stat = "rmse" ; 
var now_prof = "001" ; 

function changeMapButton(ele) {
var myid=ele.id ;
nums=myid.split("_") ;
if (nums[1] == "var" ) { 
document.getElementById("button_var_"+ now_var).style.backgroundColor = '#FCFCFC' ;
now_var=nums[2] ;
} else if (nums[1] == "stat") {
document.getElementById("button_stat_"+ now_stat).style.backgroundColor = '#FCFCFC' ;
ele.style.backgroundColor = '#A4A4A4' ;
now_stat = nums[2] ;
}
if (now_var == "tv") {
now_var_map="t";
}else {
now_var_map=now_var;
}
MapImgPath=ImgPathBase + "../sub/map/" + now_var_map + "_anal_rmse_map_d1_sub_0001.png"  ;
document.getElementById("img_map").src = MapImgPath ;
document.getElementById(now_prof).style.backgroundColor = 'blue' ;
now_prof='001';
MapImgPath=ImgPathBase + "/profile_" + now_var + "_" + now_prof +"_0001.png"  ;
document.getElementById("img_graph").src = MapImgPath ;
document.getElementById(now_prof).style.backgroundColor = 'red' ;
location.href = '?var='+ now_var ;
}
function changeMapButtonDomain(ele) {
var myid=ele.id ;
document.getElementById(now_prof).style.backgroundColor = 'blue' ;
now_prof = myid;
MapImgPath=ImgPathBase + "/profile_" + now_var + "_" + now_prof +"_0001.png"  ;
document.getElementById("img_graph").src = MapImgPath ;
ele.style.backgroundColor = 'red' ;
}

function changeMapMouseOn(ele){
ele.style.backgroundColor = 'red' ;
}
function changeMapMouseOff(ele){
if (ele.id != now_prof){
ele.style.backgroundColor = 'blue' ;
}
}


//-->
</script>
</head>

<body>

<div align="left">
  <p class="p_title">LETKF monitor (sounding) </p>
</div>


<div align="left">
<table border="0" cellpadding="0" cellspacing="0">

  <tr>
  <td valign="top" style="text-align: left">
  </td>
   <td valign="top" style="text-align: left">
<table border="0" cellpadding="0" cellspacing="0">
<tr>
<td><p>Variable:</p></td>
<td>
<input id="button_var_u" type="button" value="U" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color:<?php if($cvar=='u'){echo '#A4A4A4';} else {echo'#FCFCFC';}?>;"  onclick="changeMapButton(this)">
<input id="button_var_v" type="button" value="V" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color:<?php if($cvar=='v'){echo '#A4A4A4';} else {echo'#FCFCFC';}?>;"  onclick="changeMapButton(this)">
<input id="button_var_tv" type="button" value="Tv" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color:<?php if($cvar=='tv'){echo '#A4A4A4';} else {echo'#FCFCFC';}?> ;"  onclick="changeMapButton(this)">
<input id="button_var_q" type="button" value="Q" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color:<?php if($cvar=='q'){echo '#A4A4A4';} else {echo'#FCFCFC';}?>;"  onclick="changeMapButton(this)">
<br>
</td>
</tr>
</table>
  </td>
  </tr>

  <tr>
  <td valign="top" style="text-align: left">
  <img id="img_graph" class=\"img_show\" src="./monitor/letkf_monitor/sounding/profile_<?php echo $cvar;?>_001_0001.png" border="0" width="600px">
  </td>
   <td valign="top" style="text-align: left">
<div  align="left" style="position:relative; top:0px; left:0px;" >
  <img id="img_map" class=\"img_show\" src="./monitor/letkf_monitor/sub/map/<?php echo $cvar_map; ?>_anal_rmse_map_d1_sub_0001.png" border="0" width="800px" onclick="changeMapButton(this)" style="position:absolute; top:0px; left:0px;">

<?php

$contents=file('./json/location_prof_'.$cvar.'.txt');

$nprof = count($contents);

$vtop=0.75;
$vleft=0.10;

$fact_left=960;
$fact_top=960;

 for ($iprof = 0; $iprof < $nprof; $iprof++) {
$ele=explode("  ",$contents[$iprof]);
$cprof=$ele[0];
$locleft=intval((floatval($ele[1])-$vleft) * $fact_left) ;
$loctop=intval(($vtop-floatval($ele[2])) * $fact_top)  ;
if( $iprof == 0 ){
  $mcolor = "red" ;
 }else{
  $mcolor = "blue" ;
}

echo "<input id=\"$cprof\" type=\"button\" style=\"position:absolute; top:".$loctop."px; left:".$locleft."px; width:6px; height:6px; margin:0px; padding:0px; background-color:".$mcolor."; outline:none; border-style:solid; border-color:Transparent;\" onclick=\"changeMapButtonDomain(this)\" onmouseover=\"changeMapMouseOn(this)\" onmouseout=\"changeMapMouseOff(this)\">";

}
?>
</div>


   </td>
  </tr>
</table>
</div>

</body>
</html>
