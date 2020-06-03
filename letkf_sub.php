<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf8">
<meta http-equiv="Cache-Control" content="no-cache">
<title>LETKF monitor test (subdomain)</title>
<link href="css/css_main.css" rel="stylesheet" type="text/css">

<script language="javascript" type="text/javascript">
//<!--

var ImgPathBase = "./monitor/letkf_monitor/sub/" ;

var now_var = "u" ; 
var now_stat = "rmse" ; 
var now_dom = "000" ; 

function changeMapButton(ele) {
var myid=ele.id ;
nums=myid.split("_") ;
if (nums[1] == "var" ) { 
document.getElementById("button_var_"+ now_var).style.backgroundColor = '#FCFCFC' ;
now_var=nums[2] ;
} else if (nums[1] == "stat") {
document.getElementById("button_stat_"+ now_stat).style.backgroundColor = '#FCFCFC' ;
now_stat = nums[2] ;
}
MapImgPath=ImgPathBase + "map/" + now_var + "_anal_" + now_stat + "_map_d1_sub_0001.png"  ;
document.getElementById("img_map").src = MapImgPath ;
MapImgPath=ImgPathBase + now_dom + "/" + now_var + "_letkf_0001.png"  ;
document.getElementById("img_graph").src = MapImgPath ;
ele.style.backgroundColor = '#A4A4A4' ;
}
function changeMapButtonDomain(ele) {
var myid=ele.id ;
document.getElementById(now_dom).style.color = 'Transparent' ;
now_dom = myid;
MapImgPath=ImgPathBase + now_dom + "/" + now_var + "_letkf_0001.png"  ;
document.getElementById("img_graph").src = MapImgPath ;
ele.style.color = 'red' ;
}

function changeMapMouseOn(ele){
ele.style.color = 'red' ;
}
function changeMapMouseOff(ele){
if (ele.id != now_dom){
ele.style.color = 'Transparent' ;
}
}


//-->
</script>
</head>

<body>

<div align="left">
  <p class="p_title">LETKF monitor (subdomain) </p>
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
<input id="button_var_u" type="button" value="U" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #A4A4A4;"  onclick="changeMapButton(this)">
<input id="button_var_v" type="button" value="V" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
<input id="button_var_t" type="button" value="T" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
<input id="button_var_q" type="button" value="Q" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
<input id="button_var_ps" type="button" value="Ps" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
<br>
</td>
</tr>
<tr>
<td><p>Stats:</p></td>
<td>
<input id="button_stat_rmse" type="button" value="RMSE" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #A4A4A4;"  onclick="changeMapButton(this)">
<input id="button_stat_bias" type="button" value="BIAS" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
<input id="button_stat_nobs" type="button" value="Nobs" style="font-size:18px; width:80px; margin:10px; padding:0px; background-color: #FCFCFC;"  onclick="changeMapButton(this)">
</td>
</tr>
</table>
  </td>
  </tr>

  <tr>
  <td valign="top" style="text-align: left">
  <img id="img_graph" class=\"img_show\" src="./monitor/letkf_monitor/sub/000/u_letkf_0001.png" border="0" width="800px">
  </td>
   <td valign="top" style="text-align: left">
<div  align="left" style="position:relative; top:0px; left:0px;" >
  <img id="img_map" class=\"img_show\" src="./monitor/letkf_monitor/sub/map/u_anal_rmse_map_d1_sub_0001.png" border="0" width="800px" onclick="changeMapButton(this)" style="position:absolute; top:0px; left:0px;">

<?php

$contents=file('./json/location_d1_sub.txt');

$icount = 0 ;
$nhdom = 16 ;
$nvdom = 12 ;

$vtop=0.745;
$vleft=0.115;

$fact_left=960;
$fact_top=960;
#$bgc="#A4A4A4";
#echo "<input id=\"000\" type=\"button\" value=\"$cdom\"
#style=\"position:absolute; top:200px; left:800px; width:30px; margin:0px; padding:0px;background-color: $bgc;\"
#onclick=\"changeMapButtonMem(this)\" onmouseover=\"changeMapMouseOn(this)\" onmouseout=\"changeMapMouseOff()\">";

 for ($iv = 0; $iv < $nvdom ; $iv++) {
 for ($ih = 0; $ih < $nhdom; $ih++) {
$icount = $iv * $nhdom + $ih ;
$ele=explode("  ",$contents[$icount]);
$cdom=$ele[0];
$locleft=intval((floatval($ele[1])-$vleft) * $fact_left) ;
$loctop=intval(($vtop-floatval($ele[2])) * $fact_top)  ;
if( $icount == 0 ){
  $tcolor = "red" ;
 }else{
#  $bgc = "#FCFCFC" ;
  $tcolor = "Transparent" ;
}

$count=str_pad(strval($icount), 3, "0", STR_PAD_LEFT) ;
echo "<input id=\"$cdom\" type=\"button\" value=\"$cdom\" style=\"position:absolute; top:".$loctop."px; left:".$locleft."px; width:28px; margin:0px; padding:0px;background-color:Transparent; outline:none; color:".$tcolor."; border-style:solid; border-color:Transparent;\" onclick=\"changeMapButtonDomain(this)\" onmouseover=\"changeMapMouseOn(this)\" onmouseout=\"changeMapMouseOff(this)\">";

}
}
?>
</div>


   </td>
  </tr>
</table>
</div>

</body>
</html>
