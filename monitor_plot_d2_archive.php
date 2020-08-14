<?php
Header("Content-Type: image/png");
date_default_timezone_set('UTC');

  $basetime = $argv[1] ;
  if ( ! isset($argv[1])){
   exit("time YYYY-MM-DD HH:MM:SS is not set."); 
  };

  $TimeUinit=strtotime($basetime);

  $timeinit=date("YmdHis",$TimeUinit );
 


  $cmem = 'mdet';

  $Xcampus = 800 ;
  $Ycampus = 520  ;
  
  $image = ImageCreateTrueColor($Xcampus,$Ycampus);

  $black  = ImageColorAllocate($image,128,128,128);
  $gray   = ImageColorAllocate($image,192,192,192);
  $green  = ImageColorAllocate($image,51,255,51);
  $limegreen  = ImageColorAllocate($image,51,205,51);
  $red    = ImageColorAllocate($image,255,51,51);
  $blue   = ImageColorAllocate($image,51,51,255);
  $yellow = ImageColorAllocate($image,255,255,0);
  $white = ImageColorAllocate($image,255,255,255);

  $backgroundcolor = ImageColorAllocateAlpha($image,255,255,255,0);
  $bgc = $backgroundcolor ;

  $dashedline_blue =  array($blue,$blue,$blue,$blue,$blue,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_black = array($black,$black,$black,$black,$black,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_gray = array($gray,$gray,$gray,$gray,$gray,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_red = array($red,$red,$red,$red,$red,$bgc,$bgc,$bgc,$bgc,$bgc);
  
  $CycleSecond = 21600 ;
  $CycleH = $CycleSecond / 3600 ;  

  $NcyclePast = 7 ;
  $FcstPastLimit = 8 ;
  $NcycleFuture = 3;

  $Xleft = 140;
  $Xsize = 60; 

  $PlotH = 1 ;
  $XsizePlot = $Xsize * $PlotH / $CycleH ; 

  $Ytop = 80;
  $Ysize = 100;
  $YsizeF = 150;
  $YsizeBar = 25;

  $YsizeBarFcst = 12;
  $YsizeIntvFcst = 5;

  $Ynow=10;
  $YCD=30;
  $YCHF=50;
 
  $time_offset=file_get_contents('time_offset.txt');
 
  $TimeUinit=strtotime($basetime);

  $timeinit=date("YmdHis",$TimeUinit );
  $nowH = date("H:i",$TimeUinit );
  $nowU = date("U",$TimeUinit );;

  $latestU = $nowU - $nowU % $CycleSecond ;
  $Xoffset = ( $nowU % $CycleSecond ) / $CycleSecond * $Xsize ; 

  for ( $i = (0 - max($NcyclePast,$FcstPastLimit)) ; $i <= $NcycleFuture ; $i++ ) {
  $CycleTimes[$i] = getdate( $latestU + $i * $CycleSecond );
  $CycleHF[$i] = sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleDH[$i] = sprintf('%02d',$CycleTimes[$i]['mday']).sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleMDH[$i] = sprintf('%02d',$CycleTimes[$i]['mon']).sprintf('%02d',$CycleTimes[$i]['mday']).sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleDF[$i] = sprintf('%02d',$CycleTimes[$i]['mon'])."/".sprintf('%02d',$CycleTimes[$i]['mday']);
  if ( $CycleHF[$i] != "00" ) {
  $CycleDF[$i] = '' ;
  } ;
} ;

 $XlocMarker = array() ;
 $YlocMarker = array() ;
 $InitTimes = array() ;

  ImageFill($image,0,0,$backgroundcolor);
  
  ImageSetThickness($image, 2);

/* Prepbufr */		
/*		       
				       
  $Yrect = 90 ;
  $Yrecb = $Yrect + $YsizeBar ;

  $fh = fopen("monitor/monitor_obs.txt", 'rb');

  while ( $line = fgets($fh)) {
    $MMDDHH = substr($line,4,6) ;
    $stat = substr($line,11,5) ;
if (substr($line,13,1) == ":"){
 $stathh  = (substr($line,11,2) + 15) % 24 ;  
$stat = $stathh . substr($line,13,3) ; 
};
for  ( $i = (2 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
   if ( $MMDDHH == $CycleMDH[$i] ) {
    ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$green);
    ImageRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$gray);
    ImageString($image, 2, $Xloc-$Xsize+3, $Yrecb-13, $stat , $black);
   } ;
   } ;
} ;

  fclose($fh) ;

*/
/* GFS analysis */				       
/*				       
  $Yrect = 130 ;
  $Yrecb = $Yrect + $YsizeBar ;

  $fh = fopen("monitor/monitor_gfs.txt", 'rb');

  while ( $line = fgets($fh)) {
    $MMDDHH = substr($line,4,6) ;
    $stat = substr($line,11,5) ;
    $chec = substr($line,17,1) ;

    if (substr($line,13,1) == ":"){
 $stathh  = (substr($line,11,2) + 15) % 24 ;  
$stat = $stathh . substr($line,13,3) ; 
};
for  ( $i = (2 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
   if ( $MMDDHH == $CycleMDH[$i] ) {
    ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$green);
      if ( $chec == "C" ) {
        ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$yellow);
      };
    ImageRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$gray);
    ImageString($image, 2, $Xloc-$Xsize+3, $Yrecb-13, $stat , $black);
   } ;
   } ;
} ;

  fclose($fh) ;

*/
/* Analysis */				     
/*
  $Yrect = 170 ;
  $Yrecb = $Yrect + $YsizeBar ;

$fh = fopen("monitor/monitor_cycle_temp.txt", 'rb');

while ( $line = fgets($fh)) {
    $MMDDHH = substr($line,4,6) ;
    $stat = substr($line,11,13) ;
    if (substr($line,13,1) == ":"){
 $stathh  = (substr($line,11,2) + 15) % 24 ;  
$stat = $stathh . substr($line,13,3) ; 
};
for  ( $i = (2 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
   if ( $MMDDHH == $CycleMDH[$i] && substr($stat,0,1) != "[" ) {
    ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$green);
    ImageRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$gray);
    ImageString($image, 2, $Xloc-$Xsize+3, $Yrecb-13, $stat , $black);
   } ;
   } ;
} ;

fclose($fh) ;

$fh = fopen("monitor/monitor_cycle.txt", 'rb');
 while ( $line = fgets($fh)) {
  $MMDDHH = substr($line,4,6) ;
  $stat = substr($line,15,13) ;

for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

   if ( $MMDDHH == $CycleMDH[$i-1] ) {  // * reference time gap treatment

  if (substr($stat,0,4) == 'plot' ){
  ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,$green);
}elseif (substr($stat,2,1) == '%') {
  $istat=substr($stat,0,2);
  ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc-$Xsize+($istat / 100)
  *$Xsize,$Yrecb,$green);
}elseif (substr($stat,1,1) == '%') {
  $istat=substr($stat,0,1);
  ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc-$Xsize+($istat / 100)
  *$Xsize,$Yrecb,$green);
};

  ImageSetStyle($image,$dashedline_gray);
  ImageRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,IMG_COLOR_STYLED);
  ImageString($image, 2, $Xloc-$Xsize+3, $Yrecb-13, trim($stat) , $black);
 };
 };
};
fclose($fh) ;

*/

/* Forecast */	

  $icount = 0;  

  $Hfcst = 24 ;
  $Yrecb = 120 - $YsizeIntvFcst ;


  for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
   $MMDDHH = $CycleMDH[$i] ;
###
    $test_fsec="";
    exec ("ls -1 data/d2/$dir/$cmem/sfc_prcp/ | tail -n 1 | cut -c 11-16", $test_fsec, $ret);
    if ( $ret == 0 && file_exists( "data/d2/".substr($timeinit,0,4).$MMDDHH.'0000' ) ) {
###
    $Hfcst = intval($test_fsec[0]) / 3600 ;
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

    $Yrect = $Yrecb + $YsizeIntvFcst ;
    $Yrecb = $Yrect + $YsizeBarFcst ;

    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$gray);
###    ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3,    $Yrecb-13, $stat , $black);

   for ( $istep = 0 ; $istep <= ($Hfcst/$PlotH) ; $istep++ ) { 
     $XlocMarker[0][$icount][$istep] = $Xloc+($istep * $XsizePlot);
     $YlocMarker[0][$icount][$istep] = $Yrect; 
     $InitTimes[0][$icount] = $MMDDHH; 
    }; 
    $icount = $icount + 1 ; 
    } ;
} ;


/* MSM forecast */				       

  $icount = 0;  

  $Hfcst = 33 ;
  $Yrecb = 240 - $YsizeIntvFcst ;


for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
    $MMDDHH = $CycleMDH[$i] ;
###
    if ( file_exists( "data/d2/".substr($timeinit,0,4).$MMDDHH.'0000' ) ) {
###

    $Yrect = $Yrecb + $YsizeIntvFcst ;
    $Yrecb = $Yrect + $YsizeBarFcst ;

    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$gray);
###    ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3,    $Yrecb-13, $stat , $black);
   for ( $istep = 0 ; $istep <= ($Hfcst/$PlotH) ; $istep++ ) { 
     $XlocMarker[1][$icount][$istep] = $Xloc+($istep * $XsizePlot);
     $YlocMarker[1][$icount][$istep] = $Yrect; 
     $InitTimes[1][$icount] = $MMDDHH; 
    }; 
    $icount = $icount + 1 ; 
   } ;
} ;


/* JMA precip analysis and forecast */				       
/* 1h refresh (30min not suppoted yet ) */

/* Turned off since 2020.04.01 */
/* Radar data is instead used since 2020.07.14 */

  $Hfcst = 24 ;
  $Yrecb = 380 - $YsizeIntvFcst ;

$sout='';


/*
 exec("ls -l1 data/JMA_precip/anal_d2/realtime/anal_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d2/realtime/anal_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_latest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d2/realtime/fcst_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d2/realtime/fcst_* | tail -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_latest=$sout[0];
*/

/*
 exec("ls -l1 data/JMA_precip/nowcast_d2/realtime/radar_*0000.png | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d2/realtime/radar_*0000.png | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_latest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d2/realtime/nowcast_*0000.png | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d2/realtime/nowcast_*0000.png | tail -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_latest=$sout[0];


$TimeUAnalOldest=strtotime(substr($times_anal_oldest,0,4)."-".substr($times_anal_oldest,4,2)."-".substr($times_anal_oldest,6,2)." ".substr($times_anal_oldest,8,2).":00:00");
$TimeUAnalLatest=strtotime(substr($times_anal_latest,0,4)."-".substr($times_anal_latest,4,2)."-".substr($times_anal_latest,6,2)." ".substr($times_anal_latest,8,2).":00:00");

$TimeUFcstLatest=strtotime(substr($times_fcst_latest,0,4)."-".substr($times_fcst_latest,4,2)."-".substr($times_fcst_latest,6,2)." ".substr($times_fcst_latest,8,2).":00:00");
*/

$TimeUXleft = $latestU + (1-max($NcyclePast,$FcstPastLimit)) * $CycleSecond  +$CycleSecond ;
$TimeUXright = $latestU + $NcycleFuture * $CycleSecond  ;

$TimeUAnalLatest=$TimeUXright;
$TimeUAnalOldest=$TimeUXleft;
$TimeUFcstLatest=$TimeUAnalLatest;

$MMDDHH = date("mdH",$TimeUAnalOldest) ;

$TimeUintv=$TimeUAnalOldest-$TimeUXleft;
$TimeUlengthA=$TimeUAnalLatest-$TimeUAnalOldest;
$TimeUlengthF=$TimeUFcstLatest-$TimeUAnalLatest;
$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv / $CycleSecond * $Xsize);
$Yrect = $Yrecb + $YsizeIntvFcst ;
$Yrecb = $Yrect + $YsizeBarFcst ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);

   for ( $istep = 0 ; $istep <= ($TimeUlength/3600/$PlotH) ; $istep++ ) { 
     $XlocMarker[2][0][$istep] = $Xloc+($istep * $XsizePlot);
     $YlocMarker[2][0][$istep] = $Yrect; 
     $InitTimes[2][0] = $MMDDHH; 
    }; 

/*
   for ( $istep = 0 ; $istep <= ($TimeUlengthF/3600/$PlotH-1) ; $istep++ ) { 
     $XlocMarker[2][1][$istep] = $Xloc+($TimeUlengthA / $CycleSecond * $Xsize)+($istep * $XsizePlot);
     $YlocMarker[2][1][$istep] = $Yrect; 
     $InitTimes[2][1] = $MMDDHH; 
    }; 
*/

/* Legends */

  ImageFilledRectangle($image,1,1,$Xleft,$Ycampus,$white);
#  ImageString($image, 5, 10, 103, 'prepbufr' , $black);
#  ImageString($image, 5, 10, 143, 'GFS_analysis', $black);
#  ImageString($image, 5, 10, 183, 'Analysis' , $black);
  ImageString($image, 5, 10, 103, 'Forecast' , $black);
#  ImageString($image, 5, 10, 248, 'online_nest' , $black);
  ImageString($image, 5, 10, 263, 'MSM_Forecast' , $black);
  ImageString($image, 5, 10, 378, 'JMA_precip' , $black);

/* Vertical lines */

  $Ylinet = $Ytop ;
  $Ylineb = $Ycampus ;

  for ( $i = (1 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {

  $Xloc = $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

  ImageSetThickness($image, 1);
  ImageSetStyle($image,$dashedline_black);
  ImageLine($image,$Xloc,$Ylinet,$Xloc,$Ylineb,IMG_COLOR_STYLED);
  
/*  ImageFilledRectangle($image,5,1,995,55,$white); */

  ImageString($image, 5, $Xloc, $YCD, $CycleDF[$i] , $black);
  ImageString($image, 5, $Xloc, $YCHF, $CycleHF[$i].'Z' , $black);


  } ;


  $Xloc = $Xleft + $Xsize * ( $NcyclePast -1 ) + $Xoffset ; 

  ImageSetStyle($image,$dashedline_red);
  ImageLine($image,$Xloc,$YCD,$Xloc,$Ylineb,IMG_COLOR_STYLED);
  ImageString($image, 5, $Xloc, $Ynow, 'now'.' '.$nowH , $red);

  ImageRectangle($image,0,0,$Xcampus-1,$Ycampus-1,$black);
  
  ImageSaveAlpha($image,TRUE);
  ImagePNG($image,'./monitor/archive/monitor_plot_d2_'.$timeinit.'.png');
  ImageDestroy($image);

  $success = file_put_contents('./json/archive/xlocs_d2_'.$timeinit.'.json', json_encode($XlocMarker)); 
  $success = file_put_contents('./json/archive/ylocs_d2_'.$timeinit.'.json', json_encode($YlocMarker)); 
  $success = file_put_contents('./json/archive/itimes_d2_'.$timeinit.'.json', json_encode($InitTimes)); 

?>
