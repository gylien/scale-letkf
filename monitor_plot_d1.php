<?php
Header("Content-Type: image/png");
date_default_timezone_set('UTC');

  $Xcampus = 1500 ;
  $Ycampus = 480  ;
  
  $image = ImageCreateTrueColor($Xcampus,$Ycampus);

  $black  = ImageColorAllocate($image,128,128,128);
  $gray   = ImageColorAllocate($image,192,192,192);
  $green  = ImageColorAllocate($image,51,255,51);
  $red    = ImageColorAllocate($image,255,51,51);
  $blue   = ImageColorAllocate($image,51,51,255);
  $yellow = ImageColorAllocate($image,255,255,0);
  $white = ImageColorAllocate($image,255,255,255);

  $backgroundcolor = ImageColorAllocateAlpha($image,255,255,255,0);
  $bgc = $backgroundcolor ;

  $dashedline_blue = array($blue,$blue,$blue,$blue,$blue,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_black = array($black,$black,$black,$black,$black,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_gray = array($gray,$gray,$gray,$gray,$gray,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_red = array($red,$red,$red,$red,$red,$bgc,$bgc,$bgc,$bgc,$bgc);
  
  $CycleSecond = 21600 ;
  $CycleH = $CycleSecond / 3600 ;  

  $NcyclePast = 7 ;
  $FcstPastLimit = 8 ;
  $NcycleFuture = 20;
  $Hfcst = 120 ;

  $Xleft = 140;
  $Xsize = 45; 

  $PlotH = $CycleH ;
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

 
  $nowU = date("U") + intval($time_offset);
  $nowH = date("H:i",$nowU);
<<<<<<< HEAD
=======
  echo $time_offset.' '.$nowH; 

>>>>>>> 504b3ad... initial commit

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

/* GFS analysis */				       
				       
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




/* Analysis */				     

  $Yrect = 300 ;
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
   if ( $MMDDHH == $CycleMDH[$i] ) {
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
  ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc-$Xsize+($istat / 100)  *$Xsize,$Yrecb,$green);
}elseif (substr($stat,1,1) == '%') {
  $istat=substr($stat,0,1);
  ImageFilledRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc-$Xsize+($istat / 100)  *$Xsize,$Yrecb,$green);
};
  ImageSetStyle($image,$dashedline_gray);
  ImageRectangle($image,$Xloc-$Xsize,$Yrect,$Xloc,$Yrecb,IMG_COLOR_STYLED);
  ImageString($image, 2, $Xloc-$Xsize+3, $Yrecb-13, trim($stat) , $black);
 };
 };
};
fclose($fh) ;




/* Forecast */	

  $icount = 0;  

  $Yrecb = 340 - $YsizeIntvFcst ;

 $stimes='';
 exec("ls -l1 data/d1/ | tail -n 20 | awk '{print $9 \" \" $8}'",$stimes ,$ret);



 foreach ($stimes as $line) {
 $MMDDHH = substr($line,4,6) ;
 $dir=substr($line,0,14);
 $test_fsec='';
 exec ("ls -t1 data/d1/$dir/sfc_prcp/ | grep sfc_prcp_f | head -n 1 | cut -c 11-16", $test_fsec, $ret);
 if ( $ret == 0 && intval($test_fsec[0]) == $Hfcst * 3600 ) {
    $stat = substr($line,15) ;
    if (substr($line,17,1) == ":"){
 $stathh  = (substr($line,15,2) + 15) % 24 ;  
 $stat = $stathh . substr($line,17,3) ; 
};

for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

   if ( $MMDDHH == $CycleMDH[$i] ) {

    $Yrect = $Yrecb + $YsizeIntvFcst ;
    $Yrecb = $Yrect + $YsizeBarFcst ;

    if ( strpos($stat,":") !== false ) {
    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$gray);

   for ( $istep = 0 ; $istep <= ($Hfcst/$PlotH) ; $istep++ ) { 
     $XlocMarker[0][$icount][$istep] = $Xloc+($istep * $XsizePlot);
     $YlocMarker[0][$icount][$istep] = $Yrect; 
     $InitTimes[0][$icount] = $MMDDHH; 
    }; 
    $icount = $icount + 1 ; 

    }    elseif ( strpos($stat,"%") !== false ){ 
    $istat=preg_replace("/[^0-9]/","",$stat); 
    $stat = $istat.'%'; 
    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($istat / 100) *($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageSetStyle($image,$dashedline_gray);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+ ($Hfcst / $CycleH * $Xsize),$Yrecb,IMG_COLOR_STYLED); 
    } ;

    if ( strpos($stat,0,4) == "init" ){ 
    $stat = "init" ; 
    ImageSetStyle($image,$dashedline_gray);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+ ($Hfcst / $CycleH * $Xsize),$Yrecb,IMG_COLOR_STYLED); 
    } ;
    
    ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3, $Yrecb-13, $stat , $black);
   } ;
   } ;
} ;
} ;


$fh = fopen("monitor/monitor_fcst.txt", 'rb');
 while ( $line = fgets($fh)) {
  $MMDDHH = substr($line,4,6) ;
  $stat = substr($line,15) ;

for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
   if ( $MMDDHH == $CycleMDH[$i] ) {

  $Yrect = $Yrecb + $YsizeIntvFcst ;
  $Yrecb = $Yrect + $YsizeBarFcst ;

  if (substr($stat,0,4) == 'plot' ){
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc +($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
}elseif (substr($stat,2,1) == '%') {
  $istat=substr($stat,0,2);
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($istat / 100) *($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
}elseif (substr($stat,1,1) == '%') {
  $istat=substr($stat,0,1);
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($istat / 100) *($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
};
  ImageSetStyle($image,$dashedline_gray);
  ImageRectangle($image,$Xloc,$Yrect,$Xloc+ ($Hfcst / $CycleH * $Xsize),$Yrecb,IMG_COLOR_STYLED); 
  ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3, $Yrecb-13, trim($stat) , $black);	
 };

 };
};

fclose($fh) ;


/* GFS forecast */				       

  $icount = 0 ; 

  $Yrecb = 170 - $YsizeIntvFcst ;

$stimes='';
 exec("ls -l1 data/gfs/ | tail -n 20 | awk '{print $9 \" \" $8}'",$stimes  ,$ret);

 foreach ($stimes as $line) {
    $MMDDHH = substr($line,4,6) ;
    $stat = substr($line,15,5) ;
    if (substr($line,17,1) == ":"){
 $stathh  = (substr($line,15,2) + 15) % 24 ;  
 $stat = $stathh . substr($line,17,3) ; 
};

for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;
   if ( $MMDDHH == $CycleMDH[$i] ) {
    $Yrect = $Yrecb + $YsizeIntvFcst ;
    $Yrecb = $Yrect + $YsizeBarFcst ;

    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$gray);
    ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3, $Yrecb-13, $stat , $black);

   for ( $istep = 0 ; $istep <= ($Hfcst/$PlotH) ; $istep++ ) { 
     $XlocMarker[1][$icount][$istep] = $Xloc+($istep * $XsizePlot); 
     $YlocMarker[1][$icount][$istep] = $Yrect;
     $InitTimes[1][$icount] = $MMDDHH; 
   };  
    $icount = $icount+1 ; 
   } ;
   } ;
} ;




/* Legends */

  ImageFilledRectangle($image,1,1,$Xleft,$Ycampus,$white);
  ImageString($image, 5, 10, 103, 'prepbufr' , $black);
  ImageString($image, 5, 10, 143, 'GFS_analysis', $black);
  ImageString($image, 5, 10, 183, 'GFS_forecast', $black);
  ImageString($image, 5, 10, 313, 'Analysis' , $black);
  ImageString($image, 5, 10, 353, 'Forecast' , $black);

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
  ImagePNG($image,'./monitor/monitor_plot_d1.png');
  ImageDestroy($image);

  $success = file_put_contents('json/xlocs.json', json_encode($XlocMarker)); 
  $success = file_put_contents('json/ylocs.json', json_encode($YlocMarker)); 
  $success = file_put_contents('json/itimes.json', json_encode($InitTimes)); 

?>
