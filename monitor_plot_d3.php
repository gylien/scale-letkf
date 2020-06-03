<?php
Header("Content-Type: image/png");
date_default_timezone_set('UTC');

  $cmem = 'mdet';

  $Xcampus = 2000 ;
  $Ycampus = 520  ;
  
  $image = ImageCreateTrueColor($Xcampus,$Ycampus);

  $black  = ImageColorAllocate($image,128,128,128);
  $gray   = ImageColorAllocate($image,192,192,192);
  $green  = ImageColorAllocate($image,51,255,51);
  $limegreen  = ImageColorAllocate($image,51,205,51);
  $red    = ImageColorAllocate($image,255,51,51);
  $blue   = ImageColorAllocate($image,51,255,255);
  $yellow = ImageColorAllocate($image,255,255,0);
  $white = ImageColorAllocate($image,255,255,255);

  $backgroundcolor = ImageColorAllocateAlpha($image,255,255,255,0);
  $bgc = $backgroundcolor ;

  $dashedline_blue =  array($blue,$blue,$blue,$blue,$blue,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_black = array($black,$black,$black,$black,$black,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_gray = array($gray,$gray,$gray,$gray,$gray,$bgc,$bgc,$bgc,$bgc,$bgc);
  $dashedline_red = array($red,$red,$red,$red,$red,$bgc,$bgc,$bgc,$bgc,$bgc);
  
  $intvD3 = 600 ;
  $CycleSecond = 3600 ;
  $CycleH = $CycleSecond / 3600 ;  

  $CycleMin = $CycleSecond / 60 ;  
  $NcyclePast = 24 ;
  $FcstPastLimit = 24 ;
  $NcycleFuture = 18;

  $Xleft = 130;
  $Xsize = 40; 

  $PlotMin = 10 ;
  $XsizePlot = $Xsize * $PlotMin / $CycleMin ;

  $PlotMin_JMA = 60 ;
  $XsizePlot_JMA = $Xsize * $PlotMin_JMA / $CycleMin ;
  $PlotMin_JMA_radar = 5 ;
  $XsizePlot_JMA_radar = $Xsize * $PlotMin_JMA_radar / $CycleMin ;

  $Ytop = 80;
  $Ysize = 100;
  $YsizeF = 150;
  $YsizeBar = 30;

  $YsizeBarFcstD2 = 60;
  $YsizeBarFcstD3 = 30;

  $YsizeBarJMA = 30;
  $YsizeBarMSM = 30;

  $YsizeIntvFcst = 15;

  $Ynow=10;
  $YCD=30;
  $YCHF=50;
 
  $time_offset=file_get_contents('time_offset.txt');
 
  $nowU = date("U") + intval($time_offset);
  $nowH = date("H:i",$nowU+9*3600);



  $latestU = $nowU - $nowU % $CycleSecond ;
  $Xoffset = ( $nowU % $CycleSecond ) / $CycleSecond * $Xsize ; 

  for ( $i = (1 - max($NcyclePast,$FcstPastLimit)) ; $i <= $NcycleFuture ; $i++ ) {
  $CycleTimes_U[$i] = getdate( $latestU + $i * $CycleSecond );
  $CycleTimes[$i] = getdate( $latestU + $i * $CycleSecond + 9 * 3600 );
  $CycleHF[$i] = sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleDH[$i] = sprintf('%02d',$CycleTimes[$i]['mday']).sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleMDH[$i] = sprintf('%02d',$CycleTimes_U[$i]['mon']).sprintf('%02d',$CycleTimes_U[$i]['mday']).sprintf('%02d',$CycleTimes_U[$i]['hours']);
  $CycleDF[$i] = sprintf('%02d',$CycleTimes[$i]['mon'])."/".sprintf('%02d',$CycleTimes[$i]['mday']);
  if ( $CycleHF[$i] != "00" ) {
  $CycleDF[$i] = '' ;
  } ;
} ;

  ImageFill($image,0,0,$backgroundcolor);
  
  ImageSetThickness($image, 2);


/* Forecast D2 */	

  $HfcstD2 = 24 ;
  $Yrecb = 110 - $YsizeIntvFcst ;

  $icount = 0 ;
  $jcount=0 ;

 $stimes='' ;
 exec("ls -l1 data/d2/ | tail -n 4 | awk '{print $9 \" \" $8}' ",$stimes  ,$ret);

 foreach ($stimes as $line) {
    $dirname = substr($line,0,14) ; 
    $MMDDHHD2 = substr($line,4,6) ; 
    $stathh  = (substr($line,15,2) + 15) % 24 ;  
    $stat = $stathh . substr($line,17,3) ; 

  for  ( $i = (1 - $FcstPastLimit) ; $i <= 0 ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

   if ( $MMDDHHD2 == $CycleMDH[$i] ) {

    $Yrect = $Yrecb + $YsizeIntvFcst ;
    $Yrecb = $Yrect + $YsizeBarFcstD2 ;


    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($HfcstD2 * 60 / $CycleMin * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($HfcstD2 * 60 / $CycleMin * $Xsize),$Yrecb,$gray);
    ImageString($image, 2, $Xloc+($HfcstD2 * 60 / $CycleMin * $Xsize) -$Xsize + 3,    $Yrecb-13, $stat , $black);

   for ( $istep = 0 ; $istep <= ($HfcstD2/1) ; $istep++ ) { 
     $XlocMarker[1][$icount][$istep] = $Xloc+($istep * $Xsize/1);  ### hourly
     $YlocMarker[1][$icount][$istep] = $Yrect; 
     $InitTimes[1][$icount] = $MMDDHHD2; 
    }; 
 $icount = $icount + 1 ;


/* Forecast D3 */	

  $stimesD3='' ;
 exec("ls -l1 data/d3/ref_$dirname/ | tail -n 4 | grep -v total | awk '{print $9 \" \" $8}' ",$stimesD3  ,$ret);

 foreach ($stimesD3 as $line) {
    $MMDDHH = substr($line,4,6) ; 
    $stathh  = (substr($line,15,2) + 15) % 24 ;  
    $stat = $stathh . substr($line,17,3) ; 

    $YmdHMS= substr($line,0,14)  ;
    $fmaxsec="";
    exec("ls -l1 data/d3/ref_$dirname/$YmdHMS/$cmem/sfc_prcp | grep sfc_prcp | tail -n 1 | cut -d '_' -f 3 | cut -c 2-7", $fmaxsec ,$ret);
 
  $HfcstD3 = intval($fmaxsec[0]) / $CycleSecond ;

  for  ( $i = (1 - $FcstPastLimit) ; $i <= $NcycleFuture ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

   if ( $MMDDHH == $CycleMDH[$i] ) {

    $YrecbD3 = $Yrect + $YsizeBarFcstD3 ;

    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($HfcstD3 * 60 / $CycleMin * $Xsize),$YrecbD3,$blue);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($HfcstD3 * 60 / $CycleMin * $Xsize),$YrecbD3,$gray);
    ImageString($image, 2, $Xloc+($HfcstD3 * 60 / $CycleMin * $Xsize) -$Xsize + 3,    $YrecbD3-13, $stat , $black);

   for ( $istep = 0 ; $istep <= ($HfcstD3 * 60 / $PlotMin) ; $istep++ ) { 
     $XlocMarker[0][$jcount][$istep] = $Xloc + ($istep * $XsizePlot) ;
     $YlocMarker[0][$jcount][$istep] = $Yrect ;
     $InitTimes[0][$jcount] = $MMDDHH ;  
     $ParentTimes[0][$jcount] =$MMDDHHD2 ;
}

 $jcount = $jcount + 1 ;
    } ;
   } ;
} ;



$fh = fopen("monitor/monitor_fcst_d3.txt", 'rb');
 while ( $line = fgets($fh)) {
  $MMDDHHref = substr($line,4,6) ;
  $MMDDHH = substr($line,19,6) ;
  $stat = substr($line,30,13) ;

  $HfcstD3temp=6 ;

for  ( $i = (1 - $FcstPastLimit) ; $i <= $NcycleFuture ; $i++ ) {
    $Xloc= $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

   if ( $MMDDHH == $CycleMDH[$i] and $MMDDHHref == $MMDDHHD2 ) {
  $YrecbD3 = $Yrect + $YsizeBarFcstD3 ;

  if (substr($stat,0,4) == 'plot' ){
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc +($HfcstD3temp * 60 / $CycleMin * $Xsize),$YrecbD3,$blue);
}elseif (substr($stat,2,1) == '%') {
  $istat=substr($stat,0,2);
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($istat / 100) *($HfcstD3temp * 60 / $CycleMin * $Xsize),$YrecbD3,$blue);
}elseif (substr($stat,1,1) == '%') {
  $istat=substr($stat,0,1);
  ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($istat / 100) *($HfcstD3temp * 60 / $CycleMin * $Xsize),$YrecbD3,$green);
};

  ImageSetStyle($image,$dashedline_gray);
  ImageRectangle($image,$Xloc,$Yrect,$Xloc+ ($HfcstD3temp * 60 / $CycleMin * $Xsize),$YrecbD3,IMG_COLOR_STYLED); 
  ImageString($image, 2, $Xloc+($HfcstD3temp * 60 / $CycleMin * $Xsize) -$Xsize + 3, $YrecbD3-13, trim($stat) , $black);	
 };

 };
};
fclose($fh) ;

    } ;
   } ;
} ;



/* MSM forecast */				       

  $icount = 0;  

  $Hfcst = 33 ;
  $Yrecb = 340 - $YsizeIntvFcst ;
  $PlotH = 1 ;

 $stimes='';
 exec("ls -l1 data/msm_d3/ | tail -n 20 | awk '{print $9 \" \" $8}'",$stimes ,$ret);

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
    $Yrecb = $Yrect + $YsizeBarMSM ;

    ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$green);
    ImageRectangle($image,$Xloc,$Yrect,$Xloc+($Hfcst / $CycleH * $Xsize),$Yrecb,$gray);
    ImageString($image, 2, $Xloc+($Hfcst / $CycleH * $Xsize) -$Xsize + 3,    $Yrecb-13, $stat , $black);
   for ( $istep = 0 ; $istep <= $Hfcst ; $istep++ ) { 
     $XlocMarker[2][$icount][$istep] = $Xloc+($istep * 6 * $XsizePlot); ### TORIAEZU
     $YlocMarker[2][$icount][$istep] = $Yrect; 
     $InitTimes[2][$icount] = $MMDDHH; 
    }; 
    $icount = $icount + 1 ; 
   } ;
   } ;
} ;


/* JMA precip analysis and forecast */				       
/* 1h refresh (30min not suppoted yet ) */

/* TORI AEZU DISABLED

  $Yrecb = 415 - $YsizeIntvFcst ;

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d3/realtime/anal_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d3/realtime/anal_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_latest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d3/realtime/fcst_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/anal_d3/realtime/fcst_* | tail -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_latest=$sout[0];

$TimeUAnalOldest=strtotime(substr($times_anal_oldest,0,4)."-".substr($times_anal_oldest,4,2)."-".substr($times_anal_oldest,6,2)." ".substr($times_anal_oldest,8,2).":00:00");
$TimeUAnalLatest=strtotime(substr($times_anal_latest,0,4)."-".substr($times_anal_latest,4,2)."-".substr($times_anal_latest,6,2)." ".substr($times_anal_latest,8,2).":00:00");

$TimeUFcstLatest=strtotime(substr($times_fcst_latest,0,4)."-".substr($times_fcst_latest,4,2)."-".substr($times_fcst_latest,6,2)." ".substr($times_fcst_latest,8,2).":00:00");

$TimeUXleft = $latestU + (1-max($NcyclePast,$FcstPastLimit)) * $CycleSecond ;

$TimeUAnalOldest=max($TimeUAnalOldest,$TimeUXleft);
$TimeUintv=$TimeUAnalOldest-$TimeUXleft;
$TimeUlengthA=$TimeUAnalLatest-$TimeUAnalOldest;
$TimeUlengthF=$TimeUFcstLatest-$TimeUAnalLatest;
$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv / $CycleSecond * $Xsize);
$Yrect = $Yrecb + $YsizeIntvFcst ;
$Yrecb = $Yrect + $YsizeBarJMA ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);


$MMDDHH = date("mdH",$TimeUAnalOldest) ;

   for ( $istep = 0 ; $istep <= ($TimeUlength/60/$PlotMin_JMA) ; $istep++ ) { 
     $XlocMarker[1][0][$istep] = $Xloc+($istep * $XsizePlot_JMA);
     $YlocMarker[1][0][$istep] = $Yrect; 
     $InitTimes[1][0] = $MMDDHH; 
    }; 


*/

/* JMA precip radar and nowcast 1km */				       
/* 5min refresh  */

/* TORI AEZU DISABLED

  $Yrecb = 465 - $YsizeIntvFcst ;

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d3/realtime/radar_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d3/realtime/radar_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_anal_latest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d3/realtime/nowcast_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d3/realtime/nowcast_* | tail -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_fcst_latest=$sout[0];

$TimeUAnalOldest=strtotime(substr($times_anal_oldest,0,4)."-".substr($times_anal_oldest,4,2)."-".substr($times_anal_oldest,6,2)." ".substr($times_anal_oldest,8,2).":".substr($times_anal_oldest,10,2).":00");
$TimeUAnalLatest=strtotime(substr($times_anal_latest,0,4)."-".substr($times_anal_latest,4,2)."-".substr($times_anal_latest,6,2)." ".substr($times_anal_latest,8,2).":".substr($times_anal_latest,10,2).":00");
$TimeUFcstLatest=strtotime(substr($times_fcst_latest,0,4)."-".substr($times_fcst_latest,4,2)."-".substr($times_fcst_latest,6,2)." ".substr($times_fcst_latest,8,2).":".substr($times_fcst_latest,10,2).":00");

$TimeUXleft = $latestU + (1-max($NcyclePast,$FcstPastLimit)) * $CycleSecond ;

$TimeUAnalOldest=max($TimeUAnalOldest,$TimeUXleft);
$TimeUintv=$TimeUAnalOldest-$TimeUXleft;
$TimeUlengthA=$TimeUAnalLatest-$TimeUAnalOldest;
$TimeUlengthF=$TimeUFcstLatest-$TimeUAnalLatest;
$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv *$Xsize / $CycleSecond);
$Yrect = $Yrecb + $YsizeIntvFcst ;
$Yrecb = $Yrect + $YsizeBarJMA ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);


$MMDDHH = date("mdH",$TimeUAnalOldest) ;

   for ( $istep = 0 ; $istep <= ($TimeUlength/60/$PlotMin_JMA_radar) ; $istep++ ) { 
     $XlocMarker[2][0][$istep] = $Xloc+($istep * $XsizePlot_JMA_radar);
     $YlocMarker[2][0][$istep] = $Yrect; 
     $InitTimes[2][0] = $MMDDHH; 
    }; 


*/

/* Legends */

  ImageFilledRectangle($image,1,1,$Xleft,$Ycampus,$white);
  ImageString($image, 5, 10, 105, 'SCALE_LETKF' , $black);
  ImageString($image, 5, 10, 120, 'Forecast' , $black);
  ImageString($image, 5, 10, 135, " Domain 2" , $green);
  ImageString($image, 5, 10, 150, " Domain 3" , $blue);
  ImageString($image, 5, 10, 380, 'MSM' , $black);
#  ImageString($image, 5, 10, 380, 'JMA_precip' , $black);
#  ImageString($image, 5, 10, 425, 'anal_fcst' , $black);
#  ImageString($image, 5, 10, 465, 'radar_nowcast' , $black);



/* Vertical lines */

  $Ylinet = $Ytop ;
  $Ylineb = $Ycampus ;

  for ( $i = (1 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {

  $Xloc = $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;

  if ( ($CycleHF[$i] % 3) == 0 ) {
  ImageSetThickness($image, 2);
  }else{
  ImageSetThickness($image, 1);
  };
  ImageSetStyle($image,$dashedline_black);
  ImageLine($image,$Xloc,$Ylinet,$Xloc,$Ylineb,IMG_COLOR_STYLED);


/*  ImageFilledRectangle($image,5,1,995,55,$white); */

  if ( ($CycleHF[$i] % 3) == 0 ) {
  ImageString($image, 5, $Xloc-10, $YCD, $CycleDF[$i] , $black);
  ImageString($image, 5, $Xloc-10, $YCHF, $CycleHF[$i].":00" , $black);
  } ;

  } ;


  $Xloc = $Xleft + $Xsize * ( $NcyclePast -1 ) + $Xoffset ; 

  ImageSetStyle($image,$dashedline_red);
  ImageLine($image,$Xloc,$YCD,$Xloc,$Ylineb,IMG_COLOR_STYLED);
  ImageString($image, 5, $Xloc, $Ynow, 'now'.' '.$nowH , $red);

  ImageRectangle($image,0,0,$Xcampus-1,$Ycampus-1,$black);
  
  ImageSaveAlpha($image,TRUE);
  ImagePNG($image,'./monitor/monitor_plot_d3.png');
  ImageDestroy($image);


  $success = file_put_contents('./json/xlocs_d3.json', json_encode($XlocMarker)); 
  $success = file_put_contents('./json/ylocs_d3.json', json_encode($YlocMarker)); 
  $success = file_put_contents('./json/itimes_d3.json',  json_encode($InitTimes)); 
  $success = file_put_contents('./json/ptimes_d3.json', json_encode($ParentTimes)); 

?>
