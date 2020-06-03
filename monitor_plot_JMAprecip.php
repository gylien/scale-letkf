<?php Header("Content-Type: image/png");
date_default_timezone_set('UTC');


  $Xcampus = 1500 ;
  $Ycampus = 260  ;
  
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
  
  $CycleSecond = 3600 ;
  $CycleMin = $CycleSecond / 60 ;  
  $NcyclePast = 12 ;
  $FcstPastLimit = 12 ;
  $NcycleFuture = 8;

  $Xleft = 130;
  $Xsize = 72; 

  $PlotMin = 10 ;
  $XsizePlot = $Xsize * $PlotMin / $CycleMin ;

  $PlotMin_JMA = 60 ;
  $XsizePlot_JMA = $Xsize * $PlotMin_JMA / $CycleMin ;
  $PlotMin_JMA_radar = 5 ;
  $XsizePlot_JMA_radar = $Xsize * $PlotMin_JMA_radar / $CycleMin ;

  $Ytop = 80;
  $Ysize = 100;
  $YsizeF = 150;
  $YsizeBar = 50;

  $YsizeBarFcstD2 = 60;
  $YsizeBarFcstD3 = 40;

  $YsizeBarJMA = 50;

  $YsizeIntvFcst = 15;

  $Ynow=10;
  $YCD=30;
  $YCHF=50;
  
  $nowH = date("H:i",time()+9*3600);
#  $nowH = date("H:i");
  $nowU = date("U");


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

/* JMA precip analysis and forecast */				       
/* 1h refresh (30min not suppoted yet ) */

  $Yrecb = 100 - $YsizeIntvFcst ;

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
     $XlocMarker[0][0][$istep] = $Xloc+($istep * $XsizePlot_JMA);
     $YlocMarker[0][0][$istep] = $Yrect; 
     $InitTimes[0][0] = $MMDDHH; 
    }; 

/* JMA precip radar and nowcast 1km */				       
/* 5min refresh  */


  $Yrecb = 180 - $YsizeIntvFcst ;

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
     $XlocMarker[1][0][$istep] = $Xloc+($istep * $XsizePlot_JMA_radar);
     $YlocMarker[1][0][$istep] = $Yrect; 
     $InitTimes[1][0] = $MMDDHH; 
    }; 



/* Legends */

  ImageFilledRectangle($image,1,1,$Xleft,$Ycampus,$white);
  ImageString($image, 5, 10, 120, 'anal_fcst' , $black);
  ImageString($image, 5, 10, 200, 'radar_nowcast' , $black);



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
  ImagePNG($image,'./monitor/monitor_plot_JMAprecip.png');
  ImageDestroy($image);


  $success = file_put_contents('json/xlocs_JMA.json', json_encode($XlocMarker)); 
  $success = file_put_contents('json/ylocs_JMA.json', json_encode($YlocMarker)); 
  $success = file_put_contents('json/itimes_JMA.json',  json_encode($InitTimes)); 

?>
