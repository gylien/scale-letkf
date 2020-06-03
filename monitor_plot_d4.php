<?php
Header("Content-Type: image/png");
date_default_timezone_set('UTC');

  $Xcampus = 1000 ;
  $Ycampus = 300  ;
  
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
  
  $CycleSecond = 300 ;
  $NcyclePast = 6 ;
  $NcycleFuture = 6;

  $Xleft = 180;
  $Xsize = 60; 

  $XsizePlot = $Xsize / ($CycleSecond / 30) ;

  $Ytop = 80;
  $Ysize = 100;
  $YsizeF = 150;

  $YsizeBar = 40;
  $YsizeIntv = 15;

  $Ynow=10;
  $YCD=30;
  $YCHF=50;

  $zref='z01127m';

  $time_offset=file_get_contents('time_offset.txt');
 
  $nowU = date("U") + intval($time_offset);
  $nowH = date("H:i",$nowU+9*3600);

  $latestU = $nowU - $nowU % $CycleSecond ;
  $Xoffset = ( $nowU % $CycleSecond ) / $CycleSecond * $Xsize ; 

  for ( $i = (1 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {
  $CycleTimes[$i] = getdate( $latestU + $i * $CycleSecond + 9 * 3600 );
  $CycleMF[$i] = sprintf('%02d',$CycleTimes[$i]['minutes']);
  $CycleHF[$i] = sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleDH[$i] = sprintf('%02d',$CycleTimes[$i]['mday']).sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleMDH[$i] = sprintf('%02d',$CycleTimes[$i]['mon']).sprintf('%02d',$CycleTimes[$i]['mday']).sprintf('%02d',$CycleTimes[$i]['hours']);
  $CycleDF[$i] = sprintf('%02d',$CycleTimes[$i]['mon'])."/".sprintf('%02d',$CycleTimes[$i]['mday']);
  if ( $CycleMF[$i] != "00" ) {
  $CycleDF[$i] = '' ;
  } ;
} ;

  ImageFill($image,0,0,$backgroundcolor);
  
  ImageSetThickness($image, 2);



/* PAWR and SCALE-LETKF */				       
/* 30sec refresh */

  $Yrecb = 120 - $YsizeIntv ;

$sout='';
 exec("ls -l1 data/d4/realtime/".$zref."/anal_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/d4/realtime/".$zref."/anal_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_latest=$sout[0];

$sout='';
 exec("ls -l1 data/d4/realtime/".$zref."/fcst_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_oldest=$sout[0];
$sout_fcst = $sout[0] ;
$sout='';
 exec("ls -l1 data/d4/realtime/".$zref."/fcst_* | tail -n 1 | awk '{print $9}' | rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_latest=$sout[0];

$TimeURadarOldest=strtotime(substr($times_radar_oldest,0,4)."-".substr($times_radar_oldest,4,2)."-".substr($times_radar_oldest,6,2).
" " .substr($times_radar_oldest,8,2).":".substr($times_radar_oldest,10,2).":".substr($times_radar_oldest,12,2));
$TimeURadarLatest=strtotime(substr($times_radar_latest,0,4)."-".substr($times_radar_latest,4,2)."-".substr($times_radar_latest,6,2).
" ".substr($times_radar_latest,8,2).":".substr($times_radar_latest,10,2).":".substr($times_radar_latest,12,2));
if (isset($sout_fcst)){
$TimeUFcstLatest=strtotime(substr($times_nowcast_latest,0,4)."-".substr($times_nowcast_latest,4,2)."-".substr($times_nowcast_latest,6,2).
" ".substr($times_nowcast_latest,8,2).":".substr($times_nowcast_latest,10,2).":".substr($times_nowcast_latest,12,2));
}else{
$TimeUFcstLatest=$TimeURadarLatest ;
}

$TimeUXleft = $latestU + (1-$NcyclePast) * $CycleSecond ;

$TimeURadarOldest=max($TimeURadarOldest,$TimeUXleft);
$TimeUintv=$TimeURadarOldest-$TimeUXleft;
$TimeUlengthA=$TimeURadarLatest-$TimeURadarOldest;
$TimeUlengthF=$TimeUFcstLatest-$TimeURadarLatest;
$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv / $CycleSecond * $Xsize);
$Yrect = $Yrecb + $YsizeIntv ;
$Yrecb = $Yrect + $YsizeBar ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);

/* PAWR nowcastF */				       
/* 30sec refresh */

/* Currently not used

  $Yrecb = 170 - $YsizeIntv ;

$sout='';
 exec("ls -l1 data/d4/realtime/radar_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/d4/realtime/radar_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_latest=$sout[0];

$sout='';
 exec("ls -l1 data/d4/realtime/nowcast_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/d4/realtime/nowcast_* | tail -n 1 | awk '{print $9}' | rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_latest=$sout[0];

$TimeURadarOldest=strtotime(substr($times_radar_oldest,0,4)."-".substr($times_radar_oldest,4,2)."-".substr($times_radar_oldest,6,2)." ".substr($times_radar_oldest,8,2).":".substr($times_radar_oldest,10,2).":00");
$TimeURadarLatest=strtotime(substr($times_radar_latest,0,4)."-".substr($times_radar_latest,4,2)."-".substr($times_radar_latest,6,2)." ".substr($times_radar_latest,8,2).":".substr($times_radar_latest,10,2).":00");

$TimeUNowcastLatest=strtotime(substr($times_nowcast_latest,0,4)."-".substr($times_nowcast_latest,4,2)."-".substr($times_nowcast_latest,6,2)." ".substr($times_nowcast_latest,8,2).":".substr($times_nowcast_latest,10,2).":00");

$TimeUXleft = $latestU + (1-$NcyclePast) * $CycleSecond ;

$TimeURadarOldest=max($TimeURadarOldest,$TimeUXleft);
$TimeUintv=$TimeURadarOldest-$TimeUXleft;
$TimeUlengthA=$TimeURadarLatest-$TimeURadarOldest;
$TimeUlengthF=$TimeUNowcastLatest-$TimeURadarLatest;

$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv / $CycleSecond * $Xsize);
$Yrect = $Yrecb + $YsizeIntv ;
$Yrecb = $Yrect + $YsizeBar ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);

*/

/* JMA radar and nowcast */				       
/* 5min refresh */

  $Yrecb = 220 - $YsizeIntv ;

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d4/realtime/radar_* | head -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d4/realtime/radar_* | tail -n 1 | awk '{print $9}' |   rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_radar_latest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d4/realtime/nowcast_* | head -n 1 | awk '{print $9}' |  rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_oldest=$sout[0];

$sout='';
 exec("ls -l1 data/JMA_precip/nowcast_d4/realtime/nowcast_* | tail -n 1 | awk '{print $9}' | rev | cut -d '/' -f 1 | rev | sed -e 's/[^0-9]//g' ",$sout  ,$ret);
$times_nowcast_latest=$sout[0];

$TimeUJMARadarOldest=strtotime(substr($times_radar_oldest,0,4)."-".substr($times_radar_oldest,4,2)."-".substr($times_radar_oldest,6,2)." ".substr($times_radar_oldest,8,2).":".substr($times_radar_oldest,10,2).":00");
$TimeUJMARadarLatest=strtotime(substr($times_radar_latest,0,4)."-".substr($times_radar_latest,4,2)."-".substr($times_radar_latest,6,2)." ".substr($times_radar_latest,8,2).":".substr($times_radar_latest,10,2).":00");

$TimeUJMANowcastLatest=strtotime(substr($times_nowcast_latest,0,4)."-".substr($times_nowcast_latest,4,2)."-".substr($times_nowcast_latest,6,2)." ".substr($times_nowcast_latest,8,2).":".substr($times_nowcast_latest,10,2).":00");

$TimeUXleft = $latestU + (1-$NcyclePast) * $CycleSecond ;
$TimeUJMARadarOldest=max($TimeUJMARadarOldest,$TimeUXleft);
$TimeUJMARadarLatest=min($TimeUJMARadarOldest+1800,$TimeUJMARadarLatest);
$TimeUJMANowcastLatest=min($TimeUJMARadarOldest+2400,$TimeUJMANowcastLatest);

$TimeUintv=$TimeUJMARadarOldest-$TimeUXleft;
$TimeUlengthA=$TimeUJMARadarLatest-$TimeUJMARadarOldest;
$TimeUlengthF=$TimeUJMANowcastLatest-$TimeUJMARadarLatest;
$TimeUlength=$TimeUlengthA+$TimeUlengthF;

$Xloc=$Xleft+($TimeUintv / $CycleSecond * $Xsize);
$Yrect = $Yrecb + $YsizeIntv ;
$Yrecb = $Yrect + $YsizeBar ;
 ImageFilledRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrecb,$limegreen);
 ImageFilledRectangle($image,$Xloc+($TimeUlengthA / $CycleSecond * $Xsize),$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$green);
 ImageRectangle($image,$Xloc,$Yrect,$Xloc+($TimeUlength / $CycleSecond * $Xsize),$Yrecb,$gray);



/* Time bar output */

$TimeUOldest=min($TimeURadarOldest,$TimeUJMARadarOldest);
#$TimeULatest=max($TimeUNowcastLatest,$TimeUFcstLatest,$TimeUJMARadarLatest);
$TimeULatest=max($TimeUFcstLatest,$TimeUJMANowcastLatest);
$TimeUlength=$TimeULatest-$TimeUOldest;

   for ( $istep = 0 ; $istep <= ($TimeUlength/30) ; $istep++ ) { 
     $XlocMarker[0][0][$istep] = $Xloc+($istep * $XsizePlot);
     $YlocMarker[0][0][$istep] = $Ytop; 
     $InitTimes[0][0] = date('YmdHis',$TimeUOldest); 
    }; 


/* Legends */

  ImageFilledRectangle($image,1,1,$Xleft,$Ycampus,$white);
  ImageString($image, 5, 10, 120, 'PAWR_SCALE_LETKF' , $black);
  ImageString($image, 5, 10, 220, 'JMA_radar_nowcast' , $black);

/* Vertical lines */

  $Ylinet = $Ytop ;
  $Ylineb = $Ycampus ;

  for ( $i = (1 - $NcyclePast) ; $i <= $NcycleFuture ; $i++ ) {

  $Xloc = $Xleft + $Xsize * ( $i - (1 - $NcyclePast) ) ;


  ImageSetThickness($image, 2);
  ImageSetStyle($image,$dashedline_black);
  ImageLine($image,$Xloc,$Ylinet,$Xloc,$Ylineb,IMG_COLOR_STYLED);


/*  ImageFilledRectangle($image,5,1,995,55,$white); */

  if ( $CycleMF[$i] == "00" ){
  ImageString($image, 5, $Xloc-20, $YCHF, $CycleHF[$i].":".$CycleMF[$i] , $black);
  ImageString($image, 5, $Xloc-20, $YCHF - 25, $CycleDF[$i] , $black);
 } else{
  ImageString($image, 5, $Xloc, $YCHF, $CycleMF[$i] , $black);
}
  } ;


  $Xloc = $Xleft + $Xsize * ( $NcyclePast -1 ) + $Xoffset ; 

  ImageSetStyle($image,$dashedline_red);
  ImageLine($image,$Xloc,$YCD,$Xloc,$Ylineb,IMG_COLOR_STYLED);
  ImageString($image, 5, $Xloc, $Ynow, 'now'.' '.$nowH , $red);

  ImageRectangle($image,0,0,$Xcampus-1,$Ycampus-1,$black);
  
  ImageSaveAlpha($image,TRUE);
  ImagePNG($image,'./monitor/monitor_plot_d4.png');
  ImageDestroy($image);

  $success = file_put_contents('./json/xlocs_d4.json', json_encode($XlocMarker)); 
  $success = file_put_contents('./json/ylocs_d4.json', json_encode($YlocMarker)); 
  $success = file_put_contents('./json/itimes_d4.json', json_encode($InitTimes)); 
?>
