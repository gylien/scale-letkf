program main

implicit real(a-h,o-z)

!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331!

  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! TC vital
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! TC vital
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! TC vital
!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_ref_zero_obs=4004
  INTEGER,PARAMETER :: id_radar_vr_obs=4002
  INTEGER,PARAMETER :: id_radar_prh_obs=4003
!
! Himawari-8 (H08) observations
!
  INTEGER,PARAMETER :: id_H08IR_obs=8800

  INTEGER,PARAMETER :: nobtype=24
!
  CHARACTER(6),PARAMETER :: obtypelist(nobtype)= &
     (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR', 'PHARAD', 'H08IRB', 'TCVITL'/) ! H08


integer,parameter::nmax=100000
real(4)::vlons(nmax),vlats(nmax)
integer::ivars(nmax)
real(4)::wk(8)
!character*120::cfile='/data9/amemiya/realtime_OFP/data/ncepobs_gdas_letkf/2019031800/obs_20190318000000.dat'
character*120::cfile='/data9/amemiya/realtime_OFP/data/ncepobs_gdas/run/dec_prepbufr/outdir_2/obs_20190319000000.dat'

character*40::outdir='./outdir_2/adpsfc_tv/'

character*40::cfile_loc='./loc_upper_temp.txt'
character*7 ::cfile_out='XXX.txt'
character*4 ::cXXXX

character*40::cdummy


integer,parameter::nsta_max=9999
real(4)::vlons_sta(nsta_max)
real(4)::vlats_sta(nsta_max)
integer,parameter::nrec_max=999
integer::irecs(nsta_max)
real(4)::vtime(nsta_max,nrec_max)
real(4)::vlev(nsta_max,nrec_max)
real(4)::vvar(nsta_max,nrec_max)
real(4)::vsdv(nsta_max,nrec_max)



iobstype_target=8
!ivar_target=id_u_obs
!ivar_target=id_t_obs
ivar_target=id_tv_obs
!
open(11,file=trim(cfile_loc),form='formatted')
 ios=0
 iline=0
 do while (ios.eq.0)
  read(11,*,iostat=ios) idummy, vlon, vlat,idummy2
  iline=iline+1
  vlons_sta(iline)=vlon
  vlats_sta(iline)=vlat
 end do
close(11)

nsta=iline  !!! # of stations


irecs=0


ivars=0
ios=0
open(90,file=trim(cfile),form='unformatted')
  do while (ios.eq.0) 
   read(90,iostat=ios) wk

   do i=1,nsta
    if (int(vlons_sta(i)*100.0).eq.int(wk(2)*100.0).and.int(vlats_sta(i)*100.0).eq.int(wk(3)*100.0))then
     if (int(wk(7)).eq.iobstype_target.and.int(wk(1)).eq.ivar_target)then
      irecs(i)=irecs(i)+1
      vtime(i,irecs(i)) = wk(8)
      vlev (i,irecs(i)) = wk(4)
      vvar (i,irecs(i)) = wk(5)
      vsdv (i,irecs(i)) = wk(6)
      exit
     end if    
    end if
   end do

  end do
close(90)


do i=1,nsta
 if (irecs(i).ge.1)then
  write(cXXXX,'(I4)') 1000+i
  write(cfile_out(1:3),'(A3)')cXXXX(2:4)
  open(21,file=trim(outdir)//trim(cfile_out),form='formatted')
   do j=1,irecs(i)
    write(21,*) vtime(i,j), vlev(i,j), vvar(i,j), vsdv(i,j)
   end do
  close(21)
 end if
end do


do istation=1,nstation
 if (ivars(istation).ge.1) write(*,*) istation,vlons(istation),vlats(istation),ivars(istation)
end do
 
stop
end program main
