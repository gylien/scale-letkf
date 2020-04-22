program read_obs_dep
  implicit real(a-h,o-z)
  character(100) :: cfile='./test_data_ens/obsdep.dat'
  integer :: nobs
  real(4) :: wk(11)
  integer,parameter::nmem=50
  real(4) :: wk_ext(nmem)
  integer :: n, iunit


  vlon_ref=135.77
  vlat_ref=33.45

  ielm_ref=3074


  write(*,*) '     elm  lon      lat      lev      dat      err      typ    dif    qc    omb      oma'
   iunit=92
  open (iunit, file=cfile, form='unformatted', access='sequential', convert='big_endian')
  do n = 1, 200000
    read(iunit,iostat=ios) wk
    read(iunit,iostat=ios) wk_ext
    if (ios.ne.0) exit

    ielm=wk(1)
    vlon=wk(2)
    vlat=wk(3)
    vlev=wk(4)
    vdat=wk(5)
    verr=wk(6)
    ityp=wk(7)
    vdif=wk(8)
    iqc=wk(9)
    vomb=wk(10)
    voma=wk(11)
   if (iqc.eq.0.and.ityp.eq.1 .and. &
      abs(vlon-vlon_ref).le.0.01.and.abs(vlat-vlat_ref).le.0.01.and. &
      ielm.eq.ielm_ref ) then

    write(*,'(I5,I5,5F9.2,I4,F9.2,I4,5E10.2)') &
         n,ielm,vlon,vlat,vlev,vdat,verr,ityp,vdif,iqc,vomb,voma,wk_ext(1:3)
   end if

!   do imem=1,nmem
!    write(*,*) imem, wk_ext(imem)
!   end do 
  end do
  close (iunit)

stop
end program read_obs_dep

