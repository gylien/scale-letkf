program read_obs_dep
  implicit real(a-h,o-z)
  character(100) :: cfile='./test_data/obsdep.dat'
  integer :: nobs
  real(4) :: wk(11)
  integer :: n, iunit

  write(*,*) '     elm  lon      lat      lev      dat      err      typ    dif    qc    omb      oma'
   iunit=92
  open (iunit, file=cfile, form='unformatted', access='sequential', convert='big_endian')
  do n = 1, 30000
    read(iunit,iostat=ios) wk
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

  if (ityp.eq.8)then
    write(*,'(I5,I5,3F9.2,2E10.2,I4,F9.2,I4,2E10.2)') &
         n,ielm,vlon,vlat,vlev,vdat,verr,ityp,vdif,iqc,vomb,voma
  end if

  end do
  close (iunit)

!    wk(1) = real(obs(set(n))%elm(idx(n)), r_sngl)
!    wk(2) = real(obs(set(n))%lon(idx(n)), r_sngl)
!    wk(3) = real(obs(set(n))%lat(idx(n)), r_sngl)
!    wk(4) = real(obs(set(n))%lev(idx(n)), r_sngl)
!    wk(5) = real(obs(set(n))%dat(idx(n)), r_sngl)
!    wk(6) = real(obs(set(n))%err(idx(n)), r_sngl)
!    wk(7) = real(obs(set(n))%typ(idx(n)), r_sngl)
!    wk(8) = real(obs(set(n))%dif(idx(n)), r_sngl)
!    wk(9) = real(qc(n), r_sngl)
!    wk(10) = real(omb(n), r_sngl)
!    wk(11) = real(oma(n), r_sngl)
!    select case (nint(wk(1)))
!    case (id_u_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    case (id_v_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    case (id_t_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    case (id_tv_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    case (id_q_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!    case (id_ps_obs)
!      wk(5) = wk(5) * 0.01 ! Pa -> hPa
!      wk(6) = wk(6) * 0.01 ! Pa -> hPa
!    case (id_rh_obs)
!      wk(4) = wk(4) * 0.01 ! Pa -> hPa
!      wk(5) = wk(5) * 100.0 ! percent output
!      wk(6) = wk(6) * 100.0 ! percent output
!    case (id_tcmip_obs)
!      wk(5) = wk(5) * 0.01 ! Pa -> hPa
!      wk(6) = wk(6) * 0.01 ! Pa -> hPa
!    end select
!    write (iunit) wk
!  end do
!  close (iunit)

!  return

stop
end program read_obs_dep

