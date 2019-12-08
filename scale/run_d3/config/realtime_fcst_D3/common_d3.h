integer,parameter::nmem=50

integer,parameter::nlon=220
integer,parameter::nlat=180

real(4),parameter::axlon(nlon)=(/( 133.84+0.01*real(ilon), ilon=1,nlon)/)
real(4),parameter::axlatSN(nlat)=(/( 33.79+0.01*real(ilat), ilat=1,nlat)/)
real(4),parameter::axlatSNh(nlat+1)=(/( 33.785+0.01*real(ilat), ilat=1,nlat+1)/)

real,parameter::range_lonl=133.8
real,parameter::range_lonr=136.1

real,parameter::range_latl=33.7
real,parameter::range_latr=35.65

integer,parameter::nlonadd=20
integer,parameter::nlatadd=20
integer,parameter::nlon_ext=nlon+nlonadd
integer,parameter::nlat_ext=nlat+nlatadd
real(4),parameter::axlon_ext(nlon_ext)=(/( 133.84+0.01*real(ilon-nlonadd/2), ilon=1,nlon_ext)/)
real(4),parameter::axlatSN_ext(nlat_ext)=(/( 33.79+0.01*real(ilat-nlatadd/2), ilat=1,nlat_ext)/)

integer,parameter::nz=60

real(4),parameter::axz(nz)=(/  &
    55., 165., 275., 385., 495., 608.08, 727.49245, 853.592, 986.75315,  & 
    1127.37135, 1275.86415, 1432.67255, 1598.2623, 1773.1251, 1957.78025,  & 
    2152.77615, 2358.6919, 2576.13905, 2805.7633, 3048.24655, 3304.30895,  & 
    3574.7108, 3860.2551, 4161.7898, 4480.2102, 4816.46215, 5171.5442,  & 
    5546.51075, 5942.47535, 6360.61375, 6802.16795, 7268.4492, 7760.842,  & 
    8280.80855, 8829.89305, 9409.7261, 10022.0298, 10668.62255, 11351.4243,  & 
    12072.4629, 12842.8018, 13642.8018, 14442.8018, 15242.8018, 16042.8013,  & 
    16842.8008, 17642.8008, 18442.8008, 19242.8008, 20042.8008, 20842.8008,  & 
    21642.8008, 22442.8008, 23242.8008, 24042.8008, 24842.8008, 25642.8008,  & 
    26442.8008, 27242.8008, 28042.8008  & 
				 /)

  real(4),parameter::axzh(nz+1)=(/ 0.0, & 
            110.0000,   220.0000,   330.0000,   440.0000,   550.0000,  & 
            666.1600,   788.8249,   918.3591,  1055.1472,  1199.5955,  & 
           1352.1328,  1513.2123,  1683.3123,  1862.9379,  2052.6226,  & 
           2252.9297,  2464.4541,  2687.8240,  2923.7026,  3172.7905,  & 
           3435.8274,  3713.5942,  4006.9160,  4316.6636,  4643.7568,  & 
           4989.1675,  5353.9209,  5739.1006,  6145.8501,  6575.3774,  & 
           7028.9585,  7507.9399,  8013.7441,  8547.8730,  9111.9131,  & 
           9707.5391, 10336.5205, 11000.7246, 11702.1240, 12442.8018,  & 
          13242.8018, 14042.8018, 14842.8018, 15642.8018, 16442.8008,  & 
          17242.8008, 18042.8008, 18842.8008, 19642.8008, 20442.8008,  & 
          21242.8008, 22042.8008, 22842.8008, 23642.8008, 24442.8008,  & 
          25242.8008, 26042.8008, 26842.8008, 27642.8008, 28442.8008  & 
			      /)

integer,parameter::nlev=9
real(4),parameter::plev(nlev)=(/ 1000.0, 925.0, 850.0, 700.0, 500.0, 300.0, 200.0, 100.0, 50.0 /) 
  real(4),parameter::plevh(nlev+1)=(/ 1000.0, 965.0, 885.0, 775.0, 600.0, 400.0, 250.0, 150.0, 75.0, 25.0 /) 

real,parameter::t_ref=273.0
real,parameter::rd=287.0
real,parameter::er=6.4e6
real,parameter::cp=1004.0
real,parameter::grav=9.81
real,parameter::drad=3.141563 / 180.0

real,parameter::rmiss=-9.9999e30

character*120,parameter::cfilename='history.pe000000.nc'

!!! to be inserted by admin.sh !!!
!--cdir_base_fcst--


