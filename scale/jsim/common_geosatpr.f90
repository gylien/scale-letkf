  module common_geosatpr
    use common, only: r_size
    implicit none

    ! Array
    integer :: is,ie,js,je,ks,ke

    ! Constants
    real(r_size), parameter :: q_min_condensate = 1.0e-5  ! siginificant minimum condensate mixing ratio [g/m3]
    real(r_size), parameter :: const_g = 9.807e0 ! graviational acceleration on earth at sea level [m/s2]
    real(r_size), parameter :: const_Rd = 287.053 ! dry-air gas constant [J/K/kg] 
    real(r_size), parameter :: const_Rd_Rv = 0.622e0 ! epsilon 
    real(r_size), parameter :: const_pi = 3.14159e0 ! pai 
    real(r_size), parameter :: const_degrad = const_pi/180.d0
    real(r_size), parameter :: const_eosqv=(1.0/const_Rd_Rv-1.0)  ! conversion factor for degrees to radian
    real(r_size), parameter :: const_Kel2Cel  = 273.15e0 !conversion constant from Kelvin to Celcius
    real(r_size), parameter :: const_tmelt = 273.16e0 ! 0 deg Celcius
    real(r_size), parameter :: undefined = -999.d0

    ! Output
    real(r_size) :: z_out ! Radar reflectivity
    real(r_size) :: dbz ! Radar reflectivity [dBZ]
    !integer :: qc ! Quality control

    ! Input
    ! Atmosphere parameters at layer level
    type atmos_parameter
        real(r_size)  :: hgt      ! heigh [km]
        !real(r_size)  :: dhgt     ! thickness of layer [km]
        !real(r_size)  :: t_air    ! layer air temperature [degK]
        !real(r_size)  :: press    ! layer pressure [hPa]
        !real(r_size)  :: rh       ! relative humidity [%]
        real(r_size)  :: omega    ! verticall velocity  [m/s]
        !real(r_size)  :: exner    ! exner function [-]
    end type atmos_parameter
    type ( atmos_parameter ), allocatable, dimension(:,:,:) :: atmos
    
    ! Atmosphere parameters at staggered level
    type atmos_stagger_parameter
        real(r_size)  :: hgt      ! height at interface [km]
        real(r_size)  :: t_air    ! layer air temperature [degK]
        real(r_size)  :: omega    ! verticall velocity  [m/s]
    end type atmos_stagger_parameter
    type ( atmos_stagger_parameter ), allocatable, dimension(:,:,:) :: atmos_stag

    type surface_parameter
        real(r_size)  :: lat
        real(r_size)  :: lon
        real(r_size)  :: elev
        integer :: k1layer     ! lowest layer that has defined values above the surface
    end type surface_parameter
    type ( surface_parameter ), allocatable, dimension(:,:) :: surface

    ! 1-moment microphysics parameter
    type particle_gce ! particle type
        real(r_size)  :: cloud,rain,ice,snow,graupel,hail
    end type particle_gce
    
    type ( particle_gce ), allocatable, dimension(:,:,:) :: &
        q_gce,  & ! particle mixing ratio [g/m3]
        re_gce    ! particle effective radius [micron]
    type ( particle_gce ), allocatable, dimension(:,:) :: &
        qcol_gce  ! column integrated particle amount [kg/m2]
    type ( particle_gce ) :: & ! scalar parameters
        n0_gce, & ! intercept for exponential PSD [1/m4]
        rho_gce   ! bulk density [kg/m3]

    real(r_size),allocatable,dimension(:) :: hgt_lev, hgt_lay
    real(r_size),allocatable :: topo(:,:)     ! topography
    real(r_size),allocatable :: rho(:,:,:)    ! dry air density [kg/m3]

    ! scaler
    ! Atmosphere parameters at layer level
    type atmos_s_parameter
        real(r_size)  :: hgt      ! heigh [km]
        real(r_size)  :: dhgt     ! thickness of layer [km]
        real(r_size)  :: t_air    ! layer air temperature [degK]
        real(r_size)  :: press    ! layer pressure [hPa]
        real(r_size)  :: rh       ! relative humidity [%]
        real(r_size)  :: omega    ! verticall velocity  [m/s]
        real(r_size)  :: exner    ! exner function [-]
    end type atmos_s_parameter
    type ( atmos_s_parameter ) :: atmos_s
    type ( particle_gce ) :: q_gce_s,re_gce_s

    real(r_size) :: umu       ! Direction cosine (assigned automatically)
    !integer,parameter :: nln = ke   ! num. of total layers 
    integer :: nln    ! num. of total layers 

    integer :: numDims, numAtts
    integer :: ncid,varid
    integer :: start_2d(3),count_2d(3)
    integer :: start_2db(3),count_2db(3)
    integer :: start_3d(4),count_3d(4),count_3dh(4)
    integer :: k1, w1, w2

    ! #############################################################################
    ! ###############  Configure Single-Scattering LUTs Options ###################
    ! #############################################################################
    integer,parameter :: ice_refraction_func = 1
    ! Effective refraction functions for frozen particles for Microwave/Radar simulator (integer)
    ! 1: Oblique Maxwell-Garnett function that assumes ice inclusion within air matrix.
    ! 2: Oblique Maxwell-Garnett function that assumes air inclusion within ice matrix.
    ! 3: Effective-Medium function that assumes homogeneous mixing.

    integer,parameter :: melt_opt = 0   
    ! Effective refraction functions for melting particles for Microwave/Radar simulator (integer)
    ! 0: Does not account melting particle
    ! 1: Oblique Maxwell-Garnett function that assumes ice inclusion within water matrix.
    ! 2: Oblique Maxwell-Garnett function that assumes water inclusion within ice matrix.
    ! 3: Oblique Maxwell-Garnett function averaging option 1 and 2 --> RECOMMENDED
    ! 4: Effective-Medium function that assumes homogeneous mixing.

    ! #############################################################################
    ! ###############  Configure Geostationary Satellite  #########################
    ! #############################################################################
    integer,parameter :: beamcnv_type = 0 
    ! switch for beam convolution
    ! 0: no beam convolution (single point)
    ! 1: beam convolution (rectangular average)
    ! 2: beam convolution (geosat)
    
    real(r_size),parameter :: pi = 4.d0*atan(1.d0)
    real(r_size),parameter :: d2r = pi/180.d0
    real(r_size),parameter :: Rearth = 6371.22E+3   !< radius of earth (m)
    real(r_size),parameter :: sat_lon_e = 135.D0*d2r !< subsatellite longitude (rad)
    real(r_size),parameter :: sat_lat_e = 0.D0*d2r  !< subsatellite latitude (rad)
    real(r_size),parameter :: sat_lev_e = 35786E+3  !< height of satellite (m)

    real(r_size),parameter :: drange = 500.d0    ! range resolution (m)

    !
    !integer,parameter :: ndiv = 2 ! (ndiv*2+1)**2 
    !real(r_size),parameter :: bwidth = 0.032*d2r ! beam width (rad)
    !real(r_size),parameter :: bsigma = (bwidth**2)/(8.*log(2.)) ! sigma**2 (rad^2)
    !real(r_size),parameter :: bound = bwidth ! upper bound to integrate in gpr coordinate (rad)
    !real(r_size),parameter :: dr = bound/ndiv ! -bwidth < theta,phi < bwidth (rad)

    ! fixed dr 
    real(r_size),parameter :: bwidth_1km = 0.0016*d2r ! beam width (rad)
    real(r_size),parameter :: bwidth_3km = 0.0048*d2r ! beam width (rad)
    real(r_size),parameter :: bwidth_5km = 0.008*d2r ! beam width (rad)
    real(r_size),parameter :: bwidth_10km = 0.016*d2r ! beam width (rad)
    real(r_size),parameter :: bwidth_20km = 0.032*d2r ! beam width (rad)
    !real(r_size),parameter :: bwidth = bwidth_10km
    real(r_size),parameter :: bwidth = bwidth_20km
    real(r_size),parameter :: bsigma = (bwidth**2)/(8.*log(2.)) ! sigma**2 (rad^2)
    real(r_size),parameter :: bound = sqrt(bsigma)*2 ! upper bound in a FOV (rad)
    real(r_size),parameter :: dr = bwidth_3km ! increment (rad)
    integer,parameter :: ndiv = bound / dr

    ! #############################################################################
    ! ##########################  Configure Radar Sensor  #########################
    ! #############################################################################
    logical,parameter :: attenuation = .false.
    integer,parameter :: mxfreq_radar = 1         ! # of channels
    real(r_size),parameter :: min_echo = 20.d0    ! minimal detactable echo [dBZ]
    real(r_size),parameter :: freq_radar = 13.8   ! Channel frequencies [GHz]
    real(r_size),parameter :: k2 = 0.925          ! Radar constant |k^2| defaults
    real(r_size),parameter :: view_angle_radar = 12.13    ! viewing angle [deg] 12.13 is derived from mean of 1/mu (0 ~ 17)
    real(r_size),parameter :: lambda = 2.997925/freq_radar/10.d0 ! [m]
    real(r_size),parameter :: radar_fact = k2 * (pi**5) / (lambda**4) ! Z --> sigma

  end module common_geosatpr
