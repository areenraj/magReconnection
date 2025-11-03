module mod_usr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This experiment simulates solar flares in 2.5D.
! The experiment was last developed by Malcolm Keith Druett at KU Leuven, last update 2024-02 (Druett et al. 2023, 2024)
! It is developed from the work of Ruan et al. (2020), which, in turn was developed from modelling by Yokoyama and Shibata (2001)
!
! Reconnection is triggered via the anomalous resistivity model in four phases which are control via the usr_list: 
! t<t_imp         precusor (_pre): A precusor phase - generally with low resititivty.
!                                Can create a set of loops before triggering the flare
! t_imp<t<t_acc   implusive (_imp): Implusive phase - resistivity is increased, triggering reconnection at the chosen location
! t_acc<t<t_decay acceleration (_acc): The energetic electrons are switched on between these times
! t_decay<t       decay phase (_decay): The anomalous resistivity is switched off
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine list/map
! * Root "mod_usr routines" are listed below with capitals
! * Many of these are "pointed to" new names defined with PEARL/Loop Annotation Syntax(LASY) syntax
! * Only the calls to other user defined routines are listed
!
!   USR_INIT
!     Define coord system, usr routine names, physics module, and units
!     | -> usr_params_read  reads user parameters to define
!   USR_SET_PARAMETERS      => INITGLOBALDATA_USR
!     Sets up HS Equilib model + B-field for flare e-, resistivity params, gravity, ghost zones
!     | -> inithdstatic  (see above) exteral HS Equilib model
!        |-> init_refine_bound sets up the refined region for the flare loops
!     | -> init_Bfield   (see above) sets up the arrays and B-field line locations for accelerated particles
!   USR_GRAVITY             => GRAVITY
!     Sets values for gravity in 2D on grid. Calc: in routine getggrav. DOESNT INCLUDE GHOST CELLS.
!     | -> getggrav  (see above)
!   USR_INIT_ONE_GRID       => INITONEGRID_USR
!     Sets up initial atmosphere for experiment
!     | -> specialset_B0 (if B0field is not true (split b-field))
!     | -> mhd_to_conserved() (main subroutine, converts w from primitive to conserved vars)
!   USR_SOURCE              => SPECIAL_SOURCE
!     Implements energetics of special energy/mom/mass terms defined by user in SPECIAL_GLOBAL
!     Here: non-local Energy transport: extracted from acceleration site (ohmic heating term)
!                                       deposited down field-lines via beam method.
!     | -> getbQ : Unphysical general background heating to simulate atmos layers
!     |            returns bQgrid which is added to energy e = e + qdt * bQgrid
!     | -> getcQ : If over beam injection time, calls this to get ohmic heating
!     |            extracted from accleration site, for fast electrons
!     |            e = e - qdt * cQgrid
!     | -> getlQ : If over beam injection time, calls this to get footpoint heating
!                from the fast electrons on the grid, then adds e = e + qdt * lQgrid
!   USR_SPECIAL_BC          => SPECIALBOUND_USR
!     User defined special boundary conditions (applied to y direction only)
!     | -> mhd_to_conserved() for converting to conserved variables after calculation
!     | -> specialset_B0 (but doesnt call as B0field is true (split b-field))
!     | -> mpistop() on error
!   USR_VAR_FOR_ERREST      => P_FOR_ERREST
!     Call pressure calc, store in "var", flag if negative values
!     | -> mhd_get_pthermal (main code routine to calc thermal pressure)
!   USR_REFINE_GRID         => SPECIAL_REFINE_GRID
!     Enforce additional refinement or coarsening
!       max refinement level is enforced on any block that touches the chromosphere (bottom 3 Mm)
!       refinement elsewhere is set to default to coarsen.
!       However, then calls get_refine region user procedure for further implementation
!       later in this experiment the RF region is updated from special_global via calls to update_refine_bound
!     | -> get_refine_region
!   USR_PROCESS_GLOBAL      => SPECIAL_GLOBAL
!     Called at beginning of time steps.
!     Controls beam electrons, acceleration, beam and background heating, flare B-field info, flare file output.
!     Doesnt add these terms to the energy equation. That is done via USR_SOURCE => SPECIAL_SOURCE.
!     | -> special_global
!        | -> update_refine_bound
!        | -> get_flare_eflux
!           |-> update_Bfield   updates B-field positions/strengths before
!           |   |               the fast electron acceleration and energy deposition.
!           |   | -> trace_Bfield
!           |   | -> trace_Bfield
!           |-> locate_midpoint   locates "tops" of fieldlines to use for electron beam
!           |                     plasma depth = 0 reference points.
!           |-> update_heating_table  Creates a table of values for the formula A.12
!           |                         Ruan 2020 in terms of logarithmic plasma depth 
!           |                         scale from Nmin to Nmax. Does not include Flux_0 
!           |                         or B-field strengths which are added for each line:
!           |-> get_heating_rate   Uses the heating table to calculate the heating rate
!           |                      in each line section of each field line.
!           |-> get_Qe             Calculates beam particle heating rate on the grid
!           |-(if convert)
!           |   |-> get_spectra   Calculates (UV line emission spectra?)
!           |   |                 from fast electrons.
!           |   |-> get_HXR_line   Calculates Hard X-ray emission 
!           |   |                  from fieldline locations.
!           |   |-> interp_HXR   Converts fieldline HXR emission to cell grid emission.
!           |-> split_bfield   handles creation/removal of B-field lines
!           |                  if their separation is too high/not high enough in the 
!           |                  region of interest for fast electrons.
!        | -> get_flare_preflux   Outputs field line file before electrons are switched on
!   USR_AUX_OUTPUT          => SPECIALVAR_OUTPUT
!     This subroutine can be used in the "convert" stage of amrvac, 
!     It saves/adds additional (auxiliary) variables to the standard vtk file output
!     |-> mhd_get_pthermal              for T
!     |-> divvector                     for divB and divV
!     |-> get_normalized_divb           for divB
!     |-> get_current, cross_product    for j1234, E123
!     |-> special_eta                   for resistivity eta
!     |-> getlQ, getcQ, getbQ           for electron heating, electron accel, background heating
!     |-> getvar_cooling                for radiative cooling "rad"
!     |-> get_HXR, get_sxr_flare        for HXR and SXR
!     |-> get_refine_region             for RF
!   USR_ADD_AUX_NAMES       => SPECIALVARNAMES_OUTPUT 
!     user-added variable names saved
!   USR_SET_B0              => SPECIALSET_B0
!     B-field split into const background (wB0) + perturbation component
!     Here the time-independent background magnetic field is set
!   USR_SET_J0              => SPECIALSET_J0
!     current density is split into const background (wJ0) + perturbation
!     Here the time-independent background current density is set
!   USR_SPECIAL_RESISTIVITY => SPECIAL_ETA
!     Set the "eta" array
!     Also called from getcQ and specialvar_output
!   USRSPECIAL_CONVERT
!     Some user defined routines to calculate output values from the code
!     Calls these as routines listed after this point
!
  use mod_mhd

  implicit none
  double precision :: q_e, unit_currentdensity
  double precision :: unit_electricfield

!!!!!!!!!!!!!!!!!!! For Atmosphere !!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav
  double precision :: heatunit,gzone,SRadius,bQ0,dya
  double precision, allocatable :: pa(:),ra(:),ya(:),Ta(:)
  integer, parameter :: jmax=20000

  integer :: numxQ^D,numFL,numLP
  double precision :: xQmin^D,xQmax^D,dxQ^D,dFh
  double precision, allocatable :: xQ(:^D&,:),Qe(:^D&)
  double precision, allocatable :: xFLb(:,:),xFRb(:,:)
  double precision, allocatable :: xFL0(:,:),xFR0(:,:)
  integer,allocatable :: numRL(:),numRR(:)
  integer,allocatable :: numTurnL(:),numTurnR(:)
  double precision, allocatable :: EtotL(:),EtotR(:)
  double precision, allocatable :: vs2BL(:),vs2BR(:)
  double precision, allocatable :: vp2L(:),vp2R(:)
  double precision, allocatable :: Bxy0(:)
  double precision :: eta1,eta2,eta3,etam,heta,heta2,tar,vc,v2cr=1.28d7
  integer :: numValidL,numValidR,filenr,iApexL,iApexR
  double precision, allocatable :: mu_electrons(:^D&)
  double precision, allocatable :: electron_flux(:^D&)

  integer :: numN
  double precision :: Nmin,Nmax,dNlog
  double precision, allocatable :: HXR(:^D&)

  double precision :: t_update_Qe=-1.0e-7,dt_update_Qe=5.0e-4,vmax=40.d0
  logical :: readEtot
  logical :: opened = .false.
  double precision :: Eplus,Eminus

  !-------------- for special refinement -------------------------!
  double precision :: ybRFmin,ybRFmax,dLRF
  double precision, allocatable :: xbRFL(:),xbRFR(:),ybRF(:)
  integer :: numybRF
  !-------------- for special refinement -------------------------!
  !
  !
  !-------------- from parfile &usr_list -------------------------!
  ! CONVERSION
  ! field_convert
  !          true: outputs all field line data (line number ~ 1407, 1562)
  !          used in: specialset_B0 (line number ~3454), specialset_J0 (line number ~3487)
  logical :: field_convert=.false.
  ! INITIALISATION
  ! parb: transition distance between bipolar field patches/ current sheet width
  double precision :: parb = 5.d0/3.d0
  ! special resistivity
  ! The flare simulation has 4 phases:
  ! (1) gentle precusor, (2) impuslive, (3) particle acceleration, and (4) gentle decay
  ! Time switches t_{phasename} are used to move between these phases
  ! used in: special_eta (line numbers ~3447), get_heating_rate (line number ~2222)  
  ! PRECURSOR phase variables
  double precision :: eta0_pre=5.d-3, r_eta_pre=0.24d0, h_eta_pre=5.d0
  ! IMPULSIVE phase variables
  double precision :: t_imp=5.98d2/7.7829146115813d1, eta0_imp=3.d-2, r_eta_imp=0.24d0, h_eta_imp=5.d0
  ! PARTICLE ACCELERATION phase variables
  double precision :: t_acc=6.0d2/7.7829146115813d1, alpha_acc=2.d-4, h_eta_acc=5.d0, h_s_acc=1.d0, v_c_acc=1.d3, eta_max_acc=2.d-1
  logical :: eta_decay_acc=.false., beam_deposition=.true.
  ! DECAY phase variables
  double precision :: t_decay=2.0d3/7.7829146115813d1
  ! Coefficient for asymmetric Bfield, Basi = Busr*asy_coeff
  ! Coefficient for where active region, maintan*xmax
  double precision :: asy_coeff=1.0d0, maintan=4.d-1
  logical :: active_region=.false.
  !-------------- from parfile &usr_list -------------------------!


contains


  subroutine usr_init()
    ! Purpose:
    !   Define the coord system, usr routine names, physics module, and units
    !   The standard user routine "named" routines will be called from the main code.
    ! Parameter list:   none, only global parameters
    ! Calls:
    ! | -> set_coordinate_system
    ! | -> mhd_activate    (note: could be hd, 2fluid, magnetofriction _activate)
    ! Called from
    !   Main code (ESSENTIAL mod_usr.t routine for experiment definition)
    
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    call usr_params_read(par_files)

    if(mype .eq. 0) then
      print*, "User parameters: "
      print*, "field_convert       ", field_convert
      print*, "parb                ", parb
      print*, "beam_deposition     ", beam_deposition
      print*, "t_{imp,acc,decay}   ", t_imp, t_acc, t_decay
      print*, "eta0_{pre,imp,acc}  ", eta0_pre, eta0_imp, alpha_acc
      print*, "h_eta_{pre,imp,acc} ", h_eta_pre, h_eta_imp, h_eta_acc
      print*, "r_eta_{pre,imp,acc} ", r_eta_pre, r_eta_imp, h_s_acc
      print*, "v_c_acc             ", v_c_acc
      print*, "eta_decay_acc       ", eta_decay_acc
      print*, "eta_max_acc         ", eta_max_acc
      print*, "asy_coeff           ", asy_coeff
      print*, "active_region       ", active_region
      print*, "maintan             ", maintan
    endif

    ! define "named" routines that will be called by the main code.
    usr_set_parameters      => initglobaldata_usr
    usr_init_one_grid       => initonegrid_usr
    !usr_special_bc          => specialbound_usr
    usr_aux_output          => specialvar_output
    usr_add_aux_names       => specialvarnames_output 
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0
    usr_special_resistivity => special_eta
    usr_var_for_errest      => p_for_errest
    usr_special_convert     => usrspecial_convert
    usr_process_global      => special_global
    !usr_gravity             => gravity
    usr_source              => special_source
    usr_refine_grid         => special_refine_grid
!    usr_internal_bc         => special_boundary_sponge
    usr_set_field_w         => special_field_w

    call mhd_activate()
    ! unit of current density
    unit_currentdensity=unit_magneticfield*const_c/unit_length/4.d0/dpi
    ! unit of electric field
    unit_electricfield=unit_magneticfield*unit_length/(unit_time*const_c)
    ! unit of charge
    q_e=unit_currentdensity/unit_numberdensity/unit_velocity
    if(mype==0) print*,'unit of charge',q_e
    ! dimensionless charge of electron
    q_e=1.60217653d-19/q_e
    if(mype==0) print*,'dimensionless e',q_e

    if (mype==0) then
      print *, 'rho', unit_density
      print *, 'B', unit_magneticfield
      print *, 'p', unit_pressure
      print *, 'v', unit_velocity
      print *, 't', unit_time
      print *, 'E', unit_electricfield
      print *, 'J', unit_currentdensity
      print *, 'eta', unit_electricfield/unit_currentdensity
    endif

!    write(*,*) "end   usr_init",mype
  end subroutine usr_init


  subroutine usr_params_read(files)
    ! Purpose:  Read usr_list parameters from input parameter file
    !           These define timings and parameters of different resistivity schemes
    ! Parameter list:   "files" name of input file for experiment, default = "amrvac.par"
    ! Calls:       None
    ! Called from: USR_INIT
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ field_convert, parb, &
         eta0_pre, r_eta_pre, h_eta_pre, &
         t_imp, eta0_imp, r_eta_imp, h_eta_imp, &
         t_acc, alpha_acc, h_eta_acc, h_s_acc, v_c_acc, eta_decay_acc, eta_max_acc, &
         t_decay,asy_coeff,maintan,active_region,beam_deposition

!    write(*,*) "begin usr_params_read",mype
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

!    write(*,*) "end   usr_params_read",mype
  end subroutine usr_params_read


  subroutine initglobaldata_usr()
    ! Purpose:
    !   Sets up a Hydrostatic Equilib model for calling in usr_init_one_grid
    !   Sets up B-field parameters and resistivity for fast electron transport
    !   These tasks are subcontracted out to user defined routines
    !   Also sets up user gravity, ghost zones,
    !   dya: defines high-res cell-size in y
    ! Calls:
    !   | -> inithdstatic  (see above)
    !   |  | -> init_refine_bound
    !   | -> init_Bfield   (see above)
    ! Called from:
    !   Main code (Substitute for "usr_set_parameters")
    use mod_global_parameters

!    write(*,*) "begin initglobaldata_usr",mype
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s
    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.0d-2/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    Busr=Busr/unit_magneticfield ! magnetic field strength at the bottom
    SRadius=69.61d0 ! Solar radius
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic
    ! for fast electron heating
    call init_Bfield()
!    write(*,*) "end   initglobaldata_usr",mype
  end subroutine initglobaldata_usr


  subroutine inithdstatic
    ! Purpose:
    !   Sets up a Hydrostatic Equilib model for calling in usr_init_one_grid
    !   Based on val-c atmosphere 0Mm up to 2.543 Mm
    !   Based on hydrostatic equilibrium above this
    !   (for my future experiments, we can base on hydrostatic equilbirium below too!)
    ! Calls:
    !   | -> init_refine_bound
    ! Called from:
    !   Main code (Substitute for "usr_set_parameters")
    use mod_global_parameters
  !! initialize the table in a vertical line through the global domain
    use mod_global_parameters

    integer :: j,na,nb,ibc
    double precision, allocatable :: gg(:)
    double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa

    integer :: n_val=49, i
    double precision :: h_val(49),t_val(49)

    double precision :: ai(49),bi(49),ci(49),di(49)
    double precision :: hi

    !rpho=0.71d15/unit_numberdensity ! number density at the bottom relaxla
    rpho=0.9d15/unit_numberdensity ! number density at the bottom relaxla
    Tpho=1.d4/unit_temperature ! temperature of chromosphere
    Ttop=2.d6/unit_temperature ! estimated temperature in the top
    htra=0.2543d0 ! height of initial transition region
    wtra=0.01d0 ! width of initial transition region 
    !Ttr=1.d5/unit_temperature ! lowest temperature of upper profile
    Ttr=4.47d5/unit_temperature ! lowest temperature of upper profile
    Fc=6*2.d5/heatunit/unit_length ! constant thermal conduction flux
    kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3

   !VAL-C
    data    h_val   / 0., 50, 100, 150, 250, &
                      350, 450, 515, 555, 605, &
                      655, 705, 755, 855, 905, &
                      980, 1065, 1180, 1280, 1380, &
                      1515, 1605, 1785, 1925, 1990, &
                      2016, 2050, 2070, 2080, 2090, &
                      2104, 2107, 2109, 2113, 2115, &
                      2120, 2129, 2160, 2200, 2230, &
                      2255, 2263, 2267, 2271, 2274, &
                      2280, 2290, 2298, 2543 /

    data    t_val   / 6420, 5840, 5455, 5180, 4780, &
                      4465, 4220, 4170, 4230, 4420, &
                      4730, 5030, 5280, 5650, 5755, &
                      5925, 6040, 6150, 6220, 6280, &
                      6370, 6440, 6630, 6940, 7160, &
                      7360, 7660, 7940, 8180, 8440, &
                      9500, 10700, 12300, 18500, 21000, &
                      22500, 23000, 23500, 24000, 24200, &
                      24500, 25500, 28000, 32000, 37000, &
                      50000, 89100, 141000, 447000 /

!    write(*,*) "begin inithdstatic",mype
    h_val(1:n_val)=h_val(1:n_val)/1.d4
    t_val(1:n_val)=t_val(1:n_val)/1.d6

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    do j=1,jmax
      ya(j)=(dble(j)-0.5d0)*dya-gzone
      if(ya(j)>=htra) then
        Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
      else
        do i=1,n_val
          if (ya(j)<h_val(i+1)) then
            Ta(j)=t_val(i)+(ya(j)-h_val(i))*(t_val(i+1)-t_val(i))/(h_val(i+1)-h_val(i))
            exit
          endif
        enddo
      endif
      gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
    enddo

    !! solution of hydrostatic equation 
    nb=int(gzone/dya)
    ra(1)=rpho
    pa(1)=rpho*Ta(1)
    invT=gg(1)/Ta(1)
    invT=0.d0
    do j=2,jmax
       invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
       pa(j)=pa(1)*dexp(invT*dya)
       ra(j)=pa(j)/Ta(j)
    end do
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

    call init_refine_bound()

!    write(*,*) "end   inithdstatic",mype
  end subroutine inithdstatic


  subroutine init_refine_bound()
  ! Purpose
  !    Refinements for boundary conditions when the grid is refined nearby
  !    But this is an "initial" rountine...
  ! Calls: none
  ! Called from:
  !    inithdstatic (which is the substitute for "usr_set_parameters")

    integer :: refine_factor,ix2

!    write(*,*) "begin init_refine_bound",mype
    ybRFmin=xprobmin2
    ybRFmax=xprobmax2
    refine_factor=2**(refine_max_level-1)
    dLRF=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor)
    numybRF=domain_nx2*refine_factor

    allocate(xbRFL(numybRF),xbRFR(numybRF),ybRF(numybRF))

    do ix2=1,numybRF
      ybRF(ix2)=ybRFmin+(ix2-0.5d0)*dLRF
      xbRFL(ix2)=-0.025d0
      xbRFR(ix2)=0.025d0
    enddo

!    write(*,*) "end init_refine_bound",mype
  end subroutine init_refine_bound


  subroutine init_Bfield()
    ! Purpose:
    !   initialization of variables for fast electrons and interpolation
    !   includes anomalous resistivity parameters,
    !     fieldline coords, indices, energy quantities, perp and para e- vels 
    ! Parameter list
    !   ggrid: 1D array of gravitational values for cells in y-direction
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   x:     Cell coords
    ! Calls:        none
    ! Called from:  initglobaldata_usr
    use mod_global_parameters

    integer :: ix^D,refine_factor
    double precision :: lengthFL

!    write(*,*) "begin init_Bfield",mype
    tar=0.4d0    ! Time for switching on electrons, t+31.2s in SI units
    !tar=1.0d1/7.7829146115813d1    ! Time for switching on electrons, denominator=unit time
    !eta1=1.d-1
    !vc=1.d3
    !eta2=1.d-3
    !etam=1.d0
    eta1=5.d-2   ! eta_0 in eq 10, Ruan 2020
    eta3=1.d-2   ! What is this? not used!
    vc=1.d3      ! v_c threshold speed for resistivity in eq 11, Ruan 2020
    eta2=1.d-4   ! alpha_eta, resisitivty coefficient in eq 11, Ruan 2020
    etam=1.d-1   ! max eta, from "min[...,0.1]" condition in eq 11, Ruan 2020
    heta=5.d0    ! height of eta central location h_eta in eq 11, Ruan 2020
    heta2=1.d0   ! decay length from central eta location h_s in eq 11, Ruan 2020

    refine_factor=2**(refine_max_level-1)
    dFh=(xprobmax2-xprobmin2)/(domain_nx2*refine_factor) ! dy at max refinement
    dt_update_Qe=max(dFh/vmax,dt)

    ! table for interpolation
    xQmin1=-0.4d0                 ! B-zone: region where B-field lines are tracked (left)
    xQmax1=0.4d0                  ! B-zone: region where B-field lines are tracked (right)
    xQmin2=0.d0                  ! B-zone: region where B-field lines are tracked (bottom)
    xQmax2=0.95d0                 ! B-zone: region where B-field lines are tracked (top, 9.5Mm)
    dxQ1=dFh                     ! dx : copied from max refinement dy
    dxQ2=dFh                     ! dy : copied from max refinement dy
    numXQ1=floor((xQmax1-0.d0)/dxQ1)*2  ! number of fieldline sections initially tracked in x (initialised)
    numXQ2=floor((xQmax2-xQmin2)/dxQ2)  ! number of fieldline sections initially tracked in y direction (init)

    ! max length of a field line for allocating arrays so they dont get overfilled or allocate too much mem
    lengthFL=(xprobmax2-xprobmin2)*2.d0 ! lengthFL : max length of field-line + twice experiment size in y
    numLP=floor(lengthFL/(dFh))         ! numLP:     max number of fieldline pieces
    !numFL=numXQ1                        ! numFL:     max numb of lines tracked (array size), init: numb x steps at max ref in B-zone
    numFL=100

    allocate(xFLb(numFL,ndim),xFRb(numFL,ndim))          ! xFLb, xFRb: location of field line base (left/right)
    allocate(xQ(numXQ1,numXQ2,ndim),Qe(numXQ1,numXQ2))   ! xQ, Qe : coords of field lines, energy into accel
    allocate(numRL(numFL),numRR(numFL))                  ! numRL, numRR : number of points in each fieldline array with a physical value
    allocate(numTurnL(numFL),numTurnR(numFL))            ! numTurnL/R : array locations of field-tops x~0 (s_0, in Ruan 2020)
    allocate(HXR(numXQ1,numXQ2))                         ! HXR : hard xray signals over each line segment
    allocate(EtotL(numFL),EtotR(numFL))                  ! EtotL/R : total energy in e-s for each line segment
    allocate(vs2BL(numFL),vs2BR(numFL))                  ! vs2BL/R : perp mean velocity of electron dist
    allocate(vp2L(numFL),vp2R(numFL))                    ! vp2L/R : parallel mean velocity of electron dist
    allocate(Bxy0(numFL))                                ! Bxy0 : Magnetic field at reference point
    allocate(mu_electrons(numXQ1,numXQ2))                ! mu : pitch-angle along field lines
    allocate(electron_flux(numXQ1,numXQ2))               ! electron_flux: electron flux along field lines

    ! starting points of the field lines
    do ix1=1,numFL
      xFLb(ix1,1)=0.d0-0.5d0*dxQ1-(ix1-0.5)*dxQ1
      xFLb(ix1,2)=0.d0
      xFRb(ix1,1)=0.d0+0.5d0*dxQ1+(ix1)*dxQ1
      xFRb(ix1,2)=0.d0
    enddo

    ! coordinates for the interpolation table
    do ix1=1,numXQ1
      do ix2=1,numXQ2
        xQ(ix1,ix2,1)=0.d0+(ix1-0.5d0*numxQ1-0.5)*dxQ1
        xQ(ix1,ix2,2)=0.d0+(ix2-0.5)*dxQ2
      enddo
    enddo

    ! boundaries of the interpolation table
    xQmin1=xQ(1,1,1)-0.5d0*dxQ1
    xQmax1=xQ(numXQ1,1,1)+0.5d0*dxQ1
    xQmin2=xQ(1,1,2)-0.5d0*dxQ2
    xQmax2=xQ(1,numXQ2,2)+0.5d0*dxQ2

    ! initial fast electron energy and pitch angle
    EtotL=0.d0
    EtotR=0.d0
    Qe=0.d0
    vs2BL=0.d0
    vs2BR=0.d0
    vp2L=0.d0
    vp2R=0.d0
    Bxy0 = 0.d0

    if (iprob>=3 .and. restart_from_file /= undefined)  then
      readEtot=.true.
    else
      readEtot=.false.
    endif

!    write(*,*) "end   init_Bfield",mype
  end subroutine init_Bfield


  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    ! Purpose:
    !   Sets up values for gravity in 2D on grid cell centres
    !   returns them in gravity_field(ix indicies, 2 (x-y-vals)
    !   Calculation subcontracted out to user defined routine getggrav
    ! Parameter list
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   wCT:   Hmm, some kind of variable values... why CT?
    !   x:     Cell coords
    !   gravity_field: output, values of gravity at gridcell cell centres
    !                  DOESNT INCLUDE GHOST CELLS
    ! Calls:
    !   | -> getggrav  (see above)
    ! Called from:
    !   Main code (Substitute for "usr_gravity")
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)
!    write(*,*) "begin gravity",mype

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)
!    write(*,*) "end   gravity",mype
  end subroutine gravity


  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    ! Purpose:
    !   Provides gravity values for cell indicies back through GGRID
    !   Calculates using usr_grav value defined at photosphere, and
    !   quadratic decrease with dist from solar centre, 1D code.
    ! Parameter list
    !   ggrid: 1D array of gravitational values for cells in y-direction
    !   ixi^l: limits (min max) of input cell indices (no ghost cells)
    !   ixo^l: limits (min max) of output cell indices (inc ghost cells)
    !   x:     Cell coords
    ! Calls:         none
    ! Called from:   gravity
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine


  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! Purpose:
    !   Sets up initial atmosphere for experiment
    ! Parameter list
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   w:   variable values
    !   x:     Cell coords
    ! Calls:
    !   | -> specialset_B0 (but doesnt as B0field is true (split b-field)
    !   | -> mhd_to_conserved() (main subroutine, converts variables in w from primitive
    !                           i.e. pressure etc, into only conserved, i.e. into
    !                           the saved values that the code is run based on
    ! Called from:
    !   Main code (Substitute for "usr_init_one_grid")
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: rhoConst, tempConst

    double precision :: res,Bf(ixI^S,1:ndir)
    integer :: ix^D,na
    logical, save :: first=.true.

!    write(*,*) "begin initonegrid_usr",mype

    rhoConst=4.d15/unit_numberdensity 
    tempConst=2.d6/unit_temperature 

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
      endif
      first=.false.
    endif
    ! interpolate from "inithdstatic" model (hydrostatic+val-c) onto the grid values
    ! use this to define denisty rho_, and pressure, p_ in variable array w
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=rhoConst
        w(ix^D,p_)  =rhoConst*tempConst
    {end do\}
    w(ixO^S,mom(:))=zero
    ! splitting b-field or not? (i.e. into main + perturbation part.)
    ! if split B0field is true and the perturbation is set to zero
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else
      call specialset_B0(ixI^L,ixO^L,x,Bf)
      w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
    endif
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

!    write(*,*) "end   initonegrid_usr",mype
  end subroutine initonegrid_usr


  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    ! Purpose:
    !   Controls energetics of special energy/mom/mass defined by user
    !   Here used for energy: non-local transport
    !   energy extracted from acceleration site (ohmic heating term)
    !   energy deposited down field-lines via beam method.
    ! Parameter list:
    !   qdt :  heating "timestep" delta t_Q
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   iw^l:  limits of variable indices
    !   ?qtC : some energy terms? handed to source subroutines
    !   ?wCT : some variable terms handed to source subroutines
    !   qt :   experiment time
    !   w:     variable values
    !   x:     Cell coords
    !   iprob: remember this is a parfile integer set by the user to control the
    !          experiment version, e.g. if higher that 2 will include all the
    !          version 2 and later stuff
    ! Calls:
    !   | -> getbQ : Unphysical general background heating to simulate atmos layers
    !   |            returns bQgrid which is added to energy e = e + qdt * bQgrid
    !   | -> getcQ : If over beam injection time, calls this to get ohmic heating
    !   |            extracted from accleration site, for fast electrons
    !   |            e = e - qdt * cQgrid
    !   | -> getlQ : If over beam injection time, calls this to get footpoint heating
    !                from the fast electrons on the grid, then adds e = e + qdt * lQgrid
    ! Called from:
    !   Main code (Substitute for "usr_source")
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S),cQgrid(ixI^S)

!    write(*,*) "begin special_source",mype
    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

    !! add localized external heating lQ concentrated at feet of loops
    if(iprob >= 3)then
      !energy loss owing to particle acceleration
      cQgrid=0.d0
      if (qt>t_acc) then
        call getcQ(cQgrid,ixI^L,ixO^L,qtC,wCT,x)
        w(ixO^S,e_)=w(ixO^S,e_)-qdt*cQgrid(ixO^S)
      endif

      lQgrid=0.d0
      !heating because of fast electron deposition
      if (qt>t_acc) then
        call getlQ(lQgrid,ixI^L,ixO^L,qtC,wCT,x)
        w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
      endif
    endif

!    write(*,*) "end   special_source",mype
  end subroutine special_source


  subroutine getcQ(cQgrid,ixI^L,ixO^L,qt,w,x)
    ! Purpose:
    !   Calculate energy going into fast particle accleration.
    !     Transfers ohmic losses out from the energy of the cells
    !     and puts them into the particle acceleration model.
    ! Parameter list
    !   cQgrid : array for return of ohmic heating terms that are supplied to fast particles
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   qt : experiment time
    !   w:   variable values
    !   x:     Cell coords
    ! Calls:
    !   | special_source
    !   | -> getcQ
    !   | | -> get_current : main code procedure to calculate current, returns in  
    !   | |                  3D array current[index inner coordiante span, 3]
    !   | | -> special_eta : Special resistitivty user routine (sub for usr_special_resitivity)
    !   | |                 calculates eta values on grid, returns to eta array for cell index
    !   | | -> mhd_get_pthermal : main code uses EOS stated to calculate pressure
    ! Called from:
    !   special_source
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision :: cQgrid(ixI^S)

    integer :: idir,idirmin,ix^D  ! loop veriables
    double precision :: current(ixI^S,3),jpara(ixI^S),j2(ixI^S) ! current, parallel current, current^2 
    double precision :: Btotal(ixI^S,ndir),B2(ixI^S),Alfv(ixI^S),Ma(ixI^S) ! Btotal:B0+B_perturb, B^2, Alfven speed, ?
    double precision :: v2(ixI^S),pth(ixI^S),Te(ixI^S),eta(ixI^S) ! vy, gas pressure, kinetic temp, resititivty
    double precision :: Bn(ixI^S),v(ixI^S,ndir),Efield(ixI^S,ndir),jE3(ixI^S) ! ?unused? like Ma and Alfv

!    write(*,*) "begin getcQ",mype
    ! if after electron start then do, otherwise skip
    ! get total B field instead of splitting into w(mag) and background B0
    ! then calculate...
    ! v2cr + critical velocity to switch on particle acceleration, declared at start of mod_usr.t
    if (qt>t_acc) then
      if(B0field) then
        Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
      else
        Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
      endif
      cQgrid=zero
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      call special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
      Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
      v2(ixO^S)=w(ixO^S,mom(2))/w(ixO^S,rho_)
      j2=zero
      do idir=1,3
        j2(ixO^S)=j2(ixO^S)+current(ixO^S,idir)**2
      enddo
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (Te(ix^D)>5.d0 .and. v2(ix^D)<(-v2cr/unit_velocity)) cQgrid(ix^D)=eta(ix^D)*j2(ix^D)
      {end do\}
    else
      cQgrid=zero
    endif

!    write(*,*) "end   getcQ",mype
  end subroutine getcQ


  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    ! Purpose:
    !   Calculate background heating bQ
    !     As per eq.5 Ruan 2020
    !     H_b = (c_0[y-y_0/h_0]^(-2.7))
    !           -----------------------
    !              exp[h1/(y-y1)]-1
    !     c_0=bQ0=1.0d-2/heatunit as per paper, defined in initonegrid
    !     y_0=0.1=1Mm, h_0=0.5=5Mm, ?(2/7)=(2.7), h_1=0.3=3Mm, y1=0.201~2Mm
    !     Note the typo in Ruan 2020 where it says 2/7, 2.7 is correct
    ! Parameter list
    !   bQgrid : array for return of background heating terms 
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (inc ghost cells)
    !   qt : experiment time
    !   w:   variable values
    !   x:     Cell coords
    ! Calls: None
    !   | special_source
    !   | -> getbQ
    ! Called from:
    !   special_source
    use mod_global_parameters
    use mod_radiative_cooling

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S),ens(ixI^S)
    double precision :: res
    integer :: na,ix^D

!    write(*,*) "begin getbQ",mype
    !bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/6.d0)
    !bQgrid=0.d0   
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if (x(ix^D,2)>0.0) then
        bQgrid(ix^D)=bQ0*(abs(x(ix^D,2)-0.01)/0.05)**(-2.7)
        bQgrid(ix^D)=bQgrid(ix^D)/(exp(0.03/(x(ix^D,2)-0.0201))-1.0)
        if (bQgrid(ix^D)<0.0) bQgrid(ix^D)=0.0
      endif
    {end do\}

!    write(*,*) "end   getbQ",mype
  end subroutine getbQ


  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
    ! Purpose
    !    Calculate fast electron heating rate in the (experiment) cell grid 
    !    (H_e in Ruan 2020) from heating rates on the hyperfine grid 
    !    via weightings decribed in appendix D.
    !    output: lQgrid
    ! Parameter list
    !    lQgrid : fast electron heating rate in the (experiment) cell grid 
    !    ixI, ixO : mesh coordinate indices, input and output?
    !    O includes ghost cells
    !    qt is solar time in experiment
    !    w is array of variables at each index
    !    x: coordinates at each index
    ! Calls:         none
    ! Called from:   get_flare_eflux
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: lQgrid(ixI^S)

    integer :: ix^D,ixO^D,ixQ^L,numQ
    double precision :: xc^L,sumQ

!    write(*,*) "begin getlQ",mype
    ! remember iprob is just a user integer switch
    ! in this case >=3 give full switching on of electrons etc
    select case(iprob)    
    case(3)
      ! if the electrons have started being injected
      if ( (qt>t_acc) .and. (beam_deposition .eqv. .true.) ) then
        lQgrid(ixO^S)=zero
        {do ixO^DB=ixOmin^DB,ixOmax^DB\}
          xcmin^D=x(ixO^DD,^D)-0.5d0*dxlevel(^D);
          xcmax^D=x(ixO^DD,^D)+0.5d0*dxlevel(^D);
          ixQmin^D=floor((xcmin^D-xQmin^D)/dxQ^D)+1;
          ixQmax^D=floor((xcmax^D-xQmin^D)/dxQ^D);
          sumQ=0.d0
          numQ=0
          ! loop over hyperfine grid and add components to total heating rate
          ! for the experiment grid cell in question
          do ix1=ixQmin1,ixQmax1
            do ix2=ixQmin2,ixQmax2
              if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
                sumQ=sumQ+Qe(ix1,ix2)
                numQ=numQ+1
              endif
            enddo
          enddo
          if (numQ>0) lQgrid(ixO^D)=sumQ/numQ
        {enddo\}
      else
        lQgrid=zero
      endif
    end select

!    write(*,*) "end   getlQ",mype
  end subroutine getlQ


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    ! Purpose:
    !   User defined special boundary conditions
    !     Original param files gives values below,
    !     but MKD switched the sides to open to avoid shock reflection at boundary
    !     Boundary rho    momx   momy   momz   energy magx   magy   magz
    !     min1 (x) symm   asymm  symm   asymm  symm   asymm  symm   asymm
    !     max1 (x) symm   asymm  symm   asymm  symm   asymm  symm   asymm
    !     min2 (y) special-----------------------------------------------
    !              Photosphere
    !              rho:       intial
    !              mom->v:    antisymm condition
    !              e->p:      initial pressure
    !              mag:       initial (i.e. mag=0, B0=intial)
    !     max2 (y) special-----------------------------------------------
    !              corona
    !              rho:       same as last vertical cell value
    !              mom:       same as last vertical cell value
    !              e->T:      rho/p of last cell (dT/dy=0), thresholded back if over 5MK
    !              e->P:      const if T<5mk, otherwise thresholded back to 0 in 2 ghost cells
    !              mag1:      antisymmetric Bx (to impose field aligned heat flux perp to boundary)
    !              mag2 mag3: symmetric By and Bz
    ! Parameter list
    !   ixi^l: limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: limits (min max) of outer cell indices (includes ghost cells)
    !   NOTE
    !   ixA coordinates are introduced to handle ghost cells specifically in old code
    !   NOW:
    !   ixi^l: limits (min max) of inner cell indices (only ghost cells)... depending on the case we are in!
    !   ixo^l: limits (min max) of outer cell indices (includes all cells)
    !   qt :   experiment time
    !   iB:    Boundary index iB=[1,2,3,4]->[min1,max1,min2,max2]
    !   w:     variable values
    !   x:     Cell coords
    ! Calls: None
    !   | -> mhd_to_conserved() for converting to conserved variables after calculation
    !   | -> mpistop() on error
    ! Called from:
    !   main code, replacement for usr_special_bc
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L ! indices outer, boundary index, indices inner
    double precision, intent(in) :: qt, x(ixI^S,1:ndim) ! time, coords
    double precision, intent(inout) :: w(ixI^S,1:nw) ! variables

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),Te(ixI^S) ! gas pressure, temp, gravity, kinetic Temp
    double precision :: delydelx, Bf(ixI^S,1:ndir) ! hmm, a gradient?, Bf stores field values called from specialset_B0
    double precision :: dpdh,dh,dp,Tu ! pressure/height grad, height and pressure steps, kinetic temp in coronal (upper) ghost zone
    integer :: ix^D,idir,ixInt^L,ixIntB^L ! cell indices int loop ix1 and ix2, idir=loop for dimensions, B boundary limit var names
    integer :: ixA^L ! cell indices for ghost cells

!    write(*,*) "begin specialbound_usr",mype
    select case(iB)
    case(1)
      ! Here is the magic for transferring to ghost cell references... I guess in old code style from looking at case(3)
      ! first line creates copies of cell limits over whole span of processor block
      ! second line over-writes ixAmin and ixAmax to the ghost cell references
      ! then the ixA values are handed to the processors
      !
      ! actually seems to load in the wrong coordinates! (seems to be for ib=2)
      ixA^L=ixO^L;
      ixAmin1=ixOmax1+1;ixAmax1=ixOmax1+nghostcells;

      ! copy pressure from ghost cell regions to pth
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      !w(ixO^S,rho_)=w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,rho_)
      !w(ixO^S,p_)=pth(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2)

      ! create copy of all velocities into "momentum" array of w
      ! Nonsense BC that overwrite the whole internal domain with the the boundary cell values in v
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmax1+1^%1ixO^S,mom(1))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmax1+1^%1ixO^S,mom(2))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmax1+1^%1ixO^S,mom(3))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmax1+1^%1ixO^S)
      enddo
      ! Nonsense BC that overwrites the whole internal domain with ((-next+3nextnext)/4 in x dir) for  B
      do ix1=ixOmax1,ixOmin1,-1
        w(ix1,ixOmin2:ixOmax2,mag(:))=(1.0d0/3.0d0)* &
                   (-w(ix1+2,ixOmin2:ixOmax2,mag(:)) &
              +4.0d0*w(ix1+1,ixOmin2:ixOmax2,mag(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ! replacement coords for BCs ixA, again seems to be for iB=1 instead
      ! not used, same general nonsense as iB=1 case
      ixA^L=ixO^L;
      ixAmin1=ixOmin1-nghostcells;ixAmax1=ixOmin1-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,rho_)=w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmin1-1^%1ixO^S)
      enddo
      do ix1=ixOmin1,ixOmax1
        w(ix1,ixOmin2:ixOmax2,mag(:))=(1.0d0/3.0d0)* &
                   (-w(ix1-2,ixOmin2:ixOmax2,mag(:)) &
              +4.0d0*w(ix1-1,ixOmin2:ixOmax2,mag(:)))
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      ! Case 3: Photospheric boundary
      !---------------------------------------------------------------------------------------------
      ! fixed zero velocity (Wenzhi, what did you mean by this? isnt it a standard anti-symmetric?)
      ! 
      ! The loop below wouldnt work if things were normal here:
      ! in the second array index  we go from ixomin2:ixomax2 by expanding the w(ixo^s...)
      ! so normally this would be the whole domain.
      ! However, the right side in the second index is from ixOmax2+nghostcells:ixOmax2+1:-1
      ! So a span of nghostcells, stepping backwards.
      !---------------------------------------------------------------------------------------------
      ! Set v(ghostcells)=-v(edge of domain) antisymm reflecting boundary.
      ! Looks like top of domain (coronal boundary), but actually photospheric.
      ! Remember the "ixomax2" here is the last ghost cell below the std domain
      do idir=1,ndir
        w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      ! fixed b1 b2 b3 = B0 values
      if(B0field) then
        ! no pertubration to B-field in ghost cells
        w(ixO^S,mag(:))=0.d0
      else
        ! Or set ghost cell values to specialset_B0
        call specialset_B0(ixI^L,ixO^L,x,Bf)
        w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir)
      endif
      ! rho and gas pressure equal to initial condition values
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ! Case 4: Coronal Boundary
      ! Switch to primitive (rho, v, pressure, temp etc, mag) in ghost and internal cells
      ! ixA ununsed
      ixA^L=ixO^L;
      ixAmin2=ixOmin2-2;ixAmax2=ixOmin2-1;
      call mhd_to_primitive(ixI^L,ixA^L,w,x)
      call mhd_to_primitive(ixI^L,ixO^L,w,x)
      ! Kinetic temp= gas pressure/rho
      Te(ixO^S)=w(ixO^S,p_)/w(ixO^S,rho_)
   
      ! magnetic field
      ! manually set as anitsymmetric in x (to impose field aligned heat flux perp to boundary)
      w(ixOmin2^%2ixO^S,mag(1))=-w(ixOmin2-1^%2ixO^S,mag(1))
      w(ixOmin2+1^%2ixO^S,mag(1))=-w(ixOmin2-2^%2ixO^S,mag(1))
      ! manually set as symmetric in y/z direction
      w(ixOmin2^%2ixO^S,mag(2))=w(ixOmin2-1^%2ixO^S,mag(2))
      w(ixOmin2+1^%2ixO^S,mag(2))=w(ixOmin2-2^%2ixO^S,mag(2))
      w(ixOmin2^%2ixO^S,mag(3))=w(ixOmin2-1^%2ixO^S,mag(3))
      w(ixOmin2+1^%2ixO^S,mag(3))=w(ixOmin2-2^%2ixO^S,mag(3))

      ! momentum
      ! manually set to be constant, equal to final value inside the domain
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))
        w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))
        w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))
      enddo

      ! pressure and density
      ! pressure change is zero if final internal Temp in column is <5MK
      ! else the pressure is thresholded back down to zero in 2 ghost cells
      ! 
      do ix1=ixOmin1,ixOmax1
        dpdh=0.d0
        Tu=w(ix1,ixOmin2-1,p_)/w(ix1,ixOmin2-1,rho_)
        if (Tu>5.d0) then
          dpdh=-w(ix1,ixOmin2-1,p_)/2.d0
        endif
        do ix2=ixOmin2,ixOmax2
          dh=x(ix1,ix2,2)-x(ix1,ix2-1,2)
          dp=dpdh*dh
          w(ix1,ix2,rho_)=w(ix1,ix2-1,rho_)
          w(ix1,ix2,p_)=w(ix1,ix2-1,p_)+dp
        enddo
      enddo

      ! switch back to conserved variables: rho, mom, e, mag
      call mhd_to_conserved(ixI^L,ixA^L,w,x)
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select

!    write(*,*) "end   specialbound_usr",mype
  end subroutine specialbound_usr


  subroutine p_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    ! Purpose:
    !   Call pressure calc, store in "var", flag if negative values
    ! Parameter list
    !   ixi^l:  Limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l:  Limits (min max) of outer cell indices (inc ghost cells)
    !   iflag : Binary flag negative values
    !   w:      Variable values
    !   x:      Cell coords
    !   var:    Store pressure values
    ! Calls:
    !   | -> mhd_get_pthermal (main code routine to calc thermal pressure)
    ! Called from:
    !   Main code as substitute for "usr_var_for_errest"
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,var)
    
  end subroutine p_for_errest


  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    ! Purpose:
    !   Enforce additional refinement or coarsening
    !     iprob is a user defined switch from the parameter input file.
    !     iprob=1 option: old option used for testing refinement setups. Ignore.
    !     
    !     iprob>=2 option: 
    !         If in the bottom 3Mm at the base of the model
    !             And grid refinement is less than the max refinement level,
    !             then increase ref level and dont let it coarsen
    !         If it is already at the max level,
    !             dont allow refining or coarsening.
    !         If not in the bottom 3Mm, then only allow coarsening?!
    !         call "get_refine_region" (which identifies the flare area) to define rfgrid
    !             anything in that get refined too
    ! Parameter list
    !   igrid: main program gives -grid number and current refinement level to inspect
    !   level: main program gives grid number and -current refinement level to inspect
    !   ixi^l: Limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: Limits (min max) of outer cell indices (inc ghost cells)
    !   qt:    Experiment time
    !   w:     Variable values
    !   x:     Cell coords
    ! Calls:
    !   | -> special_refine_grid
    !        | -> get_refine_region
    ! Called from:
    !   Main code (Substitute for "usr_refine_grid")
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    logical :: rfgrid

!    write(*,*) "begin special_refine_grid",mype
    ! fix the bottom layer to the 4 level
    if (iprob==1) then
      if (any(x(ixO^S,2)<=xprobmin2+0.3d0)) then
        if (level<4) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      else
        refine=-1
        coarsen=1
      endif
    endif

    if (iprob>=2) then
!      if ( any( x(ixO^S,2)<=xprobmin2+0.3d0 ) .OR. &
!      any( abs(x(ixO^S,1))>=(xprobmax1-(xprobmax1-xprobmin1)/float(domain_nx1)) ) ) then
      if ( any( x(ixO^S,2)<=xprobmin2+0.3d0 ) ) then
        if (level<refine_max_level) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      else
        refine=-1
        coarsen=1
      endif

      call get_refine_region(ixI^L,ixO^L,w,x,rfgrid)
      if (rfgrid) then
        refine=1
        coarsen=-1
      endif

    endif

!    write(*,*) "end   special_refine_grid",mype
  end subroutine special_refine_grid


  subroutine get_refine_region(ixI^L,ixO^L,w,x,rfgrid)
    ! Purpose:
    !   Enforce additional refinement or coarsening
    !     Checking edges of the "block" being refined.
    !     Requires at least two edges to be inside the B-field tracking region,
    !     if so, then asks for refinement.
    ! Parameter list
    !   ixi^l: Limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: Limits (min max) of outer cell indices (inc ghost cells)
    !   w:     Variable values
    !   x:     Cell coords
    !   rfgrid:Boolean true/false output
    ! Calls: None
    ! Called from:
    !   | -> special_refine_grid
    !        | -> get_refine_region
    !   | -> specialvar_output to give RF map for output inspection
    !   Main code (Substitute for "usr_refine_grid")
    use mod_global_parameters

    integer :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    logical :: rfgrid

    integer :: ix^D,iyb,rfindex ! loop for dimension labels, ? , refinement index

!    write(*,*) "begin get_refine_region",mype
    rfgrid=.false.
    rfindex=0
    ! check whether inside the footpoint of the B-fields being tracked?
    ix2=ixOmin2
    iyb=floor((x(ixOmin1,ix2,2)-ybRFmin)/dLRF+0.5d0)
    ! lower-left corner
    ix1=ixOmin1
    if (iyb>=1 .and. iyb<=numybRF) then
      if (x(ix^D,1)>xbRFL(iyb) .and. x(ix^D,1)<xbRFR(iyb)) then
        rfindex=rfindex+1
      endif
    endif
    ! lower-right corner
    ix1=ixOmax1
    if (iyb>=1 .and. iyb<=numybRF) then
      if (x(ix^D,1)>xbRFL(iyb) .and. x(ix^D,1)<xbRFR(iyb)) then
        rfindex=rfindex+1
      endif
    endif

    ix2=ixOmax2
    iyb=floor((x(ixOmin1,ix2,2)-ybRFmin)/dLRF+0.5d0)
    ! upper-left corner
    ix1=ixOmin1
    if (iyb>=1 .and. iyb<=numybRF) then
      if (x(ix^D,1)>xbRFL(iyb) .and. x(ix^D,1)<xbRFR(iyb)) then
        rfindex=rfindex+1
      endif
    endif
    ! upper-right corner
    ix1=ixOmax1
    if (iyb>=1 .and. iyb<=numybRF) then
      if (x(ix^D,1)>xbRFL(iyb) .and. x(ix^D,1)<xbRFR(iyb)) then
        rfindex=rfindex+1
      endif
    endif

    if (rfindex>=2) rfgrid=.true.

!    write(*,*) "begin end_refine_region",mype
  end subroutine get_refine_region


  subroutine special_global(iit,qt)
    ! Purpose:
    !   Special global is called at the beginning of each time step.
    !     MPI communications must be specified by usr
    !     
    !     Here it is used to control the beam electrons, acceleration, heating
    !	  Also background heating and B-field info in the flare, plus flare file output
    !     (1) Updates refinement boundary: seems to mean region where the flare loops are.
    !           This is so that the refinement processes will later refine any region with flare loops.
    !     (2) If "convert" is specified, then the master process reads out files and broadcasts them.
    !     (3) If restarting experiment, reads flare data an boardcasts it.
    !     (4) Updates flare heating every 10 * max(dFh/vmax, dt) via get_flare_eflux
    !           This is the main task, where the energies are taken and deposited elsewhere,
    !           then the beam statistics are also updated, i.e. the mean pitch angle and flux.
    !     (5) If "convert" then outputs files of flare loop data.
    ! Parameter list
    !   iit:   iteration number
    !   qt:    Experiment time
    ! Calls:
    !   | -> special_global
    !        | -> update_refine_bound
    !        | -> get_flare_eflux
    !           |-> update_Bfield   updates B-field positions/strengths before
    !           |   |               the fast electron acceleration and energy deposition.
    !           |   | -> trace_Bfield
    !           |   | -> trace_Bfield
    !           |-> locate_midpoint   locates "tops" of fieldlines to use for electron beam
    !           |                     plasma depth = 0 reference points.
    !           |-> update_heating_table  Creates a table of values for the formula A.12
    !           |                         Ruan 2020 in terms of logarithmic plasma depth 
    !           |                         scale from Nmin to Nmax. Does not include Flux_0 
    !           |                         or B-field strengths which are added for each line:
    !           |-> get_heating_rate   Uses the heating table to calculate the heating rate
    !           |                      in each line section of each field line.
    !           |-> get_Qe             Calculates beam particle heating rate on the grid
    !           |-(if convert)
    !           |   |-> get_spectra   Calculates (UV line emission spectra?)
    !           |   |                 from fast electrons.
    !           |   |-> get_HXR_line   Calculates Hard X-ray emission 
    !           |   |                  from fieldline locations.
    !           |   |-> interp_HXR   Converts fieldline HXR emission to cell grid emission.
    !           |-> split_bfield   handles creation/removal of B-field lines
    !           |                  if their separation is too high/not high enough in the 
    !           |                  region of interest for fast electrons.
    ! Called from:
    !   Main code (Substitute for "usr_process_global")
    use mod_global_parameters

    integer, intent(in) :: iit
    double precision, intent(in) :: qt
    character(50) :: fname,tempc
    integer :: ixFL,iyLH,ix1,ix2
    double precision :: t1,t2,t_output
    integer :: ifile,tempi

!    write(*,*) "begin special_global",mype
    ! main processor to print out message when particle beams are activated
    if (qt>=t_acc .and. qt-dt<t_acc) then
      if (mype==0) then
        print *, '!------------------------------------------------!'
        print *, '!---------- fast electron activatived -----------!'
        print *, '!------------------------------------------------!'
      endif
    endif

    !---------- special refinement---------------------------------!
    call update_refine_bound()
    !---------- special refinement---------------------------------!

    ! for converting data, reading fast electron info
    if (iprob>=3 .and. convert .and. qt>t_acc) then
      dt_update_Qe=max(dFh/vmax,dt)

      ! read Etot
      if (mype==0) then
        filenr=snapshotini
        write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
        open(1,file=fname,action='READ')
        read(1,*) tempc
        read(1,*) tempi
        read(1,*) tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc
        do ix1=1,numFL
          read(1,*) tempi,xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                    vs2BL(ix1),vs2BR(ix1),vp2L(ix1),vp2R(ix1),Bxy0(ix1)
        enddo
        close(1)
      endif
      call MPI_BCAST(xFLb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(xFRb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vs2BL,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vs2BR,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vp2L,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vp2R,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(Bxy0,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

      !call get_flare_eflux()
    endif

    ! for restarting, reading fast electron info
    if (readEtot)  then
      dt_update_Qe=max(dFh/vmax,dt)
      readEtot=.false.

      ! read Etot
      if (mype==0) then
        filenr=snapshotini
        write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
        open(1,file=fname,action='READ')
        read(1,*) tempc
        read(1,*) tempi
        read(1,*) tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc,tempc
        do ix1=1,numFL
          read(1,*) tempi,xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                    vs2BL(ix1),vs2BR(ix1),vp2L(ix1),vp2R(ix1),Bxy0(ix1)
        enddo
        close(1)
      endif
      call MPI_BCAST(xFLb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(xFRb(:,1),numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotL    ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(EtotR    ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vs2BL    ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vs2BR    ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vp2L     ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(vp2R     ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(Bxy0     ,numFL,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

      call MPI_BCAST(filenr,1,MPI_INTEGER,0,icomm,ierrmpi)
      filenr=filenr+1

    endif

    ! calculate heating
    if (iprob>=3 .and. qt>t_acc) then
      dt_update_Qe=max(dFh/vmax,dt)
      if (qt>t_update_Qe+10*dt_update_Qe) then
        t_update_Qe=qt-1.0e-7
      endif
      if (qt>t_update_Qe) then
        t1=mpi_wtime()
        call get_flare_eflux()
        t_update_Qe=t_update_Qe+dt_update_Qe
        t2=mpi_wtime()
        !if (mype==0) print *, iit,t2-t1
      endif
    endif

    ! output fast electron info for restarting and converting data
    if (convert .eqv. .false.) then 
      if (iit==0) then
        filenr=0
      endif
      t_output=filenr*dtsave(2)
      if (qt>=t_output .and. qt<t_output+dt) then     
        if (mype==0) then
          write(fname, '(a,i4.4,a)') trim(base_filename),filenr,'_Etot.txt'
          open(1,file=fname)
          write(1,*) 'numFL'
          write(1,*) numFL
          write(1,*) 'iFL xFL(1) xFR(1) EtotL EtotR vs2BL vs2BR vp2L vp2R Bxy0'
          do ix1=1,numFL
            write(1,'(i8, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7, e15.7)') ix1,&
                  xFLb(ix1,1),xFRb(ix1,1),EtotL(ix1),EtotR(ix1),&
                  vs2BL(ix1),vs2BR(ix1),vp2L(ix1),vp2R(ix1),Bxy0(ix1)
          enddo
          close(1)
        endif
        filenr=filenr+1
      endif
    endif

    ! output field data for case that the electrons havent started yet
    if (iprob>=3 .and. convert .and. qt .LE. t_acc) then
        call get_flare_preflux()
    endif

!    write(*,*) "end   special_global",mype
  end subroutine special_global


  subroutine update_refine_bound()
  ! Purpose
  !   * calculate xbRF{L/R} : (x B-field refine left/right)
  !   * Refinements of region in which B-field lines are tracked
  !     Region is defined by 2 vertical stacks of x values stored in arrays xbRFL, xbFLR
  !   * This is done by choosing starts point on the left/right side of the recon region
  !     e.g. xfs(1,1:2) (-2.5Mm,50Mm) for the left, (2.5Mm,50Mm) for the right
  !     then trace the field "backward" and "forward" from this point, down the loop system.
  !   * We will have max refinement inside these, therefore from the max refinement stepsize
  !     we can calculate gridcell y ordinates
  !     therefore, these new values can be compared with the stored boundaries in xbRF{L,R}
  ! Calls:
  !    trace_field_single (from module mod_trace_field included in subroutine)
  ! Called from:
  !    special_global (which is the substitute for "usr_process_global")
    use mod_trace_field

    integer :: ix2, ip
    integer :: npmax=5000,numRTf=0
    !double precision :: xfs(5000,1:ndim),wfnum(5000,1),wf(5000,1:nw+ndir)
    double precision :: xfs(5000,1:ndim),wPnull(5000,1), wLnull(1)
    logical :: forwardf,wRT(nw+ndir)
    character(len=std_len) :: fieldtype, tcondit

!    write(*,*) "begin update_refine_bound",mype
    xbRFL=-0.025d0 ! set default boundaries at +-2.5Mm
    xbRFR=0.025d0

    ! Trace from left side of recon region first
    ! backwards first
    numRTf=0
    forwardf=.false.
    xfs=0.d0
    xfs(1,1)=-0.025d0
    xfs(1,2)=0.5d0
    ! numRTF 1 element matrix for number of valid field points. 1 element as trace single. initialised to 0 
    fieldtype="Bfield"
    tcondit="none"
    !call trace_field_single(xf,wPm   ,wL    ,dL  ,numP ,nwP,nwL,forwardf,fieldtype,tcondit)
    call trace_field_single(xfs,wPnull,wLnull,dLRF,npmax,1,0,forwardf,fieldtype,tcondit)
    numRTf=int(wLnull(1))
    do ip=1,numRTf
      ix2=floor((xfs(ip,2)-ybRFmin)/dLRF+0.5d0)
      ix2=max(1,ix2)
      ix2=min(numybRF,ix2)
      if (xbRFL(ix2)>xfs(ip,1)) xbRFL(ix2)=xfs(ip,1)
    enddo
    ! now forwards
    numRTf=0
    wLnull(1)=0.d0
    forwardf=.true.
    xfs=0.d0
    xfs(1,1)=-0.025d0
    xfs(1,2)=0.5d0
    call trace_field_single(xfs,wPnull,wLnull,dLRF,npmax,1,0,forwardf,fieldtype,tcondit)
    numRTf=int(wLnull(1))
    do ip=1,numRTf
      ix2=floor((xfs(ip,2)-ybRFmin)/dLRF+0.5d0)
      ix2=max(1,ix2)
      ix2=min(numybRF,ix2)
      if (xbRFL(ix2)>xfs(ip,1)) xbRFL(ix2)=xfs(ip,1)
    enddo

    ! Now trace from right side of recon region    
    ! backwards first
    forwardf=.false.
    xfs=0.d0
    xfs(1,1)=0.025d0
    xfs(1,2)=0.5d0
    call trace_field_single(xfs,wPnull,wLnull,dLRF,npmax,1,0,forwardf,fieldtype,tcondit)
    numRTf=int(wLnull(1))
    do ip=1,numRTf
      ix2=floor((xfs(ip,2)-ybRFmin)/dLRF+0.5d0)
      ix2=max(1,ix2)
      ix2=min(numybRF,ix2)
      if (xbRFR(ix2)<xfs(ip,1)) xbRFR(ix2)=xfs(ip,1)
    enddo
    ! now forwards    
    forwardf=.true.
    xfs=0.d0
    xfs(1,1)=0.025d0
    xfs(1,2)=0.5d0
    call trace_field_single(xfs,wPnull,wLnull,dLRF,npmax,1,0,forwardf,fieldtype,tcondit)
    numRTf=int(wLnull(1))
    do ip=1,numRTf
      ix2=floor((xfs(ip,2)-ybRFmin)/dLRF+0.5d0)
      ix2=max(1,ix2)
      ix2=min(numybRF,ix2)
      if (xbRFR(ix2)<xfs(ip,1)) xbRFR(ix2)=xfs(ip,1)
    enddo
!    write(*,*) "end   update_refine_bound",mype
  end subroutine update_refine_bound


  subroutine get_flare_eflux()
    ! Purpose:
    !    Calculates energy flux deposited in the atmosphere due to fast electrons
    !
    ! Paramter list:
    !    None, only uses global parameters
    !
    ! Calls:
    ! |-> update_Bfield   updates B-field positions/strengths before
    ! |   |               the fast electron acceleration and energy deposition.
    ! |   | -> trace_Bfield
    ! |   | -> trace_Bfield
    ! |-> locate_midpoint   locates "tops" of fieldlines to use for electron beam
    ! |                     plasma depth = 0 reference points.
    ! |-> update_heating_table  Creates a table of values for the formula A.12
    ! |                         Ruan 2020 in terms of logarithmic plasma depth 
    ! |                         scale from Nmin to Nmax. Does not include Flux_0 
    ! |                         or B-field strengths which are added for each line:
    ! |-> get_heating_rate   Uses the heating table to calculate the heating rate
    ! |                      in each line section of each field line.
    ! |-> get_Qe             
    ! |-(if convert)
    ! |   |-> get_spectra   Calculates (UV line emission spectra?)
    ! |   |                 from fast electrons.
    ! |   |-> get_HXR_line   Calculates Hard X-ray emission 
    ! |   |                  from fieldline locations.
    ! |   |-> interp_HXR   Converts fieldline HXR emission to cell grid emission.
    ! |   |-> for 
    ! |-> split_bfield   handles creation/removal of B-field lines
    ! |                  if their separation is too high/not high enough in the 
    ! |                  region of interest for fast electrons.
    ! 
    ! Called from
    !    special_global
    ! calculate electron deposit energy flux
    use mod_global_parameters

    integer :: ix^D,j ! integer loop param for each ^d (dimension), j?
    ! xFL/R coords of loops, QeL/R heating-> acceleration components on each 1D loop
    ! wBL/R interpolated w(variable) values on 1D-B-loops
    ! eFluxL/R total energy for particles in each 1D loop
    ! HXRL/R hard x-ray signitures (1D loops)
    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: QeL(numFL,numLP),QeR(numFL,numLP),QaccL(numFL,numLP),QaccR(numFL,numLP)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: eFluxL(numFL),eFluxR(numFL)
    double precision :: HXRL(numFL,numLP),HXRR(numFL,numLP)
    double precision :: muL(numFL,numLP),muR(numFL,numLP)
    double precision :: eFluxLP(numFL,numLP),eFluxRP(numFL,numLP)

    double precision, allocatable :: QeN(:) ! heating rate as func of plasma depth

    double precision :: delta,Ec,dEe,keV_erg ! beam power law index, cutoff, energy steps for HXR
    integer :: numEe,iEe ! number of energy bins, looper
    double precision, allocatable :: Ee(:),spectra(:,:) ! HXR spectral bins
    character(50) :: fname
    character(50) :: fname_field
    integer :: idim, filenr_field, itemp

    !write(*,*) "begin get_flare_eflux",mype
    ! L indicates that we the starting point is located at where x<0
    ! R indicates that we the starting point is located at where x>0
    ! initial fast electron spectral parameters
    delta=4.0
    Ec=20.0  ! cutoff energy [keV]
    dEe=1.0  
    numEe=280
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg

    xFL=0
    xFR=0
    wBL=0
    wBR=0
    QeR=0.d0

    ! set footpoints of field lines
    xFL(:,1,1)=xFLb(:,1)
    xFL(:,1,2)=xFLb(:,2)
    xFR(:,1,1)=xFRb(:,1)
    xFR(:,1,2)=xFRb(:,2)

    ! trace Bfield to prepare for heating
    !write(*,*) "begin update_bfield",mype
    !if (mype==0) then
    !        print*, 'Before update_Bfield: xFL, xFR, wBL, wBR = ', xFL(180,1,1), xFR(180,1,1), wBL(180,1,rho_), wBR(180,1,rho_)
    !endif
    call update_Bfield(xFL,xFR,wBL,wBR)
    !if (mype==0) then
    !        print*, 'After update_Bfield: xFL, xFR, wBL, wBR = ', xFL(180,1,1), xFR(180,1,1), wBL(180,1,rho_), wBR(180,1,rho_)
    !endif
    !if (mype==0) then
    !    write(*,*) ("update_bfield printout")
    !    do itemp = 1,nw+ndir
    !        write(*,*) "i=  ",itemp," ",maxval(wBL(:,:,itemp))
    !    enddo
    !endif
    !write(*,*) "end update_bfield",mype

    ! find the s_0 starting reference point for plas column depth, as close to x=0 symmmetry line as possible
    !write(*,*) "begin locate_midpoint",mype
    call locate_midpoint(xFL,xFR,wBR,wBL)

    ! update heating table QeN as func of column depth "N" using densities of traced loops
    ! why does this need to be done every loop?!
    numN=1001
    allocate(QeN(numN))
    Nmin=1.0e18
    Nmax=1.0e23
    dNlog=log10(Nmax/Nmin)/(numN-1.0)
    !write(*,*) "begin update_heating_table",mype
    call update_heating_table(QeN,delta,Ec)
    !if (mype==0) then
    !    write(*,*) ("update_heating_table printout")
    !    write(*,*) "QeN  ",itemp," ",maxval(QeN)
    !endif

    Eplus=0.d0
    Eminus=0.d0 
 
    ! Qet heating rate for each field line Qe"L" and Qe"R", using Qe"N" and knowledge of how N varies
    ! along the loops just traced
    !write(*,*) "begin get_heating_rate1",mype
    !if (mype==0) then
    !    write(*,*) ("get_heat_rate pre-printout")
    !    write(*,*) "xFL=  ",maxval(xFL(:,2:,:))
    !    write(*,*) "wBL=  ",maxval(wBL(:,:,:))
    !    write(*,*) "QeL=  ",maxval(QeL(:,:))
    !    write(*,*) "QeN=  ",maxval(QeN)
    !    write(*,*) "QaccL=  ",maxval(QaccL(:,:))
    !    write(*,*) "numRL=  ",maxval(numRL)
    !    write(*,*) "numValidL=  ",numValidL
    !    write(*,*) "iApexL=  ",iApexL
    !    write(*,*) "EtotL=  ",maxval(EtotL)
    !    write(*,*) "EfluxL=  ",maxval(EfluxL)
    !endif
    call get_heating_rate(xFL,wBL,QeL,QeN,QaccL,numRL,numTurnL,numValidL,iApexL,EtotL,vs2BL,vp2L,eFluxL)
    !write(*,*) "begin get_heating_rate1",mype
    call get_heating_rate(xFR,wBR,QeR,QeN,QaccR,numRR,numTurnR,numValidR,iApexR,EtotR,vs2BR,vp2R,eFluxR)
    ! total number flux of fast electrons
    eFluxL=eFluxL*heatunit*unit_length*(delta-2)/(Ec*keV_erg)**2
    eFluxR=eFluxR*heatunit*unit_length*(delta-2)/(Ec*keV_erg)**2
    !if (mype==0) then
    !    write(*,*) ("get_heat_rate post-printout")
    !    write(*,*) "xFL=  ",minval(xFL(:,2:,:)),maxval(xFL(:,2:,:))
    !    write(*,*) "wBL=  ",maxval(wBL(:,:,:))
    !    write(*,*) "QeL=  ",minval(QeL(:,:)),maxval(QeL(:,:))
    !    write(*,*) "QeN=  ",maxval(QeN)
    !    write(*,*) "QaccL=  ",maxval(QaccL(:,:))
    !    write(*,*) "numRL=  ",maxval(numRL)
    !    write(*,*) "numValidL=  ",numValidL
    !    write(*,*) "iApexL=  ",iApexL
    !    write(*,*) "EtotL=  ",maxval(EtotL)
    !    write(*,*) "EfluxL=  ",maxval(EfluxL)
    !endif
    
    ! get local "grid" heating rate via interpolation
    !write(*,*) "begin get_Qe",mype
    call get_Qe(xFL,xFR,wBL,wBR,QeL,QeR)


    ! for calculating HXR flux
    allocate(Ee(numEe),spectra(numN,numEe))
    if (convert) then
      ! get fast electron spectra table
      call get_spectra(delta,Ec,dEe,Ee,spectra,numEe)
  
      ! get HXR flux for each field line
      call get_HXR_line(Ee,spectra,numEe,dEe,xFL,wBL,numValidL,numTurnL,numRL,eFluxL,vs2BL,vp2L,HXRL,muL,eFluxLP)
      call get_HXR_line(Ee,spectra,numEe,dEe,xFR,wBR,numValidR,numTurnR,numRR,eFluxR,vs2BR,vp2R,HXRR,muR,eFluxRP)
  
      ! get local HXR flux via interpolation
      call interp_HXR(xFL,xFR,wBL,wBR,HXRL,HXRR,muL,muR,eFluxLP,eFluxRP)
      !HXR=0.d0
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! HEREHEREHERE Malcolm additional output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! issues: huge as a text file, save a header file and an dat file with binary data
      !         only doing half (L vs R)
      if (field_convert) then
        if (mype==0) write(*,*) "field convert", field_convert
        filenr_field=filenr-1
        if (filenr_field .LT. 0) filenr_field=0
        if (mype==0) then
          !write(*, '(a,i4.4,a)') trim(base_filename), filenr_field, '_fieldvalues.txt'
          write(fname_field, '(a,i4.4,a)') trim(base_filename),filenr_field,'_fieldvalues.txt'
          !write(*,*) fname_field
          open(2,file=fname_field)
          write(2,*) '# numFL, numLP'
          write(2,'(2i6.5)') numFL, numLP
          do ix1=1,numFL
            write(2,*) '# iFL efluxL efluxR'
            write(2,'(i8.4, ES17.7E3, ES17.7E3)') ix1, eFluxL(ix1), eFluxR(ix1)
            write(2,*) '# coordinates L'
            if (xFL(ix1,2,2) .LE. xFL(ix1,1,2)) then
               write(2,'(ES17.7E3)') xFL(ix1,1,1)
            else
               do idim=1,ndim
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") xFL(ix1,ix2,idim)
                  enddo
                  write(2,*)
               enddo
               do idim=1,nw+ndir
                  write(2,*) '# variable',idim
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") wBL(ix1,ix2,idim)
                  enddo
                  write(2,*)
               enddo
               if (eFluxL(ix1) .GT. 0.0) then
                  write(2,*) '# Variables: QeL, if the totals are above 0'
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") QeL(ix1,ix2)
                  enddo
                  write(2,*)
                  write(2,*) '# Variables: QaccL, if the total is above 0'
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") QaccL(ix1,ix2)
                  enddo
                  write(2,*)
               endif
            endif
            write(2,*) '# coordinates R'
            if (xFR(ix1,2,2) .LE. xFR(ix1,1,2)) then
               write(2,'(ES17.7E3)') xFR(ix1,1,1)
            else
               do idim=1,ndim
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") xFR(ix1,ix2,idim)
                  enddo
                  write(2,*)
               enddo
               do idim=1,nw+ndir
                  write(2,*) '# variable',idim
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") wBR(ix1,ix2,idim)
                  enddo
                  write(2,*)
               enddo
               if (eFluxR(ix1) .GT. 0.0) then
                  write(2,*) '# Variables: QeR, if the totals are above 0'
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") QeR(ix1,ix2)
                  enddo
                  write(2,*)
                  write(2,*) '# Variables: Qaccr, if the total is above 0'
                  do ix2=1,numLP
                     write(2,'(ES17.7E3)', advance="no") QaccR(ix1,ix2)
                  enddo
                  write(2,*)
               endif
            endif

          enddo
          close(2)
        endif
        filenr_field=filenr_field+1
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ENDENDENDENDEND Malcolm additional output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else 
      HXR=0.d0
    endif

    ! Are field lines too far apart at some point, or too close together?
    ! split_Bfield adds/removes lines as needed. 
    call split_Bfield(xFL,xFLb,wBL,EtotL,vs2BL,vp2L,numRL,numValidL)
    call split_Bfield(xFR,xFRb,wBR,EtotR,vs2BR,vp2R,numRR,numValidR)

    deallocate(Ee,spectra)
    ! This is why we need to run QeN each time...! does it not get saved?
    deallocate(QeN)

!    write(*,*) "end   get_flare_eflux",mype
  end subroutine get_flare_eflux


  subroutine get_flare_preflux()
    ! Purpose:
    !    Calculates field values along fieldlines when electron injection is turned off
    ! Paramter list:
    !    None, only uses global parameters
    ! Calls:
    ! |-> update_Bfield   updates B-field positions/strengths before
    ! |   |               the fast electron acceleration and energy deposition.
    ! |   | -> trace_field_multi
    ! |   | -> trace_field_multi
    ! |-> locate_midpoint   locates "tops" of fieldlines to use for electron beam
    ! |                     plasma depth = 0 reference points.
    ! Called from
    !    special_global
    use mod_global_parameters

    integer :: ix^D,j ! integer loop param for each ^d (dimension), j?
    ! xFL/R coords of loops, QeL/R heating-> acceleration components on each 1D loop
    ! wBL/R interpolated w(variable) values on 1D-B-loops
    ! eFluxL/R total energy for particles in each 1D loop
    ! HXRL/R hard x-ray signitures (1D loops)
    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)

    double precision, allocatable :: QeN(:) ! heating rate as func of plasma depth

    character(50) :: fname
    character(50) :: fname_field
    character(4) :: varlist(11)
    integer :: idim, filenr_field

!    write(*,*) "begin get_flare_preflux",mype
    varlist=['rho ','mom1','mom2','mom3','e   ','B1  ','B2  ','B3  ','J1  ','J2  ','J3  ']
    xFL=0
    xFR=0
    wBL=0
    wBR=0
    ! set footpoints of field lines
    xFL(:,1,1)=xFLb(:,1)
    xFL(:,1,2)=xFLb(:,2)
    xFR(:,1,1)=xFRb(:,1)
    xFR(:,1,2)=xFRb(:,2)
    ! trace Bfield to prapare for heating
    call update_Bfield(xFL,xFR,wBL,wBR)
    ! find the s_0 starting reference point for plas column depth, as close to x=0 symmmetry line as possible
    call locate_midpoint(xFL,xFR,wBL,wBR)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HEREHEREHERE Malcolm additional output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! issues: huge as a text file, save a header file and an dat file with binary data
    !         only doing half (L vs R)
    if (field_convert) then
       filenr_field=filenr-1
       if (filenr_field .LT. 0) filenr_field=0
       if (mype==0) then
          !write(*, '(a,i4.4,a)') trim(base_filename), filenr_field, '_fieldvalues.txt'
          write(fname_field, '(a,i4.4,a)') trim(base_filename),filenr_field,'_fieldvalues.txt'
          !write(*,*) fname_field
          open(2,file=fname_field)
          write(2,*) '# numFL, numLP'
          write(2,'(2i6.5)') numFL, numLP
          do ix1=1,numFL
             write(2,*) '# iFL efluxL efluxR'
             write(2,'(i8.4, ES17.7E3, ES17.7E3)') ix1, 0.0, 0.0
             write(2,*) '# coordinates L'
             if (xFL(ix1,2,2) .LE. xFL(ix1,1,2)) then
                write(2,'(ES17.7E3)') xFL(ix1,1,1)
             else
                do idim=1,ndim
                   do ix2=1,numLP
                      write(2,'(ES17.7E3)', advance="no") xFL(ix1,ix2,idim)
                   enddo
                   write(2,*)
                enddo
                do idim=1,nw+ndir
                   write(2,*) '# variable ',idim,' ',trim(varlist(idim))
                   do ix2=1,numLP
                      write(2,'(ES17.7E3)', advance="no") wBL(ix1,ix2,idim)
                   enddo
                   write(2,*)
                enddo
             endif
             write(2,*) '# coordinates R'
             if (xFR(ix1,2,2) .LE. xFR(ix1,1,2)) then
                write(2,'(ES17.7E3)') xFR(ix1,1,1)
             else
                do idim=1,ndim
                   do ix2=1,numLP
                      write(2,'(ES17.7E3)', advance="no") xFR(ix1,ix2,idim)
                   enddo
                   write(2,*)
                enddo
                do idim=1,nw+ndir
                   write(2,*) '# variable ',idim,' ',trim(varlist(idim))
                   do ix2=1,numLP
                      write(2,'(ES17.7E3)', advance="no") wBR(ix1,ix2,idim)
                   enddo
                   write(2,*)
                enddo
             endif

          enddo
          close(2)
       endif
       filenr_field=filenr_field+1
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ENDENDENDENDEND Malcolm additional output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    write(*,*) "end   get_flare_preflux",mype
  end subroutine get_flare_preflux


  subroutine update_Bfield(xFL,xFR,wBL,wBR)
    ! Purpose
    !    Trace B field and calculate column depths along each field line.
    !    Check whether any field lines are outisde the region of interest, 
    !    and remove any that are from the list of those being traced.
    !    Tracing the B is handled explicitly by "trace_Bfield",
    !    which also returns interpolated plasma parameters along the field line.
    !    including plasma depth, it would appear.
    ! 
    ! Parameter list
    !    xFL(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (left side)
    !    xFR(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (right side)
    !    wBL(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !    wBR(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !       L and R indicate left and right side of experiment again.
    !
    ! Calls
    ! | -> trace_Bfield
    ! | -> trace_Bfield
    ! 
    ! Called from
    !    get_flare_eflux
    use mod_global_parameters
    use mod_usr_methods
    use mod_trace_field
    use mod_point_searching
    !use mod_trace_Bfield

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    ! call variables
    character(len=std_len) :: ftype, tcondi
    double precision :: xF(numFL,numLP,ndim),wLm(numFL),wPm(numFL,numLP,nw+ndir) 
    double precision :: xpp(ndim),wpp(nw)
    integer :: numP, nL, nwP, nwL, iL, imype
    logical :: forwardm(numFL)
    
    character(50) :: fname

    integer :: ix^D,j
    logical :: interp
    double precision :: xmin
    !logical :: wTrans(nw+ndir)
    double precision :: t0,t1
    !double precision :: xF(numLP,ndim),wB(numLP,nw+ndir)
    !integer :: numRT

    ! numValid -- how many field lines are traced
    ! numR -- how many points in the field line are valid
    
!    write(*,*) "begin update_Bfield",mype
    t0=MPI_wtime()

    numValidL=numFL ! set initial number of field lines to track = numFL (=max,size of array) (shrunk in TRACE1 loop, if possible)
    numValidR=numFL

    ! for starting points at x<0
    xF=0.d0
    xF(1:numValidL,1,:)=xFL(1:numValidL,1,:)        ! xF(numFL,numLP,ndim): field coords
    ftype='Bfield'                                  ! ftype: field type
    tcondi='cons+j'                                 ! tcondi: user branch condition (We want conserved variables + current)
    do iL=1,numValidL                               ! forwardm: Forwards or backwards field tracing
      xpp(1:ndim)=xF(iL,1,1:ndim)
      !write(*,*) "start get_point_w",mype
      call get_point_w(xpp,wpp,'conserved')
      !write(*,*) "end get_point_w",mype
      if (wpp(mag(2))>0.d0) then
        forwardm(iL)=.true.
      else
        forwardm(iL)=.false.
      endif
    enddo
    wPm=0.d0    ! wPm: the variables defined on the fieldline
                ! This should be calculated using the helper function special_field_s and the branching condition "tcondi"
    wLm=0.d0    ! wLm(numFL): This is the variables produced as a single value from a fieldline,
                ! In this case, just 1 variable, the number of valid fieldline points, set all to zero initially
    ! dFh is from the main code from global parameters, fieldline step size= max refinement step
    ! numFL Is the max number of lines that can be traced (array size)
    ! numLP: integer, number of points on the fieldline (max from the array)
    nwP=nw+ndir ! nwP: integer number of "per point" variables, here the conserved +the current
    nwL=0       ! nwL: integer number of "per field line" variables (excluding valid points along it)
    !write(*,*) "start trace_field_multi",mype
    !if (mype==0) then
    !   !tempfloat=maxval(wLm)
    !   write(*,*) "numFL",numFL
    !   write(*,*) "numP",numLP
    !   write(*,*) "nwL",nwL
    !   write(*,*) "nwP",nwP    
    !endif
    !if (mype==0) then
    !        print*, 'LEFT - Before trace_field_multi: xF, wBRL = ', xF(180,1,1), wPm(180,1,rho_)
    !endif
    call trace_field_multi(xF,wPm,wLm,dFh,numFL,numLP,nwP,nwL,forwardm,ftype,tcondi)
    !if (mype==0) then
    !        print*, 'LEFT - After trace_field_multi: xF, wBRL = ', xF(180,1,1), wPm(180,1,rho_)
    !endif
    !if (mype==0) then
    !   !tempfloat=shape(wLm)
    !   !write(*,*) "trace_field call, shape(wLm)=",shape(wLm)
    !   write(*,*) maxval(wLm)
    !endif
    !write(*,*) "end trace_field_multi",mype
    ! shrink number of field lines to trace if possible
    numRL=int(wLm)
    TRACE1: do ix1=1,numFL
      !if field line is outside of the reconnection region, stop tracing
      xmin=abs(xF(ix1,1,1))
      do ix2=1,numRL(ix1)
        if (abs(xF(ix^D,1))<xmin) xmin=abs(xF(ix^D,1))  ! min value of |x|
      enddo
      if (xmin>0.5) then
        numValidL=numFL
        exit TRACE1
      endif
      ! Transfer valid calculated values from temporary arrays into return arrays
      xFL(ix1,1:numRL(ix1),:)=xF(ix1,1:numRL(ix1),:)
      wBL(ix1,1:numRL(ix1),:)=wPm(ix1,1:numRL(ix1),:)
    enddo TRACE1
    ! Transfer valid calculated values from temporary arrays into return arrays
    !xFL(1:numvalidL,1:numRL,:)=xF(1:numvalidL,1:numRL,:)
    !wBL(1:numvalidL,1:numRL,:)=wPm(1:numvalidL,1:numRL,:)

    ! for starting points at x>0
    xF=0.d0
    xF(1:numValidR,1,:)=xFR(1:numValidR,1,:)        ! xF(numFL,numLP,ndim): field coords
    ftype='Bfield'                                  ! ftype: field type
    tcondi='cons+j'                                 ! tcondi: user branch condition (We want conserved variables + current)
    do iL=1,numValidR                               ! forwardm: Forwards or backwards field tracing
      xpp(1:ndim)=xF(iL,1,1:ndim)
      !write(*,*) "start get_point_w",mype
      call get_point_w(xpp,wpp,'conserved')
      !write(*,*) "end get_point_w",mype
      if (wpp(mag(2))>0.d0) then
        forwardm(iL)=.true.
      else
        forwardm(iL)=.false.
      endif
    enddo
    wPm=0.d0    ! wPm: the variables defined on the fieldline
                ! This should be calculated using the helper function special_field_s and the branching condition "tcondi"
    wLm=0.d0    ! wLm(numFL): This is the variables produced as a single value from a fieldline,
                ! In this case, just 1 variable, the number of valid fieldline points, set all to zero initially
    ! dFh is from the main code from global parameters, fieldline step size= max refinement step
    ! numFL Is the max number of lines that can be traced (array size)
    ! numLP: integer, number of points on the fieldline (max from the array)
    nwP=nw+ndir ! nwP: integer number of "per point" variables, here the conserved +the current
    nwL=0       ! nwL: integer number of "per field line" variables (excluding valid points along it)
    !write(*,*) "start trace_field_multi",mype
    !if (mype==0) then
    !        print*, 'RIGHT - Before trace_field_multi: xF, wBRL = ', xF(180,1,1), wPm(180,1,rho_)
    !endif
    call trace_field_multi(xF,wPm,wLm,dFh,numFL,numLP,nwP,nwL,forwardm,ftype,tcondi)
    !if (mype==0) then
    !        print*, 'RIGHT - After trace_field_multi: xF, wBRL = ', xF(180,1,1), wPm(180,1,rho_)
    !endif
    !write(*,*) "end trace_field_multi",mype
    ! shrink number of field lines to trace if possible
    numRR=int(wLm)
    TRACE2: do ix1=1,numFL
      !if field line is outside of the reconnection region, stop tracing
      xmin=abs(xF(ix1,1,1))
      do ix2=1,numRR(ix1)
        if (abs(xF(ix^D,1))<xmin) xmin=abs(xF(ix^D,1))  ! min value of |x|
      enddo
      if (xmin>0.5) then
        numValidR=numFL
        exit TRACE2
      endif
      ! Transfer valid calculated values from temporary arrays into return arrays
      xFR(ix1,1:numRR(ix1),:)=xF(ix1,1:numRR(ix1),:)
      wBR(ix1,1:numRR(ix1),:)=wPm(ix1,1:numRR(ix1),:)
    enddo TRACE2
    ! Transfer valid calculated values from temporary arrays into return arrays
    !xFR(1:numvalidR,:,:)=xF(1:numvalidR,:,:)
    !wBR(1:numvalidR,:,:)=wPm(1:numvalidR,:,:)

    t1=MPI_wtime()

    !do imype=0,15
    !  if (mype==imype) then 
    !    write(fname,'(a,i4.4,a,i2.2,a)') trim(base_filename), filenr,'_',imype,'.txt'
    !    open(4,file=fname)
    !    write(4,*) '# LR, ix1, ix2, xFLR, wBLR'
    !    do ix1=1,numFL
    !      write(4,*) 'numRL = ', numRL(ix1)
    !      do ix2=1,numRL(ix1)
    !         write(4,*) 'Left',ix1,ix2,xFL(ix1,ix2,:),wBL(ix1,ix2,rho_)
    !      enddo
    !    enddo
    !    do ix1=1,numFL
    !      write(4,*) 'numRR = ', numRR(ix1)
    !      do ix2=1,numRR(ix1)
    !        write(4,*) 'Right',ix1,ix2,xFR(ix1,ix2,:),wBR(ix1,ix2,rho_)
    !      enddo
    !    enddo
    !  endif
    !  close(4)
    !enddo

!    write(*,*) "end   update_Bfield",mype
  end subroutine update_Bfield


  subroutine locate_midpoint(xFL,xFR,wBL,wBR)
    ! Purpose
    !    Find the reference point s_0 along each field line for the starting.
    !       point of the fast particles that descend the loops.
    !    Returns this via setting values in the variables
    !       numTurnL(numValidL), numTurnR(numValidR), iApexL, iApexR
    !    New fast electrons have a given pitch angle mu0 at this point, 
    !       (composite average of old and new particles handled elsewhere)
    ! Parameter list
    !    xFL(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (left side)
    !    xFR(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (right side)
    !    wBL(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !    wBR(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !       L and R indicate left and right side of experiment again.
    ! Calls:         none
    ! Called from:   get_flare_eflux
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)

    double precision :: sign,signNext

    integer :: iFL,iLP,numRT,ix^D
!    write(*,*) "begin locate_midpoint",mype
    numTurnL=1
    numTurnR=1
    iApexL=1
    iApexR=1
    ! Updated by Dubart M., generalised from Ruan W. 
    !    Find the apex of the field line where By changes sign
    !    Label the index of the part of the fieldline where this is as the 
    !    turning point ( numTurnL(line_index) = line_part_number )
    !    If this is at the photospheric base or not found before the end
    !    set the apex value to half of the total (?).
    ! Look for the last field line that goes through the downflow region
    !    iApexR: maximum field line ID on the right side.
    !    iApexL: maximum field line ID on the left side. 
    !    ? what about the minimum index, I guess stepping from the x=0 line this 
    !    ? is unnecessary due to the geometry of the problem.

    ! For B field lines which start from the photosphere in left part (x<0)
    do iFL=1,numValidL
      do iLP=1,numRL(iFL)-1
        sign= wBL(iFL,iLP,mag(2))/abs(wBL(iFL,iLP,mag(2)))
        signNext= wBL(iFL,iLP+1,mag(2))/abs(wBL(iFL,iLP+1,mag(2))) 
        if (signNext*sign<=0.0) then
          if (iLP<=2 .or. iLP>=numRL(iFL)-1) then
            numTurnL(iFL)=(1+numRL(iFL))/2
          else
            numTurnL(iFL)=iLP   ! index of the turning point
            iApexL=iFL          ! index of the last closed field line
          endif
        endif
      enddo
    enddo

    ! For B field lines which start from the photosphere in right part (x>0)
    !    Similar for right hand side
    do iFL=1,numValidR
      do iLP=1,numRR(iFL)-1
        sign= wBR(iFL,iLP,mag(2))/abs(wBR(iFL,iLP,mag(2)))
        signNext= wBR(iFL,iLP+1,mag(2))/abs(wBR(iFL,iLP+1,mag(2))) 
        if (signNext*sign<=0.0) then
          if (iLP<=2 .or. iLP>=numRR(iFL)-1) then
            numTurnR(iFL)=(1+numRR(iFL))/2
          else
            numTurnR(iFL)=iLP   ! index of the turning point
            iApexR=iFL          ! index of the last closed field line
          endif
        endif
      enddo
    enddo

!    write(*,*) "end locate_midpoint",mype
  end subroutine locate_midpoint


  subroutine update_heating_table(QeN,delta,Ec)
    ! Purpose
    !    Calculate heating table based on shock compression ratio
    !    These heating coefficients are calculated for the set of numN values
    !    for the plasma depth from Nmin (1e18) to Nmax(1e23) on a logarithmic
    !    plasma depth scale. They do not include the B-field or Flux_0 values
    !    These are included later for individual fieldlines in get_heating_rate()
    ! 
    !    QeN = Efactor * Beta(fract of log dist Nmin to Nmax)*(Ncol/Nc)^(-delta/2)
    !       Efactor = 2 pi e^4 Lambda_Coulomb (delta -2) / 2 E_c^2
    !       Beta function numerically intergrated over 100,000 pieces.
    !    These are the coefficients of He from Ruan 2020 formula A.12
    !       except for the B-field strengths and the intial flux
    !
    !    Therefore QeN without these factors should be the same each time, why do we
    !       recalculate them??
    !
    ! Parameter list
    !    QeN(numN)   The dE/ds values for heating along each step of plasma depth 
    !                numN is the number of pieces that the logarithmic distribution
    !                of plasma depth is split between the min and max values 
    !                (see get_flare_eflux: last checked 1e18 and 1e23).
    !                (except B and incident flux factors, which can be multiplied in later)
    !    delta   Spectral index of power law energy distribution of 
    !            fast electron energies.
    !    Ec      lower cutoff in energy for the power law of fast electron energies
    ! Calls:         None
    ! 
    ! Called from:   get_flare_eflux
    ! calculate heating table based on shock compression ratio
    use mod_global_parameters
    use mod_usr_methods

    double precision :: QeN(numN),Ncol(numN)
    double precision :: delta,Ec   

    integer :: j,iBx,numBx,iBc,iNcol
    double precision :: Atom,temp
    double precision :: Cllog,gammaN,K,Nc,Nb,Efactor
    double precision :: dtx,tx,au,al,b,maxh
    double precision, allocatable :: Bxc(:)
    double precision :: Ec_erg

!    write(*,*) "begin update_heating_table",mype
    ! parameters for electron deposition
    Cllog=25.0  ! Coulomb logarithm
    gammaN=Cllog  ! fully ionized plasma assumed, not great for chromosphere?
    K=2.0*dpi*const_e**4
    Ec_erg=Ec*1.0e3*const_ev
    Nc=Ec_erg**2/(2.0*gammaN*K)
    Efactor=K*gammaN*(delta-2.0)/(2*Ec_erg**2)

    ! incomplete beta function
    ! Creates a table for 100001 x values from x=0 to 1
    ! Beta_x( a=delta/2 , b= 1.0/ 2) = sum_(t=0)^(t=x) t^(a-1) (1-t)^(b-1) dt
    numBx=100001
    allocate(Bxc(numBx))
    dtx=1.0/(numBx-1)
    Bxc=0
    al=delta/2.0
    b=1.0/2
    do iBx=2,numBx-1
      tx=dtx*(iBx-1.0)
      Bxc(iBx)=Bxc(iBx-1) + (tx**(al-1.0))*((1.0-tx)**(b-1.0))*dtx
    enddo
    Bxc(numBx)=Bxc(numBx-1)

    ! loops over all the numN pieces of the logarithmic plasma depth scale
    ! from Nmin (1e18 when last checked) up to Nmax using dNlog pieces
    ! Efactor = 2 pi e^4 Lambda_Coulomb (delta -2) / 2 E_c^2
    ! 
    ! QeN = Efactor * Beta(fract of log dist Nmin to Nmax)*(Ncol/Nc)^(-delta/2)
    ! The formula for He (Ruan 2020, A.12) without B-field and Flux_0 values added 
    ! These values can be added for individual fieldlines as appropriate, all the 
    ! other values are consistent as long as they share the same Ec and delta beam
    ! parameters
    do iNcol=1,numN
      Ncol(iNcol)=Nmin*10**(dNlog*(iNcol-1.0))
      iBc=int((numBx-1)*Ncol(iNcol)/Nc+0.5)+1
      if (iBc>numBx)  then
        iBc=numBx
      endif
      QeN(iNcol)=Efactor*Bxc(iBc)*((Ncol(iNcol)/Nc)**(-delta/2.0))
    enddo
    deallocate(Bxc)

!    write(*,*) "begin update_heating_table",mype
  end subroutine update_heating_table


  subroutine get_heating_rate(xF,wBLR,QeLR,QeN,QaccLR,numR,numTurn,numValid,iApex,Etot,vs2B,vp2,eFluxb)
    ! Purpose
    !    Calculate local heating rates at each point along individudal field lines
    !       Calculates particle acceleration energy from joule heating term
    !       and resistivity thresholded using distance from centre specified:
    !       Ruan 2020, equ 11:   eta= a_eta (vd/v_c_acc-1) exp[-(y - h_eta)^2/h_sz^2]
    !       Then energy to accelereation Esum = eta |J|^2 * area * dt
    !    Calculate length of field lines that will contain fast electrons
    !       Uses Ruan 2020, equ A4: mu = sqrt( 1 - B/B0(apex)*(1-mu0(apex)^2) )
    !       to find where the electrons can get to using conditions of B field.
    !       where mu reaches zero (90 degree pitch angle) is where the electrons
    !       are stopped. Create a constant (1/B0(apex)*(1-mu0(apex)^2) to do this
    !       before applying it to search the length of the loop to find the limits.
    !       Find the trapping lengths for each field line on both sides of the apex.
    !    Calculate the heating applied to each section of the field line
    !       Use the column depth step calculated here. Combine this with the heating
    !       table vs column depth found in "update_heating_table", and the individual
    !       values of flux and field strength (to find area) for each field line.
    !       QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*unit_length*(eFluxb(ix1)*Bxy0(ix1)/Bxy)/mu
    !       (See Ruan 2020 equ A12 and A14)
    !       NOTE: muMin is a threshold used to prevent total loss of energy.
    !       (heating singularity that occurs in the flux conservation approach)
    ! Parameter list
    !    xF(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines
    !    wBL(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !    QeLR   heating rates calculated for each section of each field line
    !    QeN    heating rate table created in update_heating_rate
    !    QaccLR acceleration energy calculated for each section of each field line
    !    NumR   The number of parts of each individual B-field line.  
    !              NumLP is a max valid value used for array creation, 
    !              and only numR are filled with values.
    !    numTurn   Array of indices of apexes of field lines (max yvalues)
    !    numValid   Number of field lines (starting nearest x=0) that have apexes
    !    iApex    The max fieldline index involve in reconnection, so like numFL (in use)
    !    Etot     Energy budgets of each field line
    !    vs2B     vs^2 / B, V_perp^2 / B first adiabatic const? Ruan 2020 Sec 2 
    !    vp2      v_parallel^2 (in accleration region?)
    !    velocity variable names:
    !    {v}      velocity {e/p/s/e} everything/parallel/"square" i.e. perp
    !             {2} squared {N} new timestep
    !    eFluxb   electron fluxes at lower boundary
    ! Calls:         none
    ! Called from:   get_flare_eflux
    use mod_global_parameters
    use mod_usr_methods

    integer :: numR(numFL),numTurn(numFL)
    double precision :: xF(numFL,numLP,ndim),QeLR(numFL,numLP),QaccLR(numFL,numLP)
    double precision :: wBLR(numFL,numLP,nw+ndir)
    double precision :: QeN(numN),Etot(numFL),vs2B(numFL),eFluxb(numFL)
    double precision :: vp2(numFL)
    integer :: numValid,iApex
    ! Etot, vs2B, vp2 are saved each time so we can average for each fieldline with old values
    double precision :: Bv(numFL,numLP,ndim),Np(numFL,numLP)
    double precision :: Ncol(numLP,2)
    integer :: ixx^D,ix^D,iNlog,iFL,j,i
    double precision :: dl^D,dl,Bxyb,Bxy,mu0,mu,const
    logical :: Ereturn,ForReturn,BackReturn
    integer :: ireturn,iFreturn,iBreturn,iLPmin,iLPmax
    double precision :: ve,length,tau_e,dEdt
    double precision :: width,width0,widthb,area,el,dNcol,Te,v2,J2
    double precision :: eta,reta,rad,rad_acc,vd,muMin,eFlux
    double precision :: muT,vs2Bsum,Esum,vs2,ve2,vs2N,vp2N
    double precision :: elp(-1:1)

    character(50) :: fname
    integer :: filenr_field

!    write(*,*) "begin get_heating_rate",mype
    ! eFluxb -- flux at lower boundary
    ! Bxyb -- magnetic field at lower boundary
    ! widthb -- width of tube at lower boundary
    ! Bxy0(ix1) -- magnetic field at the reference point
    ! QeLR -- heating rate
    ! QaccLR -- acceleration energy rate
    muMin=0.1 ! limit for mu in calculating heating to avoid strong heating near
              ! returning point of fast electron
    ve=1.d10/unit_velocity           ! total velocity
    ve2=1.0d20/(unit_velocity**2)    ! velocity sqaured in acceleration site
    vs2=1.0d19/(unit_velocity**2)    ! v_s=v_square=v_perpendicular, "like right angle"
                                     ! perperdicular velocity squared in
                                     ! acceleration site
                                     ! gives initial mu=18.43 degrees (arcsin(sqrt(0.1))
    ! vs2=((0.00174532837)**2)*ve2   ! for mu=.1 degrees (arcsin(0.1)=0.00174532837)
    ! vs2=((0.01745240644)**2)*ve2   ! for mu=1. degrees (arcsin(1)=0.01745240644)
    ! vs2=((0.08715574275)**2)*ve2   ! for mu=5. degrees (arcsin(5)=0.08715574275)
    ! vs2=((0.17364817767)**2)*ve2   ! for mu=10 degrees (arcsin(10)=0.17364817767)
    ! vs2=((0.2588190451)**2)*ve2    ! for mu=15 degrees (arcsin(15)=0.2588190451)
    ! vs2=((0.5)**2)*ve2             ! for mu=30 degrees (arcsin(30)=0.5)
    ! vs2=((0.70710678119)**2)*ve2   ! for mu=45 degrees (arcsin(45)=0.70710678119)
    ! vs2=((0.86602540378)**2)*ve2   ! for mu=60 degrees (arcsin(60)=0.86602540378)
    ! unit_velocity is 128 km/sec according to Ruan 2020, sec 2.3
    ! density and magnetic field at the field line
    Np(:,:)=wBLR(:,:,rho_)*unit_numberdensity ! density of cgs unit
    do j=1,ndim
      Bv(:,:,j)=wBLR(:,:,mag(j))  ! magnetic field vector
    enddo


    !if (mype==0) then
    !        print*, 'In get_heating_rate: xF, wBRL = ', xF(180,1,1), wBLR(180,1,rho_)
    !endif

    ! calculate fast electron energy flux and heating rate for each field line
    eFluxb=0
    QeLR=0
    QaccLR=0
    do ix1=2,iApex
      vs2Bsum=0.d0                 ! v_perp squared times B (at top?)
      Esum=0.d0

      ! width of the flux tube element at the lower boundary 
      Bxyb=dsqrt(Bv(ix1,1,1)**2+Bv(ix1,1,2)**2)
      widthb=abs(xF(ix1+1,1,1)-xF(ix1-1,1,1))/2.d0
      ! B at reference point, numTurn(line_reference)
      Bxy0(ix1)=dsqrt(Bv(ix1,numTurn(ix1),1)**2+Bv(ix1,numTurn(ix1),2)**2)

      ! calculate particle acceleration rate
      !--------------------------------------------------------------------------!
      !if (mype==0) then
      !  write(fname,'(a,i4.4,a)') trim(base_filename), filenr, '_test.txt'
      !  open(3,file=fname)
      !  write(3,*) '# xF1, xF2, wBLR'
      !  do ixx1=1,numFL
      !    do ixx2=1,numLP
      !      write(3,'(ES17.7E3,ES17.7E3,ES17.7E3)') xF(ixx1,ixx2,1), xF(ixx1,ixx2,2), wBLR(ixx^D,rho_)
      !    enddo
      !  enddo
      !  close(3)
      !endif

      do ix2=2,numR(ix1)-1
        !if (wBLR(ix^D,rho_)==0.0) then
        !  print*, "ix^D = " , ix1, ix2
        !  print*, "xFL =  " , xF(ix1,ix2,:)
        !  print*, "mype = " , mype
        !endif
        Te=wBLR(ix^D,p_)/wBLR(ix^D,rho_)    
        v2=wBLR(ix^D,mom(2))/wBLR(ix^D,rho_)    ! vy! "2 != squared" here

        ! calculate acceleration rate at given region
        if (Te>5.d0 .and. v2<(-v2cr/unit_velocity)) then
          ! Bxy at top of loop calculated from averaging over length points
          ! (2*sqrt(|B(ix2)|^2)+sqrt(|B(ix2-1)|^2)+sqrt(|B(ix2+1)|^2))/4
          Bxy=0.d0
          Bxy=Bxy+2.d0*dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2-1,1)**2+Bv(ix1,ix2-1,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2+1,1)**2+Bv(ix1,ix2+1,2)**2)
          Bxy=Bxy/4.d0
          const=vs2/Bxy         ! first adiabatic constant at looptop Ruan 2020 Sec 2 p4
                                ! but the vs (v_perp) formulation is just a const!

          ! calc area of acceleration region for this flux tube
          dl1=xF(ix1,ix2,1)-xF(ix1,ix2-1,1)    ! deltax
          dl2=xF(ix1,ix2,2)-xF(ix1,ix2-1,2)    ! deltay
          dl=dsqrt(dl1**2+dl2**2) ! "ds" distance between two points in field line
          width=widthb*Bxyb/Bxy   ! "dw" "perp" distance, like "cross section"
          area=dl*width           ! multiply to get area of acceleration region

          ! local acceleration rate
          do i=-1,1
            J2=0.d0
            do j=1,ndir
              J2=J2+wBLR(ix1,ix2+i,nw+j)**2   ! current |J|^2
            enddo

            ! radial distance coefficient (see exponent in Ruan 2020 2.3 equ 11)
            ! acceleration decreases with distance from height h_eta_acc= 50Mm
            ! radial distance for decrease controlled by parameter h_s_acc= 10Mm
            rad_acc=exp(-((xF(ix1,ix2+i,2)-h_eta_acc)**2)/h_s_acc**2)
            vd=sqrt(J2)/wBLR(ix1,ix2+i,rho_)/q_e
                                 ! I guess Ne = rho here, no need for cgs Np var
            eta=0.d0
            if (vd>v_c_acc .and. global_time<t_decay) then
              if (eta_decay_acc) then
                eta=alpha_acc*(vd/v_c_acc-1.d0)*rad_acc   ! Ruan 2020 2.3 Eq 11 if dricer field, with spatial threshold
              else
                eta=alpha_acc*(vd/v_c_acc-1.d0)        ! Ruan 2020 2.3 Eq 11 if dricer field, without spatial threshold
              endif
            endif
            if (eta>eta_max_acc) eta=eta_max_acc   ! Threshold resistivity eta at value eta_max
            elp(i)=eta*J2            ! elp=Joule heating term, we can put this into all the electrons
          enddo
          ! local acceleration energy on B-field line 
          ! (i.e. the Joule heating term via Simpsons rule, Ruan 2020, eq 12)
          el=(elp(-1)+2.d0*elp(0)+elp(1))/4.d0  ! local acceleration rate
          !
          ! improvement point 1 here
          !          using Fleishman info, calc number of electrons
          !          energy accelerates all, total energy affects spectrum
          !
          ! energy got from local current/electric field
          if (const*Bxy0(ix1)>=1) then
            ! if 1e19/unit_velocity/Bxy(local field) *Bxy0(ix1)(apex field)>1
            ! the particle can arrive the reference point
            ! ???? Is this condition to do with particle trapping? Equ A4?
            QaccLR(ix1,ix2)=el
            Esum=Esum+el*area*dt_update_Qe              ! integration of energy
            vs2Bsum=vs2Bsum+el*area*dt_update_Qe*const  ! integration of vs^2/B
                                                        ! first adiabatic const?

          endif
        endif
      enddo
!      if(mype .eq. 0) write(*,*) "heating rate2",alpha_acc

      if (Esum>0) then          ! if there is energy on the fieldline
        vs2N=vs2Bsum*Bxy0(ix1)/Esum  ! vs2_new : perpendicular v^2 at reference point
        if (vs2N>ve2) vs2N=ve2  ! error check, if new perp vel is greater than total v, threshold it!
        vp2N=ve2-vs2N           ! vp2N : parallel v^2 at reference point
      else
        vs2N=0.d0               ! No fast electrons
        vp2N=0.d0               ! No fast electrons
      endif

      if (Etot(ix1)+Esum>0.d0) then    ! if any energy left in loop
        vs2B(ix1)=(Etot(ix1)*vs2B(ix1)+vs2Bsum)/(Etot(ix1)+Esum)  ! average vs^2 from old and new
        vp2(ix1)=(Etot(ix1)*vp2(ix1)+Esum*vp2N)/(Etot(ix1)+Esum)  ! average vp^2 from old and new
        Etot(ix1)=Etot(ix1)+Esum
      else
        Etot(ix1)=0.d0
        vs2B(ix1)=0.d0
        vp2(ix1)=0.d0
      endif
      ! end acceleration rate on this line
      !---------------------------------------------------------------------!

      !#####################################################################!
      ! particle trapping on this line
      length=0
      tau_e=0
      ForReturn=.false.   ! for trapping
      BackReturn=.false.  ! for trapping
      iFreturn=numR(ix1)  ! initial index of one returning point
      iBreturn=1          ! initial index for the other one
      ! If there is fast electron energy in the field line calculate 
      ! local instantaneous pitch angle mu0 using Ruan 2020 equ c11
      ! Perp velocity at the spatial point in question (ix1)
      ! calculated via the first adiabatic constant vs2B = v_perp^2/B=const
      ! Are we assuming that v_parallel is the same all the way down?
      if (Etot(ix1)>0.d0) then
        mu0=sqrt(vp2(ix1)/(vp2(ix1)+vs2B(ix1)*Bxy0(ix1))) ! cosine pitch angle at
                                                     ! at reference point
      ! Else there is no pitch angle
      else
        mu0=0.d0
      endif
      ! then we need to calculate general pitch angle using Ruan 2020 equ A4
      ! const to calculate pitch angle, then need sqrt(1-Blocal*const)
      const=(1.d0-mu0**2)/Bxy0(ix1)  ! const to calculate pitch angle

      ! for the field lines have fast electron, calculate
      ! trapping length for one half of the loop from the apex down.
      ! The lengths along the loop halves are then added up to get a total
      ! variable "length" for the region in which the electrons are trapped.
      if (mu0>0 .and. Etot(ix1)>0) then
        iFreturn=numTurn(ix1) ! numTurn(ix1) -- index of reference point
        iBreturn=numTurn(ix1)
        forward: do ix2=numTurn(ix1),numR(ix1)-1   ! check mu down field lines
          ! mag field av by simpsons rule
          Bxy=0.d0
          Bxy=Bxy+2.d0*dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2-1,1)**2+Bv(ix1,ix2-1,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2+1,1)**2+Bv(ix1,ix2+1,2)**2)
          Bxy=Bxy/4.d0

          ! Check for trapping of electrons owing to field loop geometry (expansion)
          ! This is checking whether mu (Ruan 2020 equ A4) would reach 90 degrees
          if (const*Bxy>=1) then
            ForReturn=.true.  
            iFreturn=ix2-1    ! index for one returning point
            exit forward
          endif

          ! scan area
          dl1=xF(ix1,ix2,1)-xF(ix1,ix2-1,1)
          dl2=xF(ix1,ix2,2)-xF(ix1,ix2-1,2)
          dl=dsqrt(dl1**2+dl2**2)
          length=length+dl  ! length for trapping region
        enddo forward

        ! now check the other side of the loop
        backward: do ix2=numTurn(ix1),2,-1
          ! pitch angle
          Bxy=0.d0
          Bxy=Bxy+2.d0*dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2-1,1)**2+Bv(ix1,ix2-1,2)**2)
          Bxy=Bxy+dsqrt(Bv(ix1,ix2+1,1)**2+Bv(ix1,ix2+1,2)**2)
          Bxy=Bxy/4.d0

          ! return of electrons owing to loop expension
          if (const*Bxy>=1) then
            BackReturn=.true.
            iBreturn=ix2+1  ! index for one returning point
            exit backward
          endif

          ! scan area
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)
          length=length+dl  ! length for trapping region
        enddo backward
      endif
      ! end trapping region
      !#################################################################!

      ! time scale for fast electron moving
      tau_e=length/ve
      if (tau_e<dt_update_Qe) tau_e=dt_update_Qe
      
      ! energy flux of fast electron for each field line
      eFluxb(ix1)=Etot(ix1)/widthb/tau_e

      ! region that has heating
      iLPmin=1
      iLPmax=numR(ix1)
      if (ForReturn) iLPmax=iFreturn
      if (BackReturn) iLPmin=iBreturn
      ! heating region iLPmin:iLPmax

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !calculate heating for field lines have fast electrons
      if (mu0>0 .and. Etot(ix1)>0) then
        ! calculate heating for one side
        Ncol=0
        dEdt=0
        if (numTurn(ix1)<iLPmax-1) then
          do ix2=numTurn(ix1)+1,iLPmax-1
            ! local magnetic field
            Bxy=0.d0
            Bxy=Bxy+2.d0*dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
            Bxy=Bxy+dsqrt(Bv(ix1,ix2-1,1)**2+Bv(ix1,ix2-1,2)**2)
            Bxy=Bxy+dsqrt(Bv(ix1,ix2+1,1)**2+Bv(ix1,ix2+1,2)**2)
            Bxy=Bxy/4.d0

            ! local pitch angle
            mu=sqrt(1-const*Bxy)
            if (mu<muMin) mu=muMin  ! limit for calculating heating

            ! column depth
            dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
            dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
            dl=dsqrt(dl1**2+dl2**2)
            dNcol=dl*unit_length*(Np(ix1,ix2)+Np(ix1,ix2-1))/(2.d0*mu)
            Ncol(ix2,1)=Ncol(ix2-1,1)+dNcol
            width=widthb*Bxyb/Bxy ! local width of the tube

            ! local heating rate
            iNlog=int(log10((Ncol(ix2,1)+1.d0)/Nmin)/dNlog+0.5)+1
            if (iNlog<1) iNlog=1
            if (iNlog>numN) iNlog=numN
            !QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*unit_length*(eFluxb(ix1)*Bxy0(ix1)/Bxy)/mu
            eFlux=eFluxb(ix1)*Bxy/Bxyb ! local energy flux 
            QeLR(ix^D)=eFlux*QeN(iNlog)*(Np(ix^D)*unit_length)/mu  ! local heating rate
            dEdt=dEdt+QeLR(ix^D)*width*dl ! energy for heating
          enddo
        endif

        ! reduce total energy in the beam.
        Etot(ix1)=Etot(ix1)-dEdt*dt_update_Qe

        ! calculate heating for the other side
        Ncol=0
        dEdt=0
        if (numTurn(ix1)>iLPmin+1) then
          do ix2=numTurn(ix1),iLPmin+1,-1
            !Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
            ! local magnetic field
            Bxy=0.d0
            Bxy=Bxy+2.d0*dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
            Bxy=Bxy+dsqrt(Bv(ix1,ix2-1,1)**2+Bv(ix1,ix2-1,2)**2)
            Bxy=Bxy+dsqrt(Bv(ix1,ix2+1,1)**2+Bv(ix1,ix2+1,2)**2)
            Bxy=Bxy/4.d0

            ! local pitch angle
            mu=sqrt(1-const*Bxy)
            if (mu<muMin) mu=muMin  !limit for calculating heating

            ! column depth
            dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
            dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
            dl=dsqrt(dl1**2+dl2**2)
            if (ix2==numTurn(ix1)) then
              dNcol=dl*unit_length*Np(ix1,ix2)/mu
              Ncol(ix2,1)=1.d0
            else
              dNcol=dl*unit_length*(Np(ix1,ix2)+Np(ix1,ix2+1))/(2.d0*mu)
              Ncol(ix2,1)=Ncol(ix2+1,1)+dNcol
            endif
            width=widthb*Bxyb/Bxy ! local width

            ! local heating rate
            iNlog=int(log10((Ncol(ix2,1)+1.d0)/Nmin)/dNlog+0.5)+1
            if (iNlog<1) iNlog=1
            if (iNlog>numN) iNlog=numN
            !QeLR(ix^D)=QeN(iNlog)*Np(ix^D)*unit_length*(eFluxb(ix1)*Bxy0(ix1)/Bxy)/mu
            eFlux=eFluxb(ix1)*Bxy/Bxyb ! local energy flux 
            QeLR(ix^D)=eFlux*QeN(iNlog)*(Np(ix^D)*unit_length)/mu  ! local heating rate
            dEdt=dEdt+QeLR(ix^D)*width*dl ! energy for heating
          enddo
        endif

        Etot(ix1)=Etot(ix1)-dEdt*dt_update_Qe
      endif

      if (Etot(ix1)<0) Etot(ix1)=0
      ! end heating rate
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    enddo

!    write(*,*) "end   get_heating_rate",mype
  end subroutine get_heating_rate


  subroutine get_Qe(xFL,xFR,wBL,wBR,QeL,QeR)
    ! Purpose
    !    Calculate electron heating rate in the (hyperfine) cell grid 
    !    (H_e in Ruan 2020) from heating rates along field lines (QeL, QeR) 
    !    via weightings decribed in appendix D.
    !    output: Qe variable saved on hyperfine grid coords
    !    see "getlQ" for values on experiment grid in the save file...
    !    (1) Paper describes that a fine-refinement grid is created 
    !    (2) Seclect all grid cell centres within dx and dy of field line sect
    !    (3) The heating will be distributed between these grid cells
    !    (4) weigting function for heating is 1/R^f, f="weightindex"=4 
    !        H_e_i= Qe(fieldline) * (1/R_i^4) / sum_{gridcell}[ 1/R_{gridcell}^4]
    ! Parameter list
    !    xFL(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (left side field lines)
    !    xFR(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines (right side field lines)
    !    wBL(numFL, numLP, number of vars + number of directions)
    !       Variable values along left field line from interpolation of cell values.
    !    wBR(numFL, numLP, number of vars + number of directions)
    !       Variable values along right field line from interpolation of cell values.
    !    QeL   heating rates calculated for left sections of each field line
    !    QeR   heating rates calculated for right sections of each field line
    ! Calls:         none
    ! Called from:   get_flare_eflux
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: QeL(numFL,numLP),QeR(numFL,numLP) 
    double precision :: sumWeights(numxQ1,numxQ2)
    double precision :: dxMax1,dxMax2,weight,dl,Bratio,Bratio1,Bratio2
    integer :: iFL,iLP,ixQ^L,ixQ^D,weightIndex
    integer :: ix^D,j
    character(50) :: fname

!    write(*,*) "begin getQe",mype
    Qe=0
    sumWeights=0
    ! get local heating rate via interpolation
    dxMax1=1*dxQ1
    dxMax2=1*dxQ2
    weightIndex=1

    ! for the left foot, loop over data-filled fieldlines, not all fieldlines.
    do iFL=1,numValidL
      ! loop over the fieldline pieces
      do iLP=1,numRL(iFL)
        ! weightindex= decay index of heating dist against dist from fieldline loc.
        ! Bratio is the ratio of field perturbation to background field str (2D).
        ! dxMax1 (1=xdir,2=ydir): spans between grid cell centres that are in the sum.
        weightIndex=4
        Bratio=Busr/sqrt(wBL(iFL,iLP,mag(1))**2+wBL(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ! set up sizes of span for searching through hyperfine gridpoints
        ! xQ(ixQ1,ixQ2,ndim) is the set of hyperfine grid coords
        !    xQ is set up in init_bfield
        ixQmin1=floor((xFL(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFL(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFL(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFL(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        ! do not go outside boundary of experiment
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        ! Loop over the hyperfine gridcells.
        ! Set up weighting table as per Ruan 2020 appendix D.
        ! However, there is a maximum weighting (first condition) to avoid inf.
        ! Store hyperfine grid heating values in Qe.
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFL(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFL(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)+weight*QeL(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo

    ! for the right foot
    do iFL=1,numValidR
      do iLP=1,numRR(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBR(iFL,iLP,mag(1))**2+wBR(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFR(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFR(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFR(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFR(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFR(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFR(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)+weight*QeR(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo

    ! divide by total weights
    do ixQ1=1,numxQ1
      do ixQ2=1,numxQ2
        if (sumWeights(ixQ1,ixQ2)>0) then
          Qe(ixQ1,ixQ2)=Qe(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
        endif
      enddo
    enddo

!    write(*,*) "end   getQe",mype
  end subroutine get_Qe


  subroutine split_Bfield(xF,xFb,wB,Etot,vs2B,vp2,numR,numValid)
    ! Purpose
    !    Handles splitting of Bfield tubes when the width is large
    !    and merging of two tubes when then tube widths are small.
    !    Adds a new field line when maximum distance between adjacent 
    !    field lines is greater than "widthMax" (4 Delta L). 
    !    Removes a field line when the maximum distance between adjacent
    !    field lines is less than "widthMin" (Delta L).
    !    L here is labelled as "Fh", since "L" is being used for left and right.
    !    NOTE: only handles one side of the flare at a time, so called twice
    ! Parameter list
    !    xF(number of field lines, number of line parts, number of dimensions)
    !       spatial coordinates of field lines
    !    xFb(nFL, ndim)
    !       Coordinates of photospheric bases of field lines
    !    wB(numFL, numLP, number of vars + number of directions)
    !       Variable values along field line from interpolation of cell values.
    !    Etot(numFL)   Total energy in fast electrons along each field line.
    !    vs2B(numFL)   v_perpendicular^2/B, first adiabatic constant, 
    !                  relates to pitch angle see Ruan 2020, sec 2.3 p4
    !    vp2(numFL)    v_parallel^2
    !    numR(numFL)   The number of parts of each individual B-field line.  
    !                  NumLP is a max valid value used for array creation, 
    !                  and only numR are filled with values.
    !    numValid      Number of fieldlines being tracked, rather than numFL, 
    !                  which is the max valid number in the allocated arrays.
    ! Calls:           None
    ! Called from      get_flare_eflux
    ! QUESTIONS:
    !    Should the subroutine update numR and numValid? Is that handled elsewhere? Yes... look for it
    !    Should the v values have had the "Etot" factor removed?
    use mod_global_parameters

    double precision :: xF(numFL,numLP,ndim),xFb(numFL,2)
    double precision :: Etot(numFL),vs2B(numFL),vp2(numFL)
    double precision :: wB(numFL,numLP,nw+ndir)
    integer :: numValid
    integer :: numR(numFL)
    
    integer :: ix^D,iFL
    double precision :: dxb,dxp,dxMax,widthMax,widthMin,Bb,Bp
    double precision :: xFnew,Enew,dxLb,dxRb,dEtot,vs2Bnew,vp2new
    logical :: splitB,mergeB

!    write(*,*) "begin split_Bfield",mype
    ! limit for the distance between two field lines
    widthMax=4.0*dFh
    widthMin=1.0*dFh

    ! start from third field line, loop over field lines
    ix1=3
    do while (ix1<numValid)
      ! dxb: B field line separation in x-dir at their base (photosphere)
      ! dxMax: Sets dxb as the max separation at the moment.
      !        Will check up fieldline in the next loop.
      dxb=abs(xF(ix1,1,1)-xF(ix1-1,1,1))
      dxMax=dxb
      ! Bb is |B| at the base of the model (photosphere)... only 2d here
      Bb=sqrt(wB(ix1,1,mag(1))**2+wB(ix1,1,mag(2))**2+wB(ix1,1,mag(3))**2)

      splitB=.FALSE.
      mergeB=.FALSE.

      ! Check the B-field values along each section of the field line
      ! and use this to derive the max separations (A) between the field lines
      ! using the flux conservation approach: A(p)B(p) = A(b)B(b).
      do ix2=1,numR(ix1)
        Bp=sqrt(wB(ix1,ix2,mag(1))**2+wB(ix1,ix2,mag(2))**2+wB(ix1,ix2,mag(3))**2)
        if (Bp<1.d0) Bp=1.d0
        dxp=Bb*dxb/Bp
        if (dxp>dxMax) dxMax=dxp  ! max distance between two fields
      enddo

      ! Conditions to add (split) or remove (merge) a fieldline.
      if (dxMax>widthMax) splitB=.TRUE.
      if (dxMax<widthMin) mergeB=.TRUE.


      ! for Bfield split, add a field line starting at location "xFnew" 
      ! mid-way between two field lines (ix1 & ix1-1)
      !   Reduce energy from old field lines/tubes on each side and add
      !   energy to new field line/tube
      !   Reduction factor Etot_left*0.5*dx_base/(dx_base+dx_base_left)
      !   Reduction factor Etot_right*0.5*dx_base/(dx_base+dx_base_right)
      !   both put into E of new tube, same factor used for vs2Bnew and vp2new
      if (splitB) then
        xFnew=0.5d0*(xF(ix1,1,1)+xF(ix1-1,1,1)) ! starting point of the new line
        Enew=0.d0
        vs2Bnew=0.d0

        ! dxLb = x separation for line ix1 and ix1+1, dxb = same for ix1-1 and ix1
        ! reduce energy from old field line/tube and add
        ! energy to new field line/tube
        dxLb=abs(xF(ix1+1,1,1)-xF(ix1,1,1)) ! distance for line ix1 and ix1+1
        dEtot=Etot(ix1)*(0.25*dxb)/(0.5*(dxLb+dxb))
        Etot(ix1)=Etot(ix1)-dEtot ! reduce energy from line ix1
        Enew=Enew+dEtot ! add energy to new line
        vs2Bnew=dEtot*vs2B(ix1) ! vs2B for new line
        vp2new=dEtot*vp2(ix1) ! vp2 for new line
        dxRb=abs(xF(ix1-1,1,1)-xF(ix1-2,1,1)) ! distance for line ix1-1 and ix1-2
        dEtot=Etot(ix1-1)*(0.25*dxb)/(0.5*(dxRb+dxb)) 
        Etot(ix1-1)=Etot(ix1-1)-dEtot ! reduce energy from line ix1-1
        Enew=Enew+dEtot ! add energy to new line
        vs2Bnew=vs2Bnew+dEtot*vs2B(ix1-1) ! vs2B for new line
        vp2new=vp2new+dEtot*vp2(ix1-1) ! vp2 for new line

        ! add new field line into table
        do iFL=numFL,ix1+1,-1
          xF(iFL,1,1)=xF(iFL-1,1,1)
          xF(iFL,1,2)=0.d0
          Etot(iFL)=Etot(iFL-1)
          vs2B(iFL)=vs2B(iFL-1)
          vp2(iFL)=vp2(iFL-1)
        enddo
        xF(ix1,1,1)=xFnew
        Etot(ix1)=Enew
        vs2B(ix1)=0.d0
        vp2(ix1)=0.d0
        if (Enew>0.d0) vs2B(ix1)=vs2Bnew/Enew
        if (Enew>0.d0) vp2(ix1)=vp2new/Enew
      endif

      ! for Bfield merge, remove field line ix1
      !   All Energy moved to other field line/tube, note factor 0.5
      !   difference from introducing new field line, so all energy included. 
      !   Reduction factor Etot_left*dx_base/(dx_base+dx_base_left)
      !   Reduction factor Etot_right*dx_base/(dx_base+dx_base_right)
      !   both put into new tube, same factor used for vs2Bnew and vp2new
      if (mergeB) then
        ! add energy to field line ix1+1
        dxLb=abs(xF(ix1+1,1,1)-xF(ix1,1,1)) ! distance for ix1 & ix1+1
        dEtot=Etot(ix1)*(0.5*dxb)/(0.5*dxb+0.5*dxLb)
        vs2Bnew=Etot(ix1+1)*vs2B(ix1+1)+dEtot*vs2B(ix1) 
        vp2new=Etot(ix1+1)*vp2(ix1+1)+dEtot*vp2(ix1) 
        Etot(ix1+1)=Etot(ix1+1)+dEtot ! move energy to ix1+1
        ! for pitch angle of ix1+1
        if (Etot(ix1+1)>0) then
          vs2B(ix1+1)=vs2Bnew/Etot(ix1+1)
          vp2(ix1+1)=vp2new/Etot(ix1+1)
        else
          vs2B(ix1+1)=0.d0
          vp2(ix1+1)=0.d0
        endif

        ! add energy to field line ix1-1
        dEtot=Etot(ix1)*(0.5*dxLb)/(0.5*dxb+0.5*dxLb)
        vs2Bnew=Etot(ix1-1)*vs2B(ix1-1)+dEtot*vs2B(ix1) 
        vp2new=Etot(ix1-1)*vp2(ix1-1)+dEtot*vp2(ix1) 
        Etot(ix1-1)=Etot(ix1-1)+dEtot ! move energy to ix1-1
        ! for pitch angle of ix1-1
        if (Etot(ix1-1)>0) then
          vs2B(ix1-1)=vs2Bnew/Etot(ix1-1)
          vp2(ix1-1)=vp2new/Etot(ix1-1)
        else
          vs2B(ix1-1)=0.d0
          vp2(ix1-1)=0.d0
        endif

        ! remove the field from table
        do iFL=ix1,numFL-1
          xF(iFL,1,1)=xF(iFL+1,1,1)
          xF(iFL,1,2)=0.d0
          Etot(iFL)=Etot(iFL+1)
          vs2B(iFL)=vs2B(iFL+1)
          vp2(iFL)=vp2(iFL+1)
        enddo
        xF(numFL,1,1)=2*xF(numFL-1,1,1)-xF(numFL-2,1,1)
        Etot(numFL)=0.d0
      endif

      ix1=ix1+1
    enddo

    xFb(:,1)=xF(:,1,1)
    xFb(:,2)=0.d0
    xF(ix1,1,1)=xFnew

!    write(*,*) "end   split_Bfield",mype
  end subroutine split_Bfield


  subroutine get_spectra(delta,Ec,dEe,Ee,spectra,numEe)
    ! Purpose
    !    Calculate fast electron energy distribution spectrum, split into E bins
    !    Use the power distribution and split over the 1D depth locations.
    !    This is for calculating the HXR signals only, not for distributing heating.
    ! Parameter list
    !    delta: spectral index of fast electrons
    !    Ec: Lower energy cutoff of fast electrons
    !    dEe: step size of energy bins
    !    Ee: energy values
    !    Spectra(numN, numEe) : the number of electrons in each energy bin
    !             numN is the number of column depth points (1001)
    !    numEe: the number of energy bins
    ! Calls:         none
    ! Called from:   get_flare_eflux (in convert mode)
    use mod_global_parameters

    double precision :: delta,Ec,dEe
    integer :: numEe
    double precision :: Ee(numEe),spectra(numN,numEe)
    integer :: iEe,iNcol
    double precision :: Ncol(numN)
    double precision :: K,lambda_,beta
    double precision :: Fe0,Fe,qt,B0,const,const2,keV_erg
    double precision :: Eph,alpha,mec2,r0,N0,c,temp
    double precision :: E0min,E0max,spec0min,spec0max,dEe0
    double precision :: sigma0
    character (50) :: fname

!    write(*,*) "begin get_spectra",mype
    mec2=511.0      ! energy of a static electron [keV]
    alpha=1.0/137.0 ! fine structure constant
    c=2.9979d10     ! light speed [cm/s]
    r0=2.8179d-13   ! radiu of electron [cm]
    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]
    K=2*dpi*const_e**4
    lambda_=25.0
    beta=2.0

    ! initial fast electron spectra
    spectra=0.d0
    const=(delta-1)/Ec
    do iEe=1,numEe
      Ee(iEe)=Ec+(iEe-1)*dEe
      spectra(1,iEe)=const*(Ee(iEe)/Ec)**(-delta)
    enddo

    ! change of electrons energy
    Ncol=0.d0
    do iNcol=2,numN
      Ncol(iNcol)=Nmin*10**(dNlog*(iNcol-1.d0))
      temp=2.d0*lambda_*K*Ncol(iNcol)/keV_erg**2
      do iEe=1,numEe
        E0min=sqrt(Ee(iEe)**2+temp)
        E0max=sqrt((Ee(iEe)+dEe)**2+temp)
        dEe0=E0max-E0min
        spec0min=const*(E0min/Ec)**(-delta)
        spec0max=const*(E0max/Ec)**(-delta)
        spectra(iNcol,iEe)=0.5d0*(spec0min+spec0max)*dEe0/dEe
      enddo
    enddo

    spectra=spectra*keV_erg   ! [electrons cm^-2 s^-1 keV^-1]

!    write(*,*) "end   get_spectra",mype
  end subroutine get_spectra


  subroutine get_HXR_line(Ee,spectra,numEe,dEe,xF,wBLR,numValid,numTurn,numR,eFluxb,vs2B,vp2,HXRLR,muLR,eFluxLRP)
    ! Purpose
    !    Calculate HXR spectra along fieldline pieces
    !    See Ruan 2020, Appendix E, equations E1, E2, E3
    ! Parameter list
    !    Ee: energy values
    !    Spectra(numN, numEe) : the number of electrons in each energy bin
    !             numN is the number of column depth points (1001)
    !    numEe: the number of energy bins
    !    dEe: step size of energy bins
    !    xF:
    !    wBLR : varaible values along left and right fieldlines
    !    numValid : number of active fieldlines
    !    numTurn : The points/locations at which the electrons are considered
    !              to originate from in the loops
    !    numR   : the number of fieldlines in the right hand set
    !    eFluxb : the electron flux at the base of the model?
    !    vs2B   : vel squared parallel or perp to fieldline
    !    vp2    : vel squared parallel or perp to fieldline
    !    HXRLR  : output array for hard x-ray contributions
    ! Calls:         none
    ! Called from:   get_flare_eflux (in convert mode)
    use mod_global_parameters

    double precision :: dEe
    integer :: numEe,numValid
    double precision :: Ee(numEe),spectra(numN,numEe)
    integer :: numR(numFL),numTurn(numFL)
    double precision :: xF(numFL,numLP,ndim),QeLR(numFL,numLP)
    double precision :: wBLR(numFL,numLP,nw+ndir)
    double precision :: eFluxb(numFL),vs2B(numFL),vp2(numFL),HXRLR(numFL,numLP)
    integer :: j,ix^D,iNlog,ireturn
    double precision :: Ncol(numLP),Np(numFL,numLP),Bv(numFL,numLP,ndim)
    double precision :: mu0,mu,Bxy,Bxyb,Bratio,dl1,dl2,dl,const,muMin
    logical :: Ereturn
    double precision :: Ephmin,Ephmax,dEph,Eph
    integer :: numEph,iEph,iEe
    double precision :: sigma0,sigmaBH,keV_erg,temp,Fe
    double precision :: HXRLRs(numFL,numLP)
    double precision :: muLR(numFL,numLP),muLRs(numFL,numLP)
    double precision :: eFluxLRP(numFL,numLP),eFluxLRPs(numFL,numLP)
    double precision :: BxyTurn,ve,ve2,vs2,eFlux

!    write(*,*) "begin get_HXR_line",mype
    HXRLR=0.d0
    HXRLRs=0.d0
    muLR=0.d0    
    eFluxLRP=0.d0

    keV_erg=1.0d3*const_ev  ! 1 keV = * erg
    sigma0=7.9e-25  ! [cm^2 keV]

    ! energy of HXR photon    
    Ephmin=25
    Ephmax=50
    dEph=1.d0
    numEph=floor((Ephmax-Ephmin)/dEph)

    muMin=0.1
    ve=1.d10/unit_velocity
    ve2=1.0d20/(unit_velocity**2)

    ! density and magnetic field
    Np(:,:)=wBLR(:,:,rho_)*unit_numberdensity
    do j=1,ndim
      Bv(:,:,j)=wBLR(:,:,mag(j))
    enddo

    ! calculate fast electron energy flux for each field line
    do ix1=1,numValid

      ! width of the flux tube element
      Bxyb=dsqrt(Bv(ix1,1,1)**2+Bv(ix1,1,2)**2) ! B at lower boundary
      Bxy0(ix1)=dsqrt(Bv(ix1,numTurn(ix1),1)**2+Bv(ix1,numTurn(ix1),2)**2) ! B at
                                                        ! reference point
      mu0=sqrt(vp2(ix1)/(vp2(ix1)+vs2B(ix1)*Bxy0(ix1))) ! pitch angle at reference point
      const=(1.d0-mu0**2)/Bxy0(ix1)

      if (mype==mod(ix1,npe) .and. eFluxb(ix1)>0.d0) then
        Ncol=0
        FIELD1: do ix2=numTurn(ix1),numR(ix1)
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          if (const*Bxy>=1) exit FIELD1
          mu=sqrt(1-const*Bxy)
          if (mu<muMin) mu=muMin
          muLRs(ix^D)=mu

          ! column depth
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
          Ncol(ix2)=Ncol(ix2-1)+dl*(Np(ix1,ix2)+Np(ix1,ix2-1))/2.0

          ! local heating rate
          iNlog=int(log10(Ncol(ix2)/Nmin)/dNlog+0.5)+1
          if (iNlog<1) iNlog=1
          if (iNlog>numN) iNlog=numN

          ! calulate HXR flux for given energy range
          eFlux=eFluxb(ix1)*Bxy/Bxyb ! local energy flux
          eFluxLRPs(ix^D)=eFlux
          do iEph=1,numEph
            Eph=Ephmin+(iEph-1)*dEph
            do iEe=1,numEe
              if (Ee(iEe)>Eph) then
                temp=sqrt(1-Eph/Ee(iEe))
                sigmaBH=(sigma0/(Eph*Ee(iEe)))*log((1+temp)/(1-temp))
                HXRLRs(ix^D)=HXRLRs(ix^D)+eFlux*spectra(iNlog,iEe)*Np(ix^D)*sigmaBH/mu
              endif
            enddo
          enddo
        enddo FIELD1

        Ncol=0
        FIELD2: do ix2=numTurn(ix1)-1,1,-1
          Bxy=dsqrt(Bv(ix1,ix2,1)**2+Bv(ix1,ix2,2)**2)
          if (const*Bxy>=1) exit FIELD2
          mu=sqrt(1-const*Bxy)
          if (mu<muMin) mu=muMin
          muLRs(ix^D)=mu

          ! column depth
          dl1=xF(ix1,ix2+1,1)-xF(ix1,ix2,1)
          dl2=xF(ix1,ix2+1,2)-xF(ix1,ix2,2)
          dl=dsqrt(dl1**2+dl2**2)*unit_length/mu
          Ncol(ix2)=Ncol(ix2+1)+dl*(Np(ix1,ix2)+Np(ix1,ix2+1))/2.0

          ! local heating rate
          iNlog=int(log10(Ncol(ix2)/Nmin)/dNlog+0.5)+1
          if (iNlog<1) iNlog=1
          if (iNlog>numN) iNlog=numN

          ! calulate HXR flux for given energy range
          eFlux=eFluxb(ix1)*Bxy/Bxyb ! local energy flux
          eFluxLRPs(ix^D)=eFlux
          do iEph=1,numEph
            Eph=Ephmin+(iEph-1)*dEph
            do iEe=1,numEe
              if (Ee(iEe)>Eph) then
                temp=sqrt(1-Eph/Ee(iEe))
                sigmaBH=(sigma0/(Eph*Ee(iEe)))*log((1+temp)/(1-temp))
                HXRLRs(ix^D)=HXRLRs(ix^D)+eFlux*spectra(iNlog,iEe)*Np(ix^D)*sigmaBH/mu
              endif
            enddo
          enddo
        enddo FIELD2

      endif
    enddo

    call MPI_ALLREDUCE(HXRLRs,HXRLR,numFL*numLP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
    
    call MPI_ALLREDUCE(muLRs,muLR,numFL*numLP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)

    call MPI_ALLREDUCE(eFluxLRPs,eFluxLRP,numFL*numLP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)

!    write(*,*) "end   get_HXR_line",mype
  end subroutine get_HXR_line


  subroutine interp_HXR(xFL,xFR,wBL,wBR,HXRL,HXRR,muL,muR,eFluxLP,eFluxRP)
    ! Purpose
    !    Calculate HXR spectra contributions in hyperfine grid
    !    using output from calculations along fieldlines (HXRLR, in get_HXR_line)
    !    Interpolated analgously to the heating rate on that hyperfine grid
    !    only difference I can spot is the change in the "weight index"... why?
    ! Parameter list
    !    xF: coordiinates along field loops
    !    wBL : varaible values along left fieldlines
    !    wBR : varaible values along right fieldlines
    !    HXRL  : output array for hard x-ray contributions along left fieldlines
    !    HXRR  : output array for hard x-ray contributions along right fieldlines
    !
    ! Calls
    !    none
    ! 
    ! Called from
    !    get_flare_eflux (in convert mode)
    use mod_global_parameters
    use mod_usr_methods

    double precision :: xFL(numFL,numLP,ndim),xFR(numFL,numLP,ndim)
    double precision :: wBL(numFL,numLP,nw+ndir),wBR(numFL,numLP,nw+ndir)
    double precision :: HXRL(numFL,numLP),HXRR(numFL,numLP)
    double precision :: muL(numFL,numLP),muR(numFL,numLP)
    double precision :: eFluxLP(numFL,numLP),eFluxRP(numFL,numLP)
 
    double precision :: sumWeights(numxQ1,numxQ2)
    double precision :: dxMax1,dxMax2,weight,dl,Bratio
    integer :: iFL,iLP,ixQ^L,ixQ^D,weightIndex

!    write(*,*) "begin interp_HXR",mype
    HXR=0
    mu_electrons=0
    electron_flux=0
    sumWeights=0
    dxMax1=8*dxQ1
    dxMax2=8*dxQ2
    weightIndex=1

    ! for the left foot
    do iFL=1,numValidL
      do iLP=1,numRL(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBL(iFL,iLP,mag(1))**2+wBL(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFL(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFL(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFL(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFL(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFL(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFL(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)+weight*HXRL(iFL,iLP)
            mu_electrons(ixQ1,ixQ2)=mu_electrons(ixQ1,ixQ2)+weight*muL(iFL,iLP)
            electron_flux(ixQ1,ixQ2)=electron_flux(ixQ1,ixQ2)+weight*eFluxLP(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo

    ! for the right foot
    do iFL=1,numValidR
      do iLP=1,numRR(iFL)
        weightIndex=4
        Bratio=Busr/sqrt(wBR(iFL,iLP,mag(1))**2+wBR(iFL,iLP,mag(2))**2)
        if (Bratio<3) then
          dxMax1=ceiling(Bratio)*dxQ1
          dxMax2=ceiling(Bratio)*dxQ2
        else
          dxMax1=3*dxQ1
          dxMax2=3*dxQ2
        endif

        ixQmin1=floor((xFR(iFL,iLP,1)-dxMax1-xQmin1)/dxQ1)+1
        ixQmin2=floor((xFR(iFL,iLP,2)-dxMax2-xQmin2)/dxQ2)+1
        ixQmax1=floor((xFR(iFL,iLP,1)+dxMax1-xQmin1)/dxQ1)+1
        ixQmax2=floor((xFR(iFL,iLP,2)+dxMax2-xQmin2)/dxQ2)+1
        if (ixQmin1<1) ixQmin1=1
        if (ixQmin2<1) ixQmin2=1
        if (ixQmax1>numxQ1) ixQmax1=numxQ1
        if (ixQmax2>numxQ2) ixQmax2=numxQ2

        !#
        do ixQ1=ixQmin1,ixQmax1
          do ixQ2=ixQmin2,ixQmax2
            dl=sqrt((xQ(ixQ1,ixQ2,1)-xFR(iFL,iLP,1))**2+&
                    (xQ(ixQ1,ixQ2,2)-xFR(iFL,iLP,2))**2)
            if (dl<1.0d-2*dxQ1) then
              weight=(1/(1.0d-2*dxQ1))**weightIndex
            else
              weight=(1/dl)**weightIndex
            endif
            sumWeights(ixQ1,ixQ2)=sumWeights(ixQ1,ixQ2)+weight
            HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)+weight*HXRR(iFL,iLP)
            mu_electrons(ixQ1,ixQ2)=mu_electrons(ixQ1,ixQ2)+weight*muR(iFL,iLP)
            electron_flux(ixQ1,ixQ2)=electron_flux(ixQ1,ixQ2)+weight*eFluxRP(iFL,iLP)
          enddo
        enddo
        !#
      enddo
    enddo


    ! divide by total weights
    do ixQ1=1,numxQ1
      do ixQ2=1,numxQ2
        if (sumWeights(ixQ1,ixQ2)>0) then
          HXR(ixQ1,ixQ2)=HXR(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
          mu_electrons(ixQ1,ixQ2)=mu_electrons(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
          electron_flux(ixQ1,ixQ2)=electron_flux(ixQ1,ixQ2)/sumWeights(ixQ1,ixQ2)
        endif
      enddo
    enddo

!    write(*,*) "end   interp_HXR",mype
  end subroutine interp_HXR


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! Purpose
  !   This subroutine can be used in the "convert" stage of amrvac, 
  !   It adds additional (auxiliary) variables to the standard vtk file output
  !     for further analysis using tecplot, paraview, etc...
  !     these auxiliary values need to be stored in the nw+1:nw+nwauxio slots  
  !   This routine should set the values of
  !     w(:,nw+1:nw+nwauxio)
  !     normconv(:,nw+1:nw+nwauxio)
  !
  !     w: auxiliary var values stored in w using (:,nw+1:nw+nwauxio) entries
  !
  !     The array normconv can be filled in the (nw+1:nw+nwauxio) range 
  !     corresponding normalization values (default value 1)
  !     normconv is an array to store the corresponding normalization values 
  !     (default value 1) to convert into the correct units
  ! Parameter list
  !     ixI^L	: Coordinate indices inside the domain
  !     ixO^L	: Coordinate indices including the outside (ghost zones)
  !     w	: Variables stored in this array
  !     x	: coordinates
  !     normconv: unit conversion for each variable in w
  ! Calls to calculate
  !    |-> mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)                 for T
  !    |-> divvector(Btotal,ixI^L,ixO^L,divb)                    for divB
  !    |-> get_normalized_divb(Btotal,ixI^L,ixO^L,divb1)         for divB
  !    |-> get_current(w,ixI^L,ixO^L,idirmin,current_o)          for j1234
  !    |-> special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,divb)   for resistivity eta
  !    |-> getlQ(lQgrid,ixI^L,ixO^L,t,w,x)                       for electron heating lQ
  !    |-> getcQ(ens,ixI^L,ixO^L,t,w,x)                          for energy acceleration cQ
  !    |-> getbQ(ens,ixI^L,ixO^L,t,w,x)                          for background heating bQ
  !    |-> getvar_cooling(ixI^L,ixO^L,w,x,ens,fl)                for radiative cooling "rad"
  !    |-> get_HXR(HXRgrid,ixI^L,ixO^L,t,w,x)                    for HXR
  !    |-> get_mu(mugrid,ixI^L,ixO^L,t,w,x)                      for pitch-angle
  !    |-> cross_product(ixI^L,ixO^L,Btotal,v,Efield)            for E123
  !    |-> get_sxr_flare(ixI^L,ixO^L,w,x,HXRgrid,6,12)           for SXR
  !    |-> get_refine_region(ixI^L,ixO^L,w,x,rfgrid)             for RF
  !    |-> divvector(Btotal,ixI^L,ixO^L,divb)                    for divV (stored using divV var names!)
  ! Call from
  !     Main code (substitution for usr_aux_output)
  !
  ! nwauxio is set in "amrvar.par"... the names are set in specialvarnames_output
  ! Auxvars: T=Pthermal/rho, C_Alfven = sqrt(|B^2|/rho), 
  !          DivB = [0.5*(divvector(B(indx,1:3)))/|B^2| ] / (1/dxlevel_of_block),
  !          Plasma Beta= 2*Pthermal/B^2,
  !          j1,2,3 from get_current(). idirmin is a minimum current value.
  !          Resistivity Eta=special_eta() careful, stored in divB for memory use
  !          Electron heating of plasma, lQ=getlQ()
  !          Energy removal (joule heating into acceleration) getcQ
  !          Background heating get bQ, radiative cooling: rad=getvar_cooling()
  !          HXR HXRgrid=get_HXR()
  !          E1,2,3= v cross B + eta*J (v= mom/rho, careful: eta saved in divB!) 
  !          SXR saved from HXRgrid=get_SXR()
  !          RF pick out points where the mesh refinement level is 
  !             greater than some value rf=rfgrid (output from get_refine_region)
  !          divV
    use mod_global_parameters
    use mod_radiative_cooling
    use mod_thermal_emission
    use mod_mhd_phys

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision :: pth(ixI^S),B2(ixI^S),divb(ixI^S),divb1(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir),current_o(ixI^S,3)
    integer :: idir,idirmin,nw_counter
    double precision :: lQgrid(ixI^S),ens(ixI^S),HXRgrid(ixI^S), mugrid(ixI^S)
    double precision :: t, El_GOES, Eu_GOES
    double precision :: Efield(ixI^S,1:ndir),v(ixI^S,1:ndir),poynting(ixI^S,1:ndir)
    logical :: rfgrid

!    write(*,*) "begin specialvar_output",mype
    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
    if(B0field) then
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))+block%B0(ixI^S,1:ndir,0)
    else
      Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    endif
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
!    w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    call get_normalized_divb(Btotal,ixI^L,ixO^L,divb1)
!    divb=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
!    print*, "size divb", size(divb(ixO^S))
!    print*, "size divb1", size(divb1(ixO^S))
!    print*, "divB diff", maxval(abs(divb1(ixO^S)-divb(ixO^S)))
    w(ixO^S,nw+3)=divb1(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output current
    call get_current(w,ixI^L,ixO^L,idirmin,current_o)
    w(ixO^S,nw+5)=current_o(ixO^S,1)
    w(ixO^S,nw+6)=current_o(ixO^S,2)
    w(ixO^S,nw+7)=current_o(ixO^S,3)
    ! output special resistivity eta
    call special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,divb)
    w(ixO^S,nw+8)=divb(ixO^S)

    ! output heating rate
    t=global_time
    call getlQ(lQgrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+9)=lQgrid(ixO^S)

    call getcQ(ens,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+10)=ens(ixO^S)

    !! background heating
    call getbQ(ens,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+11)=ens(ixO^S)

    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,w,x,ens,rc_fl)
    w(ixO^S,nw+12)=ens(ixO^S)

    call get_HXR(HXRgrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+13)=HXRgrid(ixO^S)

    ! electric field
    do idir=1,ndir
      v(ixO^S,idir)=w(ixO^S,mom(idir))/w(ixO^S,rho_)
    enddo
    call cross_product(ixI^L,ixO^L,Btotal,v,Efield)
    do idir=1,ndir
      Efield(ixO^S,idir)=Efield(ixO^S,idir)+divb(ixO^S)*current_o(ixO^S,idir)
      w(ixO^S,nw+13+idir)=Efield(ixO^S,idir)
    enddo

!    ! SXR
!    call get_sxr_flare(ixI^L,ixO^L,w,x,HXRgrid,6,12)
!    w(ixO^S,nw+17)=HXRgrid(ixO^S)
    ! SXR
    El_GOES=1.5    ! 8Å converted to keV
    Eu_GOES=12.4   ! 1Å converted to keV
    call get_sxr_flare(ixI^L,ixO^L,w,x,HXRgrid,El_GOES,Eu_GOES)
    w(ixO^S,nw+17)=HXRgrid(ixO^S)

    ! special refinement
    call get_refine_region(ixI^L,ixO^L,w,x,rfgrid)
    if (rfgrid) then
      w(ixO^S,nw+18)=1.d0
    else
      w(ixO^S,nw+18)=zero
    endif

    ! output divV
    do idir=1,ndir
       Btotal(ixI^S,idir)=w(ixI^S,mom(idir))/w(ixI^S,rho_)
    enddo
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+19)=divb(ixO^S)
!    w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)

    call get_mu(mugrid,ixI^L,ixO^L,t,w,x)
    w(ixO^S,nw+20)=mugrid(ixO^S)

    ! output Poynting Flux (Added by Malcolm D.)
    if (nwauxio > 20) then
       nw_counter=20
       nw_counter=nw_counter+1
       call get_EUV(171,ixI^L,ixO^L,w,x,te_fl_mhd,HXRgrid)
       w(ixO^S,nw+nw_counter)=HXRgrid(ixO^S)
       nw_counter=nw_counter+1
       call get_EUV(193,ixI^L,ixO^L,w,x,te_fl_mhd,HXRgrid)
       w(ixO^S,nw+nw_counter)=HXRgrid(ixO^S)
       nw_counter=nw_counter+1
       call get_EUV(211,ixI^L,ixO^L,w,x,te_fl_mhd,HXRgrid)
       w(ixO^S,nw+nw_counter)=HXRgrid(ixO^S)
       nw_counter=nw_counter+1
       call get_EUV(304,ixI^L,ixO^L,w,x,te_fl_mhd,HXRgrid)
       w(ixO^S,nw+nw_counter)=HXRgrid(ixO^S)
       ! Poynting_flux
       call cross_product(ixI^L,ixO^L,Efield,Btotal,poynting)
       do idir=1,ndir
          nw_counter=nw_counter+1
          w(ixO^S,nw+nw_counter)=poynting(ixO^S,idir)
       enddo
    endif

!    write(*,*) "end   specialvar_output",mype
  end subroutine specialvar_output


  subroutine specialvarnames_output(varnames)
  ! Purpose
  !    user-added variable names concatenated with standard list of
  !       w_names/primnames string
  ! Parameter list
  !    varnames: string, list of auxiliary variables defined by user to save
    character(len=*) :: varnames
    if (nwauxio>20) then
        varnames= &
        'Te Alfv divB beta j1 j2 j3 eta lQ cQ bQ rad HXR E1 E2 E3  & 
        &SXR RF divV mu AIA171 AIA193 AIA211 AIA304 pf1 pf2 pf3'
    else
        varnames= &
        'Te Alfv divB beta j1 j2 j3 eta lQ cQ bQ rad HXR E1 E2 E3 SXR RF divV mu'
    endif
  end subroutine specialvarnames_output


  subroutine get_HXR(HXRgrid,ixI^L,ixO^L,qt,w,x)
    ! Purpose
    !    Calculate HXR spectra contributions on experiment grid
    !    using output from interp_HXR HXR spectra on hyperfine grid
    !    Interpolated analgously to the heating rate
    !    only difference I can spot is the change in the "weight index"... why?
    !    Used from specialvar output to save to files.
    ! Parameter list
    !    HXRgrid  : output array for hard x-ray contributions on experiment grid
    !    ixI, ixO : mesh coordinate indices, input and output?
    !    O includes ghost cells
    !    qt is solar time in experiment
    !    w is array of variables at each index
    !    x: coordinates at each index
    ! Calls:         none
    ! Called from:   specialvar_output
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: HXRgrid(ixI^S)
    integer :: ix^D,ixO^D,ixQ^L,numHXR
    double precision :: xc^L,sumHXR

!    write(*,*) "begin get_HXR",mype
    select case(iprob)    

    case(3)
      HXRgrid(ixO^S)=0
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        xcmin^D=x(ixO^DD,^D)-0.5d0*dxlevel(^D);
        xcmax^D=x(ixO^DD,^D)+0.5d0*dxlevel(^D);
        ixQmin^D=floor((xcmin^D-xQmin^D)/dxQ^D)+1;
        ixQmax^D=floor((xcmax^D-xQmin^D)/dxQ^D);

        sumHXR=0.d0
        numHXR=0
        do ix1=ixQmin1,ixQmax1
          do ix2=ixQmin2,ixQmax2
            if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
              sumHXR=sumHXR+HXR(ix1,ix2)
              numHXR=numHXR+1
            endif
          enddo
        enddo

        if (numHXR>0) then
             HXRgrid(ixO^D)=sumHXR/numHXR
        endif
      {enddo\}

    end select

!    write(*,*) "end   get_HXR",mype
  end subroutine get_HXR

  subroutine get_mu(mugrid,ixI^L,ixO^L,qt,w,x)
    ! Purpose (written by Maxime Dubart)
    !    Calculate pitch-angle contributions on experiment grid
    !    using output from interp_HXR mu_electrons on hyperfine grid
    !    Interpolated analgously to the heating rate
    !    only difference I can spot is the change in the "weight index"... why?
    !    Used from specialvar output to save to files.
    ! Parameter list
    !    mugrid  : output array for pitch-angle contributions on experiment grid
    !    ixI, ixO : mesh coordinate indices, input and output?
    !    O includes ghost cells
    !    qt is solar time in experiment
    !    w is array of variables at each index
    !    x: coordinates at each index
    ! Calls:         none
    ! Called from:   specialvar_output
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: mugrid(ixI^S)
    integer :: ix^D,ixO^D,ixQ^L,nummu
    double precision :: xc^L,summu

!    write(*,*) "begin get_HXR",mype
    select case(iprob)    

    case(3)
      mugrid(ixO^S)=2.d0
      {do ixO^DB=ixOmin^DB,ixOmax^DB\}
        xcmin^D=x(ixO^DD,^D)-0.5d0*dxlevel(^D);
        xcmax^D=x(ixO^DD,^D)+0.5d0*dxlevel(^D);
        ixQmin^D=floor((xcmin^D-xQmin^D)/dxQ^D)+1;
        ixQmax^D=floor((xcmax^D-xQmin^D)/dxQ^D);

        summu=0.d0
        nummu=0
        do ix1=ixQmin1,ixQmax1
          do ix2=ixQmin2,ixQmax2
            if (ix1>=1 .and. ix1<=numXQ1 .and. ix2>=1 .and. ix2<=numXQ2) then
              if (electron_flux(ix1,ix2) > 0.d0) then
                  summu=summu+mu_electrons(ix1,ix2)
                  nummu=nummu+1
              endif
            endif
          enddo
        enddo

        if (nummu>0) then
             mugrid(ixO^D)=summu/nummu
        endif
      {enddo\}

    end select

!    write(*,*) "end   get_HXR",mype
  end subroutine get_mu

  subroutine get_SXR_flare(ixI^L,ixO^L,w,x,flux,El,Eu)
    !synthesize thermal SXR from El keV to Eu keV
    !flux (cgs): photons cm^-5 s^-1
    !flux (SI): photons m^-3 cm^-2 s^-1
    !integration of the flux is the SXR flux observed at 1AU [photons cm^-2 s^-1]
    use mod_global_parameters
    use mod_physics

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: El,Eu
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(in)  :: w(ixI^S,nw)
    double precision, intent(out) :: flux(ixI^S) 

    integer :: ix^D,ixO^D
    integer :: iE,numE
    double precision :: I0,kb,keV,dE,Ei
    double precision :: pth(ixI^S),Te(ixI^S),kbT(ixI^S)
    double precision :: Ne(ixI^S),gff(ixI^S),fi(ixI^S) !, flux1(ixI^S)
    double precision :: EM(ixI^S), undercheck(ixI^S), checkval, checklevel

!    write(*,*) "begin get_SXR_flare",mype
!    I0=1.07d-42    ! photon flux index for observed at 1AU [photon cm s^-1 keV^-1]
    checklevel=1.0d15
    I0=1.0d0
    kb=const_kb
    keV=1.0d3*const_ev
    dE=0.1
    numE=floor((Eu-El)/dE)
    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,iw_rho)*unit_temperature
    if (SI_unit) then
      Ne(ixO^S)=w(ixO^S,iw_rho)*unit_numberdensity/1.d6 ! m^-3 -> cm-3
      EM(ixO^S)=(Ne(ixO^S))**2*1.d6 ! cm^-3 m^-3
    else
      Ne(ixO^S)=w(ixO^S,iw_rho)*unit_numberdensity
      EM(ixO^S)=(Ne(ixO^S))**2
    endif
    kbT(ixO^S)=kb*Te(ixO^S)/keV
    flux(ixO^S)=0.0d0
!    flux1(ixO^S)=0.0d0
    do iE=0,numE-1
      undercheck(ixO^S)=0.0d0
      Ei=dE*iE+El*1.d0
      gff(ixO^S)=1.d0
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(kbT(ix^D)<Ei) then
          gff(ix^D)=(kbT(ix^D)/Ei)**0.4
        endif
      {enddo\}
      fi(ixO^S)=(EM(ixO^S)*gff(ixO^S))* &
                exp(-Ei/(kbT(ixO^S)))/(Ei*sqrt(kbT(ixO^S)))
!      where (fi*dE>0.d0)
!         undercheck=flux/(fi*dE)
!      elsewhere
!         undercheck=0.0d0
!      end where
!      where (undercheck < checklevel .and. undercheck > 0.d0)
!         flux=flux+fi*dE
!      elsewhere
!         flux=flux
!      end where
!      flux(ixO^S)=flux(ixO^S)+fi(ixO^S)*dE
      flux(ixO^S)=flux(ixO^S)+fi(ixO^S)*dE*Ei*keV
    enddo
!    print*, "fluxdiff", maxval(flux1-flux)
!    print*, "fluxmean", sum(flux)/(max(1,size(flux)))
!    print*, "fluxmax", maxval(flux)
    flux(ixO^S)=flux(ixO^S)*I0
!    write(*,*) "end   get_SXR_flare",mype
  end subroutine get_SXR_flare

  double precision function active_By(xActive,Bmax)
  ! Purpose 
  !     Write active_region Bfield configuration
  ! Parameter list
  !    xActive : x coordinates
  !    Bmax: Magnetic field Bmax > Busr
  !    xplus,xminus: Constant, where the value of Bfield goes down from max
  !    xshift: shift in x for tan computation
  !    bplus,bminus: shift in Bfield for tan computation
  !    Byp,By0,Bym :  intermediary arrays for computations
  !    Hp,Hm: smooth approx of Heaviside step function
  !    hside : constant for Heaviside approx
  ! Called from:
  !     specialset_B0
    double precision :: xActive
    double precision :: Bmax,xplus,xminus,xshift,bplus,bminus,hside
    double precision :: Byp,By0,Bym,Hp,Hm 
 
    xplus=maintan*xprobmax1
    xminus=-xplus
    xshift=2.d0*xplus
    bplus=-Bmax*dtanh(xplus*parb)-(Bmax-Busr)/2.d0*dtanh(xminus*parb)
    bminus=-Bmax*dtanh(xminus*parb)-(Bmax-Busr)/2.d0*dtanh(xplus*parb)
    
    Byp=(Bmax-Busr)/2.d0*dtanh((xActive-xshift)*parb)+bplus
    By0=-Bmax*dtanh(parb*xActive)
    Bym=(Bmax-Busr)/2.d0*dtanh((xActive+xshift)*parb)+bminus
    
    hside=1.d-1
    Hp=1.d0/2.d0*(1.d0+dtanh(hside*(xActive-xplus))) 
    Hm=1.d0/2.d0*(1.d0+dtanh(hside*(-xActive+xminus))) 
  
    active_By=Byp*Hp+By0*(1.d0-Hp-Hm)+Bym*Hm 
    return 
    end function active_By

  double precision function dxBy(xActive,Bmax)
  ! Purpose 
  !     Write x derivative of By  for active region configuration
  ! Parameter list
  !    xActive : x coordinates
  !    Bmax: Magnetic field Bmax > Busr
  !    xplus,xminus: Constant, where the value of Bfield goes down from max
  !    xshift: shift in x for tan computation
  !    bplus,bminus: shift in Bfield for tan computation
  !    Byp,By0,Bym,Jzp,Jz0,Jzm,wJ1,wJ2,wJ3 :  intermediary arrays for computations
  !    Hp,Hm: smooth approx of Heaviside step function
  !    dHp,dHm: derivative of Heaviside step function
  !    hside: Constant for smooth approx of Heaviside step function
  ! Called from:
  !     specialset_J0
    double precision :: xActive
    double precision :: Bmax,xplus,xminus,xshift,bplus,bminus,hside
    double precision :: Byp,By0,Bym,Hp,Hm,dHm,dHp,Jzp,Jz0,Jzm,wJ1,wJ2,wJ3 
 
    xplus=maintan*xprobmax1
    xminus=-xplus
    xshift=2.d0*xplus
    bplus=-Bmax*dtanh(xplus*parb)-(Bmax-Busr)/2.d0*dtanh(xminus*parb)
    bminus=-Bmax*dtanh(xminus*parb)-(Bmax-Busr)/2.d0*dtanh(xplus*parb)
    
    Byp=(Bmax-Busr)/2.d0*dtanh((xActive-xshift)*parb)+bplus
    By0=-Bmax*dtanh(parb*xActive)
    Bym=(Bmax-Busr)/2.d0*dtanh((xActive+xshift)*parb)+bminus
    
    Jzp=(Bmax-Busr)*parb/2.d0/dcosh((xActive-xshift)*parb)**2
    Jz0=-Bmax*parb/dcosh(xActive*parb)**2
    Jzm=(Bmax-Busr)*parb/2.d0/dcosh((xActive+xshift)*parb)**2

    hside=1.d-1
    Hp=(1.d0/2.d0)*(1.d0+dtanh(hside*(xActive-xplus))) 
    Hm=(1.d0/2.d0)*(1.d0+dtanh(hside*(-xActive+xminus))) 
    dHp=hside/2.d0/dcosh((xActive-xplus)*hside)**2
    dHm=hside/2.d0/dcosh((-xActive+xminus)*hside)**2
    
    wJ1=Jzp*Hp+Byp*dHp
    wJ2=Jz0*(1.d0-Hp-Hm)-By0*(dHp+dHm)
    wJ3=Jzm*Hm+Bym*dHm

    dxBy=wJ1+wJ2+wJ3
  
    return 
    end function dxBy

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Purpose
  !    B-field is split into const background (wB0) + perturbation component
  !    Here the time-independent background magnetic field is set
  ! Parameter list
  !    ixI : inside/input coordinate indices (internal cells, excluded ghosts)
  !    ixO : outside/output coordinate indices (including internal + ghost cells)
  !    x   : coordinates
  !    wB0 : wB0[ixO,1:3] output background magnetic field
  !          x_component is zero, y_component is approx = -B_usr, except in the 
  !          central zone with parb width parameter, where it switches into z
  !    Following added by Maxime Dubart
  !    Bmax : value of the Bfield for asymmetric configuration
  !    kcont: Constant of integration for continuity 
  !    wlike: Constant to make the function similar to the regular symmetric configuration (1/w)
  ! Calls: none
  ! Called from:
  !   Main code (substitution for usr_set_B0)
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
    double precision :: kcont,wlike,Bmax
    integer :: ix^D

!    write(*,*) "begin specialset_B0",mype
    if (iprob==1) then
      wB0(ixO^S,1)=zero
      wB0(ixO^S,2)=Busr
      wB0(ixO^S,3)=zero
    else if (iprob>=2) then
      wB0(ixO^S,1)=zero
      Bmax=Busr*asy_coeff
      if (active_region) then
        do ix1=ixOmin1,ixOmax1
          do ix2=ixOmin2,ixOmax2
            wB0(ix^D,2)=active_By(x(ix^D,1),Bmax)
          enddo
        enddo
      else
        wlike=1.d0/parb/Busr

        ! loop over coordinates 
        do ix1=ixOmin1,ixOmax1
          do ix2=ixOmin2,ixOmax2
              if (x(ix^D,1)>=0.d0) then
                wB0(ix^D,2)=-Bmax*dtanh(x(ix^D,1)/(wlike*Bmax))
              else if (x(ix^D,1)<0.d0) then
                wB0(ix^D,2)=-Busr*dtanh(x(ix^D,1)/(wlike*Busr))
              endif
          enddo
        enddo
      endif
      kcont=1.01d0*Bmax**2
      wB0(ixO^S,3)=dsqrt(kcont-wB0(ixO^S,2)**2)
    endif

!    write(*,*) "end   specialset_B0",mype
  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Purpose
  !    current density is split into const background (wJ0) + perturbation
  !    Here the time-independent background current density is set
  ! Parameter list
  !    ixI : inside/input coordinate indices (internal cells, excluded ghosts)
  !    ixO : outside/output coordinate indices (including internal + ghost cells)
  !    x   : coordinates
  !    wJ0 : wJ0[ixO,1:3] output background current density
  !          x_component is zero, y_component is approx = parb*B_usr, except in 
  !          central zone with parb width parameter, where it switches into z
  !    Following added by Maxime Dubart
  !    Bmax : value of the Bfield for asymmetric configuration, Bmax >= Busr
  !    wlike: Constant to make the function similar to the regular symmetric configuration (1/w)
  ! Calls: none
  ! Called from:
  !   Main code (substitution for usr_set_J0)
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)
    double precision :: wlike,Bmax
    integer :: ix^D

!    write(*,*) "begin specialset_J0",mype
    if (iprob==1) then
      wJ0(ixO^S,1)=zero
      wJ0(ixO^S,2)=zero
      wJ0(ixO^S,3)=zero
    else if (iprob>=2) then
      Bmax=Busr*asy_coeff
      wJ0(ixO^S,1)=zero
      if (active_region) then
        do ix1=ixOmin1,ixOmax1
          do ix2=ixOmin2,ixOmax2
            wJ0(ix^D,3)=dxBy(x(ix^D,1),Bmax)
          enddo
        enddo
      else
        wlike=1.d0/parb/Busr
        do ix1=ixOmin1,ixOmax1
          do ix2=ixOmin2,ixOmax2
              if (x(ix^D,1)>=0.d0) then
                wJ0(ix^D,3)=-1.d0/(wlike*dcosh(x(ix^D,1)/(wlike*Bmax))**2)
              else if (x(ix^D,1)<0.d0) then
                wJ0(ix^D,3)=-1.d0/(wlike*dcosh(x(ix^D,1)/(wlike*Busr))**2)
              endif
          enddo
        enddo
      endif
      wJ0(ixO^S,2)=wJ0(ixO^S,3)*block%B0(ixO^S,2,0)/block%B0(ixO^S,3,0) 
    endif

!    write(*,*) "end   specialset_J0",mype
  end subroutine specialset_J0


  subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
  ! Purpose
  !    Set the common "eta" array for resistive MHD based on w or the
  !    "current" variable which has components between idirmin and 3.
  !
  ! Parameter list
  !    w	: Variables stored in this array
  !    ixI : inside/input coordinate indices (internal cells, excluded ghosts)
  !    ixO : outside/output coordinate indices (including internal + ghost cells)
  !    idirmin : minimum current value
  !    x   : coordinates
  !    current : the current calculated elsewhere (get_current)
  !    eta : output eta format as per Ruan 2020
  !             Eq 10 if t<t_acc, Eq 11 if t>t_acc
  ! Calls: none
  !
  ! Called from:
  !   Main code (substitution for usr_special_resistivity)
    integer, intent(in) :: ixI^L, ixO^L, idirmin
    double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
    double precision :: rad(ixI^S) !,reta,reta2
    double precision :: jc,jabs(ixI^S),vd(ixI^S),rad_acc(ixI^S)

!    write(*,*) "begin special_eta",mype

    ! iprob is user integer switch, =1 in initial version of code, now is 3
    if (iprob==1) then
      eta(ixO^S)=0.d0

    ! for later versions of the code
    else if (iprob>=2) then
      !reta = 0.8d0 * 0.3d0    ! r_eta (radius of special eta: eq 10 of Ruan 2020)
      !eta1 = 0.01d0
      !tar1= 0.0d0
      ! Flare precursor phase (use Ruan 2020 eq10)
      if (global_time<t_imp) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-h_eta_pre)**2)
        where (rad(ixO^S) .lt. r_eta_pre)
          eta(ixO^S)=eta0_pre*(2.d0*(rad(ixO^S)/r_eta_pre)**3-3.d0*(rad(ixO^S)/r_eta_pre)**2+1.d0)      ! Here is eq 10 from Ruan 2020
        elsewhere
          eta(ixO^S)=zero
        endwhere
      ! Flare impuslive phase
      elseif (global_time<t_acc) then
        rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-h_eta_imp)**2)
        where (rad(ixO^S) .lt. r_eta_imp)
          eta(ixO^S)=eta0_imp*(2.d0*(rad(ixO^S)/r_eta_imp)**3-3.d0*(rad(ixO^S)/r_eta_imp)**2+1.d0)      ! Here is eq 10 from Ruan 2020
        elsewhere
          eta(ixO^S)=zero
        endwhere
      ! Flare particle acceleration phase (use Ruan 2020 eq11)
      elseif (global_time<t_decay) then
        vd(ixO^S)=dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))/w(ixO^S,rho_)/q_e
        rad_acc(ixO^S)=exp(-((x(ixO^S,2)-h_eta_acc)**2)/h_s_acc**2)
        if (eta_decay_acc) then
          where(vd(ixO^S)>v_c_acc)
            eta(ixO^S)=alpha_acc*(vd(ixO^S)/v_c_acc-1.d0)*rad_acc(ixO^S)    ! here is eq 11 with spatial threshold
          elsewhere
            eta(ixO^S)=0.d0
          endwhere
        else
          where(vd(ixO^S)>v_c_acc)
            eta(ixO^S)=alpha_acc*(vd(ixO^S)/v_c_acc-1.d0)    ! here is eq 11 without spatial threshold
          elsewhere
            eta(ixO^S)=0.d0
          endwhere
        endif
        where(eta(ixO^S)>eta_max_acc)
          eta(ixO^S)=eta_max_acc    ! apply maximum eta condition in eq 11 (applied via "minimum of" in the paper.
        endwhere
        !eta(ixO^S)=eta0
      ! decay phase
      else
        eta(ixO^S)=0.d0
      end if

    endif

!    write(*,*) "end2  special_eta",mype
  end subroutine special_eta


  subroutine usrspecial_convert(qunitconvert)
    ! Purpose:
    !   Used here to create an output text file with values specified
    !     Here only minmax values of atmospheric parameters in the flare region
    !     But can be altered for extra outputs
    ! Parameter list
    !   qunitconvert: note sure what this does, I dont use it
    ! Calls:
    !   | -> find_minmax_values()
    !   |    | -> find_extreme_value_box()
    !   |    |    | ->   find_extreme_values_grid()
    ! Called from:
    !   Main code (Substitute for "usr_special_convert")
    integer, intent(in) :: qunitconvert
!    write(*,*) "begin usrspecial_convert",mype
    call find_minmax_values()
!    write(*,*) "end   usrspecial_convert",mype

  end subroutine usrspecial_convert


  subroutine find_minmax_values()
    ! Purpose:
    !     finds minmax values of atmospheric parameters in regions specified
    !     
    ! Parameter list
    !   qunitconvert: note sure what this does, I dont use it
    !
    ! Calls:
    !   |    | -> find_extreme_value_box()
    !   |    |    | ->   find_extreme_values_grid()
    ! Called from:
    !   | -> find_minmax_values()
    double precision :: xboundmin,xboundmax,yboundmin,yboundmax
    character(5) :: filebit

!    write(*,*) "begin find_minmax_values",mype
    call find_extreme_value_box()

    xboundmin=-5.d-2
    xboundmax=5.d-2
    yboundmin=3.d-1
    yboundmax=5.d-1
    filebit="outfl"
    call find_extreme_value_rect(xboundmin,xboundmax,yboundmin,yboundmax,filebit)

    xboundmin=-3.d-1
    xboundmax=3.d-1
    yboundmin=4.d-2
    yboundmax=8.d-2
    filebit="foot2"
    call find_extreme_value_rect(xboundmin,xboundmax,yboundmin,yboundmax,filebit)
  
!    write(*,*) "end   find_minmax_values",mype
  end subroutine find_minmax_values


  subroutine find_extreme_value_rect(xboundmin,xboundmax,yboundmin,yboundmax,filebit)
    ! Purpose:
    !     finds minmax values of atmospheric parameters in a rectangular area
    !     
    ! Parameter list
    !   qunitconvert: note sure what this does, I dont use it
    !
    ! Calls:
    !   |    |    | ->   find_extreme_values_grid()
    ! Called from:
    !   | -> find_minmax_values()
    !   |    | -> find_extreme_value_box()
    use mod_global_parameters

    double precision :: dxb^D,xb^L
    integer :: iigrid,igrid,j
    integer :: ixO^L,ixI^L,ix^D
    character(5) :: filebit

    integer :: nlog=7
    double precision :: temp
    double precision, dimension (:), allocatable :: datalog, datalogtemp
    double precision, dimension (:,:), allocatable :: loclog,loclogtemp

    integer :: nrstart, loopval, filenr_field
    integer :: varno, varnote, varnorho, varnovy, varnorad
    character(100) :: fname,varnames

    double precision :: xboundmin,xboundmax,yboundmin,yboundmax, xbmin, xbmax

    ! select box for checking them
    xbmin=xprobmin1
    xbmax=xprobmax1

    allocate(datalog(nlog))
    allocate(datalogtemp(nlog))
    allocate(loclog(nlog,ndim))
    allocate(loclogtemp(nlog,ndim))

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;
    loclogtemp=huge(0.d0)
    loclog=huge(0.d0)

    ! datalogtemp(1:2) are temperature min max
    varnote=1
    datalogtemp(varnote)=huge(0.d0)
    datalogtemp(varnote+1)=-huge(0.d0)
    ! datalogtemp(3:4) are density min max
    varnorho=3
    datalogtemp(varnorho)=huge(0.d0)
    datalogtemp(varnorho+1)=-huge(0.d0)
    ! datalogtemp(5:6) are v_y min max
    varnovy=5
    datalogtemp(varnovy)=huge(0.d0)
    datalogtemp(varnovy+1)=-huge(0.d0)
    ! datalogtemp(7) is radiative losses max
    varnorad=7
    datalogtemp(varnorad)=-huge(0.d0)

    ! find the extreme values in the processor
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>ps(igrid)
      call find_extreme_values_gridrect(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,datalogtemp,loclogtemp,xboundmin,xboundmax,yboundmin,yboundmax)
    enddo
    !write(*,*) 1,mype,datalogtemp(varnorad),loclogtemp(varnorad,:)

    ! compare the values across processors and mpi communicate the min/maxvalues
    ! Temperature
    varno=varnote
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Density
    varno=varnorho
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Velocity_y
    varno=varnovy
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Radiation
    varno=varnorad
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno)=temp
    !write(*,*) 2,mype,datalog(varnorad),loclogtemp(varnorad,:)

    ! compare values with global min max. max out locations if not the extreme value
    ! Then use mpi_allreduce(min to get the locations of the global extremes.
    do loopval=1, nlog
       if(datalog(loopval) .NE. datalogtemp(loopval)) then
           loclogtemp(loopval,:)=huge(0.0d0)
       endif
    enddo
    !write(*,*) 3,mype,datalog(varnorad),loclogtemp(varnorad,:)

    ! compare the location values across processors and mpi communicate the min/maxvalues
    ! only the globals should be less than "huge"
    do varno=1,nlog
      do loopval=1, ndim
        call MPI_ALLREDUCE(loclogtemp(varno, loopval),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
        loclog(varno, loopval)=temp
      enddo
    enddo
    !write(*,*) 4,mype,datalog(varnorad),loclog(varnorad,:)

    nrstart=20
    
    filenr_field=filenr-1
    if (filenr_field .LT. 0) filenr_field=0
    write(fname, '(a,i4.4,a)') trim(base_filename),filenr_field,'_'//filebit//'_extremevalues_.txt'
    varnames='filenr qt Tmin Tmax rhomin rhomax vmin, vmax, radmax, coords'
!    write(*,*) global_time, nlog, datalog(0), loclog(0,0)
    if (mype==0) call usr_write_log(global_time,nlog,datalog,loclog,nrstart,fname,varnames)

!    deallocate(datalog,datalogtemp,loclog,loclogtemp)

  end subroutine find_extreme_value_rect


  subroutine find_extreme_value_box()
    ! Purpose:
    !     finds minmax values of atmospheric parameters in a rectangular area
    !     
    ! Parameter list
    !   qunitconvert: note sure what this does, I dont use it
    !
    ! Calls:
    !   |    |    | ->   find_extreme_values_grid()
    ! Called from:
    !   | -> find_minmax_values()
    !   |    | -> find_extreme_value_box()
    use mod_global_parameters

    character(5) :: filebit
    double precision :: dxb^D,xb^L
    integer :: iigrid,igrid,j
    integer :: ixO^L,ixI^L,ix^D

    integer :: nlog=7
    double precision :: temp
    double precision, dimension (:), allocatable :: datalog, datalogtemp
    double precision, dimension (:,:), allocatable :: loclog,loclogtemp

    integer :: nrstart, loopval, filenr_field
    integer :: varno, varnote, varnorho, varnovy, varnorad
    character(100) :: fname,varnames

    double precision :: xboundmin,xboundmax,yboundmin,yboundmax, xbmin, xbmax

    ! select box for checking them
    xbmin=xprobmin1
    xbmax=xprobmax1
    xboundmin=-3.d-1
    xboundmax=3.d-1
    yboundmin=0.04
    yboundmax=5.d-1

    allocate(datalog(nlog))
    allocate(datalogtemp(nlog))
    allocate(loclog(nlog,ndim))
    allocate(loclogtemp(nlog,ndim))

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;
    loclogtemp=huge(0.d0)
    loclog=huge(0.d0)

    ! datalogtemp(1:2) are temperature min max
    varnote=1
    datalogtemp(varnote)=huge(0.d0)
    datalogtemp(varnote+1)=-huge(0.d0)
    ! datalogtemp(3:4) are density min max
    varnorho=3
    datalogtemp(varnorho)=huge(0.d0)
    datalogtemp(varnorho+1)=-huge(0.d0)
    ! datalogtemp(5:6) are v_y min max
    varnovy=5
    datalogtemp(varnovy)=huge(0.d0)
    datalogtemp(varnovy+1)=-huge(0.d0)
    ! datalogtemp(7) is radiative losses max
    varnorad=7
    datalogtemp(varnorad)=-huge(0.d0)

    ! find the extreme values in the processor
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>ps(igrid)
      call find_extreme_values_grid(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,datalogtemp,loclogtemp)
    enddo
    !write(*,*) 1,mype,datalogtemp(varnorad),loclogtemp(varnorad,:)

    ! compare the values across processors and mpi communicate the min/maxvalues
    ! Temperature
    varno=varnote
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Density
    varno=varnorho
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Velocity_y
    varno=varnovy
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
    datalog(varno)=temp
    call MPI_ALLREDUCE(datalogtemp(varno+1),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno+1)=temp
    ! Radiation
    varno=varnorad
    call MPI_ALLREDUCE(datalogtemp(varno),temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)    
    datalog(varno)=temp
    !write(*,*) 2,mype,datalog(varnorad),loclogtemp(varnorad,:)


    ! compare values with global min max. max out locations if not the extreme value
    ! Then use mpi_allreduce(min to get the locations of the global extremes.
    do loopval=1, nlog
       if(datalog(loopval) .NE. datalogtemp(loopval)) then
           loclogtemp(loopval,:)=huge(0.0d0)
       endif
    enddo
    !write(*,*) 3,mype,datalog(varnorad),loclogtemp(varnorad,:)


    ! compare the location values across processors and mpi communicate the min/maxvalues
    ! only the globals should be less than "huge"
    do varno=1,nlog
      do loopval=1, ndim
        call MPI_ALLREDUCE(loclogtemp(varno, loopval),temp,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)    
        loclog(varno, loopval)=temp
      enddo
    enddo
    !write(*,*) 4,mype,datalog(varnorad),loclog(varnorad,:)


    nrstart=20
    
    filenr_field=filenr-1
    if (filenr_field .LT. 0) filenr_field=0
    write(fname, '(a,i4.4,a)') trim(base_filename),filenr_field,'_extremevalues_.txt'
    varnames='filenr qt Tmin Tmax rhomin rhomax vmin, vmax, radmax, coords'
!    write(*,*) global_time, nlog, datalog(0), loclog(0,0)
    if (mype==0) call usr_write_log(global_time,nlog,datalog,loclog,nrstart,fname,varnames)

!    deallocate(datalog,datalogtemp,loclog,loclogtemp)

  end subroutine find_extreme_value_box


  subroutine find_extreme_values_gridrect(ixI^L,ixO^L,w,x,datalogtemp,loclogtemp,xboundmin,xboundmax,yboundmin,yboundmax)
    ! Purpose:
    !     Block loop that searches for extreme values
    !     These values should be initialised appropriately outside the function
    !     Then you should call the function via a loop that runs over each block in the processor
    !     being careful not to overwrite any values found earlier in your loop over the proc blocks.
    !     
    ! Parameter list
    !
    !
    ! Calls: None
    !
    ! Called from:
    !   | -> find_minmax_values()
    !   |    | -> find_extreme_value_box()
    !   |    |    | ->   find_extreme_values_grid()
    use mod_global_parameters
    use mod_radiative_cooling

    integer :: nlog=7
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: pth(ixI^S),Te(ixI^S)
    double precision, intent(inout) :: datalogtemp(:)
    double precision, intent(inout) :: loclogtemp(:,:)
    integer :: varno, varnote, varnorho, varnovy, varnorad

    integer :: ix^D,ixb, loopval
    double precision :: radgrid(ixI^S)
    double precision :: xboundmin,xboundmax,yboundmin,yboundmax, xbmin, xbmax

    ! select spatial region to check the extreme values

    varnote=1
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    varnorho=3
    varnovy=5
    varnorad=7
    call getvar_cooling(ixI^L,ixO^L,w,x,radgrid,rc_fl)
    call mhd_to_primitive(ixI^L,ixO^L,w,x)

    ! Perform a loop over the grid points in the block searching for extreme values
    do ix1=ixOmin1,ixOmax1
      do ix2=ixOmin2,ixOmax2
          ! Put a spatial restriction on the extreme value regions
          if (x(ix^D,1)>xboundmin .and. x(ix^D,1)<xboundmax) then
              if (x(ix^D,2)>yboundmin .and. x(ix^D,2)<yboundmax) then
                  ! if values are outside those found in this and previous blocks, update
                  if (Te(ix^D)<datalogtemp(varnote)) then 
                      datalogtemp(varnote)=Te(ix^D)
                      loclogtemp(varnote,:)=x(ix^D,:)
                  endif
                  if (Te(ix^D)>datalogtemp(varnote+1)) then
                      datalogtemp(varnote+1)=Te(ix^D)
                      loclogtemp(varnote+1,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,rho_)<datalogtemp(varnorho)) then 
                      datalogtemp(varnorho)=w(ix^D,rho_)
                      loclogtemp(varnorho,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,rho_)>datalogtemp(varnorho+1)) then
                      datalogtemp(varnorho+1)=w(ix^D,rho_)
                      loclogtemp(varnorho+1,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,mom(2))<datalogtemp(varnovy)) then 
                      datalogtemp(varnovy)=w(ix^D,mom(2))
                      loclogtemp(varnovy,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,mom(2))>datalogtemp(varnovy+1)) then
                      datalogtemp(varnovy+1)=w(ix^D,mom(2))
                      loclogtemp(varnovy+1,:)=x(ix^D,:)
                  endif
              endif
          endif
      ! find global maximum of the cooling
      if ( radgrid(ix^D) > datalogtemp(varnorad) ) then
        datalogtemp(varnorad)=radgrid(ix^D)
        loclogtemp(varnorad,:)=x(ix^D,:)
      endif
      enddo
    enddo
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine find_extreme_values_gridrect


  subroutine find_extreme_values_grid(ixI^L,ixO^L,w,x,datalogtemp,loclogtemp)
    ! Purpose:
    !     Block loop that searches for extreme values
    !     These values should be initialised appropriately outside the function
    !     Then you should call the function via a loop that runs over each block in the processor
    !     being careful not to overwrite any values found earlier in your loop over the proc blocks.
    !     
    ! Parameter list
    !
    !
    ! Calls: None
    !
    ! Called from:
    !   | -> find_minmax_values()
    !   |    | -> find_extreme_value_box()
    !   |    |    | ->   find_extreme_values_grid()
    use mod_global_parameters
    use mod_radiative_cooling

    integer :: nlog=7
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision :: pth(ixI^S),Te(ixI^S)
    double precision, intent(inout) :: datalogtemp(:)
    double precision, intent(inout) :: loclogtemp(:,:)
    integer :: varno, varnote, varnorho, varnovy, varnorad

    integer :: ix^D,ixb, loopval
    double precision :: radgrid(ixI^S)
    double precision :: xboundmin,xboundmax,yboundmin,yboundmax, xbmin, xbmax

    ! select spatial region to check the extreme values
    xboundmin=-3.d-1
    xboundmax=3.d-1
    yboundmin=0.04
    yboundmax=5.d-1

    varnote=1
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    varnorho=3
    varnovy=5
    varnorad=7
    call getvar_cooling(ixI^L,ixO^L,w,x,radgrid,rc_fl)
    call mhd_to_primitive(ixI^L,ixO^L,w,x)

    ! Perform a loop over the grid points in the block searching for extreme values
    do ix1=ixOmin1,ixOmax1
      do ix2=ixOmin2,ixOmax2
          ! Put a spatial restriction on the extreme value regions
          if (x(ix^D,1)>xboundmin .and. x(ix^D,1)<xboundmax) then
              if (x(ix^D,2)>yboundmin .and. x(ix^D,2)<yboundmax) then
                  ! if values are outside those found in this and previous blocks, update
                  if (Te(ix^D)<datalogtemp(varnote)) then 
                      datalogtemp(varnote)=Te(ix^D)
                      loclogtemp(varnote,:)=x(ix^D,:)
                  endif
                  if (Te(ix^D)>datalogtemp(varnote+1)) then
                      datalogtemp(varnote+1)=Te(ix^D)
                      loclogtemp(varnote+1,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,rho_)<datalogtemp(varnorho)) then 
                      datalogtemp(varnorho)=w(ix^D,rho_)
                      loclogtemp(varnorho,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,rho_)>datalogtemp(varnorho+1)) then
                      datalogtemp(varnorho+1)=w(ix^D,rho_)
                      loclogtemp(varnorho+1,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,mom(2))<datalogtemp(varnovy)) then 
                      datalogtemp(varnovy)=w(ix^D,mom(2))
                      loclogtemp(varnovy,:)=x(ix^D,:)
                  endif
                  if (w(ix^D,mom(2))>datalogtemp(varnovy+1)) then
                      datalogtemp(varnovy+1)=w(ix^D,mom(2))
                      loclogtemp(varnovy+1,:)=x(ix^D,:)
                  endif
              endif
          endif
      ! find global maximum of the cooling
      if ( radgrid(ix^D) > datalogtemp(varnorad) ) then
        datalogtemp(varnorad)=radgrid(ix^D)
        loclogtemp(varnorad,:)=x(ix^D,:)
      endif
      enddo
    enddo
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine find_extreme_values_grid

  subroutine usr_write_log(qt,nlog,datalog,loclog,nrstart,fname,varnames)
    use mod_global_parameters

    double precision :: qt
    integer :: nlog,nrstart
    double precision :: datalog(nlog), loclog(nlog,ndim)
    character(*) :: fname,varnames

    integer, parameter :: outunit=35
    integer :: i,j

!    write(*,*) "begin usr_write_log",mype
    if (mype==0) then
!      if (snapshotini==nrstart) then
        open(unit=outunit,file=trim(fname),status='replace')
        write(outunit, '(a)') trim(varnames)
!      else
!        open(unit=outunit,file=trim(fname),status='old',position='append')
!      endif

!      write(outunit,'(i8)', advance='no') snapshotini
!      write(outunit,'(e15.7)', advance='no') qt
      write(outunit,'(e16.7)') qt

      if (nlog>1) then
        do j=1,nlog-1
          write(outunit,'(e16.7)', advance='no') datalog(j)
          do i=1,ndim-1
            write(outunit,'(e16.7)', advance='no') loclog(j,i)
          enddo
          write(outunit,'(e16.7)', advance='yes') loclog(j,ndim)
        enddo
        write(outunit,'(e16.7)', advance='no') datalog(nlog)
        do i=1,ndim-1
          write(outunit,'(e16.7)', advance='no') loclog(nlog,i)
        enddo
        write(outunit,'(e16.7)', advance='no') loclog(nlog,ndim)
      endif

    endif

!    write(*,*) "end   usr_write_log",mype
  end subroutine usr_write_log


  subroutine special_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
    use mod_global_parameters

    integer, intent(in)                 :: igrid,ip,numP,nwP,nwL
    double precision, intent(in)        :: xf(numP,ndim)
    double precision, intent(inout)     :: wP(numP,nwP),wL(1+nwL)
    double precision, intent(in)        :: dL
    logical, intent(in)                 :: forward
    character(len=std_len), intent(in)  :: ftype,tcondi

    integer :: ixI^L,ixO^L
    double precision :: dxb^D
    double precision :: xpp(1:ndim),wpp(1:nwP)

!    write(*,*) "begin special_field_w",mype
    if (tcondi/='none') then
      if (tcondi=='cons+j') then
        xpp(1:ndim)=xf(ip,1:ndim)
        ^D&ixImin^D=ixglo^D;
        ^D&ixImax^D=ixghi^D;
        ^D&dxb^D=rnode(rpdx^D_,igrid);
        ! indeces of nearby cells
        ^D&ixOmin^D=floor((xpp(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
        ^D&ixOmax^D=^D&ixOmin^D+1;
        if (nwP/=nw+3) then
          call MPISTOP('Wrong number of variables in tracing field')
        else
          call get_w_local(igrid,xpp,wpp,ps(igrid)%x,ps(igrid)%w,ixI^L,ixO^L,dxb^D)
          wP(ip,1:nwP)=wpp(1:nwP)
        endif
      endif
    endif

!    write(*,*) "end   special_field_w",mype
  end subroutine special_field_w


  subroutine get_w_local(igrid,xpp,wpp,x,w,ixI^L,ixO^L,dxb^D)
    use mod_geometry

    integer, intent(in) :: igrid,ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,ndim),w(ixI^S,nw),xpp(ndim)
    double precision, intent(in) :: dxb^D
    double precision, intent(inout) :: wpp(nw+ndir)

    integer          :: ix^D,j,idirmin,idirmin0
    double precision :: xd^D
    double precision :: factor(ixO^S),wBnear(ixO^S,nw+ndir),J1(ixI^S,ndir)

    ^D&xd^D=(xpp(^D)-x(ixOmin^DD,^D))/dxb^D\

    ! factors for linear interpolation
    wBnear=zero
    {do ix^D=ixOmin^D,ixOmax^D\}
      factor(ix^D)={abs(ixOmax^D-ix^D-xd^D)*}
    {enddo\}

    ! conserved variables for interpolation
    wBnear(ixO^S,1:nw)=w(ixO^S,1:nw)
    if (B0field) then
      do j=1,3
        wBnear(ixO^S,mag(j))=ps(igrid)%B0(ixO^S,j,0)
      enddo
    endif

    ! current for interpolation
    idirmin0=7-2*ndir
    call curlvector(w(ixI^S,mag(1:ndir)),ixI^L,ixO^L,J1,idirmin,idirmin0,ndir)
    wBnear(ixO^S,nw+1:nw+ndir)=J1(ixO^S,1:ndir)
    if (B0field) wBnear(ixO^S,nw+1:nw+ndir)=wBnear(ixO^S,nw+1:nw+ndir)+ps(igrid)%J0(ixO^S,1:ndir)

    ! interpolation
    wpp=0.d0
    {do ix^D=ixOmin^D,ixOmax^D\}
      do j=1,nw+ndir
        wpp(j)=wpp(j)+wBnear(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

  end subroutine get_w_local


!  subroutine special_boundary_sponge(level,qt,ixI^L,ixO^L,w,x)
    ! boundindex: array [2*ndim]
    !             specifies which boundaries should be sponge damped 
    !             first two entries are for lower and upper boundaries of first dimension 
    !             3rd&4th are for second dimension lower and upper, and so on
    !             index value: <= 0         -> no damping
    !             index value: 0< value <1  -> craziness
    !             index value: >= 1         -> standard exponential or higher exponential damping, e.g. 2 for exponent^2
    ! boundsize:  array [2*ndim]
    !             specifies sizes of sponge zones (xi_start) in you experiments length units
    ! boundcoeff: array [2*ndim]
    !             specifies coefficients (a) for the damping
    !
    ! formula     u''=u exp ( -a (xi-xi_start)/()**index )
!    use mod_global_parameters
!    integer, intent(in) :: ixI^L,ixO^L,level
!    double precision, intent(in) :: qt
!    double precision, intent(inout) :: w(ixI^S,1:nw)
!    double precision, intent(in) :: x(ixI^S,1:ndim)
!
!    double precision :: damploc
!    double precision :: boundindex(1:2*ndim),boundsize(1:2*ndim),boundcoeff(1:2*ndim)
!    integer :: iloop, ix1, ix2, na
!
!    double precision :: p0(ixO^S),rho0(ixO^S)
!    double precision :: res
!
!    ! Sponge parameters
!    boundindex=[1.0d0,0.0d0,0.0d0,0.0d0]
!    boundsize=[dble(xprobmax1-xprobmin1)/dble(block_nx1),0.0d0,0.0d0,0.0d0]
!    boundcoeff=[1.0d-4,0.0d0,0.0d0,0.0d0]
!    
!    ! splits up the etot into its components so that you can call pressure or magnetic energy or temperature
!    ! the you need to revert back to etot using mhd_to_conserved afterwards.
!    call mhd_to_primitive(ixI^L,ixO^L,w,x)
!
!    {do ix^DB=ixOmin^DB,ixOmax^DB\}
!        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
!        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
!        rho0(ix^D)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
!        p0(ix^D)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
!    {end do\}
!
!    ! xdimension min
!    iloop=1
!    if (boundindex(iloop) .ge. 1) then
!      damploc=xprobmin1+boundsize(iloop)
!      where ( x(ixO^S,1) .lt. damploc )
!      ! Bracket tracker
!      !  a         b ba  c         d dc e   f           g     g hi         j       ji          k     kh            l     l fe
!!        w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!
!!        w(ixO^S,mag(1))=w(ixO^S,mag(1))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mag(2))=w(ixO^S,mag(2))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mag(3))=w(ixO^S,mag(3))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!
!        w(ixO^S,rho_)=rho0 + (w(ixO^S,rho_)-rho0) *(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,p_)=p0 + (w(ixO^S,p_)-p0) *(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,1))/boundsize(iloop))**boundindex(iloop) ))
!
!      end where
!    endif
!
!    ! xdimension max
!    iloop=2
!    if (boundindex(iloop) .ge. 1) then
!      damploc=xprobmax1-boundsize(iloop)
!      where ( x(ixO^S,1) .gt. damploc )
!
!      ! Bracket tracker
!      !  a           b ba  c           d dc e   f           g     g hi j       j        i          k     kh            l     l fe
!!        w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!
!!        w(ixO^S,mag(1))=w(ixO^S,mag(1))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mag(2))=w(ixO^S,mag(2))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,mag(3))=w(ixO^S,mag(3))*(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!
!        w(ixO^S,rho_)=rho0 + (w(ixO^S,rho_)-rho0) *(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        w(ixO^S,p_)=p0 + (w(ixO^S,p_)-p0) *(exp(-boundcoeff(iloop)*((x(ixO^S,1)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!
!      end where
!    endif
!
!    if (ndim .gt. 1) then
!      ! ydimension min
!      iloop=3
!      if (boundindex(iloop) .ge. 1) then
!        damploc=xprobmin2+boundsize(iloop)
!        where ( x(ixO^S,2) .lt. damploc )
!
!        ! Bracket tracker
!        !  a           b ba  c           d dc e   f           g     g hi         j       ji          k     kh            l     l fe
!!          w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!
!!          w(ixO^S,mag(1))=w(ixO^S,mag(1))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mag(2))=w(ixO^S,mag(2))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mag(3))=w(ixO^S,mag(3))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
! 
!          w(ixO^S,rho_)=rho0 + (w(ixO^S,rho_)-rho0) *(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,p_)=p0 + (w(ixO^S,p_)-p0) *(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,2))/boundsize(iloop))**boundindex(iloop) ))
!
!        end where
!      endif
!    
!      ! ydimension max
!      iloop=4
!      if (boundindex(iloop) .ge. 1) then
!        damploc=xprobmax2-boundsize(iloop)
!        where ( x(ixO^S,2) .gt. damploc )
!
!        ! Bracket tracker
!        !  a           b ba  c           d dc e   f           g     g hi j       j        i          k     kh            l     l fe
!!          w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!
!!          w(ixO^S,mag(1))=w(ixO^S,mag(1))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mag(2))=w(ixO^S,mag(2))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mag(3))=w(ixO^S,mag(3))*(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
! 
!          w(ixO^S,rho_)=rho0 + (w(ixO^S,rho_)-rho0) *(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,p_)=p0 + (w(ixO^S,p_)-p0) *(exp(-boundcoeff(iloop)*((x(ixO^S,2)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!
!        end where
!      endif
!    endif
!
!!    if (ndim .gt. 2) then
!!      ! ydimension min
!!      iloop=5
!!      if (boundindex(iloop) .ge. 1) then
!!        damploc=xprobmin3+boundsize(iloop)
!!        where ( x(ixO^S,3) .lt. damploc )
!!        ! Bracket tracker
!!        !  a           b ba  c           d dc e   f           g     g hi         j       ji          k     kh            l     l fe
!!          w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,3))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,3))/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((damploc-x(ixO^S,3))/boundsize(iloop))**boundindex(iloop) ))
!!        end where
!!      endif
!!   
!!      ! ydimension max
!!      iloop=6
!!      if (boundindex(iloop) .ge. 1) then
!!        damploc=xprobmax3-boundsize(iloop)
!!        where ( x(ixO^S,3) .gt. damploc )
!!        ! Bracket tracker
!!        !  a           b ba  c           d dc e   f           g     g hi j       j        i          k     kh            l     l fe
!!          w(ixO^S,mom(1))=w(ixO^S,mom(1))*(exp(-boundcoeff(iloop)*((x(ixO^S,3)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(2))=w(ixO^S,mom(2))*(exp(-boundcoeff(iloop)*((x(ixO^S,3)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!          w(ixO^S,mom(3))=w(ixO^S,mom(3))*(exp(-boundcoeff(iloop)*((x(ixO^S,3)-damploc)/boundsize(iloop))**boundindex(iloop) ))
!!        end where
!!      endif
!!    endif
!
!    call mhd_to_conserved(ixI^L,ixO^L,w,x)
!
!  end subroutine special_boundary_sponge


end module mod_usr
