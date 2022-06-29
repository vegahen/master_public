!==============================================================================!
!==============================================================================!
module user_variables
  implicit none
  save

  logical, parameter ::                 & ! read grid from .txt instead of
    readagnfit = .false.                  ! computing it 

  integer, parameter ::                 &
    n_sets = 1,                         & ! number of MC set
    n_start = 1,                        & ! injected particles/set
    iseed_shift = 0                       ! positive shift of random seed 

  double precision, parameter ::        &
    z_initial = 3d0,                    & ! z / R_S
    z_max = 1d16,                       & ! z / cm
    electron_energythr = 1d8,           & ! threshold for keeping electrons
    photon_energythr = 1d4                ! threshold for keeping photons

  character*60 ::                       &
    filename = '_test'                    ! name in output
  
end module user_variables
!==============================================================================!
!==============================================================================!
module constants
  implicit none
  save

  double precision, parameter ::        & 
    pi = 3.1415926536d0,                &
    two_pi = 2.d0*pi,                   &
    alpha2k2on2pihbarc = 3.189408141d-9,& ! eV / cm K^2
    twome2c4onk = 6.060341847d15,       & ! eV K
    gon4pisigma = 9.366646346d-5,       & ! cm^3 s K^4 / g^2

    honk = 4.799243073d-11,             & ! K Hz^-1
    twome2c4onh = 1.262770348d26,       & ! eV Hz
    alpha2honc = 7.346078957d-30,       & ! eV cm^-1 Hz^-2
    twohonme2c4 = 3.167638522d-26,      & ! eV^-1 Hz^-1
    eightalpha2h2onthreeme2c5 = 3.102629692d-55,  & ! cm^-1 Hz^-3
    fourthgonfourpisigma = 9.837755636d-2,  & ! K (s cm^3 g^-2)^.25
    wien_constant = 5.878925757d10,     & ! Hz K^-1
    twogonc2 = 1.485232054d-28,         & ! cm g^-1
    me2c4 = 2.611199269d11,             & ! eV^2

    twoh = 8.271335392d-15,             & ! eV Hz^-1
    ame = 5.11d5                          ! electron mass

  double precision zero, infty            ! initialized in init_general

end module constants
!==============================================================================!
!==============================================================================!
module heap_agn
  implicit none
  save

  type event
    integer icq, id
    double precision z, e, w
  end type event

  integer n, s, n_injected, n_created, n_out, n_ignored
  type(event), allocatable, dimension(:), target :: events

end module heap_agn
!==============================================================================!
!==============================================================================!
module particle
  implicit none
  save
  
  integer, parameter ::                 &
    pid_max = 15,                       & ! maximal pid 
    n_stable = 11                         ! number of stable particles
    
  integer stable_pid(n_stable)
  
  data stable_pid                       &
    /-8,-7,-5,-4,-1,0,1,4,5,7,8/          !...gamma,e,nue,numu,p

end module particle
!==============================================================================!
!==============================================================================!
module internal
  use particle
  implicit none
  save

  integer iseed

  double precision, parameter ::        & ! all energies in eV
    E_min = 1d4,                        & ! minimal energy for bining
    dn = 0.05d0

  integer n_tot

  double precision d_f

  integer, parameter :: debug = 0 ! not used, needed for ELMAG functions

end module internal
!==============================================================================!
!==============================================================================!
module result
  use particle
  implicit none
  save

  integer, parameter ::                 &
    n_enbin = 140                         ! ten decades starting in E_min

  double precision                      &
    En_f(-pid_max:pid_max,n_enbin),     &
    En_f_tot(-pid_max:pid_max,n_enbin) 

end module result
!==============================================================================!
!==============================================================================!
module agn_fit
  implicit none
  save

  integer, parameter :: n_zpair = 121, n_epair = 221, zmin_pair = 12, &
    emin_pair = 4, d_zpair = 20, d_epair = 20, n_zics = 121, n_eics = 181, &
    zmin_ics = 12, emin_ics = 6, d_zics = 20, d_eics = 20

  double precision :: heights_pair(n_zpair), energies_pair(n_epair), &
    rates_pair(n_epair, n_zpair), heights_ics(n_zics), energies_ics(n_eics), &
    rates_ics(n_eics, n_zics)

  character*60 :: agn_pair_filename = 'agn_fit/v4_pair_81x221.txt'
  character*60 :: agn_ics_filename = 'agn_fit/v4_ics_81x181.txt'

end module agn_fit
!==============================================================================!
!==============================================================================!
module agn_data
  use constants, only : twogonc2, alpha2k2on2pihbarc, twome2c4onk, gon4pisigma
  implicit none
  save

  integer, parameter :: n_x = 100, n_mu = 100

  double precision, parameter ::        &
    bh_m = 2d42,                        &
    bh_mdot = 1d27,                     &
    r_s = twogonc2 * bh_m,              &
    bh_disksize = 100d0 * r_s

end module agn_data
!==============================================================================!
!==============================================================================!
module propagation_params
  use user_variables, only : z_max
  implicit none

  double precision, parameter ::        &
    incr = 10**(log10(2d0)/4d0),        & ! normal increase in dz
    decr = 2d0,                         & ! normal decrease in dz
    incr_more = 2d0,                    & ! increase in dz when prob very low
    incr_lots = 9d0,                    & ! dz / z when rate is zero
    tol_prob = 1d-2,                    & ! max prob of interaction on dz
    tol_rise = log10(1.01d0),           & ! 1/1.01 < rate(z2)/rate(z1) < 1.01
    tol_prob_highrise = 1d-3,           & ! max prob in dz when rise not satisfied
    tol_prob_verylow = 9d-1 /           &
      incr_more * tol_prob_highrise,    & ! max prob in dz when increase dz faster
    tol_prob_zero = 1d-3                  ! prob when rate(z1) or rate(z2) = 0
    
end module propagation_params
!==============================================================================!
!==============================================================================!










module test_module
  implicit none

  logical, parameter :: show_interactions = .false.                             !!!!! testing inter !!!!!
  logical, parameter :: show_propagation_wr = .false.                           !!!!! testing prop (wr) !!!!!
  logical, parameter :: show_propagation_i = .false.                            !!!!! testing prop (i) !!!!!
  logical, parameter :: show_iterations = .false.                               !!!!! testing it !!!!!
  logical, parameter :: show_integration = .false.

!  integer, parameter :: n_checks = 2000                                         !!!!! testing zvals !!!!!
!  double precision zvalueschecked(n_checks)                                     !!!!! testing zvals !!!!!

  integer pause_propagation                                                     !!!!! testing prop !!!!!
  integer pause_iterations                                                      !!!!! testing it !!!!!

end module test_module
