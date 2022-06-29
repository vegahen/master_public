!==============================================================================!
!==============================================================================!
subroutine init(myid)
  implicit none
  integer myid
  character(8) date
  character(10) time

  call init_integration
  call init_dyn_eve
  
  call init_test()                                                              !!!!! testing !!!!!

  call init_general(myid)
  if (myid==0) call init_rates_AGN

  call date_and_time(date, time)
  write(*,*) "start at: ", time(:2), ":", time(3:4), ":", time(5:6)
  
end subroutine init
!==============================================================================!
!==============================================================================!
subroutine init_general(myid)
  use internal; use user_variables
  use constants, only : zero, infty
  implicit none
  integer myid

  zero = 0d0
  infty = -log(zero)

! initialisation for random number (NumRec):
  iseed = 15321 + 2*(1+iseed_shift)*(myid+1)

  d_f = log10(E_min)-0.1d0         ! init internal variables

end subroutine init_general
!==============================================================================!
!==============================================================================!
subroutine init_rates_agn
  use user_variables
  use agn_fit
  implicit none

  integer rw, i_z, i_e

  double precision foo, bar, foobar

  double precision integrate_mu_pair, integrate_mu_ics

  rw = 11
  open(unit = rw, file = agn_pair_filename)

  if (readagnfit) then
    write(*,*) "attempting to read interaction rates, if this fails set readagnfit to .false."
    call gridfromfile(heights_pair, energies_pair, rates_pair, n_zpair, &
      n_epair, rw)
  else
    write(*,*) "calculating interaction rates"
    call init_lin_scale_nomax(heights_pair, n_zpair, zmin_pair, d_zpair, rw)
    call init_lin_scale_nomax(energies_pair, n_epair, emin_pair, d_epair, rw)
    do i_z = 1, n_zpair
      if (modulo(i_z, 10) == 0) write(*,*) i_z, n_zpair
      do i_e = 1, n_epair
        rates_pair(i_e, i_z) = log10(integrate_mu_pair(10**heights_pair(i_z), 10**energies_pair(i_e), -1d0, -1d0, foo, bar))
        write(rw, *) rates_pair(i_e, i_z)
      end do
      write(rw, *) "endrow"
    end do
    call emptyinfo(rw, 10)
  end if

  close(rw)
  open(unit = rw, file = agn_ics_filename)
  
  if (readagnfit) then
    call gridfromfile(heights_ics, energies_ics, rates_ics, n_zics, n_eics, rw)
  else
    call init_lin_scale_nomax(heights_ics, n_zics, zmin_ics, d_zics, rw)
    call init_lin_scale_nomax(energies_ics, n_eics, emin_ics, d_eics, rw)
    do i_z = 1, n_zics
      if (modulo(i_z, 10) == 0) write(*,*) i_z, n_zpair
      do i_e = 1, n_eics
        rates_ics(i_e, i_z) = log10(integrate_mu_ics(10**heights_ics(i_z), 10**energies_ics(i_e), -1d0, -1d0, foo, bar, foobar))
        write(rw, *) rates_ics(i_e, i_z)
      end do
      write(rw, *) "endrow"
    end do
    call emptyinfo(rw, 10)
  end if

  close(rw)

end subroutine init_rates_agn
!==============================================================================!
!==============================================================================!
subroutine inject
  use heap_agn, only : event; use user_variables, only : z_initial
  use agn_data, only : r_s
  implicit none
  double precision, parameter ::  E_min=1.d10,E_max=1d15,alpha=2.d0
  double precision r  
  type(event) new_event
  
  double precision ran0

  r = ran0() ! keep to ensure consistent seeding
  !if (r<0.5d0) then
  !   new_event%icq = 0
  !else
  !   new_event%icq = 1
  !end if
  new_event%icq = 0
  new_event%z = z_initial * r_s

! energy & weight:
  r = log(E_max/E_min) * ran0()

  new_event%e = E_min * exp(r)            ! sample from uniform log distribution
  new_event%e = 1d13                      ! comment if you want to use power law
  
  new_event%w = (E_min/new_event%e)**(alpha-1.d0)

  call heap_insert_events(new_event)

end subroutine inject
!==============================================================================!
!==============================================================================!










subroutine init_test()
  use test_module
  implicit none

  if (show_propagation_i) then                                                  !!!!! testing prop !!!!!
    pause_propagation = 0                                                       !!!!! testing prop !!!!!
  else                                                                          !!!!! testing prop !!!!!
    pause_propagation = -1                                                      !!!!! testing prop !!!!!
  end if                                                                        !!!!! testing prop !!!!!

  if (show_iterations) then                                                     !!!!! testing it !!!!!
    pause_iterations = 0                                                        !!!!! testing it !!!!!
  else                                                                          !!!!! testing it !!!!!
    pause_iterations = -1                                                       !!!!! testing it !!!!!
  end if                                                                        !!!!! testing it !!!!!

end subroutine init_test
