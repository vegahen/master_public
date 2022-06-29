!==============================================================================!
!==============================================================================!
program agn_cascade
  use mpi
  use result;  use user_variables, only : n_sets
  implicit none
  integer myid,n_proc,ierr,n_array
  integer set
  real secs(2), secd
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,ierr)

  call init(myid)

  secd = dtime(secs)

  open(unit=69, file="test/interactionpoints")

  do set = 1,n_sets

     call start_particle(set,myid,n_proc)

     n_array = (2*pid_max+1)*n_enbin
     call MPI_REDUCE(En_f,En_f_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                     MPI_COMM_WORLD,ierr)        ! sum individal arrays En_f
     
     if (myid==0) call output(set,n_proc)
  end do

  close(69)

  secd = dtime(secs)
  write(*,*) "run took seconds:", secd

  close(99)

  call delete_dyn_eve
  call delete_integration
  
  call output_test

end program agn_cascade
!==============================================================================!
!==============================================================================!
subroutine start_particle(set,myid,n_proc)
  use user_variables, only : n_start
  use heap_agn, only : n, n_injected
  implicit none
  integer set,myid,n_proc

  n_injected = 0

  do while (n_injected < n_start .or. n > 0)
     if (n == 0) then
        n_injected = n_injected + 1
        call inject
     end if
     call tracer
     !if (myid==0 .and. mod(n_injected*100,n_start)==0 .and. sec==0) &
        !write(*,*) set,n_injected*n_proc
  end do

end subroutine start_particle
!==============================================================================!
!==============================================================================!
subroutine tracer
  use heap_agn
  implicit none

  integer, pointer :: pid
 
  pid => events(1)%icq
 
  select case (pid)
  case (-1:1) 
    call propagate
  case default       
     write(*,*) 'pid',pid
     call error('wrong particle typ in tracer',0)
  end select

end subroutine tracer
!==============================================================================!
!==============================================================================!
subroutine propagate
  use heap_agn; use constants, only : infty; use propagation_params;
  use test_module                                                               !!!!! testing !!!!!
  implicit none
  
  logical done, move
  integer i, j
  double precision ds, ds_max, r, r1, r1log, r2, r2log, z_next, prob, rise

  double precision ran0, Rate_AGN

  integer, pointer :: icq, id
  double precision, pointer :: e, z, w

  icq => events(1)%icq
  id => events(1)%id
  e => events(1)%e
  z => events(1)%z
  w => events(1)%w

  if (.not. z < z_max) then
    if (.not. z < infty) call error ("keep: Propagating particle &
      that has escaped", 0)
    if (z > z_max) call error("test: particle is beyond z_max, &
      should only happen if injected beyond", 0)
    !if (z == z_max) call error("test: particle is at z_max, &
    !  this should happen if created here", 0)
!    zvalueschecked(1) = z                                                       !!!!! testing zvals !!!!!
    z = infty
    done = .true.
  else if (z < 0d0) then
    call error("Propagating particle that has already interacted", 0)
  else
    r1 = rate_agn(icq, z, e, r1log)
    ds_max = z * incr_lots
    if (.not. r1 > 0d0) then
      ds = ds_max
    else
      ds = 9.5d-1 * tol_prob / r1
      if (ds > ds_max) ds = ds_max
    end if
    z_next = z + ds
    done = .false.
  end if

  if (pause_propagation == 0 .and. .not. pause_iterations < 0) then             !!!!! testing it !!!!!
    pause_iterations = 0                                                        !!!!! testing it !!!!!
  end if                                                                        !!!!! testing it !!!!!
  i = 0
  do while (.not. done)
    i = i + 1
!    zvalueschecked(i) = z                                                       !!!!! testing zvals !!!!!

    if (pause_propagation == 0 .and. .not. pause_iterations < 0) then           !!!!! testing it !!!!!
      write(*,*)                                                                !!!!! testing it !!!!!
      write(*,*) id, "- id", i, "- iteration"                                   !!!!! testing it !!!!!
    end if                                                                      !!!!! testing it !!!!!

    if (z_next > z_max) then
      if (pause_propagation == 0 .and. .not. pause_iterations < 0) &            !!!!! testing it !!!!!
        write(*,*) "changed from z_next", z_next, "and ds", ds                  !!!!! testing it !!!!!
      z_next = z_max
      ds = z_next - z
    end if

    if (pause_propagation == 0 .and. .not. pause_iterations < 0) then           !!!!! testing it !!!!!
      write(*,*) z                                                              !!!!! testing it !!!!!
      write(*,*) z_next                                                         !!!!! testing it !!!!!
      write(*,*) ds                                                             !!!!! testing it !!!!!
    end if                                                                      !!!!! testing it !!!!!

    r2 = rate_agn(icq, z_next, e, r2log)

    prob = ds * (r1 + r2) / 2d0
    if (pause_propagation == 0 .and. .not. pause_iterations < 0) &              !!!!! testing it !!!!!
      write(*,*) "prob", prob                                                   !!!!! testing it !!!!!

    if (r1 > 0 .and. r2 > 0) then
      if (prob < tol_prob) then
        rise = abs(r1log - r2log)
        if (pause_propagation == 0 .and. .not. pause_iterations < 0) &          !!!!! testing it !!!!!
          write(*,*) "rise", rise                                               !!!!! testing it !!!!!
        if (rise < tol_rise .or. prob < tol_prob_highrise) then
          if (pause_propagation == 0 .and. .not. pause_iterations < 0) &        !!!!! testing it !!!!!
            write(*,*) 1                                                        !!!!! testing it !!!!!
          move = .true.
        else
          if (pause_propagation == 0 .and. .not. pause_iterations < 0) &        !!!!! testing it !!!!!
            write(*,*) 2                                                        !!!!! testing it !!!!!
          move = .false.
        end if
      else
        if (pause_propagation == 0 .and. .not. pause_iterations < 0) &          !!!!! testing it !!!!!
          write(*,*) 3                                                          !!!!! testing it !!!!!
        move = .false.
      end if
    else
      if (prob < tol_prob_zero) then
        if (pause_propagation == 0 .and. .not. pause_iterations < 0) &          !!!!! testing it !!!!!
          write(*,*) 4                                                          !!!!! testing it !!!!!
        move = .true.
      else
        if (pause_propagation == 0 .and. .not. pause_iterations < 0) &          !!!!! testing it !!!!!
          write(*,*) 5                                                          !!!!! testing it !!!!!
        move = .false.
      end if
    end if 
    
    if (move) then
      if (pause_propagation == 0 .and. .not. pause_iterations < 0) &            !!!!! testing it !!!!!
        write(*,*) "we move"                                                    !!!!! testing it !!!!!
      if (prob > 0) then
        r = ran0()
        if (pause_propagation == 0 .and. .not. pause_iterations < 0) &          !!!!! testing it !!!!!
          write(*,*) "ran", r                                                   !!!!! testing it !!!!!
        if (r < prob) then
          if (pause_propagation == 0 .and. .not. pause_iterations < 0) &        !!!!! testing it !!!!!
            write(*,*) "success at i =", i                                      !!!!! testing it !!!!!
          done = .true.
        end if
      end if
      
      if (.not. prob > 0d0) then
        ds = z_next * incr_lots
      else if (prob < tol_prob_verylow) then
        ds = ds * incr_more
      else
        ds = ds * incr
      end if

      z = z_next
      z_next = z + ds
      r1log = r2log
      r1 = r2

      if (.not. done .and. .not. z < z_max) then
        if (show_propagation_wr) write(*,*) "particle  ", id, &                 !!!!! testing prop (wr) !!!!!
          "escaped to  ", z, "after iterations:", i, e                          !!!!! testing prop (wr) !!!!!
!        zvalueschecked(i+1) = z                                                 !!!!! testing zvals !!!!!
        z = infty
        done = .true.
      end if
    else
      ds = ds / decr
      z_next = z + ds
    end if
    if (pause_propagation == 0) then                                            !!!!! testing it !!!!!
      if (pause_iterations == 0) then                                           !!!!! testing it !!!!!
        write(*,*) "enter next iteration to stop at"                            !!!!! testing it !!!!!
        read(*,*) pause_iterations                                              !!!!! testing it !!!!!
      end if                                                                    !!!!! testing it !!!!!
      pause_iterations = pause_iterations - 1                                   !!!!! testing it !!!!!
    end if                                                                      !!!!! testing it !!!!!
  end do
  
  if (.not. z > z_max) then
    if (show_propagation_wr) write(*,*) "particle  ", id, &                     !!!!! testing prop (wr) !!!!!
      "interacts at", z, "after iterations:", i, r1                             !!!!! testing prop (wr) !!!!!
!    zvalueschecked(i+1) = z                                                     !!!!! testing zvals !!!!!
    call interaction_AGN(E,z,w,icq)
  else
    call store(icq,E,w)
  end if

!  zvalueschecked(i+2) = z                                                       !!!!! testing zvals !!!!!
!  do j = i+3, 2000                                                              !!!!! testing zvals !!!!!
!    !write(*,*) j                                                               !!!!! testing zvals !!!!!
!    zvalueschecked(j) = 0                                                       !!!!! testing zvals !!!!!
!  end do                                                                        !!!!! testing zvals !!!!!

  if (pause_propagation == 0) then                                              !!!!! testing prop (i) !!!!!
    write(*,*) "enter next propagation to stop at"                              !!!!! testing prop (i) !!!!!
    read(*,*) pause_propagation                                                 !!!!! testing prop (i) !!!!!
  else if (show_propagation_wr) then
    write(*,*)
  end if                                                                        !!!!! testing prop (i) !!!!!
  pause_propagation = pause_propagation - 1                                     !!!!! testing prop (i) !!!!!

end subroutine propagate
!==============================================================================!
!==============================================================================!
subroutine interaction_AGN(E,z,w,icq)
  use heap_agn, only : event, n_created, n_ignored, n
  use constants, only : twoh, me2c4
  use user_variables, only : photon_energythr, electron_energythr
  use test_module, only : show_interactions, show_propagation_wr
  implicit none
  integer icq, i
  double precision E, z, w, frac, r_mu, r_nu, foo, mu_int, nu_int, s, beta
  type(event) new_events(2)
  character(8) date
  character(10) time

  double precision ran0, integrate_mu_pair, integrate_mu_ics, zpair, zics
  
  new_events(1)%z = z
  new_events(2)%z = z
  
  r_mu = ran0()
  r_nu = ran0()

  if (icq == 0) then
    new_events(1)%icq = 1
    new_events(2)%icq = -1
    foo = integrate_mu_pair(z, e, r_mu, r_nu, mu_int, nu_int)
    beta = -1d0
    
    s = twoh * e * nu_int * (1 - mu_int)
    frac = zpair(s)
  else if (icq == 1 .or. icq == -1) then
    new_events(1)%icq = icq
    new_events(2)%icq = 0
    foo = integrate_mu_ics(z, e, r_mu, r_nu, mu_int, nu_int, beta)
    
    s = me2c4 + twoh * e * nu_int * (1 - beta * mu_int)
    frac = zics(s)
  else
    call error('wrong particle typ in interaction', 0)
  end if

  if (show_propagation_wr) write(*,*) "r_mu", r_mu, "r_nu", r_nu, "foo", foo, &
    "mu_int", mu_int, "nu_int", nu_int, "s", s, "sqrt(s)", sqrt(s), "beta", &
    beta, "frac", frac, "1-frac", 1-frac

  new_events(1)%w = w
  new_events(2)%w = w

  new_events(1)%e = frac * e     ! make sure you pick the right particle here
  new_events(2)%e = (1d0-frac) * e   ! pick

  write(69,*) z
  z = -1d0
  call heap_extract_events

  do i = 1, 2
    if ((new_events(i)%icq == 0 .and. new_events(i)%e < photon_energythr) &
      .or. ((new_events(i)%icq == 1 .or. new_events(i)%icq == -1) &
      .and. new_events(i)%e < electron_energythr)) then
      n_created = n_created + 1
      n_ignored = n_ignored + 1
      if (show_interactions) then                                               !!!!! testing inter !!!!!
        write(*,*) "ignoring  ", n_created, "type", new_events(i)%icq, &        !!!!! testing inter !!!!!
          "with energy", new_events(i)%e, "at z =", new_events(i)%z             !!!!! testing inter !!!!!
      end if                                                                    !!!!! testing inter !!!!!
    else
      call heap_insert_events(new_events(i))
    end if
    if (modulo(n_created, 1000000) == 0) then
      call date_and_time(date, time)
      write(*,*) n_created, "- created", n, "- heap size           ", &
        "time: ", time(:2), ":", time(3:4), ":", time(5:6)
    end if
  end do

end subroutine interaction_AGN
!==============================================================================!
!==============================================================================!
!                     energy partition for pair production                     !
!==============================================================================!
double precision function zpair(sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for gamma-gamma interaction
!-----------------------------------------------------------------------------
  use constants, only : ame
  use internal, only : debug
  implicit none
  double precision sgam,bet,zmin,z,gb
  double precision psran

  bet=sqrt(max(0.d0,1.d0-4.d0*ame**2/sgam))
  zmin=(1.d0-bet)/2.d0
  do 
     z=.5d0*(2.d0*zmin)**psran()
     gb=(z**2/(1.d0-z)+1.d0-z+(1.d0-bet**2)/(1.d0-z)-(1.d0-bet**2)**2 &
      &/4.d0/z/(1.d0-z)**2)/(1.d0+2.d0*bet**2*(1.d0-bet**2))
     if (debug.gt.0.and.gb.gt.1.d0) write(*,*)'zpair: gb=',gb
     if (psran().lt.gb) exit
  end do
  if (psran().gt.0.5d0) z=1.d0-z
  zpair=z

end function zpair
!==============================================================================!
!==============================================================================!
!        energy partition for inverse Compton (E-fraction taken by e+-)        !
!==============================================================================!
double precision function zics(sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
  use constants, only : ame
  use internal, only : debug
  implicit none
  double precision sgam,zmin,zmax,z,gb
  double precision psran

  zmax=1.d0
  zmin=ame**2/sgam

  if (zmin.ge.zmax) then
     if (debug.gt.0) call error('zmin>zmax in zics',0)
     zics=zmin
  else
     do
        z=zmin*(zmax/zmin)**psran()
        gb=(1.d0+z*z)/2.d0-2.d0*zmin/z*(z-zmin)*(1.d0-z)/(1.d0-zmin)**2
        !if (debug>0.and.gb.gt.1.d0) write(*,*)'zics: gb=',gb,z,zmin,zmax
        if (psran().lt.gb) exit
     end do
     zics=z
  endif

end function zics
!==============================================================================!
!==============================================================================!
double precision function psran()
  implicit none
  double precision ran0
  psran=ran0()
end function psran
!==============================================================================!
!==============================================================================!
double precision function Rate_AGN(icq, z, E, rlog) ! rate/cm
  implicit none
  integer icq
  double precision z, E, rlog

  double precision interpolate_rate
  
  Rate_AGN = interpolate_rate(z, e, icq, rlog)

end function Rate_AGN
!==============================================================================!
!==============================================================================!
subroutine store(pid,En,w)
  use internal; use result
  use heap_agn, only : n_out
  implicit none
  integer pid,i
  double precision En,w,l

  l = log10(En)                         ! energy bin
  i = int( (l-d_f)/dn )
  
  if (i<=0) then 
     call error('stor_esc, E<=Emin',11)
     i=1
  end if
  if (i>n_enbin) then
     i=n_enbin
  end if

  En_f(pid,i) = En_f(pid,i) + w*En

  call heap_extract_events
  n_out = n_out + 1

end subroutine store
!==============================================================================!
!==============================================================================!
