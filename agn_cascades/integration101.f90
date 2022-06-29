!==============================================================================!
!==============================================================================!
module mu_integration
  implicit none
  save

  integer, parameter ::                 &
  nmu_estimate_pair = 51,               &
  consecutivenonzeros_tol_pair = 26,    &

  nmu_estimate_ics = 51,                &
  consecutivenonzeros_tol_ics = 26

  integer smu, nmu

  double precision, allocatable, dimension(:,:) :: mu_est

  double precision p_left, h_p_left, p_right, h_p_right, mu_curr, h_curr, &
    mu_prev, h_prev

end module mu_integration
!==============================================================================!
!==============================================================================!
module nu_integration
  implicit none
  save

  integer, parameter ::                 &
  entrynu_searches = 15,                &
  epsilon_loops = 10,                   &
  nnu_factor_pair = 1,                  &

  nnu_factor_ics = 1

  double precision, parameter ::        &
  entrynu_scale_factor = 1d0 / (1d0 + 2d0**(-21d0)),  &
  epsilon = 5d-324,                     &
  peak_tol_nu_pair = 9.99d-1,           &
  foot_tol_nu_pair = 1d-3,              &
  estimate_tol_nu_pair = 1d-3,          &

  thomson_tol = 9.99d-1,                &
  peak_tol_nu_ics = 9.99d-1,            &
  foot_tol_nu_ics = 1d-3,               &
  estimate_tol_nu_ics = 1d-3

  integer snu, nnu

  double precision, allocatable, dimension(:,:) :: nu_est

  double precision nu_a, g_a, nu_b, g_b, a_left, g_a_left, nu_prev, g_prev, &
  nu_curr, g_curr, a_right, g_a_right, b_left, g_b_left, b_midleft, &
  g_b_midleft, b_midright, g_b_midright, b_right, g_b_right 

end module nu_integration
!==============================================================================!
!==============================================================================!
subroutine init_integration
  use mu_integration; use nu_integration
  use test_module, only : show_integration
  implicit none

  if (show_integration) write(*,*) "initializing integration lists"

  smu = 1
  snu = 1
  nmu = 0
  nnu = 0

  allocate(mu_est(2,smu))
  allocate(nu_est(2,snu))

end subroutine init_integration
!==============================================================================!
!==============================================================================!
subroutine delete_integration
  use mu_integration; use nu_integration
  use test_module, only : show_integration
  implicit none
  
  if (show_integration) write(*,*) "deleting integration lists"

  smu = 0
  snu = 0
  nmu = 0
  nnu = 0
  
  deallocate(mu_est)
  deallocate(nu_est)

end subroutine delete_integration
!==============================================================================!
!==============================================================================!
subroutine increase_mu
  use mu_integration
  use test_module, only : show_integration
  implicit none

  integer i
  double precision, allocatable, dimension(:,:) :: temp

  if (smu < 1) call error("attempting to increase size of empty mu array", 0)
  smu = smu*2

  allocate(temp(2,nmu))
  do i = 1, nmu
    temp(1,i) = mu_est(1,i)
    temp(2,i) = mu_est(2,i)
  end do

  deallocate(mu_est)
  allocate(mu_est(2,smu))

  do i = 1, nmu
    mu_est(1,i) = temp(1,i)
    mu_est(2,i) = temp(2,i)
  end do
  deallocate(temp)
  
  if (show_integration) write(*,*) "increased mu_est size to", smu

end subroutine increase_mu
!==============================================================================!
!==============================================================================!
subroutine increase_nu
  use nu_integration
  use test_module, only : show_integration
  implicit none

  integer i
  double precision, allocatable, dimension(:,:) :: temp

  if (snu < 1) call error("attempting to increase size of empty nu array", 0)
  snu = snu*2

  allocate(temp(2,nnu))
  do i = 1, nnu
    temp(1,i) = nu_est(1,i)
    temp(2,i) = nu_est(2,i)
  end do

  deallocate(nu_est)
  allocate(nu_est(2,snu))

  do i = 1, nnu
    nu_est(1,i) = temp(1,i)
    nu_est(2,i) = temp(2,i)
  end do
  deallocate(temp)
  
  if (show_integration) write(*,*) "increased nu_est size to", snu

end subroutine increase_nu
!==============================================================================!
!==============================================================================!
subroutine add_mu(mu, h)
  use mu_integration
  use test_module, only : show_integration
  implicit none

  double precision, intent(in) :: mu, h
  
  if (.not. smu > nmu) call increase_mu

  nmu = nmu + 1
  mu_est(1,nmu) = mu
  mu_est(2,nmu) = h    

  if (show_integration) write(*,*) "added to mu", mu, h

end subroutine add_mu
!==============================================================================!
!==============================================================================!
subroutine add_nu(nu, g)
  use nu_integration
  use test_module, only : show_integration
  implicit none

  double precision, intent(in) :: nu, g
  
  if (.not. snu > nnu) call increase_nu
  
  nnu = nnu + 1
  nu_est(1,nnu) = nu
  nu_est(2,nnu) = g    

  if (show_integration) write(*,*) "added to nu", nu, g

end subroutine add_nu
!==============================================================================!
!==============================================================================!
double precision function integrand_pair(t, mu, nu, e)
use constants, only : twome2c4onh, honk, alpha2honc
implicit none

double precision, intent(in) :: t, mu, nu, e

double precision expn, s, beta, fp, p, res

s = e * nu * (1 - mu)
if (s > twome2c4onh) then
  expn = honk * nu / t
  if (expn < 709) then
    beta = sqrt(1 - twome2c4onh / s)
    fp = (3 - beta**4) * log((1+beta)/(1-beta)) - 2*beta*(2-beta**2)
    p = nu / (exp(expn) - 1)

    res = alpha2honc / e * p * fp
  else
    res = 0d0
  end if
else
  res = 0d0
end if

integrand_pair = res

end function integrand_pair
!==============================================================================!
!==============================================================================!
double precision function integrand_ics(t, beta, mu, nu, e)
use constants, only : honk, twohonme2c4, eightalpha2h2onthreeme2c5
use nu_integration, only : thomson_tol
implicit none

double precision, intent(in) :: t, beta, mu, nu, e

double precision expn, p, y_min, fc, res

expn = honk * nu / t
if (expn < 709) then
  p = nu**2 / (exp(expn) - 1)
  y_min = 1d0 / (1d0 + twohonme2c4*e*nu * (1d0 - beta*mu))

  if (y_min < thomson_tol) then
    fc = y_min * 7.5d-1 * (-log(y_min) / (1d0 - y_min) &
      * (1 - 4d0*y_min * (1d0 + y_min) / (1d0 - y_min)**2d0) &
      + 8d0*y_min / (1d0 - y_min)**2d0 + 5d-1 * (1 + y_min))
  else
    fc = y_min
  end if

  res = eightalpha2h2onthreeme2c5 * p * (1 - beta*mu) * fc / beta
else
  res = 0d0
end if

integrand_ics = res

end function integrand_ics
!==============================================================================!
!==============================================================================!
double precision function integrate_nu_pair(mu, z, e, r_nu, nu_int)
use nu_integration; use agn_data
use constants, only : fourthgonfourpisigma, twome2c4onh
implicit none

double precision, intent(in) :: mu, z, e, r_nu
double precision, intent(inout) :: nu_int

logical failed
integer i, epsilon_counter, leftidx, rightidx
double precision r, t, ds, peak_proximity, low, high, foot_proximity, &
  estimate_proximity, de, top, bottom, res, cum_target_nu

double precision integrand_pair

r = z * sqrt(1/mu**2 - 1)
t = fourthgonfourpisigma * (bh_m * bh_mdot / r**3d0)**2.5d-1

failed = .true.
nu_prev = entrynu_scale_factor * twome2c4onh / (e * (1-mu))
g_prev = integrand_pair(t, mu, nu_prev, e)
if (.not. g_prev == 0d0) call error("entrynu is not zero", 0)
do i = 0, entrynu_searches
  nu_curr = nu_prev * (1d0 + 2d0**(-i))
  g_curr = integrand_pair(t, mu, nu_curr, e)
  if (g_curr > 0) then
    failed = .false.
    exit
  end if
end do

if (.not. failed) then
  a_left = nu_prev
  g_a_left = g_prev
  a_right = nu_curr
  g_a_right = g_curr
  
  b_left = a_left
  g_b_left = g_a_left
  b_right = a_right
  g_b_right = g_a_right

  ds = (b_right - b_left) / 3d0
  b_midleft = b_left + ds
  g_b_midleft = integrand_pair(t, mu, b_midleft, e)
  b_midright = b_left + 2d0*ds
  g_b_midright = integrand_pair(t, mu, b_midright, e)

  do while (.not. g_b_midright > g_b_right)
    b_midright = b_right
    g_b_midright = g_b_right
    
    ds = (b_midright - b_left) / 2d0
    b_midleft = b_left + ds
    g_b_midleft = integrand_pair(t, mu, b_midleft, e)
    b_right = b_left + 3d0*ds
    g_b_right = integrand_pair(t, mu, b_right, e)
  end do

  peak_proximity = 0d0
  do while (peak_proximity < peak_tol_nu_pair)
    if (g_b_midleft < g_b_midright) then
      b_left = b_midleft
      g_b_left = g_b_midleft
    else if (g_b_midleft > g_b_midright) then
      b_right = b_midright
      g_b_right = g_b_midright
    else
      b_left = b_midleft
      g_b_left = g_b_midleft
      b_right = b_midright
      g_b_right = g_b_midright
    end if

    ds = (b_right - b_left) / 3d0
    b_midleft = b_left + ds
    g_b_midleft = integrand_pair(t, mu, b_midleft, e)
    b_midright = b_left + 2d0*ds
    g_b_midright = integrand_pair(t, mu, b_midright, e)
    if (.not. g_b_left > 0d0 .or. .not. g_b_midleft > 0d0 &
      .or. .not. g_b_midright > 0d0 .or. .not. g_b_right > 0d0) then
      peak_proximity = 0d0
    else
      if (g_b_left < g_b_right) then
        low = g_b_left
      else
        low = g_b_right
      end if
      if (g_b_midleft < g_b_midright) then
        high = g_b_midright
      else
        high = g_b_midleft
      end if
      peak_proximity = low / high
    end if
  end do

  if (g_b_midleft < g_b_midright) then
    nu_b = b_midright
    g_b = g_b_midright
  else
    nu_b = b_midleft
    g_b = g_b_midleft
  end if

  a_right = a_left + (nu_b - a_left) / 1.5d3
  g_a_right = integrand_pair(t, mu, a_right, e)
  if (.not. g_a_right > 0d0) then
    a_right = b_left
    g_a_right = g_b_left
  end if

  foot_proximity = g_a_right / g_b
  epsilon_counter = 0
  do while (foot_proximity > foot_tol_nu_pair)
    nu_curr = a_left + (a_right - a_left) / 2d0
    g_curr = integrand_pair(t, mu, nu_curr, e)
    if (g_curr > 0) then
      foot_proximity = g_curr / g_b
      a_right = nu_curr
      g_a_right = g_curr
    else
      a_left = nu_curr
    end if

    if (.not. g_a_right > epsilon) then
      epsilon_counter = epsilon_counter + 1
      if (epsilon_counter > epsilon_loops - 1) then
        exit
      end if
    end if
  end do
  nu_a = a_right
  g_a = g_a_right
  
  nnu = 0
  estimate_proximity = 1d0
  ds = (nu_b - nu_a) / (15d0 * nnu_factor_pair)

  call add_nu(nu_curr, 0d0)
  nu_prev = nu_curr
  g_prev = g_curr

  do while (estimate_proximity > estimate_tol_nu_pair)
    nu_curr = nu_curr + ds
    g_curr = integrand_pair(t, mu, nu_curr, e)

    de = 5e-1 * (nu_curr - nu_prev) * (g_prev + g_curr)
    call add_nu(nu_curr, nu_est(2,nnu) + de)

    top = (de/(ds))
    bottom = (nu_est(2,nnu)/(nu_curr-nu_a))
    if (.not. bottom > 0d0 .or. .not. top > 0d0) then
      if (de > 0d0) then
        estimate_proximity = 1d0
      else
        estimate_proximity = 0d0
      end if
    else
      estimate_proximity = top / bottom
    end if
    
    if (modulo(nnu-1, nnu_factor_pair) == 0d0) then
      if (nnu-1 > 7 * nnu_factor_pair) then
        ds = (nu_curr - nu_a) / (5 * nnu_factor_pair)
      else if (nnu-1 < 3 * nnu_factor_pair) then
        ds = 2d0*ds
      end if
    end if

    nu_prev = nu_curr
    g_prev = g_curr
    
  end do

  res = nu_est(2,nnu)

else
  res = 0d0
end if

if (.not. r_nu < 0) then
  if (res == 0) call error("ran_nu_pair - h is zero", 0)
  if (.not. r_nu > 0 .or. .not. r_nu < 1) then
    call error("ran_nu_pair - invalid r", 0)
  end if

  cum_target_nu = r_nu * nu_est(2,nnu)
  do i = 2, nnu
    if (nu_est(2,i) > cum_target_nu) then
      rightidx = i
      leftidx = i-1
      exit
    end if
  end do

  nu_int = nu_est(1,leftidx) + (cum_target_nu - nu_est(2,leftidx)) &
    * (nu_est(1,rightidx) - nu_est(1,leftidx)) &
     / (nu_est(2,rightidx) - nu_est(2,leftidx))
     
  if (.not. nu_int > nu_est(1,1) .or. .not. nu_int < nu_est(1,nnu)) then
    call error("invalid nu_int_pair", 0)
  end if

else if (.not. r_nu == -1d0) then
  call error("invalid un-ran_nu_pair", 0)
end if


integrate_nu_pair = res

end function integrate_nu_pair
!==============================================================================!
!==============================================================================!
double precision function integrate_nu_ics(beta, mu, z, e, r_nu, nu_int)
use nu_integration; use agn_data
use constants, only : fourthgonfourpisigma, wien_constant
implicit none

double precision, intent(in) :: beta, mu, z, e, r_nu
double precision, intent(inout) :: nu_int

integer i, leftidx, rightidx
double precision r, t, ds, peak_proximity, low, high, foot_proximity, &
  estimate_proximity, de, res, cum_target_nu

double precision integrand_ics

r = z * sqrt(1/mu**2 - 1)
t = fourthgonfourpisigma * (bh_m * bh_mdot / r**3d0)**2.5d-1

nu_curr = wien_constant * t
g_curr = integrand_ics(t, beta, mu, nu_curr, e)
nu_prev = 1d-1 * nu_curr
g_prev = integrand_ics(t, beta, mu, nu_prev, e)
if (.not. g_curr > 0) call error("g_ics zero", 0)

ds = (nu_curr - nu_prev) / 3d0
b_left = nu_prev
g_b_left = g_prev
b_midleft = b_left + ds
g_b_midleft = integrand_ics(t, beta, mu, b_midleft, e)
b_midright = b_left + 2d0*ds
g_b_midright = integrand_ics(t, beta, mu, b_midright, e)
b_right = nu_curr
g_b_right = g_curr

do while (.not. g_b_left < g_b_midleft)
  b_left = b_left / 1d1
  g_b_left = integrand_ics(t, beta, mu, b_left, e)
  b_right = b_midleft
  g_b_right = g_b_midleft
  
  ds = (b_right - b_left) / 3d0
  b_midleft = b_left + ds
  g_b_midleft = integrand_ics(t, beta, mu, b_midleft, e)
  b_midright = b_left + 2d0*ds
  g_b_midright = integrand_ics(t, beta, mu, b_midright, e)
end do

nu_prev = b_left
g_prev = g_b_left

peak_proximity = 0d0
do while (peak_proximity < peak_tol_nu_ics)
  if (g_b_midleft < g_b_midright) then
    b_left = b_midleft
    g_b_left = g_b_midleft
  else if (g_b_midleft > g_b_midright) then
    b_right = b_midright
    g_b_right = g_b_midright
  else
    call error("going central nu_ics", 0)
  end if

  ds = (b_right - b_left) / 3
  b_midleft = b_left + ds
  g_b_midleft = integrand_ics(t, beta, mu, b_midleft, e)
  b_midright = b_left + 2d0*ds
  g_b_midright = integrand_ics(t, beta, mu, b_midright, e)

  if (g_b_left < g_b_right) then
    low = g_b_left
  else
    low = g_b_right
  end if
  if (g_b_midleft < g_b_midright) then
    high = g_b_midright
  else
    high = g_b_midleft
  end if
  if (.not. high > 0d0) call error("high is zero", 0)

  peak_proximity = low / high

end do

if (g_b_midleft < g_b_midright) then
  nu_b = b_midright
  g_b = g_b_midright
else
  nu_b = b_midleft
  g_b = g_b_midleft
end if

a_right = nu_prev
g_a_right = g_prev
foot_proximity = g_a_right / g_b
do while (foot_proximity > foot_tol_nu_ics)
  a_right = a_right / 1d1
  g_a_right = integrand_ics(t, beta, mu, a_right, e)
  foot_proximity = g_a_right / g_b
end do

nu_a = a_right
g_a = g_a_right
nu_curr = nu_a
g_curr = g_a

nnu = 0
estimate_proximity = 1d0
ds = (nu_b - nu_a) / (19d0 * nnu_factor_ics)

call add_nu(nu_curr, 0d0)
nu_prev = nu_curr
g_prev = g_curr

do while (estimate_proximity > estimate_tol_nu_ics)
  nu_curr = nu_curr + ds
  g_curr = integrand_ics(t, beta, mu, nu_curr, e)

  de = 5e-1 * (nu_curr - nu_prev) * (g_prev + g_curr)
  call add_nu(nu_curr, nu_est(2,nnu) + de)
  estimate_proximity = (de/ds) / (nu_est(2,nnu)/(nu_curr-nu_a))

  if (modulo(nnu-1, nnu_factor_ics) == 0) then
    if (nnu-1 > 8 * nnu_factor_ics) then
      ds = (nu_curr - nu_a) / (6d0 * nnu_factor_ics)
    else if (nnu-1 < 3 * nnu_factor_ics) then
      ds = ds * 2
    end if
  end if

  nu_prev = nu_curr
  g_prev = g_curr
end do

res = nu_est(2,nnu)

if (.not. r_nu < 0) then
  if (res == 0) call error("ran_nu_ics - h is zero", 0)
  if (.not. r_nu > 0 .or. .not. r_nu < 1) then 
    call error("ran_nu_ics - invalid r", 0)
  end if

  cum_target_nu = r_nu * nu_est(2,nnu)
  do i = 2, nnu
    if (nu_est(2,i) > cum_target_nu) then
      rightidx = i
      leftidx = i-1
      exit
    end if
  end do

  nu_int = nu_est(1,leftidx) + (cum_target_nu - nu_est(2,leftidx)) &
    * (nu_est(1,rightidx) - nu_est(1,leftidx)) &
     / (nu_est(2,rightidx) - nu_est(2,leftidx))
     
  if (.not. nu_int > nu_est(1,1) .or. .not. nu_int < nu_est(1,nnu)) then
    call error("invalid nu_int_ics", 0)
  end if

else if (.not. r_nu == -1d0) then
  call error("invalid un-ran_nu_ics", 0)
end if

integrate_nu_ics = res

end function integrate_nu_ics
!==============================================================================!
!==============================================================================!
double precision function integrate_mu_pair(z, e, r_mu, r_nu, mu_int, nu_int)
use mu_integration; use agn_data
implicit none

double precision, intent(in) :: z, e, r_mu, r_nu
double precision, intent(inout) :: mu_int, nu_int

logical found_left_limit
integer passes, consecutivenonzeros, nmu_pass, i, rightidx, leftidx
double precision mu_min, mu_max, ds, de, res, cum_target_mu, foo

double precision integrate_nu_pair

mu_min = 1d0 / sqrt(1d0 + bh_disksize**2 / z**2)
mu_max = 1d0 / sqrt(1d0 + (3d0*r_s)**2 / z**2)
if (mu_max == 1) call error("mu_pair - mu_max = 1", 0)

p_left = mu_min
h_p_left = integrate_nu_pair(p_left, z, e, -1d0, nu_int)
p_right = mu_max
h_p_right = integrate_nu_pair(p_right, z, e, -1d0, nu_int)

passes = 0
consecutivenonzeros = 0
nmu_pass = nmu_estimate_pair

do while (consecutivenonzeros < consecutivenonzeros_tol_pair)
  if (passes > 0) then
    if (passes > 2) then
      call error("pair - more than 3 passes", 0)
    end if
  end if
  nmu = 0
  consecutivenonzeros = 0
  found_left_limit = .false.

  ds = (p_right - p_left) / (nmu_pass - 1)
  mu_curr = p_left
  h_curr = h_p_left
  mu_prev = mu_curr
  h_prev = h_curr

  if (h_curr > 0d0) then
    consecutivenonzeros = consecutivenonzeros + 1
    found_left_limit = .true.
  end if
  call add_mu(mu_curr, 0d0)

  do i = 1, nmu_pass - 2
    mu_curr = mu_curr + ds
    h_curr = integrate_nu_pair(mu_curr, z, e, -1d0, nu_int)
    if (mu_curr == mu_prev) then
      if (passes > 1 .and. z > 8d21) then
        write(*,*) "mu_pair - mu_curr = mu_prev in non-relevant area, &
          rate is set to zero, z:", z, "e:", e 
        integrate_mu_pair = 0d0
        return
      else
        call error("mu_pair - mu_curr = mu_prev", 0)
      end if
    end if
    if (.not. found_left_limit) then
      if (h_curr > 0) then
        p_left = mu_prev
        h_p_left = h_prev
        found_left_limit = .true.
        consecutivenonzeros = consecutivenonzeros + 1
      end if
    else
      if (.not. h_curr > 0) then
        p_right = mu_curr
        h_p_right = h_curr
        exit
      else
        consecutivenonzeros = consecutivenonzeros + 1
      end if
    end if

    de = 5d-1 * (mu_curr - mu_prev) * (h_prev + h_curr)
    call add_mu(mu_curr, mu_est(2,nmu) + de)
    mu_prev = mu_curr
    h_prev = h_curr
  end do

  if (h_p_right > 0) consecutivenonzeros = consecutivenonzeros + 1
  de = 5d-1 * (p_right - mu_prev) * (h_prev + h_p_right)
  call add_mu(p_right, mu_est(2,nmu) + de)

  passes = passes + 1
  if (.not. found_left_limit) exit
end do

res = mu_est(2, nmu)

if (.not. r_mu < 0) then
  if (.not. res > 0) then
    call error("ran_mu_pair - r is zero", 0)
  end if
  if (.not. r_mu > 0 .or. .not. r_mu < 1) then 
    call error("ran_mu_pair - invalid r", 0)
  end if

  cum_target_mu = r_mu * mu_est(2, nmu)
  do i = 2, nmu
    if (mu_est(2,i) > cum_target_mu) then
      rightidx = i
      leftidx = i-1
      exit
    end if
  end do

  mu_int = mu_est(1,leftidx) + (cum_target_mu - mu_est(2,leftidx)) &
    * (mu_est(1,rightidx) - mu_est(1,leftidx)) &
    / (mu_est(2,rightidx) - mu_est(2,leftidx))
  
  if (.not. mu_int > mu_est(1,1) .or. .not. mu_int < mu_est(1,nmu)) then
    call error("invalid mu_int_pair", 0)
  end if

  foo = integrate_nu_pair(mu_int, z, e, r_nu, nu_int)

else
  if (.not. r_mu == -1d0) call error("invalid un_ran_mu_pair", 0)
end if

integrate_mu_pair = res

end function integrate_mu_pair !perfect
!==============================================================================!
!==============================================================================!
double precision function integrate_mu_ics(z, e, r_mu, r_nu, mu_int, nu_int, &
  beta)
use mu_integration; use agn_data
use constants, only : me2c4
implicit none

double precision, intent(in) :: z, e, r_mu, r_nu
double precision, intent(inout) :: mu_int, nu_int, beta

logical found_left_limit
integer passes, consecutivenonzeros, nmu_pass, i, rightidx, leftidx
double precision mu_min, mu_max, ds, de, res, cum_target_mu, foo

double precision integrate_nu_ics

mu_min = 1d0 / sqrt(1d0 + bh_disksize**2 / z**2)
mu_max = 1d0 / sqrt(1d0 + (3d0*r_s)**2 / z**2)
if (mu_max == 1) call error("mu_ics - mu_max = 1", 0)
if (e**2 < me2c4) call error("ics - e < mec2", 0)
beta = sqrt(1 - me2c4/e**2)

p_left = mu_min
h_p_left = integrate_nu_ics(beta, p_left, z, e, -1d0, nu_int)
p_right = mu_max
h_p_right = integrate_nu_ics(beta, p_right, z, e, -1d0, nu_int)

passes = 0
consecutivenonzeros = 0
nmu_pass = nmu_estimate_ics

do while (consecutivenonzeros < consecutivenonzeros_tol_ics)
  if (passes > 0) call error("ics - more than one pass", 0)
  nmu = 0
  consecutivenonzeros = 0
  found_left_limit = .false.

  ds = (p_right - p_left) / (nmu_pass - 1)
  mu_curr = p_left
  h_curr = h_p_left
  mu_prev = mu_curr
  h_prev = h_curr

  if (h_curr > 0d0) then
    consecutivenonzeros = consecutivenonzeros + 1
    found_left_limit = .true.
  end if
  call add_mu(mu_curr, 0d0)

  do i = 1, nmu_pass - 2
    mu_curr = mu_curr + ds
    h_curr = integrate_nu_ics(beta, mu_curr, z, e, -1d0, nu_int)
    if (mu_curr == mu_prev) call error("mu_ics - mu_curr = mu_prev", 0)
    if (.not. found_left_limit) then
      if (h_curr > 0) then
        p_left = mu_prev
        h_p_left = h_prev
        found_left_limit = .true.
        consecutivenonzeros = consecutivenonzeros + 1
      end if
    else
      if (.not. h_curr > 0) then
        p_right = mu_curr
        h_p_right = h_curr
        exit
      else
        consecutivenonzeros = consecutivenonzeros + 1
      end if
    end if

    de = 5d-1 * (mu_curr - mu_prev) * (h_prev + h_curr)
    call add_mu(mu_curr, mu_est(2,nmu) + de)
    mu_prev = mu_curr
    h_prev = h_curr
  end do

  if (h_p_right > 0) consecutivenonzeros = consecutivenonzeros + 1
  de = 5d-1 * (p_right - mu_prev) * (h_prev + h_p_right)
  call add_mu(p_right, mu_est(2,nmu) + de)

  passes = passes + 1
end do

res = mu_est(2, nmu)

if (.not. r_mu < 0) then
  if (.not. res > 0) then
    call error("ran_mu_ics - r is zero", 0)
  end if
  if (.not. r_mu > 0 .or. .not. r_mu < 1) then 
    call error("ran_mu_ics - invalid r", 0)
  end if

  cum_target_mu = r_mu * mu_est(2, nmu)
  do i = 2, nmu
    if (mu_est(2,i) > cum_target_mu) then
      rightidx = i
      leftidx = i-1
      exit
    end if
  end do

  mu_int = mu_est(1,leftidx) + (cum_target_mu - mu_est(2,leftidx)) &
    * (mu_est(1,rightidx) - mu_est(1,leftidx)) &
    / (mu_est(2,rightidx) - mu_est(2,leftidx))
  
  if (.not. mu_int > mu_est(1,1) .or. .not. mu_int < mu_est(1,nmu)) then
    call error("invalid mu_int_ics", 0)
  end if

  foo = integrate_nu_ics(beta, mu_int, z, e, r_nu, nu_int)

else
  if (.not. r_mu == -1d0) call error("invalid un_ran_mu_ics", 0)
end if

integrate_mu_ics = res

end function integrate_mu_ics
!==============================================================================!
!==============================================================================!
