!==============================================================================!
!==============================================================================!
double precision function interpolate_rate(z, e, icq, rlog)
  use agn_fit; use constants, only : infty
  implicit none

  double precision, intent(in) :: z, e
  integer, intent(in) :: icq
  double precision, intent(inout) :: rlog

  double precision zlog, elog, g1, g2, g
  integer ij(2), i, j

  zlog = log10(z)
  elog = log10(e)

  if (icq == 0) then
    call interpolationsite(zlog, elog, heights_pair(1), energies_pair(1), &
      n_zpair, n_epair, d_zpair, d_epair, ij)
  else if (icq == 1 .or. icq == -1) then
    call interpolationsite(zlog, elog, heights_ics(1), energies_ics(1), &
      n_zics, n_eics, d_zics, d_eics, ij)
  end if

  i = ij(1)
  j = ij(2)

  if (icq == 0) then
    if (.not. rates_pair(j,i) > -infty &
      .or. .not. rates_pair(j+1,i) > -infty &
      .or. .not. rates_pair(j,i+1) > -infty &
      .or. .not. rates_pair(j+1,i+1) > -infty) then
      g = -infty
    else
      g1 = rates_pair(j,i) + (elog-energies_pair(j))/(energies_pair(j+1)-energies_pair(j)) &
        * (rates_pair(j+1,i)-rates_pair(j,i))
      g2 = rates_pair(j,i+1) + (elog-energies_pair(j))/(energies_pair(j+1)-energies_pair(j)) &
        * (rates_pair(j+1,i+1)-rates_pair(j,i+1))
      g = g1 + (zlog-heights_pair(i))/(heights_pair(i+1)-heights_pair(i)) * (g2-g1)
    end if
  else if (icq == 1 .or. icq == -1) then
    if (.not. rates_ics(j,i) > -infty &
      .or. .not. rates_ics(j+1,i) > -infty &
      .or. .not. rates_ics(j,i+1) > -infty &
      .or. .not. rates_ics(j+1,i+1) > -infty) then
      g = -infty
    else
      g1 = rates_ics(j,i) + (elog-energies_ics(j))/(energies_ics(j+1)-energies_ics(j)) &
        * (rates_ics(j+1,i)-rates_ics(j,i))
      g2 = rates_ics(j,i+1) + (elog-energies_ics(j))/(energies_ics(j+1)-energies_ics(j)) &
        * (rates_ics(j+1,i+1)-rates_ics(j,i+1))
      g = g1 + (zlog-heights_ics(i))/(heights_ics(i+1)-heights_ics(i)) * (g2-g1)
    end if
  else
    call error("interpolate rate of unknown particle", 0)
  end if

  rlog = g
  interpolate_rate = 10**g

end function interpolate_rate
!==============================================================================!
!==============================================================================!
subroutine interpolationsite(x, y, x0, y0, n_x, n_y, d_x, d_y, ij)
  implicit none

  double precision, intent(in) :: x, y, x0, y0
  integer, intent(in) :: n_x, n_y, d_x, d_y
  integer, dimension(2), intent(inout) :: ij

  double precision i, j

  i = d_x * (x - x0)
  j = d_y * (y - y0)
  if (i < 0) call error("x-coordinate is below grid", 0)
  if (j < 0) call error("y-coordinate is below grid", 0)
  if (i > n_x - 1) call error("x-coordinate is above grid", 0)
  if (j > n_y - 1) call error("y-coordinate is above grid", 0)
  if (.not. i < n_x - 1) i = n_x - 2
  if (.not. j < n_y - 1) j = n_y - 2
  ij(1) = int(i+1)
  ij(2) = int(j+1)

end subroutine interpolationsite
!==============================================================================!
!==============================================================================!
