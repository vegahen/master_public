
!==============================================================================!
!==============================================================================!
!                                error handling                                !
!==============================================================================!
subroutine error(string,s)
  character (len=*), intent(in) :: string
  integer, intent(in) :: s
  integer, save :: n_warn  
  integer :: n_warn_max = 100

  if (s==1.or.s==11) then                   ! warning message
     write(*,*) 
     write(*,*) 'Warning:'
     write(*,*) string
     write(99,*) string
     if (s==1) then                         ! warning 
        n_warn = n_warn+1
        if (n_warn>n_warn_max) then
           write(*,*)
           write(*,*) 'more than',n_warn_max,' warnings!'
           write(*,*)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!  ELMAG 3.01 stops program excecution  !'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           stop
        endif
     endif
  endif
  if (s==0) then                    ! error
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!   a serious error:                    !'
     write(99,*)'!   a serious error:                    !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)
     write(*,*) string
     write(99,*)string
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!  ELMAG 3.01 stops program excecution  !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop
  endif

end subroutine error
!==============================================================================!
!==============================================================================!
!          random number generator from Numerical Recipes (Fortran90)          !
!==============================================================================!
function ran0()
  use internal, only : iseed
  implicit none
  integer, parameter :: K4B=selected_int_kind(9)
!  integer(K4B), intent(inout) :: iseed
  double precision ran0
  integer(K4B),parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  real, save :: am
  integer(K4B), save :: ix=-1, iy=-1,k

  if (iseed <= 0 .or. iy < 0) then
     am = nearest(1.e0,-1.e0)/IM
     iy = ior(ieor(888889999,abs(iseed)),1)
     ix = ieor(777755555,abs(iseed))
     iseed = abs(iseed)+1
  end if
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy-k*IQ)-IR*k
  if (iy<0) iy = iy+IM
  ran0 = am*ior(iand(IM,ieor(ix,iy)),1)
end function ran0
!==============================================================================!
!==============================================================================!
subroutine init_lin_scale_nomax(x, n, min, d, w)
  implicit none

  double precision, dimension(n), intent(inout)       :: x
  integer, intent(in)                                 :: n
  integer, intent(in)                                 :: min
  integer, intent(in)                                 :: d
  integer, intent(in)                                 :: w

  integer                                             :: i
  double precision                                    :: idbl

  do i = 1, n
    idbl = i
    x(i) = min + (idbl - 1) / d
  end do
  if (w /= 0) then
    do i = 1, n
      write(w,*) x(i)
    end do
    write(w,*) "end"
  end if

end subroutine init_lin_scale_nomax
!==============================================================================!
!==============================================================================!
subroutine emptyinfo(w, n)
  implicit none

  integer, intent(in) :: w, n

  integer i

  do i = 1, n
    write(w, *) 0
  end do
  write(w, *) "end"

end subroutine emptyinfo
!==============================================================================!
!==============================================================================!
subroutine gridfromfile(x, y, grid, nx, ny, r)
  implicit none

  double precision, intent(inout) :: x(nx), y(ny), grid(ny,nx)
  integer, intent(in) :: nx, ny, r

  double precision foo
  character*60 :: bar
  integer i, j

  do i = 1, nx
    read(r,*) x(i)
  end do
  read(r,*) bar
  if (.not. bar == "end") call error("x-axis size does not match file", 0)
  do j = 1, ny
    read(r,*) y(j)
  end do
  read(r,*) bar
  if (.not. bar == "end") call error("y-axis size does not match file", 0)
  do i = 1, nx
    do j = 1, ny
      read(r,*) grid(j,i)
    end do
    read(r,*) bar
    if (.not. bar == "endrow") call error("end of row not found", 0)
  end do
  do i = 1, 10
    read(r,*) foo
  end do
  read(r,*) bar
  if (.not. bar == "end") call error("end of file not reached", 0)

end subroutine gridfromfile
!==============================================================================!
!==============================================================================!
