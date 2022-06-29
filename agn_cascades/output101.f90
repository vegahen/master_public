!==============================================================================!
!==============================================================================!
subroutine output(set,n_proc)
  use result; use internal; use user_variables
  implicit none
  integer set,n_proc,i,j,pid
  double precision l,E,m,nu_tot

  call banner(n_proc,set)

  n_tot = n_tot + n_start*n_proc 
  write(*,*) 'set,n_tot',set,n_tot

  open(20,file='Data/spec_aneut'//filename)    ! -8
  open(21,file='Data/spec_aprot'//filename)    ! -7
  open(22,file='Data/spec_anum'//filename)     ! -5
  open(23,file='Data/spec_anue'//filename)     ! -4
  open(24,file='Data/spec_pos'//filename)      ! -1
  open(25,file='Data/spec_gam'//filename)      !  0
  open(26,file='Data/spec_ele'//filename)
  open(27,file='Data/spec_nue'//filename)
  open(28,file='Data/spec_num'//filename)
  open(29,file='Data/spec_prot'//filename)
  open(30,file='Data/spec_neut'//filename)
  open(50,file='Data/spec_nu'//filename)
  do j=1,n_enbin
     nu_tot = 0.d0
     do i=1,n_stable
        pid = stable_pid(i)
        l = dble(j)*dn+d_f
        E = 10.d0**l
        m=En_f_tot(pid,j)/(dble(n_tot)*log(10.d0)*dn)         
        if (m>0.d0) write(19+i,23) E,m,log10(E),log10(m)
        !if (m>0.d0) write(20+i,23) E,m,log10(E),log10(m)
        if (abs(pid)==4.or.abs(pid)==5) nu_tot=nu_tot+m
     end do
     write(50,23) E,nu_tot
  end do
  close(20); close(21); close(22); close(23); close(24); close(25);
  close(26); close(27); close(28); close(29); close(30); close(50);
23 format(24E16.6)

  En_f_tot = 0.d0

end subroutine output
!==============================================================================!
!==============================================================================!
subroutine banner(n_proc,i)
  use user_variables
  implicit none
  integer n_proc,i

  write(*,*) 
  write(*,*) ' files saved as ',filename
  write(*,*) 
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  if (i==0) write(*,*) '!  SNR -- start of the main program      !'
  if (i==n_sets) write(*,*) '!  SNR - end of the main program         !'
  write(*,*) '!  # processes,           ',n_proc,'  !'
  write(*,*) '!  injected particles/set,   ',n_start*n_proc,'  !'
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

end subroutine banner
!==============================================================================!
!==============================================================================!










subroutine output_test
  use user_variables
  use heap_agn, only : n_injected, n_created, n_out, n_ignored;
!  use test_module                                                               !!!!! testing zvals !!!!!
  implicit none

  integer w, i
  character(8) date
  character(10) time

  w = 11

!  open(unit=w, file="test/zvalueschecked.txt")                                  !!!!! testing zvals !!!!!
!  write(*,*) "now writing checked to file"                                      !!!!! testing zvals !!!!!
!  do i = 1, n_checks                                                            !!!!! testing zvals !!!!!
!    write(w,*) zvalueschecked(i)                                                !!!!! testing zvals !!!!!
!  end do                                                                        !!!!! testing zvals !!!!!
!  write(w,*) "end"                                                              !!!!! testing zvals !!!!!
!  call emptyinfo(w, 2000)                                                       !!!!! testing zvals !!!!!
!  call emptyinfo(w, 10)                                                         !!!!! testing zvals !!!!!
!  close(w)                                                                      !!!!! testing zvals !!!!!
  
  call date_and_time(date, time)
  write(*,*) "done at: ", time(:2), ":", time(3:4), ":", time(5:6)
  write(*,*) "================================= output_test =================================="
  write(*,*) n_sets, "- n_sets    ", n_start, "- n_start  ", iseed_shift, "- iseed_shift"
  write(*,*) n_injected, "- n_injected", n_created,  "- n_created", n_out, "- n_out      ", n_ignored, "- n_ignored"
  write(*,*) z_max, "- z_max"
  write(*,*) "================================================================================"

end subroutine output_test
