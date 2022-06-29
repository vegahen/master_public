!==============================================================================!
!==============================================================================!
subroutine init_dyn_eve
  use heap_agn
  use test_module, only : show_interactions                                     !!!!! testing inter !!!!!
  implicit none

  if (show_interactions) write(*,*) "initializing heap"                         !!!!! testing inter !!!!!

  n = 0
  s = 1
  n_created = 0
  n_out = 0
  n_ignored = 0
  allocate(events(s))

end subroutine init_dyn_eve
!==============================================================================!
!==============================================================================!
subroutine delete_dyn_eve
  use heap_agn
  use test_module, only : show_interactions                                     !!!!! testing inter !!!!!
  implicit none

  if (show_interactions) write(*,*) "deleting heap"                             !!!!! testing inter !!!!!

  n = 0
  s = 0
  deallocate(events)

end subroutine delete_dyn_eve
!==============================================================================!
!==============================================================================!
subroutine increase_dyn_eve
  use heap_agn
  use test_module, only : show_interactions                                     !!!!! testing inter !!!!!
  implicit none

  integer i
  type(event), allocatable, dimension(:) :: tmp

  if (s < 1) call error("attempting to increase size of empty array", 0)
  s = s*2

  allocate(tmp(n))
  do i = 1, n
    tmp(i) = events(i)
  end do

  deallocate(events)
  allocate(events(s))

  do i = 1, n
    events(i) = tmp(i)
  end do
  deallocate(tmp)

  if (show_interactions) write(*,*) "increased heap size to", s                 !!!!! testing inter !!!!!

end subroutine increase_dyn_eve
!==============================================================================!
!==============================================================================!
recursive subroutine heapify_events(i)
  use heap_agn
  implicit none

  integer, intent(in) :: i

  integer l, r, m
  type(event) temp

  l = 2*i
  r = 2*i + 1
  m = i
  if (.not. l > n .and. events(l)%e < events(i)%e) m = l
  if (.not. r > n .and. events(r)%e < events(m)%e) m = r
  if (.not. m == i) then
    temp = events(m)
    events(m) = events(i)
    events(i) = temp
    call heapify_events(m)
  end if

end subroutine heapify_events
!==============================================================================!
!==============================================================================!
subroutine heap_insert_events(elem)
  use heap_agn
  use test_module, only : show_interactions                                     !!!!! testing inter !!!!!
  implicit none

  type(event), intent(in) :: elem

  integer i, p

  if (show_interactions) then                                                   !!!!! testing inter !!!!!
    write(*,*) "inserting ", n_created+1, "type", elem%icq, &                   !!!!! testing inter !!!!!
    "with energy", elem%e, "at z =", elem%z, n+1                                !!!!! testing inter !!!!!
  end if                                                                        !!!!! testing inter !!!!!

  if (n+1 > s) call increase_dyn_eve
  n = n+1

  i = n
  do while (i > 1)
    p = i/2
    if (elem%e < events(p)%e) then
      events(i) = events(p)
      i = p
    else
      exit
    end if
  end do
  events(i) = elem

  n_created = n_created+1
  events(i)%id = n_created

end subroutine heap_insert_events
!==============================================================================!
!==============================================================================!
subroutine heap_extract_events
  use heap_agn
  use test_module, only : show_interactions                                     !!!!! testing inter !!!!!
  implicit none

  if (show_interactions) then                                                   !!!!! testing inter !!!!!
    write(*,*) "extracting", events(1)%id, "type", events(1)%icq, &             !!!!! testing inter !!!!!
    "with energy", events(1)%e, "at z =", events(1)%z, n-1                      !!!!! testing inter !!!!!
  end if                                                                        !!!!! testing inter !!!!!

  if (n < 1) call error("attempting to extract from empty heap", 0)

  events(1) = events(n)
  n = n-1
  call heapify_events(1)

end subroutine heap_extract_events
!==============================================================================!
!==============================================================================!
