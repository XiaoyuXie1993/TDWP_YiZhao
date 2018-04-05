subroutine initial()

  use spectral_density
  use time_evolution
  
  character*12 :: ch

  open(11, file = 'input')
! paramters in spectral_density
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) alpha
  read(11, '(A)', advance = 'no') ch
  read(11, *) omega_c
  read(11, '(A)', advance = 'no') ch
  read(11, *) beta
  read(11, '(A)', advance = 'no') ch
  read(11, *) check_quantum
!  write(*, '(3f10.5, l5)') alpha, omega_c, beta, check_quantum
!  stop
! paramters in time_evolution
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) total_time
  interval_time = total_time / time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_statistic
!  write(*, '(i7, f10.5, i7)') time_steps, total_time, N_statistic
  close(11)

end subroutine

subroutine initialphi()

  use constants
  use spectral_density
  
  do i = 1, N_basis
    do j = 1, N_omega
      call random_number(phi(i, j))
      phi(i, j) = phi(i, j) * 2.0d0 * pi
    end do
  end do
  
end subroutine
  
!! initial random seed from fortran manual
subroutine init_random_seed(pid)

  use iso_fortran_env, only: int64
  
  implicit none
  
  integer, intent(in) :: pid
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8)
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
   end if
  call random_seed(put=seed)
  
  contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
  
end subroutine init_random_seed