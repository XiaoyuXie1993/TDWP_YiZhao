subroutine initial()

  use spectral_density
  use time_evolution
  
  character*20 :: ch

  open(11, file = 'input')
! Hamiltonian
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_basis
  allocate(H0(N_basis, N_basis))
  do i = 1, N_basis
    read(11, *) H0(i, :)
  end do
! paramters in spectral_density
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_omega
  read(11, '(A)', advance = 'no') ch
  read(11, *) interval_omega
  read(11, '(A)', advance = 'no') ch
  read(11, *) eta
  read(11, '(A)', advance = 'no') ch
  read(11, *) omega_c
  read(11, '(A)', advance = 'no') ch
  read(11, *) beta
  allocate(n_therm(N_omega), h(N_omega))
  allocate(phi(N_omega, 2))
! paramters in time_evolution
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) total_time
  interval_time = total_time / time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_statistic
  allocate(psi0(N_basis))
  read(11, '(A)', advance = 'no') ch
  read(11, *) psi0
  close(11)

end subroutine

subroutine initialphi()

  use constants
  use spectral_density
  
  double precision, allocatable :: U1(:), U2(:), tmp1(:), tmp2(:)
  
  allocate(U1(N_omega), U2(N_omega), tmp1(N_omega), tmp2(N_omega))
  
  call random_number(U1)
  call random_number(U2)
  tmp1 = dsqrt(-2.0d0 * dlog(U1))
  tmp2 = 2.0d0 * pi * U2
  phi(:, 1) = tmp1 * dcos(tmp2)
  phi(:, 2) = tmp1 * dsin(tmp2)

  deallocate(U1, U2, tmp1, tmp2)
  
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