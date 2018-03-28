subroutine initial()

  use constants
  use spectral_density
  use time_evolution
  
  double precision :: omega
  character*20 :: ch

  open(11, file = 'input')
! electronic Hamiltonian
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
  allocate(SP(N_basis, N_basis, N_omega), S(N_basis, N_basis, N_omega))
  allocate(phi(N_basis, N_basis, N_omega))
  read(11, '(A)', advance = 'no') ch
  read(11, *) interval_omega
! cm^{-1} to fs^{-1}
  interval_omega = interval_omega * 2.7992458 * 2.0d0 * pi * 1.0e-5
  read(11, '(A)', advance = 'no') ch
  read(11, *) temperature
  read(11, '(A)', advance = 'no') ch
  read(11, *) check_quantum
  beta = 1.0d0 / temperature * kB
! paramters in time_evolution
  allocate(psi0(N_basis))
  read(11, *)
  read(11, '(A)', advance = 'no') ch
  read(11, *) time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) total_time
  interval_time = total_time / time_steps
  read(11, '(A)', advance = 'no') ch
  read(11, *) N_statistic
  read(11, '(A)', advance = 'no') ch
  read(11, *) psi0
!  write(*, '(i7, f10.5, i7)') time_steps, total_time, N_statistic
  close(11)

! read spectral densities from files
  open(11, file = 'density.dat')
  read(11, *)
  do i = 1, N_omega
    read(11, '(f14.7)', advance = 'no') omega 
    do j1 = 1, N_basis; do j2 = 1, j1
      read(11, '(ES16.6)', advance = 'no') SP(j1, j2, i)
      if(j1 /= j2) SP(j2, j1, i) = SP(j1, j2, i)
    end do; end do
    read(11, *)
  end do
  close(11)
  
end subroutine

subroutine initialphi()

  use constants
  use spectral_density
  
  do i1 = 1, N_basis; do i2 = 1, i1
    do j = 1, N_omega
      call random_number(phi(i1, i2, j))
      phi(i1, i2, j) = phi(i1, i2, j) * 2.0d0 * pi
      if(i1 /= i2) phi(i2, i1, j) = phi(i1, i2, j)
    end do
  end do; end do
  
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
