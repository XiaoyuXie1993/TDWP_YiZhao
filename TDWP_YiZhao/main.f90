!! time evolution of a spin-boson system using parallel programing via mpi (J. Chem. Phys. 2011, 135, 134110)
program spin_boson
  
  use time_evolution
  use mpi

  integer :: ierr, num_procs, my_id, n_traj_per_para
  integer, allocatable :: stat(:)
  double precision :: density
  double complex, allocatable :: psi(:, :, :)
  double precision, allocatable :: diff_density(:)
  
  call initial()
  allocate(stat(N_basis))
  allocate(diff_density(time_steps))
  call discretization()

!  write(*, '(i4, 2f14.7)') time_steps, interval_time, total_time

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  n_traj_per_para = (N_statistic - 1) / num_procs + 1
  allocate(psi(n_traj_per_para, time_steps, N_basis))
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  diff_density = 0.0d0
  do i = 1, n_traj_per_para
!  write(*, *) i, my_id
    call init_random_seed(my_id)
    call initialphi()
!    call Runge_Kutta(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    call Chebyshev(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    if(my_id /= 0) then
      call MPI_SEND(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, 0, i, MPI_COMM_WORLD, ierr)
    else
      do k = 1, time_steps
        diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) - density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
      end do
      do j = 1, num_procs - 1
        call MPI_RECV(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, j, i, MPI_COMM_WORLD, stat, ierr)
        do k = 1, time_steps
          diff_density(k) = diff_density(k) + (density(psi(i, k, 1)) - density(psi(i, k, 2))) / (n_traj_per_para * num_procs)
        end do
      end do
    end if
  end do
  
  call MPI_FINALIZE(ierr)

  if(my_id == 0) then
    open(22, file = 'result.dat')
    write(22, '(2f14.7)') 0.0d0, density(psi0(1)) - density(psi0(2))
    time = interval_time
    do i = 1, time_steps
      write(22, '(2f14.7)') time, diff_density(i)
      time = time + interval_time
    end do
  end if
  close(11)

end program 

! density (population) of a electronic state  
double precision function density(coeff)
  
  double complex, intent(in) :: coeff

  density = real(coeff) ** 2 + imag(coeff) ** 2

end function