!! time evolution of a spin-boson system using parallel programing via mpi (J. Chem. Phys. 2011, 135, 134110)
program dynamics
  
  use time_evolution
  use mpi

  integer :: ierr, num_procs, my_id, n_traj_per_para
  integer, allocatable :: stat(:)
  double precision :: density
  double complex, allocatable :: psi(:, :, :)
  double precision, allocatable :: population(:, :)
  
  call initial()
  allocate(stat(N_basis))
  allocate(population(time_steps, N_basis))
  call discretization()

!  write(*, '(i4, 2f14.7)') time_steps, interval_time, total_time

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  n_traj_per_para = (N_statistic - 1) / num_procs + 1
  allocate(psi(n_traj_per_para, time_steps, N_basis))
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  population = 0.0d0
  do i = 1, n_traj_per_para
!  write(*, *) i, my_id
    call init_random_seed(my_id)
    call initialphi()
!    call Runge_Kutta(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    call Chebyshev(time_steps, total_time, N_basis, psi0, psi(i, :, :))
    if(my_id /= 0) then
      call MPI_SEND(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, 0, i, MPI_COMM_WORLD, ierr)
    else
      do k = 1, time_steps; do l = 1, N_basis
        population(k, l) = population(k, l) + density(psi(i, k, l)) / (n_traj_per_para * num_procs)
      end do; end do
      do j = 1, num_procs - 1
        call MPI_RECV(psi(i, :, :), time_steps * N_basis, MPI_DOUBLE_COMPLEX, j, i, MPI_COMM_WORLD, stat, ierr)
        do k = 1, time_steps; do l = 1, N_basis
          population(k, l) = population(k, l) + density(psi(i, k, l)) / (n_traj_per_para * num_procs)
        end do; end do
      end do
    end if
  end do
  
  call MPI_FINALIZE(ierr)

  if(my_id == 0) then
    open(22, file = 'result.dat')
    write(22, '(f14.7)', advance = 'no') 0.0d0
    do i = 1, N_basis
      write(22, '(f14.7)', advance = 'no') density(psi0(i))
    end do
    write(22, *)
    do i = 1, time_steps
      time = i * interval_time
      write(22, '(f14.7)', advance = 'no') time
      do j = 1, N_basis
        write(22, '(f14.7)', advance = 'no') population(i, j)
      end do
      write(22, *)
    end do
  end if
  close(11)

end program 

! density (population) of a electronic state  
double precision function density(coeff)
  
  double complex, intent(in) :: coeff

  density = real(coeff) ** 2 + imag(coeff) ** 2

end function