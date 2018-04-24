!! Runge-Kutta for differential equations
!program Runge_Kutta

!  double precision :: y
!  double precision :: interval
!  double precision :: t

!  y = 1.0d0
!  t = 0.0d0
!  interval = 0.10d0

!  do i = 1, 11
!    write(*, '(2f12.8)') t, y
!    call update(y, t, interval)
!    t = t + interval
!  end do

!end program

subroutine Runge_Kutta(n_time, t_total, n_eq, y0, y)

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: y0(n_eq)
  double complex, intent(out) :: y(n_time, n_eq)
  double precision :: time, interval

  interval = t_total / n_time
  
  do i = 1, n_time
    time = i * interval
    if(i == 1) then
      call solve(n_eq, y0, time, interval, y(i, :))
    else
      call solve(n_eq, y(i - 1, :), time, interval, y(i, :))
    end if
  end do

end subroutine

!!using fourth-order Runge-Kutta methods
subroutine solve(n_eq, y, t, interval, newy)

  integer, intent(in) :: n_eq
  double complex, intent(in) :: y(n_eq)
  double precision, intent(in) :: t, interval
  double complex, intent(out) :: newy(n_eq)
  double complex, allocatable :: K1(:), K2(:), K3(:), K4(:)

  allocate(K1(n_eq), K2(n_eq), K3(n_eq), K4(n_eq))
  call diff_func(n_eq, K1, y, t)
  call diff_func(n_eq, K2, y + 0.5d0 * interval * K1, t + 0.5d0 * interval)
  call diff_func(n_eq, K3, y + 0.5d0 * interval * K2, t + 0.5d0 * interval)
  call diff_func(n_eq, K4, y + interval * K3, t + interval)

  newy = y + interval / 6.0d0 * (K1 + 2.0d0 * K2 + 2.0d0 * K3 + K4)
  
  deallocate(K1, K2, K3, K4)

end subroutine

! time-dependent Schrodinger equation
subroutine diff_func(n_eq, dy, y, t)
  
  use constants

  integer, intent(in) :: n_eq
  double complex, intent(in) :: y(n_eq)
  double precision, intent(in) :: t
  double complex, intent(out) :: dy(n_eq)
  double precision, allocatable :: Hamiltonian(:, :)
  double complex :: alpha0, beta0

  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = 0.0d0

  allocate(Hamiltonian(n_eq, n_eq))

  call get_Hamiltonian(t, Hamiltonian)
  
  call dzgemm('N', 'N', n_eq, n_eq, n_eq, alpha0, Hamiltonian, n_eq, y, n_eq, beta0, dy, n_eq)

  deallocate(Hamiltonian)

  dy = dcmplx(0.0d0, -1.0d0) * dy / hbar

end subroutine