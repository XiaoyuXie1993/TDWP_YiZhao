
!! solve TDSE using Chebyshev polynomials methods
subroutine Chebyshev(n_time, t_total, n_eq, y0, y)

  use constants

  integer, intent(in) :: n_time, n_eq
  double precision, intent(in) :: t_total
  double complex, intent(in) :: y0(n_eq)
  double complex, intent(out) :: y(n_time, n_eq)
  double precision :: time, interval, a, b
  double precision :: tx(n_eq, n_eq)
  double precision :: jt
  double precision :: Hamiltonian(n_eq, n_eq)
  double complex :: U(n_eq, n_eq)

  interval = t_total / n_time
  
  do i = 1, n_time
    time = i * interval
    call get_Hamiltonian(time, Hamiltonian)
    call diagonal(n_eq, Hamiltonian, a, b)
    jt = - b * interval / hbar
    tx = Hamiltonian
    do j = 1, n_eq
      tx(j, j) = tx(j, j) - a
    end do
    tx = tx / b
    call expansion(n_eq, tx, jt, U)
    U = U * cdexp(dcmplx(0.0d0, -1.0d0) * interval / hbar)
    if(i == 1) then
      call zgemm('N', 'N', n_eq, 1, n_eq, 1.0d0, U, n_eq, y0, n_eq, 0.0d0, y(i, :), n_eq)
    else
      call zgemm('N', 'N', n_eq, 1, n_eq, 1.0d0, U, n_eq, y(i - 1, :), n_eq, 0.0d0, y(i, :), n_eq)
    end if
  end do
      
  
end subroutine
  
!! calculate exp(i * x * time) using Chebyshev polynomials (exp(ixt) = sum_n h * i ^ n * Jn(x) * Tn(t)£©
subroutine expansion(n_matrix, x, time, expixt)

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: x(n_matrix, n_matrix)
  double precision, intent(in) :: time
  double complex, intent(out) :: expixt(n_matrix, n_matrix)
  integer :: truncation0, truncation
  parameter(truncation0 = 100000)
  double precision :: J(truncation0)
  double precision, allocatable :: T(:, :, :)
  double precision :: h
  double complex :: pexpansion(n_matrix, n_matrix)

  call Bessel_function(time, truncation0, truncation, J)
  allocate(T(truncation, n_matrix, n_matrix))
  call Cheyshev_polynomial(n_matrix, x, truncation, T)
  
!  do i = 1, truncation
!    write(*, '(i5, 4f10.5)') i, J(i), T(i), (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
!  end do 
  
  expixt = 0.0d0
  do i = 1, truncation
    if(i == 1) then
      h = 1.0d0
    else
      h = 2.0d0
    end if
    pexpansion = h * J(i) * T(i, :, :) * (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
    expixt = expixt + pexpansion
  end do
!  write(*, *) expansion
  deallocate(T)

end subroutine

!! Cheyshev polynomials Tn(x)
subroutine Cheyshev_polynomial(n_matrix, x, n, results)

  integer, intent(in) :: n_matrix, n
  double precision, intent(in) :: x(n_matrix, n_matrix)
  double precision, intent(out) :: results(n, n_matrix, n_matrix)

  results = 0.0d0
  do i = 1, n_matrix
    results(1, i, i) = 1.0d0
  end do
  results(2, :, :) = x
  
  do i = 3, n
    call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, x, n_matrix, results(i - 1, :, :), n_matrix, 0.0d0, results(i, :, :), n_matrix)
    results(i, :, :) = 2.0d0 * results(i, :, :) - results(i - 2, :, :)
  end do

end subroutine

!! Bessel functions Jn(x) with a truncation (output)
subroutine Bessel_function(x, n, truncation, results)

  integer, intent(in) :: n
  double precision, intent(in) :: x
  integer, intent(out) :: truncation
  double precision, intent(out) :: results(n)
  integer :: truncation0
  parameter(truncation0 = 100000)
  double precision :: factorial, threshold, pBessel
  parameter(threshold = 1e-35)
  
  truncation = n
  results = 0.0d0
  do i = 1, n
    do j = 1, truncation0
      if(dabs(pBessel(x, i - 1, j - 1)) < threshold) exit
      results(i) = results(i) + pBessel(x, i - 1, j - 1)
    end do
    if(dabs(results(i)) < threshold) exit
  end do

  truncation = i

end subroutine

subroutine diagonal(n_Matrix, Matrix, ea, eb)

  integer, intent(in) :: n_Matrix
  double precision, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double precision, intent(out) :: ea, eb
  double precision :: eigenvalue(n_Matrix)
  integer :: info, lwork, liwork, lwmax
  parameter(lwmax = 100000)
  integer :: iwork(lwmax)
  double precision :: work(lwmax)
  
  lwork = -1
  liwork = -1
  
  call dsyevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  liwork = min(lwmax, iwork(1))
  
  call dsyevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  ea = 0.5d0 * (eigenvalue(n_Matrix) + eigenvalue(1))
  eb = 0.5d0 * (eigenvalue(n_Matrix) - eigenvalue(1))

end subroutine

double precision function pBessel(x, n, m)

  double precision, intent(in) :: x
  integer, intent(in) :: n, m

  pBessel = 1.0d0

  do i = 1, n + m
    pBessel = pBessel * (x * 0.50d0) / i
  end do
  
  do i = 1, m
    pBessel = pBessel * (-x * 0.50d0) / i
  end do
  
end function