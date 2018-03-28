subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double precision, intent(out) :: Hamiltonian(N_basis, N_basis)
  double precision :: delta_H(N_basis, N_basis)

  delta_H = 0.0d0
  
  do i1 = 1, N_basis; do i2 = 1, N_basis
    do j = 1, N_omega
      delta_H(i1, i2) = delta_H(i1, i2) + 2.0d0 * dsqrt(S(i1, i2, j) * hbar * interval_omega) * dcos(j * interval_omega * t + phi(i1, i2, j))
    end do
  end do; end do
  
  Hamiltonian = H0 + delta_H

end subroutine

!!Ohmic spectral density
subroutine discretization()

  use constants
  use spectral_density

  double precision :: omega

  do i = 1, N_omega
    omega = i * interval_omega
!quantum
    if(check_quantum) then
      S(:, :, i) = SP(:, :, i) / (pi * ( 1 - dexp(-beta * hbar * omega)))
!classic
    else
      S(:, :, i) = SP(:, :, i) / (pi * beta * hbar * omega)
    end if
!    do j = 1, N_basis
!      call random_number(phi(j, i))
!      phi(j, i) = phi(j, i) * 2.0d0 * pi
!    end do
  end do
!stop

end subroutine