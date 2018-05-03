subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double precision, intent(out) :: Hamiltonian(N_basis, N_basis)
  double precision, allocatable :: x_operator(:, :)
  double precision :: force, omega

  allocate(x_operator(N_basis, N_basis))
  x_operator = 0.0d0
  x_operator(1, 1) = 1.0d0
  x_operator(2, 2) = -1.0d0
  
  force = 0.0d0
  do i = 1, N_omega
    omega = i * interval_omega
    force = force + h(i) * dsqrt(2.0d0 * n_therm(i) + 1.0d0) * (phi(i, 1) * dcos(omega * t) + phi(i, 2) * dsin(omega * t))
  end do
  
  Hamiltonian = H0 + force * x_operator

  deallocate(x_operator)
  
end subroutine

!!Debye-Drude spectral density
subroutine discretization() 

  use constants
  use spectral_density

  double precision :: SP, omega

  do i = 1, N_omega
    omega = i * interval_omega
    SP = eta * omega * omega_c / (omega ** 2.0d0 + omega_c ** 2.0d0)
    n_therm(i) = 1.0d0 / (dexp(beta * hbar * omega) - 1.0d0)
    h(i) = dsqrt(SP * hbar * interval_omega / pi)
  end do

end subroutine