subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double precision, intent(out) :: Hamiltonian(N_basis, N_basis)
  double precision :: omega, force
  
  Hamiltonian = H0
  
  do i = 1, N_basis
    force = 0.0d0
    do j = 1, N_omega
      omega = j * interval_omega
      force = force + h(j) * dsqrt(2.0d0 * n_therm(j) + 1.0d0) * (phi(i, j, 1) * dcos(omega * t) + phi(i, j, 2) * dsin(omega * t))
    end do
    Hamiltonian(i, i) = Hamiltonian(i, i) + force
  end do
  
end subroutine

!!Ohmic spectral density
subroutine discretization() 

  use constants
  use spectral_density

  double precision :: SP, omega

  do i = 1, N_omega
    omega = i * interval_omega
    SP = 2.0d0 * eta * omega * omega_c / (omega ** 2.0d0 + omega_c ** 2.0d0)
    n_therm(i) = 1.0d0 / (dexp(beta * hbar * omega) - 1.0d0)
    h(i) = dsqrt(SP * hbar * interval_omega / pi)
  end do

end subroutine