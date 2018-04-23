subroutine get_Hamiltonian(t, Hamiltonian)

  use constants
  use spectral_density

  double precision, intent(in) :: t
  double precision, intent(out) :: Hamiltonian(n_basis, n_basis)
  double precision :: delta_H(n_basis, n_basis)
  double precision :: omega

  delta_H = 0.0d0
  
  do i = 1, N_basis
    do j = 1, N_omega
      omega = j * interval_omega
      delta_H(i, i) = delta_H(i, i) + S(j) * (phi(i, j, 1) * dcos(omega * t) + phi(i, j, 2) * dsin(omega * t))
!      delta_H(i, i) = delta_H(i, i) + dsqrt(2.0d0 * S(j) * hbar * interval_omega) * dcos(omega * t + phi(i, j, 1))
    end do
  end do
  
  Hamiltonian = H0 + delta_H

  Hamiltonian(2, 2) = -Hamiltonian(1, 1)
  
!  write(33, '(3f10.5)') t, Hamiltonian(1, 1), Hamiltonian(2, 2)
  
end subroutine

!!Ohmic spectral density
subroutine discretization() 

  use constants
  use spectral_density

  double precision :: omega
  
  interval_omega = omega_max / N_omega

  do i = 1, N_omega
    omega = i * interval_omega
    SP(i) = pi * 0.5d0 * alpha * omega * dexp(- omega / omega_c)
!quantum
    if(check_quantum) then
!      S(i) = SP(i) * (dexp(beta * hbar * omega) +  1.0d0) / (pi * (dexp(beta * hbar * omega) - 1.0d0))
    S(i) = dsqrt(SP(i) * hbar * interval_omega * (dexp(beta * hbar * omega) +  1.0d0) / (pi * (dexp(beta * hbar * omega) - 1.0d0 )))
!classic
    else
      S(i) = SP(i) * 2.0d0 / (pi * beta * hbar * omega)
    end if
!    write(22, '(2f14.7)') omega, S(i)
!    do j = 1, N_basis
!      call random_number(phi(j, i))
!      phi(j, i) = phi(j, i) * 2.0d0 * pi
!    end do
  end do
!stop

end subroutine