module constants

!! constants
  double precision :: pi, hbar
  parameter(pi = 3.1415926535897932)
  parameter(hbar = 1.0d0)

end module

module Hamiltonian_electronic

!! parameters of electronic Hamiltonian (two level system with spin-boson model)
! basis set of electronic states
  integer :: N_basis
  parameter(N_basis = 2)
! Hamiltonian elements
  double precision :: H0(n_basis, n_basis)
  parameter(H0 = (/0.0d0, 1.0d0, 1.0d0, 0.0d0/))

end module

module spectral_density

  use Hamiltonian_electronic

!! parameters of Ohmic spectral density J(omega) = pi / 2 * alpha * omega * exp(-omega / omega_c)
  double precision :: alpha, omega_c, beta
!  parameter(alpha = 1.2d0)
!  parameter(omega_c = 2.5d0)
!  parameter(beta = 0.2d0)
!  parameter(beta = 5.0d0)
!  parameter(alpha = 0.1d0)
!  parameter(omega_c = 7.5d0)
!! parameters for discretization of spectral density
  integer :: N_omega
  parameter(N_omega = 10000)
  double precision :: interval_omega, omega_max
  parameter(omega_max = 100.0d0)
  double precision :: SP(N_omega), S(N_omega)
  double precision :: phi(N_basis, N_omega)
  logical :: check_quantum

end module

module time_evolution

  use Hamiltonian_electronic

!! parameters of initial electronic states and time-dependent simulation
  double complex :: psi0(N_basis)
  parameter(psi0 = (/1.0d0, 0.0d0/))
  integer :: time_steps
!  parameter(time_steps = 1000)
  double precision :: interval_time, total_time
!  parameter(total_time = 15.0d0)
!! parameters of statistic average
  integer :: N_statistic
!  parameter(N_statistic = 500)

end module