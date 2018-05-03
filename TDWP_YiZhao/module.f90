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
! Hamiltonian elements
  double precision, allocatable :: H0(:, :)

end module

module spectral_density

  use Hamiltonian_electronic
  
!! parameters of Debye-Drude spectral density J(omega) = eta * omega * omega_c / (omega ** 2.0d0 + omega_c ** 2.0d0)
  double precision :: eta, omega_c
  double precision :: beta
!! parameters for discretization of spectral density
  integer :: N_omega
  double precision :: interval_omega
  double precision, allocatable :: n_therm(:), h(:)
  double precision, allocatable :: phi(:, :)

end module

module time_evolution

  use Hamiltonian_electronic

!! parameters of time-dependent simulation and initial electronic states
  integer :: time_steps
  double precision :: interval_time, total_time
  double complex, allocatable :: psi0(:)
!! parameters of statistic average
  integer :: N_statistic

end module