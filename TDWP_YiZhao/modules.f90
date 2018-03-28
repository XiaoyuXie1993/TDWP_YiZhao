module constants

!! constants
  double precision :: pi, hbar, kB
  parameter(pi = 3.1415926535897932)
! reduced Planck constant with eV*fs 
  parameter(hbar = 0.6582119514)
! Boltzmann constant with eV/K
  parameter(kB = 8.6173303e-5)

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
  
  integer :: N_omega
! frequency with eV
  double precision :: interval_omega
! temperature of the system with K
  double precision :: temperature
! beta = 1/ (kB * T) with eV^{-1}
  double precision :: beta
  logical :: check_quantum
! spectral density with eV
  double precision, allocatable :: SP(:, :, :), S(:, :, :)
  double precision, allocatable :: phi(:, :, :)

end module

module time_evolution

  use Hamiltonian_electronic

!! parameters of initial electronic states and time-dependent simulation
  double complex, allocatable :: psi0(:)
!! time evolution with fs
  integer :: time_steps
  double precision :: interval_time, total_time
!! parameters of statistic average
  integer :: N_statistic

end module