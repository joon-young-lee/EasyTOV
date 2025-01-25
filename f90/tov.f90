program tov
 use tov_rk4
 implicit none
  ! tov solver that returns M-R relation of neutron stars
  ! Constants
  ! integer, parameter :: dp = selected_real_kind(15)
  integer(dp), parameter :: NN = 10 ** 5
  integer(dp), parameter :: coe = 4 * (NN-1)
  integer, parameter :: max_iterations =  10 ** 6 
  integer(8) :: n, k
  character(len=100) :: file_name
  real(dp), allocatable :: e_data(:), p_data(:), Coefficients(:, :) ! Coefficient of cubic-poly
  real(dp) :: y, p, d_p, d_r2, d_m, del_h, p_c, M, r, p_c1

  ! Read File
  ! print *, MeV_fm_to_km
  file_name ="eft_pnm32_000001_ldm_eos_s.dat" ! "eft_pnm32_000002_ldm_eos_s.dat" ! 'empirical_reduced.dat' 
  print *, 'Attempting to open file: ', file_name
  CALL read_eos(file_name, n, e_data, p_data) ! Defines n, e_data, p_data, in km 
  CALL CubicSpline_Coefficients(n, p_data, e_data, Coefficients)

  ! p_c = 500_dp
  p_c1 = 112_dp
  del_h = -1.e-5
  ! CALL rk4(n, p_data, Coefficients, del_h, p_c, max_iterations) !, M, r)
  CALL rk4(n, p_data, Coefficients, del_h, p_c1, max_iterations) !, M, r)

  do k = 1, 20, 1
    CALL rk4(n, p_data, Coefficients, del_h, (500 - 100) / 100 * DBLE(k) + 100, max_iterations) !, M, r)
  enddo
end program tov
