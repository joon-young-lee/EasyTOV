program tov
 use tov_rk4
 implicit none
  ! tov solver that returns M-R relation of neutron stars
  ! Constants
  ! integer, parameter :: dp = selected_real_kind(15)
  integer(dp), parameter :: NN = 10 ** 5
  integer(dp), parameter :: coe = 4 * (NN-1)
  integer(8), parameter :: max_iterations =  10 ** 6
  integer(8) :: k, LengthOfEoS, num, unit_number
  character(len=100) :: file_name
  real(dp), allocatable :: e_array(:), p_array(:), &
  Coefficients(:, :), M(:), R(:) ! Coefficient of cubic-poly
  real(dp) :: y, p_float, d_p, d_r2, d_m, del_h, p_c, &
  p_c1, EoS, p_start, p_final, &
  start_time, end_time, elapsed_time, M_0, R_0, diff

  ! Read File
  ! print *, MeV_fm_to_km
  file_name = "../eos/eft_pnm32_000098_ldm_eos_s.dat" 
  ! "eft_pnm32_000002_ldm_eos_s.dat" ! 'empirical_reduced.dat' 
  print *, 'Attempting to open file: ', file_name
  CALL read_eos(file_name, LengthOfEoS, p_array, e_array) 
  ! Defines n, e_data, p_data, in km 
  CALL CubicSpline_Coefficients(LengthOfEoS, p_array, e_array, Coefficients)

  p_c = 500.0_dp
  p_c1 = 200.0_dp
  del_h = -1.e-5

  CALL rk4(LengthOfEoS, p_array, Coefficients, del_h, 500.0_dp,&
   max_iterations, M_0, R_0)
  print *, M_0
  print *, R_0
  print *, LengthOfEoS
  p_start = 100.0_dp
  p_final = 1000.0_dp
  num = 100
  diff = (p_final - p_start) / DBLE(num)
  print *, diff
  allocate(M(num+1), R(num+1))
  ! Get the start time
  CALL cpu_time(start_time)
  do k = 1, num+1, 1 ! concurrent (k  = 1:num + 1)! concurrent (k = 1:num)
    CALL rk4(LengthOfEoS, p_array, Coefficients, del_h, &
    (p_final - p_start) / DBLE(num) * (DBLE(k) - 1.0_dp) + 100.0_dp, &
     max_iterations, M_0, R_0)
    M(k) = M_0
    R(k) = R_0
    ! print *, ' '
    print *, "--------------------"
    print *, M_0, "Mass in M_0"
    print *, R_0, "Radius in km"
    print *, "--------------------"
    print *, ' '
  enddo

  CALL cpu_time(end_time)
  elapsed_time = end_time - start_time
  ! Print results
    print *, "Elapsed CPU time (seconds): ", elapsed_time
    print *, "Central Pressure", p_start, "MeV/fm^3 ~ ", p_final, "MeV/fm^3"
  ! Open the file for writing (unit number 20, replace "output.txt" with your desired file name)
    unit_number = 20
    open(unit_number, file="output.txt", status="unknown", action="write")

    ! Write header to the file
    write(unit_number, '(A)') "Mass (M_0)    Radius (km)"

    ! Write data into two columns
    do k = 1, num+1
        write(unit_number, '(F10.5, F15.5)') M(k), R(k)
    end do

    ! Close the file
    close(unit_number)

    print *, "Data written to output.txt"
    print *, "Total ", num+1, " stars"
deallocate(M, R)
end program tov
