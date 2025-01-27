module tov_rk4
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: c = 2.99792458e10_dp
    real(dp), parameter :: c2 = c**2
    real(dp), parameter :: G = 6.67259e-8_dp
    real(dp), parameter :: pi = 3.141592_dp
    real(dp), parameter :: MeV = 1.602179e-6_dp
    real(dp), parameter :: fm = 1.0e-13_dp
    real(dp), parameter :: MeV_fm_to_cgs = 1.602179e+33
    real(dp), parameter :: erg_cm_to_MeV_fm = 6.2414999e-34_dp
    real(dp), parameter :: M0 = 1.9884e33
    real(dp), parameter :: g_to_km = G/c2 * 1.e-5
    real(dp), parameter :: km_to_M0 = 1 / (M0 * g_to_km)
    real(dp), parameter :: MeV_to_km = 1.60218 * 1.e-6 * g_to_km / c2
    real(dp), parameter :: fm_to_km = 1.e-18
    real(dp), parameter :: MeV_fm_to_km = MeV_to_km / fm_to_km ** 3
contains       
subroutine  CubicSpline_Coefficients(n, x_array, y_array, Coefficients)
                implicit none
                integer, parameter :: dp = selected_real_kind(15) 
                integer(8), intent(in) :: n ! input is size of x_array
                integer(8) :: nn, i, j, nrhs, info, lda, ldb
                integer, dimension(4 * (n - 1)) :: ipiv        ! Pivot indices
                real(dp), dimension(n) :: x_array, y_array ! Input Arrays
                real(dp), allocatable :: A(:, :)
                real(dp), allocatable :: Coefficients(:, :)! , CC(:, :)
!------------------------------------------------------------------------------------------                
!------------------------------------------------------------------------------------------                
! Check if x_array is strictlly increasing
do i = 1, n-1, 1
  if (x_array(i+1) < x_array(i)) then
    print *,"x_array is not strictly increasing! " 
    exit
  else
    cycle
  end if 
end do
!------------------------------------------------------------------------------------------                
!------------------------------------------------------------------------------------------                
                nrhs = 1                
                nn = 4 * (n - 1) ! number of coefficients that needed to be decided
                allocate(A(nn, nn), Coefficients(nn, 1))!, CC(nn, 1))
                do concurrent (i = 1:nn, j =1:nn)
                                A(i, j) = 0.0_dp
                end do

                do concurrent (i = 1:nn)
                        Coefficients(i, 1) = 0.0_dp
                end do 
                
                ! print *, shape(A)
                ! print *, size(B)
        
                ! Matching Left Function Value, 1 to n - 1 (total: n - 1)
                do i = 1, n - 1, 1
                        do j = 1 + 4 * (i - 1), 4 + 4 * (i - 1), 1
                                if (mod(j, 4) == 1) then
                                        A(i, j) = x_array(i) ** 3
                                elseif (mod(j, 4) == 2) then
                                        A(i, j) = x_array(i) ** 2 
                                elseif (mod(j, 4) == 3) then
                                        A(i, j) = x_array(i) ** 1 
                                elseif (mod(j, 4) == 0) then
                                        A(i, j) = x_array(i) ** 0
                                endif
                       end do
                end do
! print *, A(1, :)
        
                ! Matching Right Function Value, n to 2n - 2 (total: n - 1)
                do i = n, 2 * n - 2, 1
                        do j = 1 + 4 * (i - n), 4 + 4 * (i - n), 1
                                if (mod(j, 4) == 1) then
                                        A(i, j) = x_array(i + 1 - (n - 1)) ** 3 
                                elseif (mod(j, 4) == 2) then 
                                        A(i, j) = x_array(i + 1 - (n - 1)) ** 2 
                                elseif (mod(j, 4) == 3) then
                                        A(i, j) = x_array(i + 1 - (n - 1)) ** 1 
                                elseif (mod(j, 4) == 0) then
                                        A(i, j) = x_array(i + 1 - (n - 1)) ** 0
                                endif
                       end do
                end do
! print *, A
                ! Matching Derivative, 2n - 1 to 3n - 4 (total: n - 2)
        do i = 2 * n - 1, 3 * n - 4, 1
                do j = 1 + 4 * (i - (2 * n - 1)), 8 + 4 * (i - (2 * n - 1)), 1 ! 4 spaces
                        if (mod(j, 8) == 1) then
                                A(i, j) = 3 * x_array(i + 1 - (2 * n - 2)) ** 2
                        elseif (mod(j, 8) == 2) then
                                A(i, j) = 2 * x_array(i + 1 - (2 * n - 2)) ** 1
                        elseif (mod(j, 8) == 3) then
                                A(i, j) = 1.0_dp
                        elseif (mod(j, 8) == 5) then
                                A(i, j) = -3 * x_array(i + 1 - (2 * n - 2)) ** 2
                        elseif (mod(j, 8) == 6) then
                                A(i, j) = -2 * x_array(i + 1 - (2 * n - 2)) ** 1
                        elseif (mod(j, 8) == 7) then
                                A(i, j) = -1.0_dp
                        else
                                A(i, j) = 0.0_dp
                        endif
                end do
        end do
! print *, A
                ! Matching Second Derivative, 3n - 3 to 4n - 6
        do i = 3 * n - 3, 4 * n - 6, 1
                do j = 1 + 4 * (i - (3 * n - 3)), 8 + 4 * (i - (3 * n - 3)), 1 
                        if (mod(j, 8) == 1) then
                                A(i, j) = 6 * x_array(i + 1 - (3 * n - 4)) ** 1
                        elseif (mod(j, 8) == 2) then
                                A(i, j) = 2.0_dp
                        elseif (mod(j, 8) == 5) then
                                A(i, j) = -6 * x_array(i + 1 - (3 * n - 4)) ** 1
                        elseif (mod(j, 8) == 6) then
                                A(i, j) = -2.0_dp
                        else
                                A(i, j) = 0.0_dp
                        endif
                end do
        end do      

                ! Matching Boundary
                A(nn - 1, 1) = 6 * x_array(1) ! A(nn - 1, 1) = 6 * x_array(1)
                A(nn - 1, 2) = 2.0_dp
                A(nn, nn - 3) = 6 * x_array(n)
                A(nn, nn - 2) = 2.0_dp 
                ! do i = 3, nn, 1
                     ! A(nn - 1, i) = 0.0_dp
                ! end do
                ! do i = 1, nn - 4, 1
                     ! A(nn, i) = 0.0_dp
                ! end do
! print *, A(nn, nn-2) 
                ! A(nn, nn-1) = 0.0_dp
                ! A(nn, nn) = 0.0_dp
                ! print *, A
                do i = 1, n - 1, 1
                      Coefficients(i, 1) = y_array(i)
                end do
                
                do i = n, 2 * n - 1, 1
                      Coefficients(i, 1) = y_array(i + 1 - (n - 1))
                end do
        ! CC = Coefficients
        ! n: Order of the matrix.
        ! nrhs: Number of right-hand sides (1 for Ax=b).
        ! A: The input matrix; overwritten with L and U factors.
        ! ipiv: Pivot indices array.
        ! b: Overwritten with the solution vector x
        ! info: Output status (0 for success).
        ! print *, A(1, :)
        ! print *, Coefficients
        ! print *, "n: ", n
        ! print *, "nn: " , nn
        ! print *, "NRHS: ", nrhs
        ! print *, A
        ! print *, "nn: ", nn
        ! print *, ipiv
        ! print *, "info: ", info
        call dgesv(nn, nrhs, A, nn, ipiv, Coefficients, nn, info)
        ! Check the results
        ! print *, A
        ! if (info == 0) then
        !   print *, "The Coefficients were properly obtained."
        !   ! print *, "Solution"
        !   ! print *, Coefficients
        !   ! print *, info      
        ! else if (info > 0) then
        !   print *, "Factorization Failed"
        !   print *, info
        
        ! else
        !   print *, "Invalid Arguments"
        
        ! end if
 ! print *, Coefficients
        
        deallocate(A)
end subroutine

! subroutine polynomial_array(x_array, coefficients, x_input, y_output)
          
       
! end subroutine

pure subroutine polynomial_float(n, x_array, Coefficients, x, y)
           use ieee_arithmetic
           implicit none
           integer, parameter :: dp = selected_real_kind(15)
           ! Double precision (15 digits, large exponent range) 
           integer(8), intent(in) :: n
           integer ::  i
           real(dp), dimension(n), intent(in) :: x_array
           real(dp), dimension(4 * (n - 1), 1), intent(in) :: Coefficients
           real(dp) :: a, b, c, d
           real(dp), intent(in) :: x 
           real(dp), intent(out) :: y
           logical :: place, NOT_CAL
           NOT_CAL = (MAXVAL(x_array) < x) .OR. (x < MINVAL(x_array))
              
      
           do i = 1, n - 1, 1
             place = (x_array(i) < x) .AND. (x <= x_array(i+1)) 
              if (NOT_CAL) then
                !   print *, "Max Pressure is: ", MAXVAL(x_array)
                !   print *, "Min Pressure is: ", MINVAL(x_array)
                !   print *, "But", x, "is inputed!"
                  y = ieee_value(y, ieee_quiet_nan) 
                  exit
                  ! print *, "Invalid x inputed to interpolation!!"
              else if (place) then
                  a = Coefficients(4 * (i - 1) + 1, 1)
                  b = Coefficients(4 * (i - 1) + 2, 1)
                  c = Coefficients(4 * (i - 1) + 3, 1)
                  d = Coefficients(4 * (i - 1) + 4, 1)
                  
                  y = a * x ** 3 + b * x ** 2 + c * x + d
                  exit
              else 
                  cycle
              end if    
            end do
end subroutine
subroutine read_eos_NoCrust(file_name, LengthOfEoS, p_data, e_data)
    implicit none
    character(len=100), intent(in) :: file_name
    integer :: i, ios, k, m, n
    integer(8), intent(out) :: LengthOfEoS
    real(dp), allocatable :: e_data(:), p_data(:)
    ! real(dp), allocatable, intent(out):: e_data(:), p_data(:)
     
    open(unit=10, file=file_name, status='old', action='read')
    ! count the number of data
    n = 0
    do 
      read(10, *, iostat=ios) ! Try to read a line
      if (ios /= 0) exit
      n = n + 1 ! Count the number of lines
    end do
    
    close(10)
   
    open(unit=10, file=file_name, status='old', action='read')
   
    allocate(e_data(n), p_data(n))
    
    ! Read the data into the arrays
    do i = 1, n
       read(10, *) p_data(n-i+1), e_data(n-i+1)
    end do
    
    e_data = e_data * MeV_fm_to_km
    p_data = p_data * MeV_fm_to_km
    close(10) 

end subroutine

subroutine read_eos(file_name, LengthOfEoS, p_array, e_array)
    implicit none
    character(len=100), intent(in) :: file_name
    integer :: i, ios, k, m, n
    integer(8), intent(out) :: LengthOfEoS
    real(dp), allocatable :: e_data(:), p_data(:), BPS_P(:), BPS_E(:)
    real(dp), allocatable, intent(out) :: e_array(:), p_array(:)
    ! real(dp), allocatable, intent(out):: e_data(:), p_data(:)
     
    open(unit=10, file=file_name, status='old', action='read')
    ! count the number of data
    n = 0
    do 
      read(10, *, iostat=ios) ! Try to read a line
      if (ios /= 0) exit
      n = n + 1 ! Count the number of lines
    end do
    
    close(10)
   
    open(unit=10, file=file_name, status='old', action='read')
   
    allocate(e_data(n), p_data(n))
    
    ! Read the data into the arrays
    do i = 1, n
       read(10, *) p_data(n-i+1), e_data(n-i+1)
    end do
    
    ! e_data = e_data * MeV_fm_to_km
    ! p_data = p_data * MeV_fm_to_km
    close(10) 

        ! Read BPS EoS and include BPS
        open(unit=10, file="BPS.txt", status='old', action='read')
        ! count the number of data
        k = 0
        do 
        read(10, *, iostat=ios) ! Try to read a line
        if (ios /= 0) exit
        k = k + 1 ! Count the number of lines
        end do
        
        close(10)

        open(unit=10, file="BPS.txt", status='old', action='read')
        ! Skip the first line
        read(10, *)
        k = k-1
        allocate(BPS_P(k), BPS_E(k))
        
        ! Read the data into the arrays
        do i = 1, k, 1
        read(10, *) BPS_E(i), BPS_P(i)
        end do
        
        close(10)

        BPS_E = BPS_E * c2 * erg_cm_to_MeV_fm
        BPS_P = BPS_P * erg_cm_to_MeV_fm

        i = 0
        do 
                if (i > SIZE(e_data)) then
                  print *, "Data EoS is not low enough! \n"
                  exit

                else if (e_data(i) >= MAXVAL(BPS_E)) then
                 m = i
                 ! print(m)
                 exit
                else if (e_data(i) < MAXVAL(BPS_E)) then
                 i = i + 1

                endif
        enddo
        print *, BPS_E(k)!  / (c2 * erg_cm_to_MeV_fm)
        print *, e_data(m)
        
        print *, 'm:', m
        print *, 'n: ', n
        print *, 'k: ', k
        LengthOfEoS = n+k-m
        print *, 'Length of EoS: ', LengthOfEoS
        allocate(p_array(LengthOfEoS), e_array(LengthOfEoS))
        do i = 1, k-1, 1
                e_array(i) = BPS_E(i)
                p_array(i) = BPS_P(i)
        enddo
        do i = k, LengthOfEoS, 1
                e_array(i) = e_data(m + i - k)
                p_array(i) = p_data(m + i - k)
        enddo
        print *, p_array(SIZE(p_array))
        print *, p_array(SIZE(e_array))
        ! print *, "Minimum Pressure in EoS: ", MINVAL(p), "MeV/fm^3"
        ! print *, "Maximum Pressure in EoS: ", MAXVAL(p), "MeV/fm^3"
        e_array = e_array * MeV_fm_to_km
        p_array = p_array * MeV_fm_to_km
        do i = 1, SIZE(p_array) - 1, 1
          if (p_array(i) < p_array(i+1)) then
            if (i == SIZE(p_array) - 1) then
              ! print *, "p_array is increasing"
            endif
            cycle
          else if (p_array(i) > p_array(i+1)) then
            ! print *, "p_array is not strictly increasing"
            exit
          else
            ! print *, "p_array is increasing"
          endif

        end do
        
        ! print *, MAXVAL(p)
        ! print *, "Minimum Pressure in EoS: ", MINVAL(p), "km"
        ! print *, "Maximum Pressure in EoS: ", MAXVAL(p), "km"
        deallocate(e_data, p_data, BPS_E, BPS_P)
end subroutine

! subroutine rk4(p_0, dh, nn)


! end subroutine rk4

pure subroutine dp_dh(LengthOfEoS, x_array, Coefficients, p, d_p)
    use ieee_arithmetic
    implicit none
    integer(8), intent(in) :: LengthOfEoS
    real(dp), dimension(LengthOfEoS), intent(in) :: x_array
    real(dp), intent(in) :: p
    real(dp), dimension(4 * (LengthOfEoS-1), 1), intent(in) :: Coefficients
    real(dp) :: EoS
    real(dp), intent(out) :: d_p
 call  polynomial_float(LengthOfEoS, x_array, Coefficients, p, EoS)
    d_p = p + EoS
    ! print *, "d_p: ", d_p
!     if (d_p == ieee_value(d_p, ieee_quiet_nan)) then
!       print *, "dp_dh NAN"
!     endif

end subroutine 

pure subroutine dr2_dh(LengthOfEoS, x_array, Coefficients, r2, p, m, d_r2)
    use ieee_arithmetic
    implicit none
    integer :: i, ios
    integer(8), intent(in) :: LengthOfEoS
    real(dp), dimension(LengthOfEoS), intent(in) :: x_array
    real(dp), intent(in) :: r2, p, m
    real(dp), dimension(4 * (LengthOfEoS-1), 1), intent(in) :: Coefficients
    real(dp), intent(out) :: d_r2
    d_r2 = -2 * r2 * (sqrt(r2) - 2 * m) / (m + 4 * pi * p * r2 ** (3.0/2.0))
    ! print *, "d_r2: ", d_r2
    
!     if (ieee_is_nan(d_r2)) then
!       print *, "dr2_dh NAN"
!     endif

end subroutine


pure subroutine dm_dh(LengthOfEoS, x_array, Coefficients, r2, p, m, d_m)
    use ieee_arithmetic
    implicit none
    integer(8), intent(in) :: LengthOfEoS
    real(dp), dimension(LengthOfEoS), intent(in) :: x_array
    real(dp), intent(in) :: r2, p, m
    real(dp), dimension(4 * (LengthOfEoS-1), 1), intent(in) :: Coefficients
    real(dp) :: EoS
    real(dp), intent(out) :: d_m
 call  polynomial_float(LengthOfEoS, x_array, Coefficients, p, EoS)
    d_m = -4 * pi * EoS * (r2 ** (3.0/2.0)) * (sqrt(r2) - 2 * m) /&
     (m + 4 * pi * p * r2 ** (3.0_dp/2.0_dp))
    ! print *, "d_m: ", d_m
!     if (ieee_is_nan(d_m)) then
!       print *, "dm_dh NAN"
!     endif
end subroutine

pure subroutine  rk4(LengthOfEoS, x_array, Coefficients, del_h, p_c,&
 max_iterations, Mass, Radius) ! Coefficients specifies EoS, p_c in MeV/fm^3
  USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  integer(8), intent(in) :: LengthOfEoS, max_iterations
  real(dp), dimension(LengthOfEoS), intent(in) :: x_array
  real(dp), dimension(4 * (LengthOfEoS-1), 1), intent(in) :: Coefficients
  integer :: i
  real(dp), intent(in) :: del_h, p_c
  real(dp), intent(out) :: Mass, Radius
  real(dp) :: p, m, r2, EoS, d_p, d_m, d_r2, &
                    k1_p, k1_m, k1_r2, &
                    p_2, m_2, r_2, &
                    k2_p, k2_m, k2_r2, &
                    p_3, m_3, r_3, &
                    k3_p, k3_m, k3_r2, &
                    p_4, m_4, r_4, &
                    k4_p, k4_m, k4_r2
  ! nan_value = IEEE_VALUE(0.0, IEEE_QUIET_NAN)                   
  ! print *, "MeV_fm_to_km: ", MeV_fm_to_km
  ! print *, "Central Pressure: ", p_c, "MeV/fm^3"
  p = p_c * MeV_fm_to_km
  ! print *, p, "Central Pressure in km"
  ! print *, x_array  
  CALL  polynomial_float(LengthOfEoS, x_array, Coefficients, p, EoS)

  r2 = 0.0_dp 
  m = 0.0_dp
  d_m = 0.0_dp
  CALL dp_dh(LengthOfEoS, x_array, Coefficients, p, k1_p)  
  k1_r2 = -3 / (2 * pi * (3 * p + EoS))
  
  p_2 = p + del_h * k1_p / 2
  r_2 = r2 + del_h * k1_r2 / 2
  CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_2, k2_p) 
  CALL polynomial_float(LengthOfEoS, x_array, Coefficients, p_2, EoS)
  k2_r2 = -3 / (2 * pi * (3 * p_2 + EoS))
  ! print *, k2_r2
  p_3 = p + del_h * k2_p / 2
  r_3 = r2 + del_h * k2_r2 / 2
  CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_3, k3_p)
  CALL polynomial_float(LengthOfEoS, x_array, Coefficients, p_3, EoS)
  k3_r2 = -3 / (2 * pi * (3 * p_3 + EoS))
  ! print *, k3_r2  
  p_4 = p + del_h * k3_p
  r_4 = r2 + del_h * k3_r2
  CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_4, k4_p)
  CALL polynomial_float(LengthOfEoS, x_array, Coefficients, p_4, EoS)
  k4_r2 = -3 / (2 * pi * (3 * p_4 + EoS))
  ! print *, k4_r2
  d_p = (k1_p + 2 * k2_p + 2 * k3_p + k4_p) * del_h / 6
  d_r2 = (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) * del_h / 6
  ! print *, d_r2,"d_r2"
  p = p + d_p
  r2 = r2 + d_r2
  
  CALL polynomial_float(LengthOfEoS, x_array, Coefficients, p, EoS)
!   print *, r2, "r2"
!   print *, m, "m"
!   print *, p, 'p'
  do i = 1, max_iterations, 1
    ! print *, p, "input pressure"
    CALL dp_dh(LengthOfEoS, x_array, Coefficients, p, k1_p)
    CALL dm_dh(LengthOfEoS, x_array, Coefficients, r2, p, m, k1_m)
    CALL dr2_dh(LengthOfEoS, x_array, Coefficients, r2, p, m, k1_r2)
    
    p_2 = p + del_h * k1_p / 2
    m_2 = m + del_h * k1_m / 2
    r_2 = r2 + del_h * k1_r2 / 2
    
    CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_2, k2_p)
    CALL dm_dh(LengthOfEoS, x_array, Coefficients, r_2, p_2, m_2, k2_m)
    CALL dr2_dh(LengthOfEoS, x_array, Coefficients, r_2, p_2, m_2, k2_r2)
    
        
    p_3 = p + del_h * k2_p / 2
    m_3 = m + del_h * k2_m / 2
    r_3 = r2 + del_h * k2_r2 / 2
    

    CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_3, k3_p)
    CALL dm_dh(LengthOfEoS, x_array, Coefficients, r_3, p_3, m_3, k3_m)
    CALL dr2_dh(LengthOfEoS, x_array, Coefficients, r_3, p_3, m_3, k3_r2)
    
        
    p_4 = p + del_h * k3_p
    m_4 = m + del_h * k3_m
    r_4 = r2 + del_h * k3_r2
    
    CALL dp_dh(LengthOfEoS, x_array, Coefficients, p_4, k4_p)
    CALL dm_dh(LengthOfEoS, x_array, Coefficients, r_4, p_4, m_4, k4_m)
    CALL dr2_dh(LengthOfEoS, x_array, Coefficients, r_4, p_4, m_4, k4_r2)
    
        
    d_p = (k1_p + 2 * k2_p + 2 * k3_p + k4_p) * del_h / 6
    d_m = (k1_m + 2 * k2_m + 2 * k3_m + k4_m) * del_h / 6
    d_r2 = (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) * del_h / 6


      if ((ieee_is_nan(d_p)) .OR. (ieee_is_nan(d_m))&
       .OR. (ieee_is_nan(d_r2))) then
       ! print *, d_p, "d_p"
       ! print *, d_m, "d_m"
       ! print *, d_r2, "d r2"
       
       Mass = m * km_to_M0
       Radius = sqrt(r2)
!        print *, "Pressure (MeV/fm^3): ", p / MeV_fm_to_km
!        print *, "Radius (km): ", Radius
!        print *, "Mass (M0): ", Mass
!        print *, i, "Iterations"
!        print *, "Ended with NAN"
!        print *, "do loop end-------------------------------------------------- \n"
       exit
       
      else if (abs(d_m) / (m+1.e-15) < 1.0e-15_dp) then
        Mass = m * km_to_M0
        Radius = sqrt(r2)
        ! print *, "Pressure (MeV/fm^3): ", p / MeV_fm_to_km
        ! print *, "Radius (km): ", Radius
        ! print *, "Mass (M0): ", Mass
        ! print *, i, "Iterations"
        ! print *, "Ended with dm"
        
        ! print *, "do loop end--------------------------------------------------\n"
        exit
      else if (p / MeV_fm_to_km < 1.e-13_dp) then
        Mass = m * km_to_M0
        Radius = sqrt(r2)
        ! print *, "Pressure (MeV/fm^3): ", p / MeV_fm_to_km
        ! print *, "Radius (km): ", Radius
        ! print *, "Mass (M0): ", Mass
        ! print *, i, "Iterations"
        ! print *, "Ended with dm"
        ! print *, "do loop end--------------------------------------------------\n"
        exit
      else
        p = p + d_p
        m = m + d_m
        r2 = r2 + d_r2
        ! print *, p, "next p"
        ! print *, m, "next m"
        ! print *, r2, "next r2"
      endif
     end do
  
end subroutine
end module tov_rk4
