module statistics
  implicit none
  ! Define double precision kind
  integer, parameter :: dp = selected_real_kind(15, 307)
  ! Calculate pi accurately
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

contains

  subroutine axial_mean_std(angles_deg, mean_deg, std_deg)
    real(dp), intent(in)  :: angles_deg(:)
    real(dp), intent(out) :: mean_deg, std_deg
    real(dp) :: x_sum, y_sum, r, std_rad_doubled
    real(dp) :: deg2rad
    integer  :: n

    n = size(angles_deg)
    
    ! Handle empty array
    if (n == 0) then
       mean_deg = 0.0_dp
       std_deg = 0.0_dp
       return
    end if

    deg2rad = pi / 180.0_dp

    ! 1 & 2. Double angles, convert to Cartesian, and sum.
    x_sum = sum(cos(2.0_dp * angles_deg * deg2rad))
    y_sum = sum(sin(2.0_dp * angles_deg * deg2rad))

    ! ==========================================
    ! MEAN ANGLE
    ! ==========================================
    if (abs(x_sum) < epsilon(1.0_dp) .and. abs(y_sum) < epsilon(1.0_dp)) then
       mean_deg = -999.0_dp ! Flag for undefined mean
    else
       mean_deg = (atan2(y_sum, x_sum) / 2.0_dp) / deg2rad
       !mean_deg = modulo(mean_deg, 180.0_dp)
    end if

    ! ==========================================
    ! STANDARD DEVIATION
    ! ==========================================
    ! Calculate Mean Resultant Vector Length (R)
    r = sqrt(x_sum**2 + y_sum**2) / real(n, dp)

    ! Handle edge cases for R
    if (r < epsilon(1.0_dp)) then
       ! Vectors perfectly cancel out; variance is theoretically maximum
       std_deg = -999.0_dp 
    else if (r >= 1.0_dp) then
       ! All angles are identical (handle potential float precision > 1.0)
       std_deg = 0.0_dp 
    else
       ! Mardia's formula for circular standard deviation (in radians)
       std_rad_doubled = sqrt(-2.0_dp * log(r))
       
       ! Halve the deviation to map back to 180-degree domain, convert to degrees
       std_deg = (std_rad_doubled / 2.0_dp) / deg2rad
    end if

  end subroutine axial_mean_std


  ! =======================================================
  ! Calculates mean and std dev for standard linear data
  ! =======================================================
  subroutine linear_mean_std(data_array, mean_val, std_val)
    real(dp), intent(in)  :: data_array(:)
    real(dp), intent(out) :: mean_val, std_val
    integer :: n

    n = size(data_array)

    ! Handle empty array
    if (n == 0) then
       mean_val = 0.0_dp
       std_val  = 0.0_dp
       return
    end if

    ! 1. Calculate Arithmetic Mean
    mean_val = sum(data_array) / real(n, dp)

    ! 2. Calculate Sample Standard Deviation
    if (n == 1) then
       ! Std dev is 0 for a single point (prevents divide by zero)
       std_val = 0.0_dp
    else
       ! Fortran array syntax applies the subtraction and squaring to every element
       std_val = sqrt( sum((data_array - mean_val)**2) / real(n - 1, dp) )
    end if

  end subroutine linear_mean_std

end module statistics