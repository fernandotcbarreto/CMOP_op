Module delft2oil
implicit none

 double precision, dimension(:, :, :), allocatable:: u_model1, u_model2, v_model1, v_model2, w_model1, w_model2
 double precision, dimension(:, :, :), allocatable:: t_model1, t_model2, s_model1, s_model2
 double precision, dimension(:, :, :), allocatable:: kz_model1, kz_model2
 double precision, dimension(:, :), allocatable:: lon_model, lat_model, lat_model_sum, lon_model_sum, coord_sum
double precision, dimension(:, :), allocatable:: lon_modelm, lat_modelm,lat_model_summ, lon_model_summ, coord_summ
double precision, dimension(:), allocatable:: lat_part_sum, lon_part_sum, coord_sum_f2
 double precision::  ui1,ui2,vi1,vi2,ti1,ti2,si1,si2, counttimeh , intertime(2), vel_interp(2), ui_vec(1), ts_vec(1),&
 di, dt_h, wi1, wi2, wi, dt_hf2, coastad
 double precision:: kz1, kz2, kz, counttimeh_r
 double precision, dimension(:, :), allocatable:: lat_part, lon_part, depth, area_grid, zf1, coastvalue,coastvalueb
 double precision, dimension(:), allocatable :: time_vec, vec_tdelf
 integer:: numlat, numlon, numlatm, numlonm, numtime, numz, tdelf1(1), tdelf2(1), lati, lonj
 integer:: lat_in(2), lon_in(2), pinf2(1), lat_in_f2(2), lon_in_f2(2), lat_inm(2), lon_inm(2)
 character(len=1024)::  path, bas, bas2, pre_end, end_prob, apist

contains 


subroutine init_delft
 implicit none 
 
 path='/mnt/c/Users/fernando.barreto/Documents/CMOP_OIL_SPILL-master/hycom_remo/'

end subroutine 

subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_VALUE_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
end subroutine 

end module
