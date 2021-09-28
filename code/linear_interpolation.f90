Module linear_interpolation
use delft2oil
use era5_2_oil

Implicit none 

  integer:: linear_interp
  integer, parameter:: slic=1
  double precision, parameter:: power = 1.0D+00
  integer, parameter:: siz_1d = (slic*2+1)*(slic*2+1)
  double precision, dimension(:,:):: slic_lat(slic*2+1, slic*2+1), slic_lon(slic*2+1, slic*2+1), slic_u1(slic*2+1, slic*2+1) &
     , slic_u2(slic*2+1, slic*2+1), slic_v1(slic*2+1, slic*2+1), slic_v2(slic*2+1, slic*2+1), slic_w1(slic*2+1, slic*2+1), &
slic_w2(slic*2+1, slic*2+1), &
    slic_kz1(slic*2+1, slic*2+1), slic_kz2(slic*2+1, slic*2+1),  slic_t1(slic*2+1, slic*2+1), slic_t2(slic*2+1, slic*2+1), &
     slic_s1(slic*2+1, slic*2+1), slic_s2(slic*2+1, slic*2+1)

  double precision, dimension(:):: lat_1d(siz_1d), lon_1d(siz_1d), u1_1d(siz_1d), u2_1d(siz_1d), v1_1d(siz_1d), v2_1d(siz_1d), &
 w1_1d(siz_1d), w2_1d(siz_1d), &
  kz1_1d(siz_1d), kz2_1d(siz_1d), t1_1d(siz_1d), t2_1d(siz_1d), s1_1d(siz_1d), s2_1d(siz_1d)

  double precision, dimension(:):: u1_1d_2(siz_1d), u2_1d_2(siz_1d), v1_1d_2(siz_1d), v2_1d_2(siz_1d), w1_1d_2(siz_1d), &
w2_1d_2(siz_1d), kz1_1d_2(siz_1d), kz2_1d_2(siz_1d), t1_1d_2(siz_1d), t2_1d_2(siz_1d), s1_1d_2(siz_1d), s2_1d_2(siz_1d)

  double precision, dimension(:,:):: slic_u1_2(slic*2+1, slic*2+1), slic_u2_2(slic*2+1, slic*2+1), slic_v1_2(slic*2+1, &
 slic*2+1), slic_v2_2(slic*2+1, slic*2+1), &
  slic_w1_2(slic*2+1, slic*2+1), slic_w2_2(slic*2+1, slic*2+1), slic_kz1_2(slic*2+1, slic*2+1), slic_kz2_2(slic*2+1, slic*2+1), &
  slic_t1_2(slic*2+1, slic*2+1), slic_t2_2(slic*2+1, slic*2+1), slic_s1_2(slic*2+1, slic*2+1), slic_s2_2(slic*2+1, slic*2+1)
  
  double precision::  out_int(1), out_int_2(1), inter_depth(2), par_dep(1)
  integer:: vert_index1, vert_index2

contains


subroutine init_interpolation(twol, lat_in_r, lon_in_r, varin1, varin2)

      integer, intent(in):: twol,lat_in_r, lon_in_r, varin1, varin2

!print*, '44444444444444444444444', varin1, varin2

        lat_1d = reshape(slic_lat, (/siz_1d/) )
        lon_1d = reshape(slic_lon, (/siz_1d/) )


        slic_u1 = u_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_u2 = u_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_v1 = v_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_v2 = v_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_t1 = t_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_t2 = t_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_s1 = s_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_s2 = s_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_w1 = w_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_w2 = w_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_kz1 = kz_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_kz2 = kz_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  

        u1_1d  = reshape(slic_u1, (/siz_1d/) )
        u2_1d  = reshape(slic_u2, (/siz_1d/) )
        v1_1d  = reshape(slic_v1, (/siz_1d/) )
        v2_1d  = reshape(slic_v2, (/siz_1d/) )
        t1_1d  = reshape(slic_t1, (/siz_1d/) )
        t2_1d  = reshape(slic_t2, (/siz_1d/) )
        s1_1d  = reshape(slic_s1, (/siz_1d/) )
        s2_1d  = reshape(slic_s2, (/siz_1d/) )
        w1_1d  = reshape(slic_w1, (/siz_1d/) )
        w2_1d  = reshape(slic_w2, (/siz_1d/) )
        kz1_1d  = reshape(slic_kz1, (/siz_1d/) )
        kz2_1d  = reshape(slic_kz2, (/siz_1d/) )



      if(twol.eq.2) then


! print*, 'bhna', varin1, varin2

        slic_u1_2 = u_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)
        slic_u2_2 = u_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)
        slic_v1_2 = v_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)
        slic_v2_2 = v_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)  
        slic_t1_2 = t_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_t2_2 = t_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_s1_2 = s_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)
        slic_s2_2 = s_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin1)  
        slic_w1_2 = w_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)
        slic_w2_2 = w_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)  
        slic_kz1_2 = kz_model1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)
        slic_kz2_2 = kz_model2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic, varin2)

        u1_1d_2  = reshape(slic_u1_2, (/siz_1d/) )
        u2_1d_2  = reshape(slic_u2_2, (/siz_1d/) )
        v1_1d_2  = reshape(slic_v1_2, (/siz_1d/) )
        v2_1d_2  = reshape(slic_v2_2, (/siz_1d/) )
        t1_1d_2  = reshape(slic_t1_2, (/siz_1d/) )
        t2_1d_2  = reshape(slic_t2_2, (/siz_1d/) )
        s1_1d_2  = reshape(slic_s1_2, (/siz_1d/) )
        s2_1d_2  = reshape(slic_s2_2, (/siz_1d/) )
        w1_1d_2  = reshape(slic_w1_2, (/siz_1d/) )
        w2_1d_2  = reshape(slic_w2_2, (/siz_1d/) )
        kz1_1d_2  = reshape(slic_kz1_2, (/siz_1d/) )
        kz2_1d_2  = reshape(slic_kz2_2, (/siz_1d/) )

      endif





end subroutine



subroutine init_interpolation_wind(lat_in_r, lon_in_r)

      integer, intent(in):: lat_in_r, lon_in_r

      lat_1d = reshape(slic_lat, (/siz_1d/) )
      lon_1d = reshape(slic_lon, (/siz_1d/) )

      slic_u1 = u10_era1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic)
      slic_u2 = u10_era2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic)
      slic_v1 = v10_era1 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic)
      slic_v2 = v10_era2 (lat_in_r-slic : lat_in_r+slic, lon_in_r - slic : lon_in_r + slic)  

      u1_1d  = reshape(slic_u1, (/siz_1d/) )
      u2_1d  = reshape(slic_u2, (/siz_1d/) )
      v1_1d  = reshape(slic_v1, (/siz_1d/) )
      v2_1d  = reshape(slic_v2, (/siz_1d/) )


end subroutine init_interpolation_wind



subroutine shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi )

!*****************************************************************************80
!
!! SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
!
!  Discussion:
!
!    This code should be vectorized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Shepard,
!    A two-dimensional interpolation function for irregularly spaced data,
!    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
!    ACM, pages 517-524, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data points.
!
!    Input, real ( kind = 8 ) ZD(ND), the data values.
!
!    Input, real ( kind = 8 ) P, the power.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) s
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yi(ni)
  integer ( kind = 4 ) z
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = 8 )

    else

      z = -1
      do j = 1, nd
        w(j) = sqrt ( ( xi(i) - xd(j) ) ** 2 + ( yi(i) - yd(j) ) ** 2 )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= -1 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
        s = sum ( w )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    zi(i) = dot_product ( w, zd )

  end do

  return
end subroutine




End module 


!b=shape(lat_model)
!print*, b(1), b(2)
!print*, shape(lat_model(10-2:10+2,100))
!stop
!a= b(1)*b(2)
!allocate(res_ar(a))
!print*, shape(res_ar)
!res_ar = reshape(lat_model, (/a/) )
!!print*, shape(reshape(lat_model, (/a/) ))
!lat_model = reshape(res_ar, (/int(b(1)), int(b(2))/) )
!print*, lat_model(10+2:30,100)
!stop
