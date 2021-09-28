Module era5_2_oil

implicit none


  double precision, dimension(:,:), allocatable:: u10_era1, u10_era2, v10_era1, v10_era2, lat_era
  double precision, dimension(:,:), allocatable:: lon_era, lat_era_sum, lon_era_sum, coord_sum_era
  integer:: lat_in_era(2), lon_in_era(2), lati_era, latj_era
  integer:: numlat_era, numlon_era, num_time_era, tera1(1), tera2(1)
  double precision, dimension(:), allocatable :: time_vec_era, vec_era
  character(len=1024)::  bas_era, bas2_era
  double precision:: uwind1, uwind2, vwind1, vwind2, intertime_era(2), vel_interp_era(2), counttimeh_era, ts_vec_w(1)

end module
