Module dissolved_fase_mod

 implicit none 

 double precision, dimension (:,:), allocatable:: xf3, yf3,zf3, lon_partf3, lat_partf3, massaf3

 integer, parameter :: numparcels_dis = 500

 integer:: parcel_dis_cont, m1_f3

 integer::  lat_in_f3(2), lon_in_f3(2)

 double precision::  dif3, uif3, vif3, wif3,  kzf3

 integer:: dissolved_fase

 double precision:: uwdf3, vwdf3


end module
