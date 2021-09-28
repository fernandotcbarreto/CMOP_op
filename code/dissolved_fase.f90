Module dissolved_fase_mod

 implicit none 

 double precision, dimension(:,:), allocatable:: xf3, yf3,zf3, lon_partf3, lat_partf3, massaf3

 integer, parameter :: numparcels_dis = 5

 integer:: parcel_dis_cont

 integer::  lat_in_f3(2), lon_in_f3(2)

 integer:: dissolved_fase


end module
