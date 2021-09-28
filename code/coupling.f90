Module coupling

 implicit none 
 integer:: coupling_ind, start_index, inf_time, tc_in, frac_in, coup_ct
 character(len=1024)::  path_ini_coup, path_frac
 integer:: numb_lines_c, num_res_par
 real:: dummy_val
 double precision, dimension(:), allocatable :: time_res, lat_ref_res, lon_ref_res, vol_res, res_par_in
 double precision:: time_ini_spread, counttimeh_coup
 double precision, dimension(:,:), allocatable:: x_height, y_height, lat_height, lon_height, height_map, &
lat_model_sum_height, lon_model_sum_height, coord_sum_height
 integer:: num_sp, lat_in_height(2), lon_in_height(2)
 double precision:: dx_h, dy_h
 

 contains

  subroutine init_coupl
   implicit none 
 
    inf_time = 1
    num_res_par = 0
    path_ini_coup='/home/valdir/Documentos/oil_model/run_summer/abrolhos_summer/surface_data_abrolhosj_jan_2.txt'
    path_frac='/home/valdir/Documentos/oil_model/run_summer/abrolhos_summer/frac_mass_abrolhos_jan_2.txt'
    num_sp = 50
    dx_h = 1000
    dy_h = 1000
!    num_sp = 120 deep spill
!    dx_h = 100 deep spill
!    dy_h = 100 deep spill

  end subroutine 

End module
