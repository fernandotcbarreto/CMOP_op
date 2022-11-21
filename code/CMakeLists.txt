CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(CMOP Fortran)

SET(VERSION 0.0.1)

SET(CMOP_exe cmop)
SET(CMOP_src coupling.f90
             datetime_module.f90
             delft2oil.f90
             dissolve_part.f90
             dissolved_fase_mod.f90
             era5_2_oil.f90
             evaporation.f90
             lagrange_continuo_new_time_loop.f90
             linear_interpolation.f90
             oil_fractions.f90
             processes.f90
             r8lib.f90
             random_variables.f90
             vertical_dispersion.f90
             wind_variables.f90
             white_space.f90
)

ADD_EXECUTABLE(${CMOP_exe} ${CMOP_src})

INSTALL(TARGETS ${CMOP_exe} DESTINATION bin)