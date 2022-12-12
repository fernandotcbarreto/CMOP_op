Module vertical_dispersion
use processes 
implicit none 

double precision, dimension(:,:), allocatable :: dropdiam, numdrops, voldrops

double precision, dimension(:), allocatable:: levelsd, deplevel, tsevol, tsevol2

double precision:: dropdiamax, dropdiamin, rn4, numparcel, massparcel, ohnumb, wenumb

double precision::  turbed, cmax, cmin   ! turbulent energy dissipation rate in breaking waves (j/m3/s)  reed et al (1995)
  
double precision:: gd, bb(2), oil_entrainment_probability   !entrainment rate  kg/(m2*s)

double precision:: deltad, viscst, qd, zini, viscstem, seafrac

double precision, dimension (:,:), allocatable:: dropdiamf2, numdropsf2, voldropsf2, areaf2, volf2, massaf2, massadropsf2

double precision, dimension (:,:), allocatable:: xf2, yf2,zf2, lon_partf2, lat_partf2, diamf2, spmtf2


integer :: numparcels

double precision, parameter:: limen = 1.  ! time interval entrainment in hours

double precision, parameter:: limassem = 0.0  ! limit mass activate entrainment in kg

integer:: parcelcont, m2, zl, indexz(1), contind, numtot, m1, hh

integer:: total_particles, contprob

double precision, dimension(:, :, :), allocatable:: MASSCOMPf2, masscompdropf2

integer, dimension(:, :), allocatable:: partcontmap,partcontmap_2

double precision, dimension(:, :), allocatable:: probmap, probmap_2, time_sum, prob_time_sum

double precision, dimension(:, :), allocatable:: mass_dist, conc

double precision::  dif2, uif2, vif2,tif2,sif2, wi1f2,wi2f2, wif2, kz1f2, kz2f2, kzf2, &
                    temp_outf2, SAL_Af2

double precision:: dt_random, xrandom, yrandom, zrandom

double precision:: step_random

contains


subroutine entrainment_param
implicit none 

 turbed = 10000  !berry et al 
 cmax =3400    !default
 cmin = 500
! cmax = 600
! cmin = 500
! print*, gravity

end subroutine




subroutine vert_disp (viscst, ro_A, windspms, deltad, dropdiam, wvp, qd, zini, kzm, seafrac)
 implicit none

  double precision:: cstar, viscst, ro_a, Hsig , ustar2  , windspms, Uth, tw, dv, dd, hb
  double precision, intent (out):: seafrac  ! fraction of sea surface hit by breaking waves
  double precision:: S, zi
  double precision:: deltad, dropdiam, wvp, kzm

  double precision, intent (out) :: qd, zini

  Uth = 2.5                                                           !threshold 10 m wind speed for onset of breaking waves assumed 5 m/s
  tw = 8.13 * (windspms / gravity)                                   !significant wave period and assuming a fully-developed sea
  s = 1                                                           !fraction of sea covered by oil
  dv = 0.0015 * windspms                     !vertical diffusion coefficient (m**2/s)

!  print*, 'compare', dv, kzm

  if (viscst.lt.132) then
      cstar = exp(-0.1023 * log(viscst) + 7.575)
  else
      cstar = exp(-1.8927 * log(viscst) + 16.313)
  endif

  ustar2 = 0.71 * (windspms**1.23)
 
  hsig = (0.243 * ustar2) / gravity
  
  print*, 'hsig', hsig

  hb= 1.5 * hsig
 
  dd = 0.0034 * gravity * ro_a * (hb**2)

  if (windspms .lt. Uth) then
     seafrac = 0.000003 * ( (windspms)**3.5 / tw)   
  else
     seafrac  = 0.032 * ( (windspms - Uth)/ tw)

     seafrac  = (1.92*0.032) * ( (windspms - Uth)/ windspms)   !li et all 2017
  endif
  

  qd = cstar * (dd**0.57) * s* seafrac * (dropdiam** 0.7) * deltad

  zi = 1.5*hb

  zini = max (kzm/wvp, zi)

  zini = zi

!print*, zi, hb, kzm/wvp, wvp, kzm
!stop
!print*, dropdiam, seafrac, dd, s, deltad, windspms, qd, viscst, cstar, hsig
!  print*, 'yyuyuddddddddddd', tw,dropdiam, seafrac, qd
!stop
!   print*, 'gggg'

end subroutine 



subroutine vert_disp_li_2007 (visc, interfacial_tension, ro_a, ro_oil, windspms, qd, zini, seafrac&
, kzm, wvp)
 implicit none

  double precision:: cstar, viscst, ro_a, Hsig , ustar2  , windspms, Uth, tw, dv, dd, hb
  double precision,intent (out):: seafrac  ! fraction of sea surface hit by breaking waves
  double precision:: S, zi,d0, delta_ro
  double precision:: deltad, dropdiam, wvp, kzm
  double precision:: visc, interfacial_tension, ro_oil
  double precision, intent (out) :: qd, zini

  
  delta_ro= abs(ro_a - ro_oil)

  d0 = 4 * ( (interfacial_tension / (delta_ro*gravity) )**(0.5) )

  Uth = 5                                                           !threshold 10 m wind speed for onset of breaking waves assumed 5 m/s
  tw = 8.13 * (windspms / gravity)                                   !significant wave period and assuming a fully-developed sea
  !tw=2*pi/(0.877*9.81/(1.17*windspms))
  s = 1                                                           !fraction of sea covered by oil
  dv = 0.0015 * windspms                     !vertical diffusion coefficient (m**2/s)

!  print*, 'compare', dv, kzm

  ustar2 = 0.71 * (windspms**1.23)
 
  hsig = (0.243 * ustar2) / gravity

  hsig = hsig 

  hb= 1.5 * hsig
 
  dd = 0.0034 * gravity * ro_a * (hb**2)

  if (windspms .lt. Uth) then
     seafrac = 0.000003 * ( (windspms)**3.5 / tw)   
  else
     seafrac  = 0.032 * ( (windspms - Uth)/ tw)

     !seafrac  = (1.92*0.032) * ( (windspms - Uth)/ windspms)   !li et all 2017
  endif
  

  wenumb = ro_a * gravity * hsig * d0  / interfacial_tension

  ohnumb =  visc / ( (ro_oil * interfacial_tension * d0)**(0.5) )

 !print*, '44',  tw, pi

 ! qd =  4.604e-10 * ( wenumb**(1.805) ) * ( ohnumb**(-1.023) ) * seafrac * ro_oil * 0.0002
  qd =  4.604e-10 * ( wenumb**(1.805) ) * ( ohnumb**(-1.023) ) * seafrac  
  
  zi = 1.5*hb

  zini = max (kzm/wvp, zi)

 ! zini = zi

! print*, zi, kzm/wvp
 ! print*, 'yyuyuuuuuuuuuuuuuu',wenumb, ohnumb
 ! print*, 'naba', qd,  wenumb**(1.805) ,  ohnumb**(-1.023), seafrac,  ( wenumb**(1.805) ) * ( ohnumb**(-1.023) ), interfacial_tension
 ! print*, 'iii', qd, ro_a ,gravity ,hsig , d0 ,  interfacial_tension
 ! print*, 'qd', qd

end subroutine 




End Module
