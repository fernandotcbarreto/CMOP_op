Module vertical_dispersion
use processes 
use ieee_arithmetic, only: ieee_is_finite
use random_variables

implicit none 

double precision, dimension(:,:), allocatable :: dropdiam, numdrops, voldrops

double precision, dimension(:), allocatable:: levelsd, deplevel, tsevol, tsevol2

double precision:: dropdiamax, dropdiamin,  numparcel, massparcel, ohnumb, wenumb

double precision::  turbed, cmax, cmin   ! turbulent energy dissipation rate in breaking waves (j/m3/s)  reed et al (1995)
  
double precision:: gd, bb(2), oil_entrainment_probability   !entrainment rate  kg/(m2*s)

double precision:: deltad, viscst, qd, zini, viscstem, seafrac

double precision, dimension (:,:), allocatable:: dropdiamf2, numdropsf2, voldropsf2, areaf2, volf2, massaf2, massadropsf2

double precision, dimension (:,:), allocatable:: xf2, yf2,zf2, lon_partf2, lat_partf2, diamf2, spmtf2

double precision, dimension(1):: rand_samples

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
 cmax = 2000
 cmin = 500
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

  Uth = 5                                                           !threshold 10 m wind speed for onset of breaking waves assumed 5 m/s   Delvigne and Sweeney(1988)
  Uth = 2.5                                                         !threshold 10 m wind speed for onset of breaking waves assumed 2.5 m/s  Lehr and Simecek-Beatty (2000) and Lehr et al. (2002)
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


  
 ! print*,"li",  ohnumb, wenumb, visc, ro_oil, interfacial_tension, d0
 ! stop
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



subroutine size_distr_li_2007 (visc, interfacial_tension, ro_a, ro_oil, windspms, qd, zini, seafrac&
, kzm, wvp)
 implicit none

  double precision:: cstar, viscst, ro_a, Hsig , ustar2  , windspms, Uth, tw, dv, dd, hb
  double precision,intent (out):: seafrac  ! fraction of sea surface hit by breaking waves
  double precision:: S, zi,d0, delta_ro
  double precision:: deltad, dropdiam, wvp, kzm
  double precision:: visc, interfacial_tension, ro_oil, d_ref
  double precision, intent (out) :: qd, zini
  double precision, dimension(:), allocatable :: droplet_spectrum_diameter, droplet_spectrum_pdf
  double precision :: d_o, r, p, q, dV_50, sd, we, oh, min_val, max_val, step
  double precision, dimension(:), allocatable :: spectrum
  integer :: iyg, nhdh
!  integer, parameter :: num_elements = 1000000
  integer, parameter :: num_elements =  10000
  double precision :: dN_50
  logical :: is_finite_flag

  allocate(droplet_spectrum_diameter(num_elements))
!  droplet_spectrum_diameter = [(iyg * (3.0e-3 - 1.0e-6) / (num_elements - 1.0d0), iyg = 1, num_elements)]
 ! droplet_spectrum_diameter = [(iyg * (3.0e-4 - 1.0e-6) / (num_elements - 1.0d0), iyg = 1, num_elements)] !mais correto
 min_val = 1.0e-6
! max_val = 3.0e-4 !mais correto
! max_val = 1.0e-4 !mais correto
! max_val = 3.0e-5 !mais oleo na coluna de agua
 max_val = 3.0e-3 ! da velocidades negativas de wvp
 
  if (windspms.eq.0) then
    windspms=0.01
 endif 
 
!num_elements = 1000000
step = (max_val - min_val) / (num_elements - 1)

do iyg = 1, num_elements
    droplet_spectrum_diameter(iyg) = min_val + (iyg - 1) * step
end do
  
  delta_ro= abs(ro_a - ro_oil)

  ustar2 = 0.71 * (windspms**1.23)
 
  hsig = (0.243 * ustar2) / gravity

  hsig = hsig      !Mar não desenvolvido, input dinâmico do vento 
  
 !print*, hsig
 ! hsig=(0.0246*windspms**2)/gravity   !Mar desenvolvido  opendrift
 ! print*, (0.0246*windspms**2)/gravity, (0.243 * ustar2) / gravity
  ! Compute parameters
  d_o = 4.0d0 * sqrt(interfacial_tension / (delta_ro * gravity))
 ! d0 = 4 * ( (interfacial_tension / (delta_ro*gravity) )**(0.5) )
  
  we = (ro_a * gravity * hsig * d_o) / interfacial_tension
  !we = (ro_a * gravity * hsig) / interfacial_tension
 ! wenumb = ro_a * gravity * hsig * d0  / interfacial_tension

  oh =  visc / ( (ro_oil * interfacial_tension * d_o)**(0.5) )
  !d_ref = 1.0d-4   ! por exemplo, 100 μm
 
  !oh = visc / sqrt(ro_oil * interfacial_tension * d_ref)
  
  !print*,"particle",  oh, we, visc, ro_oil, interfacial_tension, d_o
  !stop
  !oh = self.elements.viscosity * self.elements.density * (
   !             self.elements.density * interfacial_tension *
    !            d_o)**-0.5  
 ! ohnumb =  visc / ( (ro_oil * interfacial_tension * d0)**(0.5) )
  
  r = 1.791d0
  p = 0.460d0
  q = -0.518d0
  
  ! Median droplet diameter in volume distribution
 ! dV_50 = d_o * r * (1.0d0 + 10.0d0 * oh)**p * we**q
  dV_50 = d_o * r * oh**p * we**q

  ! Log standard deviation in log10 units
  sd = 0.4d0
  
  ! Log standard deviation in natural log units
  Sd = log(10.0d0) * sd
  
  ! Convert number distribution to volume distribution
  !dV_50 = mean(dV_50)
 ! print*, dV_50, d_o , r ,  1.0d0 + 10.0d0 * oh, oh,p, we, q
 ! stop
!dV_50 = 0.0d0
!do iyg = 1, num_elements
!  dV_50 = dV_50 + droplet_spectrum_diameter(iyg)
!end do
!dV_50 = dV_50 / num_elements
!  print*, dV_50

!  print*, 'wind =', windspms, 'we =', we, 'oh =', oh, 'dV_50 =', dV_50, "inter", interfacial_tension, "visc", visc , "ro",ro_oil

 !stop
 
  dN_50 = exp(log(dV_50) - 3.0d0 * Sd**2)  
  
  allocate(spectrum(num_elements))
  do iyg = 1, num_elements
     spectrum(iyg) = exp(-((log(droplet_spectrum_diameter(iyg)) - log(dV_50))**2) / (2.0d0 * Sd**2)) / &
                    (droplet_spectrum_diameter(iyg) * Sd * sqrt(2.0d0 * pi))
  end do  
  
!            spectrum = (np.exp(
!                -(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 /
!                (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd *
!                                 np.sqrt(2 * np.pi))  

  ! Normalize the spectrum
  droplet_spectrum_pdf = spectrum / sum(spectrum)  
  
  ! Check for validity
!  if (.not. all(isfinite(droplet_spectrum_pdf)) .or. abs(sum(droplet_spectrum_pdf) - 1.0d0) > 1.0e-6) then
 
  is_finite_flag = isfinite_array(droplet_spectrum_pdf)
 
  if (.not. is_finite_flag .or. abs(sum(droplet_spectrum_pdf) - 1.0d0) > 1.0e-6) then
  
     print *, 'Could not update droplet diameters.'
     ! Return an array or value based on your requirements
	 rand_samples(1)=888898
  else
     ! Generate random samples based on the droplet spectrum PDF
     call random_sample(droplet_spectrum_diameter, droplet_spectrum_pdf, rand_samples)
 !    print *, 'Random samples:', rand_samples(1)
	! stop
     ! Return the generated random samples
  end if  
 
   !  print *, 'Random samples:', rand_samples(1), dV_50,d_o
   ! print*, interfacial_tension

	! stop
end subroutine 


  subroutine random_sample(choices, probabilities, samples)
    double precision, dimension(:), intent(in) :: choices, probabilities
    double precision, dimension(:), intent(out) :: samples
    double precision :: ujkh
    integer :: itr, jgh
    
    do itr = 1, size(samples)
 !      ujkh = random_number()
	   CALL random_number(ujkh)
       jgh = 1
       do while (ujkh > sum(probabilities(1:jgh)))
          jgh = jgh + 1
       end do
       samples(itr) = choices(jgh)
    end do
    
  end subroutine random_sample
  
  logical function isfinite_array(arr)
    real(8), dimension(:), intent(in) :: arr
    integer :: i

    do i = 1, size(arr)
      if (.not. ieee_is_finite(arr(i))) then
        isfinite_array = .false.
        return
      end if
    end do

    isfinite_array = .true.
  end function isfinite_array



End Module
