Module random_variables
 implicit none

 double precision vapp, r, temp, evap1, evap2, CDIF_HOR, CDIF_VER, U_ALEA, V_ALEA, W_ALEA, WP, SALin
 double precision rn4
 double precision rn5
 integer count_prob_2

 contains

 subroutine rand_var(vapp, r, temp)
  double precision:: vapp,r , temp
  vapp=0.01
  r=8.208*(10.**(-5.)) 
  temp=44 !celsius
!  CDIF_HOR       = 20.0D0				
!  CDIF_HOR       = 6.0D0	 baia guanabara			
  CDIF_VER       = 1.D-3
  SALin          =  35  
  
 end subroutine rand_var


 subroutine celstokel(temp)
   double precision:: temp
   temp=temp+273.15
 end subroutine celstokel

end module random_variables
