Module processes 
  use oil_fractions
  implicit none 

  double precision:: stab_index, ucomp2,vcomp2, vel2,dir2
  double precision:: XA !fraction and mass of asphaltenes
  double precision:: XW  !fraction and mass of wax
  double precision, dimension(:), allocatable:: XAM, XWM
  double precision:: kao, kaw
  double precision:: alfa, alfa0, alfa67
  double precision, dimension(:,:), allocatable:: watf   ! w
  double precision:: maxwf  !max water fraction
  double precision, dimension(:,:), allocatable:: watcont  ! water content
  double precision:: watup  !water uptake
  double precision:: watout !water release
  double precision:: kem
  double precision:: k1N ! a rate constant for water incorporation

  double precision:: seawat_rho, viscinA20, RO_A , VIS_DIN_A , CP_A , TS_A, Z , PROF_REF , PROFPOS, T_A , SAL_A

  DOUBLE PRECISION    ::  PM_H2O           	! peso molecular agua 
  DOUBLE PRECISION    :: GRAVITY , PI , R_TERRA  ! constantes, gravidade, pi e raio da terra

  double precision:: LON_REF, LAT_REF, LONOUT, LATOUT, PROFOUT, XPOS, YPOS, ZPOS
  integer::  OPT
  double precision::  circle_radius, time_spread


  contains

  subroutine emulsify_parameter

     kao=3.3  !at 293 K
  
     kaw=200 ! at 293 k

     xa=0.4

     xw=0.4

     watf(:,:)=0

     maxwf=0.8   !Xie et al. 2007 Modeling emulsification after an oil spill in the sea 

     watcont(:,:)=0

     XAM(:)=massa(1,1)*xa

     XWM(:)=massa(1,1)*xw

     massae=massa

     watout=0
  
     watup=0

     !print*, kem
	 !stop
      kem =  0.0000005   !VALUE USED IN EMULFICATION VALIDATION
     ! kem =  0.000001   !Value used to validate VISCOSITY

     !kem =  0.000002  !

     k1N = 0.8

     PM_H2O  = 18.D0     	! peso molecular agua g/mol

     seawat_rho = 1025   !!kg/m3

     viscinA20=0.00000105   !kinematic visc água m2/s
     Z = 0
     PROF_REF = 0 
     PROFPOS=0
     T_A  = 273   ! KELVIN
     SAL_A = 33



     GRAVITY  = 9.8D0              		! gravidade
     PI       = 4.D0*DATAN(1.D0)      		! pi
     R_TERRA  = 6371000.D0			! raio da terra

  end subroutine emulsify_parameter


SUBROUTINE vel_conv
 ! IMPLICIT NONE

!  DOUBLE PRECISION  , INTENT(IN)   :: vel2, dir2
!  DOUBLE PRECISION  , INTENT(OUT)   :: ucomp2,vcomp2
  PI       = 4.D0*DATAN(1.D0)   
  if (dir2 .le. 90) then
    ucomp2 = vel2*sin(dir2*PI/180.D0)
    vcomp2 = vel2*cos(dir2*PI/180.D0)
	print*, PI
  elseif (dir2 .gt. 90 .and. dir2 .le. 180) then
    dir2=dir2-90
    ucomp2 = vel2*cos(dir2*PI/180.D0)
    vcomp2 = -vel2*sin(dir2*PI/180.D0)
  elseif (dir2 .gt. 180 .and. dir2 .le. 270) then
    dir2=dir2-180
    ucomp2 = -vel2*sin(dir2*PI/180.D0)
    vcomp2 = -vel2*cos(dir2*PI/180.D0)
  elseif (dir2 .gt. 270 .and. dir2 .le. 360) then
    dir2=dir2-270
    ucomp2 = -vel2*cos(dir2*PI/180.D0)
    vcomp2 = vel2*sin(dir2*PI/180.D0)
  endif
END SUBROUTINE

  subroutine emulsify_parameter_coupling(numtot, num_res_par)

     integer:: numtot, num_res_par

     kao=3.3  !at 293 K
  
     kaw=200 ! at 293 k

     xa=0.4

     xw=0.4

     watf(:,:)=0

     maxwf=0.8   !Xie et al. 2007 Modeling emulsification after an oil spill in the sea 

     watcont(:,:)=0

     XAM(num_res_par+1:numtot)=massa(num_res_par+1:numtot,1)*xa

     XWM(num_res_par+1:numtot)=massa(num_res_par+1:numtot,1)*xw

     massae(num_res_par+1:numtot, :)=massa(num_res_par+1:numtot, :)

     watout=0
  
     watup=0

     ! kem =  0.0000005   !VALUE USED IN EMULFICATION VALIDATION

     ! kem =  0.000001   !Value used to validate VISCOSITY

     kem =  0.000002  !

     k1N = 0.8

     PM_H2O  = 18.D0     	! peso molecular agua g/mol

     seawat_rho = 1025   !!kg/m3

     viscinA20=0.00000105   !kinematic visc água m2/s
     Z = 0
     PROF_REF = 0 
     PROFPOS=0
     T_A  = 273   ! KELVIN
     SAL_A = 33



     GRAVITY  = 9.8D0              		! gravidade
     PI       = 4.D0*DATAN(1.D0)      		! pi
     R_TERRA  = 6371000.D0			! raio da terra

  end subroutine emulsify_parameter_coupling





   
  subroutine emulsify(dt, windspms, xam_r, xwm_r)
   
    implicit none
    double precision:: dt
    double precision:: windspms
    double precision:: xam_r, xwm_r

!    print*, windspms
    

    xa=xam_r/massae(jj,ii-1)
    xw=xwm_r/massae(jj,ii-1)

!    watf=watcont(jj,ii-1)/massae(jj,ii-1)

 !   print*, 'massae',   massae(jj,ii-1)

    stab_index = xa*exp(kao*((1-xa-xw)**2) + kaw*(xw**2))*exp(-0.04*(temp_out-293))


    alfa0= log(maxwf/0.1)/3600
    alfa67=log(maxwf/0.1)/(3600*24)


    if (stab_index .ge. 1.22 ) then
       alfa=0


    else if (stab_index .le. 1.22 .and. stab_index .ge. 0.67) then
       alfa=alfa67*((1.22-stab_index)/(1.22-0.67))

!print*, 'MESO'

    else 
       alfa=alfa0-(((alfa0-alfa67)*stab_index)/0.67)

!PRINT*, 'UNST'

    endif

    
!   PRINT*, STAB_INDEX


    watup=(kem*((windspms+1)**2)*((maxwf-watf(jj,ii-1))/maxwf))*dt  !fraction by time

    watf(jj,ii)=watup+watf(jj,ii-1)
   
    watout=-alfa*watf(jj,ii-1)*dt   !is fraction by time

    watf(jj,ii)=watout+watf(jj,ii)


!    watout=watout*massae(jj, ii-1) ERRADO

!    watup=watup*massae(jj, ii-1)   ERRADO
    

!    print*, 'watf', watf

 !   print*, 'watup', watup

 !   print*, 'watout', watout
 
!    print*, watf, xa, xw
   

!    print*, 'stab', stab_index

  end subroutine



  subroutine rho_emulsion(spmt)
   double precision:: spmt

   rho_e(jj,ii)=watf(jj,ii)*seawat_rho + (1-watf(jj,ii))*spmt
   


 end subroutine rho_emulsion

  

  subroutine visc_emulsion(VIS_DIN_OIL_OUT)
     implicit none 
     double precision:: VIS_DIN_OIL_OUT

!     print*, VIS_DIN_OIL_OUT

     visc_e(jj,ii)=VIS_DIN_OIL_OUT*exp((2.5*watf(jj,ii))/(1-k1N*watf(jj,ii)))    !viscosity is in cP

!    print*, 'CALCCC', VIS_DIN_OIL_OUT

     
!     print*, visc_e
   
  end subroutine visc_emulsion






 SUBROUTINE PROP_AMBIENTE (  Z , PROF_REF , T_A , SAL_A , &
			     RO_A , VIS_DIN_A , CP_A , TS_A) 

  IMPLICIT NONE

! determinar as propriedades densidade, viscosidade, calor especifico e tensao superficial
! do ambiente e da agua dentro de vc

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN)  :: T_A , SAL_A		! temperatura e salinidade do ambiente 
   DOUBLE PRECISION , INTENT(IN)  :: Z , PROF_REF		! posicao de vc

! -------------------------------------------------------  parametros out
   DOUBLE PRECISION , INTENT(OUT) ::  RO_A  
   DOUBLE PRECISION , INTENT(OUT) ::  VIS_DIN_A , CP_A					! variaveis ambientes
   DOUBLE PRECISION , INTENT(OUT) ::  TS_A 						! tensao superficial

! ------------------------------------------------------- parametros locais
   DOUBLE PRECISION :: A1 , A2 , A3 , A4 , A
   DOUBLE PRECISION :: K1 , K2 , K3 , K0 , K
   DOUBLE PRECISION :: B1 , B2 , B
   DOUBLE PRECISION :: RO_0					! densidade ambiente de referencia 
   DOUBLE PRECISION :: VIS_0  					! viscosidade ambiente de referencia
   DOUBLE PRECISION :: T_AMB, P_AMB , S_AMB , S_AMB1
   DOUBLE PRECISION :: TEMP , SAL
   DOUBLE PRECISION :: RO , CP ,  VIS_DIN , TS , TS_0
   DOUBLE PRECISION :: P_A, GRAVITY

!-----------------------------------------------------------------------------------

! --- a densidade da agua pode ser calculada usando ou nao a pressao. Logo, temos duas opções:

! -------------------------------------------------------
! --- usar uma equacao para calcular a densidade independente da pressao -> segundo Sharqawy et al (2010) 
! -> salinidade em g/kg e temperatura em °C

   GRAVITY=10.
   T_AMB = T_A - 0.00025D0*(T_A - 273.15D0) - 273.15D0	! temperatura na International Temperature Scale of 1990


   S_AMB = 1.00472D0*SAL_A				! salinidade na escala reference-composition salinity
   
   B  = (2.D0*S_AMB - 150.D0)/150.D0
   A1 = 4.032D0*0.5D0  + 0.115D0*B + 3.26D-4*(2.D0*B**2.D0 - 1.D0)
   A2 = -0.108D0*0.5D0 + 1.571D-3*B - 4.23D-4*(2.D0*B**2.D0 - 1.D0)
   A3 = -0.012D0*0.5D0 + 1.74D-3*B - 9.D-6*(2.D0*B**2.D0 - 1.D0)
   A4 = 6.92D-4*0.5D0  - 8.7D-5*B - 5.3D-5*(2.D0*B**2.D0 - 1.D0)

   K1  = (2.D0*T_AMB - 200.D0)/160.D0
   K2  = 2.D0*K1**2.D0 - 1.D0
   K3  = 4.D0*K1**3.D0 - 3.D0*K1

   RO_A = 1.D3*(A1*0.5D0 + A2*K1 + A3*K2 + A4*K3) 

! pressao ambiente em Pascal -> quando usar a pressao no calculo da densidade

   P_A = 101325.D0 + RO_A*GRAVITY*(PROF_REF - Z)

! -------------------------------------------------------


! -------------------------------------------------------
!- conversao de P_A e T_A

   P_AMB = P_A*1.D-5		  ! pressao ambiente em bar

   T_AMB = T_A - 273.15D0         ! temperatura ambiente em oC

   S_AMB = SAL_A

! -------------------------------------------------------


! -------------------------------------------------------

!### densidade calculada segundo UNESCO(1981) -> salinidade em PSU , temperatura em oC e pressao em bar

!- calcular a densidade da agua em atmosfera padrao
   A1 = 999.841594D0 + 6.793952D-2*T_AMB - 9.09529D-3*T_AMB**2.D0 + 1.001685D-4*T_AMB**3.D0 - &
 1.120083D-6*T_AMB**4.D0 + 6.536332D-9*T_AMB**5.D0
   A2 = 8.24493D-1 - 4.0899D-3*T_AMB + 7.6438D-5*T_AMB**2.D0 - 8.2467D-7*T_AMB**3.D0 + 5.3875D-9*T_AMB**4.D0 
   A3 = -5.72466D-3 + 1.0227D-4*T_AMB - 1.6546D-6*T_AMB**2.D0 
   A4 = 4.8314D-4

   RO_0 = A1 + A2*S_AMB + A3*S_AMB**(3.D0/2.D0) + A4*S_AMB**2.D0

!- calcular a densidade da agua em P_AMB
   K1 = 19652.21D0 + 148.4206D0*T_AMB - 2.327105D0*T_AMB**2.D0 + 1.360477D-2*T_AMB**3.D0 - 5.155288D-5*T_AMB**4.D0
   K2 = 54.6746D0 - 0.603459D0*T_AMB + 1.09987D-2*T_AMB**2.D0 - 6.167D-5*T_AMB**3.D0
   K3 = 7.944D-2 + 1.6483D-2*T_AMB - 5.3009D-4*T_AMB**2.D0

   K0 = K1 + K2*S_AMB + K3*S_AMB**(3.D0/2.D0)

   A1 = 3.239908D0 + 1.43713D-3*T_AMB + 1.16092D-4*T_AMB**2.D0 - 5.77905D-7*T_AMB**3.D0
   A2 = 2.2838D-3  - 1.0981D-5*T_AMB  - 1.6078D-6*T_AMB**2.D0
   A3 = 1.91075D-4

   A = A1 + A2*S_AMB + A3*S_AMB**(3.D0/2.D0)

   B1 = 8.50935D-5 - 6.12293D-6*T_AMB + 5.2787D-8*T_AMB**2.D0
   B2 = -9.9348D-7 + 2.0816D-8*T_AMB + 9.1697D-10*T_AMB**2.D0

   B = B1 + B2*S_AMB

   K = K0 + A*P_AMB + B*P_AMB**2.D0

!   RO   = RO_0 / ( 1. - P_AMB/K )				! densidade da agua em P_AMB , T_AMB e S_AMB em kg/m3
    RO   = RO_0							! densidade da agua em T_AMB e S_AMB em kg/m3

!-----------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------
!### viscosidade dinamica calculada segundo Sharqawy et al (2010) -> salinidade em kg/kg e temperatura em o0C

!- calcular a viscosidade dinamica da agua - nao varia com a pressao
   S_AMB1 = S_AMB*1.D-3							! em kg/kg - salinidade de vc


   VIS_0 = 4.2844D-5 + ( 0.157D0*(T_AMB + 64.993D0)**2.D0 - 91.296D0 )**(-1.D0)

   A = 1.541D0 + 1.998D-2*T_AMB - 9.52D-5*T_AMB**2.D0

   B = 7.974D0 - 7.561D-2*T_AMB + 4.724D-4*T_AMB**2.D0

   VIS_DIN = VIS_0 * ( 1.D0 + A*S_AMB1 + B*S_AMB1**2.D0)		! viscosidade dinamica do ambiente em kg/m.s

!-----------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------
!### tensao interfacial calculada segundo Sharqawy et al (2010) -> salinidade em g/kg e temperatura em Celsius

!- calcular a tensao interfacial da agua
   TS_0 = 0.2358D0 * ( 1.D0 - (T_AMB + 273.15D0)/647.096D0)**1.256D0 * (1.D0 - 0.625D0*(1.D0 - (T_AMB + 273.15D0)/647.096D0))    ! em N/m


   TS = TS_0 * (1.D0 + (0.000226D0*T_AMB + 0.00946D0) * DLOG(1.D0 + 0.0331D0 * S_AMB))         ! tensao superficial da agua em N/m

!-----------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------
!### calor especifico calculado segundo Sharqawy et al (2010) -> salinidade em g/kg e temperatura em K

!- calcular o calor especifico da agua do mar 
   S_AMB = 1.00472D0*S_AMB 					! em salinidade absoluta em g/kg
   
   T_AMB = T_AMB + 273.15D0					! temperatura em kelvin

   T_AMB = (T_AMB - 0.00025D0 * (T_AMB - 273.15D0) )		! temperatura na International Temperature Scale of 1990


   A1 = 5.328D0 - 9.76D-2*S_AMB + 4.04D-4*S_AMB**2.D0
   A2 = -6.913D-3 + 7.351D-4*S_AMB - 3.15D-6*S_AMB**2.D0
   A3 = 9.6D-6 - 1.927D-6*S_AMB + 8.23D-9*S_AMB**2.D0
   A4 = 2.5D-9 + 1.666D-9*S_AMB - 7.125D-12*S_AMB**2.D0

   CP  = A1 + A2*T_AMB + A3*T_AMB**2.D0 + A4*T_AMB**3.D0		! em kJ/kg.K

   CP  = CP*1.D+3							! capacidade calorifica do ambiente em J/kg.K

!-----------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------

         RO_A = RO
         CP_A = CP
         VIS_DIN_A = VIS_DIN
         TS_A = TS


!-----------------------------------------------------------------------------------


 !----------------------
  END SUBROUTINE PROP_AMBIENTE



  SUBROUTINE  VEL_ASCENSAO_BG ( DIAM , DIAM_HYD , RO_A , RO_P , VIS_A  , TS  , WP  )

  IMPLICIT NONE

!-------- velocidade de ascensao de goticulas (m/s)

          DOUBLE PRECISION, INTENT(IN)  :: DIAM                 ! diametro da bolha do gas
          DOUBLE PRECISION, INTENT(IN)  :: DIAM_HYD             ! diametro da bolha do gas
          DOUBLE PRECISION, INTENT(IN)  :: RO_A , RO_P        	! densidade ambiente, gas (ambiente pode ser gas) 
          DOUBLE PRECISION, INTENT(IN)  :: VIS_A               	! coeficiente de viscosidade dinamica do ambiente
          DOUBLE PRECISION, INTENT(IN)  :: TS			! tensao interfacial do ambiente

          DOUBLE PRECISION, INTENT(OUT) :: WP                   ! velocidade de ascensao

          DOUBLE PRECISION :: W1 , W2 
          DOUBLE PRECISION :: RE , A , ND
          DOUBLE PRECISION :: EO , MO , H , J , CD
          DOUBLE PRECISION :: DC , H1, H2, J1, J2, UT1, UT2, D1, D2, X1, X2, Y1, Y2, A1, A2, B1, B2, EO1, EO2 
          DOUBLE PRECISION :: FAC , WE , WP0 , DEF , WP1

!---------------------------------------------------------------------------------------------------
!-------- calculo de reynolds -> Clift et al (1978)

	  ND = 4.D0*RO_A*GRAVITY*DIAM**(3.D0)*(RO_A - RO_P)/(3.D0*VIS_A**2.D0)
	  A = DLOG10(ND)

	  IF ( ND .LE. 73.D0 ) THEN
     		RE = ND/24.D0 - 1.7569D-4*ND**2.D0 + 6.9252D-7*ND**3.D0 - 2.3027D-10*ND**4.D0

	  ELSE IF ( ND .GT. 73.D0 .AND. ND .LE. 580.D0 ) THEN
     		RE = 10.D0**(-1.7095D0 + 1.33438D0*A - 0.11591D0*A**2.D0)

	  ELSE IF ( ND .GT. 580.D0 .AND. ND .LE. 1.55D7 ) THEN
    	 	RE = 10.D0**(-1.81391D0 + 1.34671D0*A - 0.12427D0*A**2.D0 + 0.006344D0*A**3.D0)

	  ELSE
     		RE = 10.D0**(5.33283D0 - 1.21728D0*A + 0.19007D0*A**2.D0 - 0.007005D0*A**3.D0)

	  END IF



!---------------------------------------------------------------------------------------------------
!-------- para achar DC , metodo Zheng and Yapa (2000) 


		MO = GRAVITY * VIS_A**(4.D0) * (RO_A - RO_P)/(RO_A**(2.D0) * TS**(3.D0))
	        H1 = 59.3D0
	        J1 = 0.94D0  *H1**0.757D0
		EO1 = 3.D0*H1/(4.D0 * MO**(-0.149D0) *(VIS_A/9.D-4)**(-0.14D0))
                D1 = (EO1*TS/(GRAVITY * (RO_A - RO_P)))**(0.5D0)
		UT1 = (VIS_A) * MO**(-0.149D0) * (J1 - 0.857D0) /(RO_A * D1)
		X1 = DLOG10(D1)
                Y1 = DLOG10(UT1)
                A1 = 0.5D0
		B1 = DLOG10(0.711D0*(GRAVITY*(RO_A - RO_P)/RO_A)**0.5D0)

		D2 = 0.015D0
	        EO2 = GRAVITY*(RO_A - RO_P)*D2**(2.D0)/TS      
		H2 = 4.D0/3.D0 * EO2 * MO**(-0.149D0) *(VIS_A/9.D-4)**(-0.14D0)
		IF ( H2 .GT. 2.D0 .AND. H2 .LE. 59.3D0) THEN
                        J2 = 0.94D0  * H2**0.757D0
		ELSE
			J2 = 3.42D0 * H2**0.441D0
		END IF
		UT2 = (VIS_A) * MO**(-0.149D0) * (J2 - 0.857D0) /(RO_A * D2)
		X2 = DLOG10(D2)
                Y2 = DLOG10(UT2)
		A2 = (Y2 - Y1)/(X2 - X1)
		B2 = Y1 - A2*X1

		DC = 10.D0**((B2 - B1)/(A1 - A2))


!-------- parametros auxiliares -> Clift et al (1978)
		MO = GRAVITY * VIS_A**(4.D0) * (RO_A - RO_P)/(RO_A**(2.D0) * TS**(3.D0))
	        EO = GRAVITY * (RO_A - RO_P) * DIAM**(2.D0) /TS      
		H = 4.D0/3.D0 * EO * MO**(-0.149D0) * (VIS_A/9.D-4)**(-0.14D0)



!-------- calculo de wp para goticulas sem hidrato
!-------- Clift et al (1978)
	  IF (DIAM .EQ. 0.D0) THEN
		WP = 0.D0

	  ELSE  IF (DIAM .LT. 0.001D0) THEN
		WP = RE*VIS_A/(RO_A*DIAM)

	  ELSE  IF ( DIAM .GE. 0.001D0 .AND. DIAM .LE. DC) THEN
		IF ( H .LE. 59.3D0) THEN
                        J = 0.94D0  * H**(0.757D0)
		ELSE
			J = 3.42D0 * H**(0.441D0)
		END IF

		WP = (VIS_A) * MO**(-0.149D0) * (J - 0.857D0) /(RO_A * DIAM)


           ELSE
		WP = 0.711D0*(GRAVITY*DIAM*(RO_A - RO_P)/RO_A)**0.5D0

	   END IF




!-------- calculo de wp para goticulas com hidrato
!-------- Bigalke et al (2010)
	IF ( DIAM_HYD .NE. 0.D0) THEN
		WP0 = WP
			
		700  WP1 = WP

 		WE = DIAM*RO_A*WP1**(2.d0) /TS
 		DEF = 2.D0/( 3.974D-3*(WE - 12.62D0)**(2.D0) - 7.186D-4*(EO - 17.87D0)**(2.D0) &
 + 3.28D-5*EO*WE*((EO - 27.77D0)*(WE - 8.405D0) + 67.08D0) + 1.130D0 )
		FAC = 9.D0/RE**(0.5D0) + 0.9D0*(0.75D0*EO**(2.D0))/(0.75d0*EO**(2.D0) + 4.5D0) 
 		CD = DEF * FAC 
 
		WP = (4.D0 * DIAM * GRAVITY *(RO_A - RO_P) /(3.D0 * CD * RO_A) )**0.5D0


  		IF (DABS(WP - WP1)/WP1 .GT. 1.D-8) then
			RE = WP*(RO_A*DIAM)/VIS_A
			GO TO 700
		END IF

		IF ( WP .GT. WP0 ) WP = WP0

	END IF


!-------------------------
  END SUBROUTINE VEL_ASCENSAO_BG








   SUBROUTINE INI_AMB( LON_REF, LAT_REF, PROF_REF, LONOUT, LATOUT, PROFOUT, XPOS, YPOS, ZPOS, OPT)  


! carrega os arquivo de onde virao os dados do ambiente

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
    DOUBLE PRECISION  , INTENT(INOUT)   :: LONOUT , LATOUT  , PROFOUT

    DOUBLE PRECISION  , INTENT(INOUT)   ::  XPOS   , YPOS   , ZPOS  

    DOUBLE PRECISION  , INTENT(IN)      ::  LAT_REF , LON_REF , PROF_REF   
    INTEGER           , INTENT(IN)      ::  OPT  ! se OPT = 1 , entra LAT , LON , PROF , sai     X_POS, Y_POS, Z_POS             
                                                 ! se OPT = 2 , sai   LAT , LON , PROF , sentra  X_POS, Y_POS, Z_POS

! ------------------------------------------------------- 
      IF(OPT .EQ. 1 ) THEN

! --- opcao 1 serve quando precisa transformar lat e lon em x e y em relacao a um referencial;
!     por exemplo, quando faz continuacao de uma rodada, calcula e guarda um lag inicial em 
!     metros entre o lat lon inicial e o referencial. 
!     Isso por que a integracao sempre comeca com x e y zerados.
!     Em z tambem, guarda um lag z, que seria o lag inicial com o referencial 
!     
         YPOS  = -( LAT_REF -  LATOUT ) * PI * R_TERRA /180.D0
         XPOS  = -( LON_REF -  LONOUT ) * PI * R_TERRA*COS(LATOUT*PI/180.D0)/180.D0
         ZPOS = -( PROFOUT - PROF_REF)


      ELSE 

! --- opcao 2 transforma x e y em lat e lon. Isso eh feito para dar continuidade a uma proxima
!     rodada;
!     por exemplo, quando para uma rodada, para continuar eu tenho que informar a lat lon e 
!     profundidade. Note que eu somei os lags antes de entrar nesta subroutina
!    
!     Aqui deve-se ter cuidado com a prof, metros entre o lat lon inicial e o referencial. 
!     Isso por que a integracao sempre comeca com x e y zerados.
!     Em z tambem, guarda um lag z, que seria o lag inicial com o referencial 
! 

         LATOUT  = LAT_REF + YPOS / ( (PI/180.D0)* R_TERRA)
         LONOUT  = LON_REF + XPOS / ( (PI/180.D0)* R_TERRA*COS(LATOUT*PI/180.D0) )   !!fix this lat_ref to particle lat
         PROFOUT = -ZPOS + PROF_REF


      ENDIF


 !----------------------
  END SUBROUTINE INI_AMB
 !----------------------



















end module 

