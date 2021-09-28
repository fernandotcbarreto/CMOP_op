program calc_api

  implicit none

  INTEGER :: NCOMP_OIL
  INTEGER ::  NCLAS_OIL				! numero de classe de tamanho da gota de oleo 
  integer :: ii,jj  !indexes to be used along all model
  double precision:: counttime, mas_ref, mass_test, vol_test
  double precision, dimension(:, :), allocatable:: height, massa1, massa2,  massa3, massa4, massa5, massa6
  double precision, dimension(:, :), allocatable:: massa7, massa8, massa9, massa10, massa11, massa12
  double precision, dimension(:, :), allocatable:: massa13, massa14, massa15, massa16, massa17, massa18
  double precision, dimension(:, :), allocatable:: massa19, massa20, massa21, massa22, massa23, massa24, massa25
  DOUBLE PRECISION , DIMENSION (:, :) , ALLOCATABLE :: rho_e, visc_e, mas_evap, porc_evap



  double precision, dimension(:, :), allocatable:: area1, area2,  area3, area4, area5, area6
  double precision, dimension(:, :), allocatable:: area7, area8, area9, area10, area11, area12
  double precision, dimension(:, :), allocatable:: area13, area14, area15, area16, area17, area18
  double precision, dimension(:, :), allocatable:: area19, area20, area21, area22, area23, area24, area25

  double precision, dimension(:, :), allocatable:: diam1, diam2,  diam3, diam4, diam5, diam6
  double precision, dimension(:, :), allocatable:: diam7, diam8, diam9, diam10, diam11, diam12
  double precision, dimension(:, :), allocatable:: diam13, diam14, diam15, diam16, diam17, diam18
  double precision, dimension(:, :), allocatable:: diam19, diam20, diam21, diam22, diam23, diam24, diam25


  
  double precision:: massdis1, massdis2,  massdis3, massdis4, massdis5, massdis6
  double precision:: massdis7, massdis8, massdis9, massdis10, massdis11, massdis12
  double precision:: massdis13, massdis14, massdis15, massdis16, massdis17, massdis18
  double precision:: massdis19, massdis20, massdis21, massdis22, massdis23, massdis24, massdis25

  double precision:: masstest

!!!!!!!!!!!!!!!!!!3D ARRAYS
double precision, dimension(:, :, :), allocatable:: MASSCOMP   !THIS ARRAY IS (PARTICLE ID, TIME, OIL COMPONENT)

!!!!!!


  

  double precision::scp1, scp2   !specific mass of fractions
  double precision::porc1,porc2  !% of volume of each fraction
  double precision:: mwf1, mwf2, vp1, vp2, P_A, aux, DIAM_HYD


     DOUBLE PRECISION :: SG , MOIL 							! gravidade espeficica e massa do oleo
     DOUBLE PRECISION    :: PM_OIL_OUT				! peso molecular do oleo
     DOUBLE PRECISION    :: RO_OIL_OUT				! densidade do oleo na boca do tubo 
     DOUBLE PRECISION    :: TS_OIL_OUT
     DOUBLE PRECISION    :: VIS_DIN_OIL_OUT
     DOUBLE PRECISION   :: spmt   !density of particle

     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE ::  DIAM_OUT_OIL   		! diametro inicial de cada classe
     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE ::  VOLU_OUT_OIL  		! proporcao inicial do volume de cada classe
     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE ::  FRAC_MASS_OUT		! fracao massica do componente do oleo na boca do tubo

     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE :: PM_COMP_OIL , SOLU_COMP_OIL 	! peso molecular e solubilidade do compomente do oleo
     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE :: RO_COMP_OIL , RO_COMP_OIL_15	! densidade do componente do oleo
     DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE :: CP_COMP_OIL , TEB_COMP_OIL	! calor especifico e temperatura de ebulicao 


								! do componente do oleo
     DOUBLE PRECISION  :: PCRI_O , VCRI_O, TCRI_O , FAC_O

    DOUBLE PRECISION , DIMENSION (:) , ALLOCATABLE :: V_COMP , MOL_COMP	! volume e mol do componente do oleo
    DOUBLE PRECISION :: MOL_TOT  , V_TOT , FRAC_TOT					! mol , volume e fracao massica de oleo em vc

     DOUBLE PRECISION  :: CP_OIL  , VIS_DIN_OIL

     double precision, dimension(:,:), allocatable:: massa, vol, massae, vole  !mass oil, vol oil, mass emulsi, vol emulsi

     double precision:: API , VAZAO_OIL_OUT,TEMP_OUT
   DOUBLE PRECISION :: TS_VC						! tensao superficial de vc
     integer comps, DISSOLVE, EMULSI, EVAP_TURN


!=========================================================

  SUBROUTINE COMPONENTS ( API , VAZAO , DT , TEMP_OUT , &
	         	     RO_OIL , PM_OIL , TS_OIL , VIS_DIN_OIL )

  IMPLICIT NONE

!### calcular o peso molecular do oleo com base nas fracoes dos componetes

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(INOUT) :: API				! escala que caracteriza tipo de oleo
   DOUBLE PRECISION , INTENT(IN)    :: VAZAO	 			! vazao de saida do oleo
   DOUBLE PRECISION , INTENT(IN)    :: TEMP_OUT
   DOUBLE PRECISION , INTENT(IN)    :: DT

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(INOUT) :: RO_OIL				! densidade do oleo
   DOUBLE PRECISION , INTENT(OUT)   :: PM_OIL				! peso molecular do oleo
   DOUBLE PRECISION , INTENT(OUT)   :: TS_OIL
   DOUBLE PRECISION , INTENT(OUT)   :: VIS_DIN_OIL

! ------------------------------------------------------- parametros locais

   DOUBLE PRECISION :: MOIL_COMP  							! massa do componente do oleo
   DOUBLE PRECISION :: API_GOLF , API_JUB , API_LULA , API_FRADE			! fator que caracteriza o tipo de oleo
   DOUBLE PRECISION :: TB
   DOUBLE PRECISION :: MOILC7 , VOILC7 , MOLINT , MOLVOL , SGC7 , TF , VOLINT , MOLC7 , PMC7 , FRAC7 , BPP
   DOUBLE PRECISION :: VABP , FRAC_VOL

   DOUBLE PRECISION :: FRAC_MOL , MABP , CABP_1 , CABP , MeABP 

   DOUBLE PRECISION :: KUOP , T_PC , P_PC , P1  , PR1 , TR 
   DOUBLE PRECISION :: A , A1 , A2 , A3 , A4 , A5 , SP
   DOUBLE PRECISION :: VIS_REF , VIS_COR 
   DOUBLE PRECISION :: VIS_CIN_1							! viscosidade cinematica a 100F em cSt
   DOUBLE PRECISION :: VIS_CIN_2							! viscosidade cinematica a 210F em cSt
   DOUBLE PRECISION :: VIS_CIN_OIL						! viscosidade cinematica de vc 
   DOUBLE PRECISION :: VIS_DIN_SP 						! viscosidade dinamica do oleo em T_VC e P_A

   DOUBLE PRECISION :: T_RANK , TPC
   DOUBLE PRECISION :: T1 , T2
   INTEGER :: I, comps


  ALLOCATE ( V_COMP(NCOMP_OIL)    ) 
  ALLOCATE ( MOL_COMP(NCOMP_OIL)  )
  ALLOCATE ( FRAC_MASS_OUT(NCOMP_OIL) )
  ALLOCATE ( PM_COMP_OIL  (NCOMP_OIL) )
  ALLOCATE ( SOLU_COMP_OIL(NCOMP_OIL) )
  ALLOCATE ( RO_COMP_OIL_15(NCOMP_OIL))
  ALLOCATE ( RO_COMP_OIL  (NCOMP_OIL) )
  ALLOCATE ( TEB_COMP_OIL (NCOMP_OIL) )     
  ALLOCATE ( CP_COMP_OIL  (NCOMP_OIL) )



!-----------------------------------------------------------------------------------
!## cada tipo de oleo tem um conjunto de fracoes diferentes - esses sao para os oleos dos campos de Golfinho, Jubarte e Lula,
!com base na composicao media do oleo brasileiro dado por Zilio, Pinho (2002) [Cenpes].

!##### Campo do Golfinho
 API_GOLF = 28.80D0

!#### Campo de Jubarte
 API_JUB = 19.30D0

!#### Campo de Lula
 API_LULA = 31.D0

!#### Campo de FRADE
 API_FRADE = 19.60D0


! Composicao oléo brasileiro -> Zilio, Pinho (2002) Cenpes
! 21 < API <= 29 -> 55% saturado ; 23% aromático ; 22% polares
!      API <= 21 -> 45% saturado ; 29% aromático ; 26% polares



! NCOMP -> Johansen(2003) 
 NCOMP_OIL = 25				   ! numero de componentes do oleo


! PM_COMP_OIL = peso molecular de cada componente do oleo em G/MOL -> Johansen(2003) 
! SOLU_COMP_OIL = solubilidade em agua de cada componente do oleo em KG/M3 -> Johansen(2003)
! RO_COMP_OIL = densidade de cada componente do oleo em KG/M3 -> Johansen(2003)
! TEB_COMP_OIL = temperatura de ebulicao de cada componente do oleo em °F -> API Handbook(1997)
! CP_COMP_OIL = calor especifico de cada componente do oleo em Btu/lb°F -> API Handbook(1997)




   I = 1  ;  PM_COMP_OIL(I)  = 37.D0     ;  SOLU_COMP_OIL(I) = 40.D-03   ;  RO_COMP_OIL_15(I) = 615.D0  ;  
             TEB_COMP_OIL(I) = 17.38D0   ;  CP_COMP_OIL(I)   = 0.3431D0
   I = 2  ;  PM_COMP_OIL(I)  = 66.D0     ;  SOLU_COMP_OIL(I) = 95.D-03   ;  RO_COMP_OIL_15(I) = 673.D0  ;  
             TEB_COMP_OIL(I) = 90.53D0   ;  CP_COMP_OIL(I)   = 0.4137D0
   I = 3  ;  PM_COMP_OIL(I)  = 80.5D0    ;  SOLU_COMP_OIL(I) = 32.5D-03  ;  RO_COMP_OIL_15(I) = 697.D0  ;  
             TEB_COMP_OIL(I) = 149.70D0  ;  CP_COMP_OIL(I)   = 0.4776D0
   I = 4  ;  PM_COMP_OIL(I)  = 99.D0     ;  SOLU_COMP_OIL(I) = 9.D-03    ;  RO_COMP_OIL_15(I) = 712.D0  ;  
	     TEB_COMP_OIL(I) = 185.98D0  ;  CP_COMP_OIL(I)   = 0.4854D0
   I = 5  ;  PM_COMP_OIL(I)  = 113.D0    ;  SOLU_COMP_OIL(I) = 4.35D-03  ;  RO_COMP_OIL_15(I) = 753.D0  ;  
	     TEB_COMP_OIL(I) = 237.95D0  ;  CP_COMP_OIL(I)   = 0.4407D0
   I = 6  ;  PM_COMP_OIL(I)  = 127.D0    ;  SOLU_COMP_OIL(I) = 2.05D-04  ;  RO_COMP_OIL_15(I) = 764.D0  ;  
	     TEB_COMP_OIL(I) = 284.88D0  ;  CP_COMP_OIL(I)   = 0.4999D0
   I = 7  ;  PM_COMP_OIL(I)  = 78.D0     ;  SOLU_COMP_OIL(I) = 1.780D0   ;  RO_COMP_OIL_15(I) = 884.D0  ; 
 	     TEB_COMP_OIL(I) = 176.16D0  ;  CP_COMP_OIL(I)   = 0.4113D0
   I = 8  ;  PM_COMP_OIL(I)  = 92.D0     ;  SOLU_COMP_OIL(I) = 5.15D-01  ;  RO_COMP_OIL_15(I) = 880.D0  ;  
	     TEB_COMP_OIL(I) = 231.13D0  ;  CP_COMP_OIL(I)   = 0.3995D0
   I = 9  ;  PM_COMP_OIL(I)  = 106.D0    ;  SOLU_COMP_OIL(I) = 1.75D-01  ;  RO_COMP_OIL_15(I) = 875.D0  ;  
	     TEB_COMP_OIL(I) = 285.15D0  ;  CP_COMP_OIL(I)   = 0.4057D0
   I = 10 ;  PM_COMP_OIL(I)  = 120.D0    ;  SOLU_COMP_OIL(I) = 5.75D-02  ;  RO_COMP_OIL_15(I) = 875.D0  ;  
	     TEB_COMP_OIL(I) = 326.84D0  ;  CP_COMP_OIL(I)   = 0.4131D0
   I = 11 ;  PM_COMP_OIL(I)  = 141.5D0   ;  SOLU_COMP_OIL(I) = 1.25D-02  ;  RO_COMP_OIL_15(I) = 880.D0  ;  
	     TEB_COMP_OIL(I) = 391.48D0  ;  CP_COMP_OIL(I)   = 0.4312D0
   I = 12 ;  PM_COMP_OIL(I)  = 140.5D0   ;  SOLU_COMP_OIL(I) = 5.75D-05  ;  RO_COMP_OIL_15(I) = 773.D0  ;  
	     TEB_COMP_OIL(I) = 325.39D0  ;  CP_COMP_OIL(I)   = 0.4549D0
   I = 13 ;  PM_COMP_OIL(I)  = 156.5D0   ;  SOLU_COMP_OIL(I) = 2.2D-05   ;  RO_COMP_OIL_15(I) = 810.D0  ;  
	     TEB_COMP_OIL(I) = 418.52D0  ;  CP_COMP_OIL(I)   = 0.4378D0 
   I = 14 ;  PM_COMP_OIL(I)  = 185.5D0   ;  SOLU_COMP_OIL(I) = 2.5D-06   ;  RO_COMP_OIL_15(I) = 816.D0  ;  
	     TEB_COMP_OIL(I) = 490.08D0  ;  CP_COMP_OIL(I)   = 0.4461D0
   I = 15 ;  PM_COMP_OIL(I)  = 215.5D0   ;  SOLU_COMP_OIL(I) = 1.D-05    ;  RO_COMP_OIL_15(I) = 823.D0  ;  
	     TEB_COMP_OIL(I) = 550.62D0  ;  CP_COMP_OIL(I)   = 0.4813D0
   I = 16 ;  PM_COMP_OIL(I)  = 238.D0    ;  SOLU_COMP_OIL(I) = 1.D-06    ;  RO_COMP_OIL_15(I) = 828.D0  ;  
	     TEB_COMP_OIL(I) = 606.45D0  ;  CP_COMP_OIL(I)   = 0.4592D0
   I = 17 ;  PM_COMP_OIL(I)  = 273.D0    ;  SOLU_COMP_OIL(I) = 1.D-06    ;  RO_COMP_OIL_15(I) = 818.D0  ;  
	     TEB_COMP_OIL(I) = 646.34D0  ;  CP_COMP_OIL(I)   = 0.4541D0
   I = 18 ;  PM_COMP_OIL(I)  = 317.5D0   ;  SOLU_COMP_OIL(I) = 1.D-06    ;  RO_COMP_OIL_15(I) = 823.D0  ;  
	     TEB_COMP_OIL(I) = 715.91D0  ;  CP_COMP_OIL(I)   = 0.4430D0
   I = 19 ;  PM_COMP_OIL(I)  = 465.D0    ;  SOLU_COMP_OIL(I) = 1.D-06    ;  RO_COMP_OIL_15(I) = 950.D0  ;  
	     TEB_COMP_OIL(I) = 812.93D0  ;  CP_COMP_OIL(I)   = 0.3607D0
   I = 20 ;  PM_COMP_OIL(I)  = 135.D0    ;  SOLU_COMP_OIL(I) = 2.75D-02  ;  RO_COMP_OIL_15(I) = 1015.D0 ;  
	     TEB_COMP_OIL(I) = 454.27D0  ;  CP_COMP_OIL(I)   = 0.3707D0
   I = 21 ;  PM_COMP_OIL(I)  = 163.D0    ;  SOLU_COMP_OIL(I) = 5.5D-03   ;  RO_COMP_OIL_15(I) = 1016.D0 ;  
	     TEB_COMP_OIL(I) = 510.12D0  ;  CP_COMP_OIL(I)   = 0.3722D0
   I = 22 ;  PM_COMP_OIL(I)  = 177.D0    ;  SOLU_COMP_OIL(I) = 3.65D-03  ;  RO_COMP_OIL_15(I) = 980.D0  ;  
	     TEB_COMP_OIL(I) = 643.D0    ;  CP_COMP_OIL(I)   = 0.2399D0
   I = 23 ;  PM_COMP_OIL(I)  = 222.5D0   ;  SOLU_COMP_OIL(I) = 1.005D-04 ;  RO_COMP_OIL_15(I) = 980.D0  ;  
	     TEB_COMP_OIL(I) = 797.9D0   ;  CP_COMP_OIL(I)   = 0.2379D0
   I = 24 ;  PM_COMP_OIL(I)  = 130.D0    ;  SOLU_COMP_OIL(I) = 51.D0     ;  RO_COMP_OIL_15(I) = 986.D0  ;  
	     TEB_COMP_OIL(I) = 381.7D0   ;  CP_COMP_OIL(I)   = 0.4872D0
   I = 25 ;  PM_COMP_OIL(I)  = 215.D0    ;  SOLU_COMP_OIL(I) = 1.5D-01   ;  RO_COMP_OIL_15(I) = 1015.D0 ;  
	     TEB_COMP_OIL(I) = 662.D0    ;  CP_COMP_OIL(I)   = 0.2248D0


! FRAC_MASS_OUT -> % do total da massa de oleo de cada componente - partindo dos valores de Johansen(2003), as fracoes foram 
!                  reajustadas para obter um valor de API proximo do inicial



 IF ( API .EQ. API_GOLF) THEN

      I = 1  ;  FRAC_MASS_OUT(I) = 1.48D-02 
      I = 2  ;  FRAC_MASS_OUT(I) = 1.19D-02 
      I = 3  ;  FRAC_MASS_OUT(I) = 1.66D-02 
      I = 4  ;  FRAC_MASS_OUT(I) = 1.19D-02   
      I = 5  ;  FRAC_MASS_OUT(I) = 3.71D-02  
      I = 6  ;  FRAC_MASS_OUT(I) = 1.18D-02 
      I = 7  ;  FRAC_MASS_OUT(I) = 6.8D-03
      I = 8  ;  FRAC_MASS_OUT(I) = 1.43D-02  
      I = 9  ;  FRAC_MASS_OUT(I) = 1.85D-02
      I = 10 ;  FRAC_MASS_OUT(I) = 1.36D-02   
      I = 11 ;  FRAC_MASS_OUT(I) = 1.2D-03 
      I = 12 ;  FRAC_MASS_OUT(I) = 2.04D-02 
      I = 13 ;  FRAC_MASS_OUT(I) = 4.24D-02 
      I = 14 ;  FRAC_MASS_OUT(I) = 6.73D-02  
      I = 15 ;  FRAC_MASS_OUT(I) = 4.86D-02
      I = 16 ;  FRAC_MASS_OUT(I) = 3.84D-02  
      I = 17 ;  FRAC_MASS_OUT(I) = 4.49D-02 
      I = 18 ;  FRAC_MASS_OUT(I) = 5.46D-02  
      I = 19 ;  FRAC_MASS_OUT(I) = 2.774D-01 
      I = 20 ;  FRAC_MASS_OUT(I) = 10.9D-03 
      I = 21 ;  FRAC_MASS_OUT(I) = 11.8D-03 
      I = 22 ;  FRAC_MASS_OUT(I) = 2.4D-03
      I = 23 ;  FRAC_MASS_OUT(I) = 2.4D-03
      I = 24 ;  FRAC_MASS_OUT(I) = 0.D0
      I = 25 ;  FRAC_MASS_OUT(I) = 2.2D-01
      
!      KUOP = 11.8

   ELSE IF (API .EQ. API_LULA) THEN

      I = 1  ;  FRAC_MASS_OUT(I) = 2.79D-02 
      I = 2  ;  FRAC_MASS_OUT(I) = 2.71D-02 
      I = 3  ;  FRAC_MASS_OUT(I) = 2.22D-02 
      I = 4  ;  FRAC_MASS_OUT(I) = 2.21D-02   
      I = 5  ;  FRAC_MASS_OUT(I) = 3.90D-02  
      I = 6  ;  FRAC_MASS_OUT(I) = 1.58D-02 
      I = 7  ;  FRAC_MASS_OUT(I) = 7.8D-03
      I = 8  ;  FRAC_MASS_OUT(I) = 1.53D-02  
      I = 9  ;  FRAC_MASS_OUT(I) = 1.75D-02
      I = 10 ;  FRAC_MASS_OUT(I) = 1.86D-02   
      I = 11 ;  FRAC_MASS_OUT(I) = 1.2D-03 
      I = 12 ;  FRAC_MASS_OUT(I) = 1.84D-02 
      I = 13 ;  FRAC_MASS_OUT(I) = 4.14D-02 
      I = 14 ;  FRAC_MASS_OUT(I) = 5.33D-02  
      I = 15 ;  FRAC_MASS_OUT(I) = 4.76D-02
      I = 16 ;  FRAC_MASS_OUT(I) = 3.74D-02  
      I = 17 ;  FRAC_MASS_OUT(I) = 4.39D-02 
      I = 18 ;  FRAC_MASS_OUT(I) = 4.06D-02  
      I = 19 ;  FRAC_MASS_OUT(I) = 2.654D-01 
      I = 20 ;  FRAC_MASS_OUT(I) = 10.9D-03 
      I = 21 ;  FRAC_MASS_OUT(I) = 11.8D-03 
      I = 22 ;  FRAC_MASS_OUT(I) = 2.4D-03
      I = 23 ;  FRAC_MASS_OUT(I) = 2.4D-03
      I = 24 ;  FRAC_MASS_OUT(I) = 0.D0
      I = 25 ;  FRAC_MASS_OUT(I) = 2.1D-01

!      KUOP = 11.8

   ELSE IF (API .EQ. API_JUB) THEN
       
      I = 1  ;  FRAC_MASS_OUT(I) = 1.02D-02 
      I = 2  ;  FRAC_MASS_OUT(I) = 0.78D-02 
      I = 3  ;  FRAC_MASS_OUT(I) = 1.06D-02 
      I = 4  ;  FRAC_MASS_OUT(I) = 0.88D-02   
      I = 5  ;  FRAC_MASS_OUT(I) = 1.84D-02  
      I = 6  ;  FRAC_MASS_OUT(I) = 0.64D-02 
      I = 7  ;  FRAC_MASS_OUT(I) = 0.56D-02  
      I = 8  ;  FRAC_MASS_OUT(I) = 1.11D-02  
      I = 9  ;  FRAC_MASS_OUT(I) = 1.03D-02 
      I = 10 ;  FRAC_MASS_OUT(I) = 1.02D-02    
      I = 11 ;  FRAC_MASS_OUT(I) = 0.16D-02  
      I = 12 ;  FRAC_MASS_OUT(I) = 1.04D-02 
      I = 13 ;  FRAC_MASS_OUT(I) = 1.32D-02 
      I = 14 ;  FRAC_MASS_OUT(I) = 1.39D-02 
      I = 15 ;  FRAC_MASS_OUT(I) = 1.37D-02 
      I = 16 ;  FRAC_MASS_OUT(I) = 1.49D-02 
      I = 17 ;  FRAC_MASS_OUT(I) = 1.38D-02 
      I = 18 ;  FRAC_MASS_OUT(I) = 1.41D-02 
      I = 19 ;  FRAC_MASS_OUT(I) = 33.56D-02 
      I = 20 ;  FRAC_MASS_OUT(I) = 7.74D-02 
      I = 21 ;  FRAC_MASS_OUT(I) = 7.80D-02 
      I = 22 ;  FRAC_MASS_OUT(I) = 2.8D-02   
      I = 23 ;  FRAC_MASS_OUT(I) = 2.6D-02 
      I = 24 ;  FRAC_MASS_OUT(I) = 0.D0  
      I = 25 ;  FRAC_MASS_OUT(I) = 2.6D-01

!      KUOP = 11.6


  ELSE IF (API .EQ. API_FRADE) THEN
       
      I = 1  ;  FRAC_MASS_OUT(I) = 1.02D-02 
      I = 2  ;  FRAC_MASS_OUT(I) = 1.80D-02 
      I = 3  ;  FRAC_MASS_OUT(I) = 1.56D-02 
      I = 4  ;  FRAC_MASS_OUT(I) = 1.38D-02   
      I = 5  ;  FRAC_MASS_OUT(I) = 1.84D-02  
      I = 6  ;  FRAC_MASS_OUT(I) = 1.64D-02 
      I = 7  ;  FRAC_MASS_OUT(I) = 1.66D-02  
      I = 8  ;  FRAC_MASS_OUT(I) = 1.61D-02  
      I = 9  ;  FRAC_MASS_OUT(I) = 1.53D-02 
      I = 10 ;  FRAC_MASS_OUT(I) = 1.52D-02    
      I = 11 ;  FRAC_MASS_OUT(I) = 0.56D-02  
      I = 12 ;  FRAC_MASS_OUT(I) = 1.52D-02 
      I = 13 ;  FRAC_MASS_OUT(I) = 1.32D-02 
      I = 14 ;  FRAC_MASS_OUT(I) = 1.49D-02 
      I = 15 ;  FRAC_MASS_OUT(I) = 1.37D-02 
      I = 16 ;  FRAC_MASS_OUT(I) = 1.49D-02 
      I = 17 ;  FRAC_MASS_OUT(I) = 1.18D-02 
      I = 18 ;  FRAC_MASS_OUT(I) = 1.01D-02 
      I = 19 ;  FRAC_MASS_OUT(I) = 31.06D-02 
      I = 20 ;  FRAC_MASS_OUT(I) = 7.04D-02 
      I = 21 ;  FRAC_MASS_OUT(I) = 7.00D-02 
      I = 22 ;  FRAC_MASS_OUT(I) = 2.0D-02   
      I = 23 ;  FRAC_MASS_OUT(I) = 2.1D-02 
      I = 24 ;  FRAC_MASS_OUT(I) = 0.D0  
      I = 25 ;  FRAC_MASS_OUT(I) = 25.3D-02


! experimento DeepBlow (johansen)
   ELSE IF (API .EQ. 36.5523D0) THEN    

      I = 1  ;  FRAC_MASS_OUT(I) = 2.70D-02 
      I = 2  ;  FRAC_MASS_OUT(I) = 2.30D-02 
      I = 3  ;  FRAC_MASS_OUT(I) = 2.95D-02 
      I = 4  ;  FRAC_MASS_OUT(I) = 2.30D-02   
      I = 5  ;  FRAC_MASS_OUT(I) = 5.73D-02  
      I = 6  ;  FRAC_MASS_OUT(I) = 1.89D-02 
      I = 7  ;  FRAC_MASS_OUT(I) = 0.75D-02  
      I = 8  ;  FRAC_MASS_OUT(I) = 1.57D-02  
      I = 9  ;  FRAC_MASS_OUT(I) = 2.02D-02 
      I = 10 ;  FRAC_MASS_OUT(I) = 1.49D-02    
      I = 11 ;  FRAC_MASS_OUT(I) = 0.14D-02  
      I = 12 ;  FRAC_MASS_OUT(I) = 3.05D-02 
      I = 13 ;  FRAC_MASS_OUT(I) = 5.37D-02 
      I = 14 ;  FRAC_MASS_OUT(I) = 8.41D-02 
      I = 15 ;  FRAC_MASS_OUT(I) = 6.12D-02 
      I = 16 ;  FRAC_MASS_OUT(I) = 4.90D-02 
      I = 17 ;  FRAC_MASS_OUT(I) = 5.67D-02 
      I = 18 ;  FRAC_MASS_OUT(I) = 6.85D-02 
      I = 19 ;  FRAC_MASS_OUT(I) = 34.70D-02 
      I = 20 ;  FRAC_MASS_OUT(I) = 0.23D-02 
      I = 21 ;  FRAC_MASS_OUT(I) = 0.34D-02 
      I = 22 ;  FRAC_MASS_OUT(I) = 0.26D-02   
      I = 23 ;  FRAC_MASS_OUT(I) = 0.26D-02   
      I = 24 ;  FRAC_MASS_OUT(I) = 0.D0
      I = 25 ;  FRAC_MASS_OUT(I) = 0.D0



   ELSE IF (API .EQ. 45) THEN    

      I = 1  ;  FRAC_MASS_OUT(I) = 1
      print*, 'caozao'


   ENDIF



!!!!Calculation API
  
mas_ref = 100.  !kg

 DO I = 1 , NCOMP_OIL
  mass_test= mass_test + FRAC_MASS_OUT(I) * mas_ref

  vol_test = vol_test + ((FRAC_MASS_OUT(I) * mas_ref) / RO_COMP_OIL_15(I))
 enddo




   MOL_TOT = 0.D0
   V_TOT 	= 0.D0
   FRAC_TOT = 0.D0
   MABP 	= 0.D0
   CABP_1 	= 0.D0
   VABP 	= 0.D0


! --------------
   SG = 141.5D0 / (API + 131.5D0)                ! gravidade especifica do oleo inicial 

   RO_OIL = SG * 999.041D0                       ! densidade inicial do oleo em kg/m3 com base no api inicial

!   MOIL = VAZAO * DT * RO_OIL                    ! massa total do oleo em VC em kg   





   MOIL = vol(1,1) * RO_OIL         !all particles at first step are equal


!   massa(:,:) = VOL(:,:) * RO_OIL 

! --------------



   DO I = 1 , NCOMP_OIL

!---- fazer o PM do componente com base na densidade e na TEB
      TB = TEB_COMP_OIL(I) + 459.67D0				! ponto de ebulicao normal em graus Rankine (1°Rankine = 1F + 459.67)
      SG = RO_COMP_OIL_15(I) / 999.041D0 			! gravidade especifica 60F/60F                       

      PM_COMP_OIL(I) = 20.486D0 * DEXP(1.165D-4*TB -7.78712D0*SG + 1.1582D-3*TB*SG) * TB**(1.26007D0) * SG**(4.98308D0)  

   ENDDO


! --------------

  CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )

! --------------



!-----------------------------------------------------------------------------------
   DO I = 1 , NCOMP_OIL

!---- massa de cada classe do oleo      
      MOIL_COMP = MOIL * FRAC_MASS_OUT(I)                 		! massa de cada componente do oleo em kg

!---- numero de mols de cada classe do oleo
      MOL_COMP(I) = MOIL_COMP * 1000.D0/ PM_COMP_OIL(I) 		! numero de mols de cada componente do oleo em mol

      MOL_TOT = MOL_TOT + MOL_COMP(I)                     		! numero de mols total do oleo no VC em mol

!---- volume de cada classe do oleo
      V_COMP(I) = MOIL_COMP / RO_COMP_OIL(I)           		! volume de cada componente do oleo em m3

      
      V_TOT = V_TOT + V_COMP(I)                        		! volume total do oleo no VC em m3

!---- verificar se a fracao total dá 1
      FRAC_TOT = FRAC_TOT + FRAC_MASS_OUT(I)                    ! verificando se a fracao da 1

   ENDDO

!!!!!!!!!!!!!!!!!!!!1 mass fractions for particles
massa1(:,:)=moil * FRAC_MASS_OUT(1)
massa2(:,:)=moil * FRAC_MASS_OUT(2)
massa3(:,:)=moil * FRAC_MASS_OUT(3)
massa4(:,:)=moil * FRAC_MASS_OUT(4)
massa5(:,:)=moil * FRAC_MASS_OUT(5)
massa6(:,:)=moil * FRAC_MASS_OUT(6)
massa7(:,:)=moil * FRAC_MASS_OUT(7)
massa8(:,:)=moil * FRAC_MASS_OUT(8)
massa9(:,:)=moil * FRAC_MASS_OUT(9)
massa10(:,:)=moil * FRAC_MASS_OUT(10)
massa11(:,:)=moil * FRAC_MASS_OUT(11)
massa12(:,:)=moil * FRAC_MASS_OUT(12)
massa13(:,:)=moil * FRAC_MASS_OUT(13)
massa14(:,:)=moil * FRAC_MASS_OUT(14)
massa15(:,:)=moil * FRAC_MASS_OUT(15)
massa16(:,:)=moil * FRAC_MASS_OUT(16)
massa17(:,:)=moil * FRAC_MASS_OUT(17)
massa18(:,:)=moil * FRAC_MASS_OUT(18)
massa19(:,:)=moil * FRAC_MASS_OUT(19)
massa20(:,:)=moil * FRAC_MASS_OUT(20)
massa21(:,:)=moil * FRAC_MASS_OUT(21)
massa22(:,:)=moil * FRAC_MASS_OUT(22)
massa23(:,:)=moil * FRAC_MASS_OUT(23)
massa24(:,:)=moil * FRAC_MASS_OUT(24)
massa25(:,:)=moil * FRAC_MASS_OUT(25)


do comps=1, NCOMP_OIL
masscomp(:,:,comps)= moil * FRAC_MASS_OUT(comps)
enddo




!!!!!!!!!!!!!!!!!!!!!!



massa(:,:)= massa1(:,:)  + massa2(:,:)+massa3(:,:)+massa4(:,:)+massa5(:,:)+massa6(:,:)&
   +massa7(:,:)+massa8(:,:)+massa9(:,:)+massa10(:,:)+massa11(:,:)+massa12(:,:)+massa13(:,:)&
  +massa14(:,:)+massa15(:,:)+massa16(:,:)+massa17(:,:)+massa18(:,:)+massa19(:,:)+massa20(:,:)+massa21(:,:)&  
  +massa22(:,:)+massa23(:,:)+massa24(:,:)+massa25



massa(:,:)=sum(masscomp(1,1,:))


!massa(:,:)=sum(masscomp(:,:,:))



moil=massa(1,1)   !mass of any first step particle

spmt=(massa(1,1))/((massa1(1,1)/RO_COMP_OIL(1))+(massa2(1,1)/RO_COMP_OIL(2)) + (massa3(1,1)/RO_COMP_OIL(3)) + &    !corrigada para a temperatura ambiente
(massa4(1,1)/RO_COMP_OIL(4))+ (massa5(1,1)/RO_COMP_OIL(5)) + (massa6(1,1)/RO_COMP_OIL(6)) + (massa7(1,1)/RO_COMP_OIL(7))+&
(massa8(1,1)/RO_COMP_OIL(8))+(massa9(1,1)/RO_COMP_OIL(9))+(massa10(1,1)/RO_COMP_OIL(10))+(massa11(1,1)/RO_COMP_OIL(11))+&
(massa12(1,1)/RO_COMP_OIL(12))+(massa13(1,1)/RO_COMP_OIL(13))+(massa14(1,1)/RO_COMP_OIL(14))+(massa15(1,1)/RO_COMP_OIL(15))+&
(massa16(1,1)/RO_COMP_OIL(16)) + (massa17(1,1)/RO_COMP_OIL(17))+ (massa18(1,1)/RO_COMP_OIL(18))+(massa19(1,1)/RO_COMP_OIL(19))+&
(massa20(1,1)/RO_COMP_OIL(20))+(massa21(1,1)/RO_COMP_OIL(21))+(massa22(1,1)/RO_COMP_OIL(22)) + (massa23(1,1)/RO_COMP_OIL(23)) + &
(massa24(1,1)/RO_COMP_OIL(24)) + (massa25(1,1)/RO_COMP_OIL(25)))   


 aux=0
do comps=1, NCOMP_OIL
  aux=aux + masscomp(1,1,comps)/RO_COMP_OIL(comps)
enddo

spmt=massa(1,1)/aux



print*, spmt
print*, 'massa1', massa1(1, 1)
!print*, 'ooooo'
!print*,  V_TOT*spmt
!print*,  massa(1,1)




!!!!

!  print*, 'casa'
!  print*, v_tot
!  print*, vol(1,1)
!  print*, RO_OIL                    !these prints check input information
!  print*, spmt
!  print*, massa(1,1)
!  print*, vol(1,1)*spmt
!  print*, moil



!----------------------------------------------------------------------------------- 
   DO I = 1 , NCOMP_OIL

!---- calculo da fracao molar de cada componente
      FRAC_MOL = MOL_COMP(I)/MOL_TOT 					! fracao molar do componente do oleo

!---- calculo da fracao volumetrica de cada componente
      FRAC_VOL = V_COMP(I)/V_TOT					! fracao volumetrica do componente do oleo

!---- calculo do ponto de ebulicao medio molar em °F
      MABP = MABP + FRAC_MOL*TEB_COMP_OIL(I)				! temperatura de ebulicao molar do oleo em °F

!---- usar no calculo do ponto de ebulicao medio cubico em °F
      CABP_1 = CABP_1 + FRAC_VOL*TEB_COMP_OIL(I)**(1.D0/3.D0)		! temperatura de ebulicao cubica do oleo em °F


   END DO


!- calculo do ponto de ebulicao medio cubico em °F
   CABP = CABP_1**3.D0

!- calculo do ponto de ebulicao medio do oleo em °F
   MeABP = (MABP + CABP)/2.D0 

! ponto de ebulicao normal em graus Rankine (1°Rankine = 1K*9/5 = 1F + 459.67)
   TB = MeABP + 459.67D0

!----------------------------------------------------------------------------------- 


!- calculo do peso molecular do oleo 
   PM_OIL = MOIL * 1.D3/ MOL_TOT                         	! peso molecular do oleo em g/mol   (DA PARTICULA NO MEU CASO)

   PM_OIL_OUT = PM_OIL
!- calculo da densidade com base nas fracoes dos componentes
   RO_OIL = MOIL / V_TOT
!   print*, ro_oil
   				  	! densidade do oleo calculada com as fracoes dos componentes em kg/m3

!- calculo do novo API
!   API = 141.5D0 / (RO_OIL / 999.041D0) - 131.5D0                 ! api calculado com a nova densidade



!-----------------------------------------------------------------------------------
!#### calcular temperatura e volume criticas do oleo seguindo o API Handbook(1997), Spencer and Daubert (1973), Ahmed (2016)

   DO I = 1 , NCOMP_OIL

!---- calculo da fracao volumetrica de cada componente
      FRAC_VOL = V_COMP(I)/V_TOT					! fracao volumetrica do componente do oleo

!---- usar no calculo do ponto de ebulicao medio cubico em °F
      VABP = VABP + FRAC_VOL*TEB_COMP_OIL(I)				! temperatura de ebulicao media volumetrica do oleo em °F

   END DO


! Kesler and Lee Correlations:(Ahmed 2016)
   TCRI_O = 341.7D0 + 811.1D0*SG + (0.4244D0 + 0.1174D0*SG)*TB + (0.4669D0 - 3.26238D0*SG)*1.D5 / TB		! temperatura critica em Rankine
   TCRI_O = TCRI_O / (9.D0/5.D0)										! temperatura critica em K
										
   PCRI_O = DEXP(8.3634D0 - 0.0566D0/SG -(0.24244D0 + 2.2898D0/SG + 0.11857D0/SG**2.D0)*1.D-3 * TB  +  	&
		 	+ (1.4685D0 + 3.648D0/SG + 0.47227D0/SG**2.D0)*1.D-7 * TB**2.D0 		& 
			- (0.42019D0 + 1.6977D0/SG**2.D0)*1.D-10 * TB**3.D0 )					! pressao critica em psi
   PCRI_O = PCRI_O / 1.4504D-4											! pressao critica em Pa


   KUOP = TB**(1.D0/3.D0)/SG					! fator de Watson   
   TR = TB/(TCRI_O * (9.D0/5.D0))


   IF ( TR .LE. 0.8D0 ) THEN
   	P1 = 14.696D0						! pressao do ponto de ebuliçao normal em psia
   	PR1 = P1/(PCRI_O * 1.4504D-4)

	FAC_O = ( DLOG(PR1) - 5.92714D0 + 6.09648D0/TR + 1.28862D0*DLOG(TR) - 0.169347D0*TR**6.D0 ) / (15.2518D0 - 15.6875D0/TR - 13.4721D0*DLOG(TR) + 0.43577D0*TR**6.D0)
   ELSE
	FAC_O = -7.904D0 + 0.1352D0*KUOP - 0.007456D0*(KUOP)**2.D0 + 8.359D0*TR + (1.408D0 - 0.01063*KUOP)/TR
   END IF


! Hall and Yarborough Correlations (Ahmed 2016)
   VCRI_O = 0.025D0 * (PM_OIL/(SG**0.69D0))**(1.15D0)		! volume critico do oleo , ft3/lbmol
   VCRI_O = VCRI_O / PM_OIL					! volume critico do oleo , ft3/lb
   VCRI_O = VCRI_O * 2.8317D-2 / 0.45359D0			! volume critico do oleo , m3/kg


! Handbook of Petroleum             
   T_PC = 10.6443D0 * DEXP(-5.1747D-4 * TB - 0.54444D0*SG + 3.5995D-4*TB*SG) * TB**(0.81067D0) * SG**(0.53691D0)	! temperatura pseudocritica em Rankine
   T_PC = T_PC / (9.D0/5.D0)												! temperatura pseudocritica em K

   P_PC = 6.162D6 * DEXP(-4.725D-3*TB -4.8014*SG + 3.1939D-3*TB*SG) * TB**(-0.4844D0) * SG**(4.0846D0)  		! pressao pseudocritica em lb/in2 (psi)
   P_PC = P_PC / 1.4504D-4												! pressao pseudocritica em Pa



!-----------------------------------------------------------------------------------
!## as goticulas de oleo terao diametros divididos em 8 classes conforme visto no experimento DeepSpill 
!## (citado tb por Yapaetal 2012)

      NCLAS_OIL = 8			! numero de classes de gota de oleo (!DEEPBLOW )


      ALLOCATE ( DIAM_OUT_OIL(NCLAS_OIL)  )
      ALLOCATE ( VOLU_OUT_OIL(NCLAS_OIL)  ) 


!---- classe, diametro da classe (m) e porcentagem volumetrica da classe (!DEEPBLOW )
    I = 1 ; DIAM_OUT_OIL(I) =  1.D0*0.001D0 ; VOLU_OUT_OIL(I)  =   5.D0
    I = 2 ; DIAM_OUT_OIL(I) =  2.D0*0.001D0 ; VOLU_OUT_OIL(I)  =  29.D0 
    I = 3 ; DIAM_OUT_OIL(I) =  3.D0*0.001D0 ; VOLU_OUT_OIL(I)  =  31.D0
    I = 4 ; DIAM_OUT_OIL(I) =  4.D0*0.001D0 ; VOLU_OUT_OIL(I)  =  15.D0 
    I = 5 ; DIAM_OUT_OIL(I) =  5.D0*0.001D0 ; VOLU_OUT_OIL(I)  =   6.D0
    I = 6 ; DIAM_OUT_OIL(I) =  6.D0*0.001D0 ; VOLU_OUT_OIL(I)  =   2.D0 
    I = 7 ; DIAM_OUT_OIL(I) =  7.D0*0.001D0 ; VOLU_OUT_OIL(I)  =   7.D0 
    I = 8 ; DIAM_OUT_OIL(I) =  8.D0*0.001D0 ; VOLU_OUT_OIL(I)  =   5.D0




!-----------------------------------------------------------------------------------
!#### calcular a tensao interfacial do oleo seguindo o API Handbook(1997)
!- obter a temperatura do ponto de liberacao em °Rankine
   T_RANK = TEMP_OUT * 9.D0/5.D0						 ! temperatura de vc em °R
!- obter a temperatura pseudocritica do oleo em °Rankine
   TPC = T_PC * (9.D0/5.D0)	    

!- obter a tensao superficial do oleo
   TS_OIL = 673.7 * ((TPC - T_RANK)/TPC)**1.232D0 /KUOP                 ! tensao superficial em dyn/cm

   TS_OIL = TS_OIL * 1.D-3				                ! tensao superficial em N/m


!-----------------------------------------------------------------------------------
!### calcular a viscosidade dinamica do oleo

!- calcular viscosidade cinematica a 100F 
   T1 = 100.D0 + 459.67D0                             			! temperatura de 100F em rankine

   A = -1.35579D0 + 8.16059D-4*TB + 8.38505D-7*TB**2.D0
   VIS_REF = 10.D0**(A)
   
   A1 = 3.49310D+1 - 8.84387D-2*TB + 6.73513D-5*TB**2.D0 - 1.01394D-8*TB**3.D0
   A2 = -2.92649D0 + 6.98405D-3*TB - 5.09947D-6*TB**2.D0 + 7.49378D-10*TB**3.D0 
   A = A1 + A2 * KUOP
   VIS_COR = 10.D0**(A) 

   VIS_CIN_1 = VIS_REF + VIS_COR                  			! viscosidade cinematica a 100F em cSt -> 1cSt = D-6 m2/s  


!- calcular viscosidade cinematica a 210F 
   T2 = 210.D0 + 459.67D0                             			! temperatura de 210F em rankine

   A = -1.92353D0 + 2.41071D-4*TB + 0.5113D0*DLOG10(TB*VIS_CIN_1)

   VIS_CIN_2 = 10.D0**(A)                           			! viscosidade cinematica a 210F em cSt -> 1cSt = D-6 m2/s
   



!- calcular viscosidade cinematica na temperatura T_VC  a 1atm em cSt
   T_RANK = TEMP_OUT * 9.D0/5.D0						 ! temperatura de vc em °R

   A1 = VIS_CIN_1 + 0.7D0 + DEXP(-1.47D0 - 1.84D0*VIS_CIN_1 - 0.51D0*VIS_CIN_1**2.D0)
   A2 = VIS_CIN_2 + 0.7D0 + DEXP(-1.47D0 - 1.84D0*VIS_CIN_2 - 0.51D0*VIS_CIN_2**2.D0)
   A3 = ( DLOG10(LOG10(A1)) - DLOG10(DLOG10(A2)) )/( DLOG10(T1) - DLOG10(T2) )

   A4 = DLOG10(DLOG10(A1)) + A3*(DLOG10(T_RANK) - DLOG10(T1))
   A5 = 10.D0**(A4)
   A  = 10.D0**(A5)

   VIS_CIN_OIL = A - 0.7D0 - DEXP(-0.7487D0 - 3.295D0*(A - 0.7D0) + 0.6119D0*(A - 0.7D0)**2.D0 - 0.3193D0*(A - 0.7D0)**3.D0)    ! em cSt
   VIS_CIN_OIL = VIS_CIN_OIL*1.D-6             				! viscosidade cinematica do oleo em m2/s


!  print*, A, temp_out, a1, a2, a3, a4, a5

!- obter a viscosidade dinamica na pressao ambiente em cP -> 1cP = D-3 Pa.s  
   VIS_DIN_SP = VIS_CIN_OIL*RO_OIL          				! viscosidade dinamica em 1atm em Pa.s
   VIS_DIN_SP = VIS_DIN_SP/(1.D-3)          				! em cP 

  
   P_A = 101325.D0
   SP = P_A/(6894.757D0)                      				! em lbf/in2 -> 1lbf/in2 = 6894.757 Pa 

   A = SP*(-0.0102D0 + 0.04042D0*VIS_DIN_SP**0.181D0)/1000.D0
   VIS_DIN_OIL = 10.D0**A*VIS_DIN_SP           				! em cP

   VIS_DIN_OIL = VIS_DIN_OIL*1.D-3          				! viscosidade dinamica do oleo em P_A em Pa.s


 
!-----------------------------------------------------------------------------------

!- caso não haja pluma liquida
!   IF ( VAZAO .EQ. 0.D0 ) THEN
!        NCLAS_OIL = 1
!        NCOMP_OIL = 1
!        RO_OIL_OUT = 1.D0
!!        PM_OIL_OUT = 1.D0
!    END IF


 !-----------------------------
  END SUBROUTINE COMPONENTS


end program


