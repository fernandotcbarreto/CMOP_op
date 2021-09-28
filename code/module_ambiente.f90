MODULE module_ambiente
use processes 
! DECLARATION SECTION

   IMPLICIT NONE

     INTEGER, PARAMETER :: NX_OCEAN = 169, NY_OCEAN = 121 ,  NZ_OCEAN = 50 , NBIT = 1
     CHARACTER(LEN=20) :: FILE_AMBIENTE , COORD_Z , COORD_Y , COORD_X ! nome dos arquivos do ambiente 
                         
     TYPE LINE_FILE
        INTEGER(KIND=8)           :: YMDH
        INTEGER(KIND=8)           :: TIME_S
        INTEGER(KIND=8)           :: REC_N
        CHARACTER (LEN=152) :: FILE
     END TYPE LINE_FILE

    TYPE (LINE_FILE) :: LINE_ANT
    TYPE (LINE_FILE) :: LINE_POS
  
    INTEGER(KIND=8)     ::  I_CENTRAL , J_CENTRAL   ! indices centrais do cubo ambiente

    DOUBLE PRECISION    :: UNDEF

    DOUBLE PRECISION    ::  PM_H2O           	! peso molecular agua 
    DOUBLE PRECISION    ::  DMA			! diamtero da molecula de agua (cm)

    DOUBLE PRECISION    ::  U_A   ,  V_A   , W_A   , RO_A , P_A , T_A ,  SAL_A 	! variaveis ambientes
    DOUBLE PRECISION    ::  VIS_DIN_A , CP_A , TS_A
    DOUBLE PRECISION    ::  CDIF_HOR , CDIF_VER

    INTEGER(KIND=8)     ::  IF_LOG    		! se voltou Ok, IF_LOG = 0, se com problemas IF_LOG = 1 

    DOUBLE PRECISION    :: HORA
    INTEGER :: I_PA , OPT

    double precision:: LON_REF, LAT_REF, LONOUT, LATOUT, PROFOUT, XPOS, YPOS, ZPOS
!========================================================= 

  CONTAINS



!========================================================= 

   SUBROUTINE INI_AMB( LON_REF, LAT_REF, PROF_REF, LONOUT, LATOUT, PROFOUT, XPOS, YPOS, ZPOS, OPT)  !!!! I AM USING THE ONE FROM processes.f90


! carrega os arquivo de onde virao os dados do ambiente

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
    DOUBLE PRECISION  , INTENT(INOUT)   :: LONOUT , LATOUT  , PROFOUT

    DOUBLE PRECISION  , INTENT(INOUT)   ::  XPOS   , YPOS   , ZPOS  

    DOUBLE PRECISION  , INTENT(IN)      ::  LAT_REF , LON_REF , PROF_REF   
    INTEGER           , INTENT(IN)      ::  OPT  ! se OPT = 1 , entra LAT , LON , PROF , sai     X_POS, Y_POS, Z_POS             
                                                 ! se OPT = 2 , sai   LAT , LON , PROF , sentra  X_POS, Y_POS, Z_POS

! ------------------------------------------------------- 

      FILE_AMBIENTE = 'd-in/file_1_tmp'		! arquivo com ymdh , time em seg , rec e caminho do dat do ambiente 

      COORD_Y = 'd-in/ydef_lst'			! arquivo com as coordenadas y do dat
      COORD_X = 'd-in/xdef_lst'			! arquivo com as coordenadas x do dat
      COORD_Z = 'd-in/zdef_lst'			! arquivo com as coordenadas x do dat
      UNDEF    = -999.D0

      GRAVITY  = 9.8D0              		! gravidade
      PI       = 4.D0*DATAN(1.D0)      		! pi
      R_TERRA  = 6371000.D0			! raio da terra

      PM_H2O  = 18.D0     	! peso molecular agua g/mol
      DMA = 2.75D-8		! diametro da molecula de agua em cm

      IF(OPT .EQ. 1 ) THEN

! --- opcao 1 serve quando precisa transformar lat e lon em x e y em relacao a um referencial;
!     por exemplo, quando faz continuacao de uma rodada, calcula e guarda um lag inicial em 
!     metros entre o lat lon inicial e o referencial. 
!     Isso por que a integracao sempre comeca com x e y zerados.
!     Em z tambem, guarda um lag z, que seria o lag inicial com o referencial 
!     

         XPOS  = ( LAT_REF -  LATOUT ) * PI * R_TERRA*COS(LON_REF*PI/180.D0)/180.D0
         YPOS  = ( LON_REF -  LONOUT ) * PI * R_TERRA /180.D0
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

         LATOUT  = LAT_REF - XPOS / ( (PI/180.D0)* R_TERRA*COS(LON_REF*PI/180.D0)  ) 
         LONOUT  = LON_REF - YPOS / ( (PI/180.D0)* R_TERRA )
         PROFOUT = -ZPOS + PROF_REF


      ENDIF


 !----------------------
  END SUBROUTINE INI_AMB
 !----------------------

!=========================================================

!========================================================= 

  SUBROUTINE AMBIENTE_ALLPROP (  LONPOS , LATPOS , PROFPOS , PROF_REF 			, &
			I_CENTRAL , J_CENTRAL  				, &
			TIME_ANT_SEC , TIME_POS_SEC , TIME_I , IF_LOG 		, &
			U_A  , V_A  , W_A  					, &
			RO_A , P_A , T_A , SAL_A				, &
			VIS_DIN_A , CP_A , TS_A , CDIF_HOR , CDIF_VER		)

  IMPLICIT NONE

! define as variaveis do ambiente para o local da parcela I_VC
! por enquanto estamos com valores fixos

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN)  :: LONPOS , LATPOS , PROFPOS 
   DOUBLE PRECISION , INTENT(IN)  :: PROF_REF		! posicao de vc 

   INTEGER(KIND=8)  , INTENT(IN)  ::  I_CENTRAL , J_CENTRAL
   INTEGER(KIND=8)  , INTENT(IN)  ::  TIME_ANT_SEC , TIME_POS_SEC
   DOUBLE PRECISION , INTENT(IN)  ::  TIME_I 

! -------------------------------------------------------  parametros out
   DOUBLE PRECISION , INTENT(OUT) ::  U_A  , V_A  , W_A 
   DOUBLE PRECISION , INTENT(OUT) ::  RO_A , P_A , T_A , SAL_A 
   DOUBLE PRECISION , INTENT(OUT) ::  VIS_DIN_A , CP_A					! variaveis ambientes
   DOUBLE PRECISION , INTENT(OUT) ::  TS_A 						! tensao superficial

   DOUBLE PRECISION , INTENT(OUT) ::  CDIF_HOR , CDIF_VER 

   INTEGER(KIND=8) , INTENT(OUT)  ::  IF_LOG

! -------------------------------------------------------  parametros locais
   DOUBLE PRECISION  :: Z
   INTEGER(KIND=8)   :: TIME

!-----------------------------------------------------------------------------------

! determinar as velocidades, temperatura e salçinidade do ambiente

  TIME = NINT(TIME_I)

  CALL PROP_AMB_HIDTER ( LONPOS , LATPOS , PROFPOS                     , &
                         I_CENTRAL , J_CENTRAL                         , & 
                         TIME_ANT_SEC , TIME_POS_SEC , TIME            , & 
		   U_A , V_A  , W_A , T_A , SAL_A  , IF_LOG      )



!-----------------------------------------------------------------------------------

! determinar as propriedades da agua

   Z = PROF_REF - PROFPOS 

    CALL PROP_AMBIENTE ( Z , PROF_REF , T_A , SAL_A , &
			RO_A , VIS_DIN_A , CP_A , TS_A)



!-----------------------------------------------------------------------------------
! pressao ambiente em Pascal 

   P_A = 101325.D0 + RO_A*GRAVITY*(PROF_REF - Z)

!-----------------------------------------------------------------------------------


!------------ coeficiente difusao
!GRANDE ESCALA
    CDIF_HOR       = 1.D-1				
    CDIF_VER       = 1.D-3  


 !----------------------
  END SUBROUTINE AMBIENTE_ALLPROP
 !----------------------

!=========================================================

!=========================================================  
 
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
   DOUBLE PRECISION :: P_A

!-----------------------------------------------------------------------------------

! --- a densidade da agua pode ser calculada usando ou nao a pressao. Logo, temos duas opções:

! -------------------------------------------------------
! --- usar uma equacao para calcular a densidade independente da pressao -> segundo Sharqawy et al (2010) 
! -> salinidade em g/kg e temperatura em °C

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
   A1 = 999.841594D0 + 6.793952D-2*T_AMB - 9.09529D-3*T_AMB**2.D0 + 1.001685D-4*T_AMB**3.D0 - 1.120083D-6*T_AMB**4.D0 + 6.536332D-9*T_AMB**5.D0
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
 !----------------------

!=========================================================

!=========================================================   

  SUBROUTINE CORRECT_SAL ( ID_OG , PROF_REF , PROFPOS , TEMP , SAL , RO_AVC , E1_K , DM1 , POL1 , DMA , RO_COMP ,  PM_COMP ,CORRECAO )


  IMPLICIT NONE

! faz correcao sa solubilidade do componente devido a salinidade do ambiente 
! seguindo Masterton (1975) e Millero (2000)

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN)  :: TEMP 			! temp de vc em kelvin
   DOUBLE PRECISION , INTENT(IN)  :: SAL 			! salinidade de vc
   DOUBLE PRECISION , INTENT(IN)  :: RO_AVC 			! densidade da agua de vc
   DOUBLE PRECISION , INTENT(IN)  :: E1_K 			! parametro de energia do nao eletrolito
   DOUBLE PRECISION , INTENT(IN)  :: DM1 			! diametro do nao eletrolito
   DOUBLE PRECISION , INTENT(IN)  :: POL1 			! polaridade do nao eletrolito
   DOUBLE PRECISION , INTENT(IN)  :: DMA 			! diametro da molecula de agua
   DOUBLE PRECISION , INTENT(IN)  :: RO_COMP 		! densidade do componente
   DOUBLE PRECISION , INTENT(IN)  :: PM_COMP 		! peso molecular do componente
   DOUBLE PRECISION , INTENT(IN)  :: PROF_REF
   DOUBLE PRECISION , INTENT(IN)  :: PROFPOS

   INTEGER , INTENT(IN) 	  :: ID_OG			! identificador oleo/gas	

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(OUT)  :: CORRECAO			! correcao de solibilidade 

! -------------------------------------------------------  parametros locais
   INTEGER :: I , NCLA_ION

   DOUBLE PRECISION :: A1 , A2, A3 , A20 , CTE
   DOUBLE PRECISION :: CL					! clorinidade
   DOUBLE PRECISION :: KA , KB , KS				! coef. de salinidade
   DOUBLE PRECISION :: T_CEL					! temp em celsius
   DOUBLE PRECISION :: G_BA_CL , N_BA_CL , PM_BA		! parametros do acido borico
   DOUBLE PRECISION :: NT , K1					! mol equivalente de ion total e coef de NT 
   DOUBLE PRECISION :: MT , K2					! peso equivalente de ion total e coef 
   DOUBLE PRECISION :: VEA					! volume molar equivalente aparente
   DOUBLE PRECISION :: K3 
   DOUBLE PRECISION :: SAL1 , RO_0				! salin e densid da agua pura
   DOUBLE PRECISION :: AV					! cte que relaciona a interacao ion-agua 
   DOUBLE PRECISION :: DM12					! media entre diametro do nao eletrolito e da agua
   DOUBLE PRECISION :: DM1I					! media entre diametro do nao eletrolito e do ion	
   DOUBLE PRECISION :: EI_K					! parametro de energia do ion
   DOUBLE PRECISION :: CI , RO_A				! concentracao do ion
   DOUBLE PRECISION :: IV , II					! forca ionica da agua 
   DOUBLE PRECISION :: V0 , VOIL				! volume molar da agua pura e do gas
   DOUBLE PRECISION :: A , B , C , D , E		
   DOUBLE PRECISION :: T_AMB , COEF_OSM	, S_AMB			! temp em Celsius e coef. osmotico da agua , sal em kg/kg
   DOUBLE PRECISION :: K4 , CS					! forca ionica da agua 
   DOUBLE PRECISION :: Z  
		 

   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: GI_CL	! razao gramas de ion por kg de agua do mar e clorinidade
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: NI_CL	! razao mol de ion por kg de agua do mar e clorinidade
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: PM_ION	! peso equivalente do ion
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: ZI		! numero de eletrons do ion
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: POLI	! polaridade do ion
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: DI		! diametro do ion
   DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: EI		! carga do ion


! -------------------------------------------------------  via parameter
! DOUBLE PRECISION :: DMA

! -------------------------------------------------------

! --- composicao da agua do mar natural (Millero, 1973; 1974)
! na agua do mar, ha 10 ions com concentracao majoritaria: Na+, Mg2+, Ca2+, K+, Sr2+, Cl-, SO42-, HCO3-, Br- e F-
	NCLA_ION = 10

	ALLOCATE ( GI_CL(NCLA_ION) )
	ALLOCATE ( NI_CL(NCLA_ION) )
	ALLOCATE ( PM_ION(NCLA_ION) )
	ALLOCATE ( ZI(NCLA_ION) )
	ALLOCATE ( POLI(NCLA_ION) )
	ALLOCATE ( DI(NCLA_ION) )
	ALLOCATE ( EI(NCLA_ION) )


	I = 1   ;   GI_CL(I) = 0.55556D0   ;    NI_CL(I) = 0.0241653D0   ;   PM_ION(I) = 22.989D0  ;			!Na+
		    ZI(I)    = 10.D0  	   ;    POLI(I)  = 0.21D-24	 ;   DI(I)     = 2.02D-8   ;   EI(I) = +1.D0   
	I = 2   ;   GI_CL(I) = 0.06680D0   ;    NI_CL(I) = 0.0054970D0   ;   PM_ION(I) = 24.305D0			!Mg2+
		    ZI(I)    = 10.D0  	   ;    POLI(I)  = 0.12D-24	 ;   DI(I)     = 1.40D-8   ;   EI(I) = +2.D0
	I = 3   ;   GI_CL(I) = 0.02125D0   ;    NI_CL(I) = 0.0010604D0   ;   PM_ION(I) = 40.078D0			!Ca2+
		    ZI(I)    = 18.D0  	   ;    POLI(I)  = 0.53D-24	 ;   DI(I)     = 1.92D-8   ;   EI(I) = +2.D0
	I = 4   ;   GI_CL(I) = 0.02060D0   ;    NI_CL(I) = 0.0005269D0   ;   PM_ION(I) = 39.098D0			!K+
		    ZI(I)    = 18.D0  	   ;    POLI(I)  = 0.87D-24	 ;   DI(I)     = 2.64D-8   ;   EI(I) = +1.D0
	I = 5   ;   GI_CL(I) = 0.00041D0   ;    NI_CL(I) = 0.0000092D0   ;   PM_ION(I) = 87.620D0			!Sr2+
		    ZI(I)    = 36.D0  	   ;    POLI(I)  = 1.42D-24	 ;   DI(I)     = 2.22D-8   ;   EI(I) = +2.D0
	I = 6   ;   GI_CL(I) = 0.99894D0   ;    NI_CL(I) = 0.0282764D0   ;   PM_ION(I) = 35.450D0			!Cl-
		    ZI(I)    = 18.D0  	   ;    POLI(I)  = 3.02D-24	 ;   DI(I)     = 3.64D-8   ;   EI(I) = -1.D0
	I = 7   ;   GI_CL(I) = 0.14000D0   ;    NI_CL(I) = 0.0029148D0   ;   PM_ION(I) = 96.056D0			!SO42-
		    ZI(I)    = 50.D0  	   ;    POLI(I)  = 5.7D-24	 ;   DI(I)     = 5.80D-8   ;   EI(I) = -2.D0
	I = 8   ;   GI_CL(I) = 0.00735D0   ;    NI_CL(I) = 0.0001205D0   ;   PM_ION(I) = 61.016D0			!HCO3-
		    ZI(I)    = 32.D0  	   ;    POLI(I)  = 4.0D-24	 ;   DI(I)     = 2.66D-8   ;   EI(I) = -1.D0
	I = 9   ;   GI_CL(I) = 0.00348D0   ;    NI_CL(I) = 0.0000436D0   ;   PM_ION(I) = 79.904D0			!Br-
		    ZI(I)    = 36.D0  	   ;    POLI(I)  = 4.17D-24	 ;   DI(I)     = 3.96D-8   ;   EI(I) = -1.D0
	I = 10  ;   GI_CL(I) = 0.000067D0  ;    NI_CL(I) = 0.0000035D0   ;   PM_ION(I) = 18.998D0			!F-
		    ZI(I)    = 10.D0  	   ;    POLI(I)  = 0.99D-24	 ;   DI(I)     = 2.64D-8   ;   EI(I) = -1.D0


! o acido borico tambem esta presente na agua do mar, mas ele nao eh um ion	 
	G_BA_CL = 0.00132D0
	N_BA_CL = 0.0000218D0
	PM_BA = 60.550D0


! --- determinar a clorinidade de vc
	CL = SAL/1.80655D0
	

! - determinar coef de NT 
	K1 = 0.D0

	K3 = 0.D0
	K4 = 0.D0
	DO I = 1 , NCLA_ION
	   K1 = K1 + 1.D0/2.D0 * NI_CL(I)
	   K3 = K3 + 1.D0/2.D0 * GI_CL(I)*EI(I)**(2.d0) /PM_ION(I)
	   K4 = K4 + GI_CL(I) /PM_ION(I)
	END DO	
	
	K1 = K1 + N_BA_CL


! - determinar MT em g/eq em mol 
	K2 = 0.D0
	DO I = 1 , NCLA_ION
	   K2 = K2 + NI_CL(I)/K1 * ( PM_ION(I)/DABS(EI(I)) )
	END DO	
	
	MT  = K2 + N_BA_CL/K1 * PM_BA


!---- para determinar KS

IF ( ID_OG .EQ. 2 ) THEN 		! para o gas

! --- parametro KA 

	T_CEL = TEMP - 273.15D0 
	A1 = 0.00888D0  - 3.8303D-5*T_CEL  + 4.24242D-7*T_CEL**2.D0
	A2 = 1.26612D6  - 6.42316D3*T_CEL  + 60.9957D0 *T_CEL**2.D0
	A3 = 1.20214D14 - 6.32104D11*T_CEL + 5.4026D9 *T_CEL**2.D0

	KA = A1 + A2*DM1 + A3*DM1**2.D0


! --- parametro KB

! - determinar VEA ml/eq em mol-> Millero (1973)
	VEA = 10.464D0 + 0.2276D0*T_CEL - 4.949686D-3*T_CEL + 6.276458D-5*T_CEL - 3.920747D-7*T_CEL

! - determinar densidade da água pura - sal = 0
	SAL1 = 0.D0
	CTE = 1.D0			! cte qualquer pq outras variaveis ambiente nao sao uteis aqui
   	Z = PROF_REF - PROFPOS 
   	CALL PROP_AMBIENTE ( Z , PROF_REF , TEMP , SAL1 , &
		          RO_0 , CTE , CTE , CTE) 

	
	RO_0 = RO_0 * 1.D-3		! em g/ml

! - determinar AV -> Millero, Lepple (1973)
	AV = 1.D-3 * (MT - RO_0*VEA) * K1

! - determinar KB
	DM12 = (DM1 + DMA)/2.D0
	A1 = -1.498D24 * (AV/0.03600D0 - 0.05040D0) * ( E1_K**(0.5D0)*DM12**3.D0 + (0.996d-21*POL1)/DM12**3.D0 )

	A20 = 0.D0
	DO I = 1 , NCLA_ION
	   DM1I = (DM1 + DI(I))/2.D0
	   EI_K = 2.28D-8 * POLI(I)**(1.5D0) * ZI(I)**(0.5D0) / DI(I)**6.D0
	   CI = GI_CL(I) / (K3 * PM_ION(I))

	   A20 = A20 + CI * EI_K**(0.5D0) * DM1I**3.D0
	END DO 

	A2 = -2.922D21 * E1_K**(0.5) * A20


	KB = A1/TEMP + A2/TEMP 

! - determinar KS
	KS = KA + KB

! - determinar CORRECAO
   	Z = PROF_REF - PROFPOS 
   	CALL PROP_AMBIENTE ( Z , PROF_REF , TEMP , SAL , &
		          RO_A , CTE , CTE , CTE) 


	RO_A = RO_A * 1.D-3		! em kg/L

	IV = K3 *CL * RO_A

	CORRECAO = DEXP(KS * IV)
 
	CORRECAO = 1.D0/CORRECAO



ELSE 				! ks para oleo

! - determinar KS para o oleo
	SAL1 = 0.D0
	CTE = 1.D0			! cte qualquer pq outras variaveis ambiente nao sao uteis aqui
   	Z = PROF_REF - PROFPOS 
   	CALL PROP_AMBIENTE ( Z , PROF_REF , TEMP , SAL1 , &
		          RO_0 , CTE , CTE , CTE) 

	V0 = PM_H2O * 1.D-3 / RO_0

	VOIL = PM_COMP * 1.D-3 / RO_COMP


      ! coeficiente osmotico da agua do mar segundo Sharqawy et al (2010)
           T_AMB = TEMP        					
	S_AMB = SAL 

	II = 19.9201D0*(S_AMB) / (1.D3 - 1.00488D0*(S_AMB))

	A = 20.660673D0 - 432.578663D0/T_AMB - 3.711944D0*DLOG(T_AMB) + 8.637627D-3*T_AMB
	B = 2.303D0/II**(3.D0/2.D0) * ( (1.D0 + II**(1.D0/2.D0)) - 1.D0/(1.D0 + II**(1.D0/2.D0)) - 2.D0*DLOG(1.D0 + II**(1.D0/2.D0)) )
	C = -831.658611D0 + 17022.3989D0/T_AMB + 157.65271D0*DLOG(T_AMB) - 0.493D0*T_AMB + 2.595060D-4*T_AMB**2.D0
	D = 553.905988D0 - 11200.445D0/T_AMB - 105.239035D0*DLOG(T_AMB) + 0.333219D0*T_AMB - 1.773514D-4*T_AMB**2.D0
	E = -0.151125D0

	COEF_OSM = 1.D0 - ( A*B*II**(1.D0/2.D0) + C*II + D*II**(3.D0/2.D0) + E*II**2.D0 ) 

	KS = 3.86D0 * V0**(1.D0/3.D0) * VOIL**(2.D0/3.D0) * COEF_OSM		! m3/mol
	KS = KS * 1.D3								! L/mol


! - determinar CORRECAO
	RO_A = RO_AVC * 1.D-3		! em kg/L

	CS = K4 *CL * RO_A		! mol/L

	CORRECAO = DEXP(KS * CS)
 
	CORRECAO = 1.D0/CORRECAO


END IF


 !----------------------
  END SUBROUTINE CORRECT_SAL
 !----------------------

!========================================================= 

!=========================================================

     SUBROUTINE FILES_BETWEEN( YMDHM_OUT , TIME_I , LINE_ANT , LINE_POS , TIME_ANT_SEC , TIME_POS_SEC) 
!determina os arquivos do mercator antes e depois da hora do vazamento

     IMPLICIT NONE

! ------------------------------------------------------- parametros in
     INTEGER(KIND=8) , INTENT(IN)             ::  YMDHM_OUT  	! data do blowout
     DOUBLE PRECISION  , INTENT(IN)           ::  TIME_I       ! tempo de rodada em segundos

! ------------------------------------------------------- parametros out
     TYPE (LINE_FILE) , INTENT(OUT) :: LINE_ANT      ! informacoes do arquivo ambiente anterior
     TYPE (LINE_FILE) , INTENT(OUT) :: LINE_POS      ! informacoes do arquivo ambiente posterior


     INTEGER(KIND=8)  , INTENT(INOUT) ::  TIME_ANT_SEC , TIME_POS_SEC ! tempo 

! ------------------------------------------------------- parametros locais
     INTEGER(KIND=8)             ::  IIFF                           ! local
     INTEGER(KIND=8)             ::  YYYYMMDD , JDAY , JHOUR , JMIN ! local
     INTEGER(KIND=8)             ::  TIME_OUT_SEC 
     INTEGER(KIND=8) 		 ::  TIME

   !  ---------------------------------------------------------------------------

     TIME = NINT(TIME_I)


         CALL JULDAY ( YMDHM_OUT / 10000, JDAY )

         JMIN  =  YMDHM_OUT - (YMDHM_OUT / 100   )*100

         JHOUR =  YMDHM_OUT/100 - (YMDHM_OUT  / 10000 )*100
         TIME_OUT_SEC    = (JDAY)*24*60*60 + (JHOUR)*60*60  +  (JMIN)*60  +  TIME  ! YMDHM_OUT em segundos


         OPEN( 50 , FILE= FILE_AMBIENTE , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
  
           READ (50,*,IOSTAT = IIFF) LINE_ANT ; IF( IIFF .NE. 0 ) GO TO 600
           CALL JULDAY ( LINE_ANT%YMDH/100, JDAY )
           JHOUR           =  LINE_ANT%YMDH - (LINE_ANT%YMDH / 100)*100
           TIME_ANT_SEC    = (JDAY)*24*60*60 + (JHOUR)*60*60         ! TIME_ANT_SEC em segundos
           TIME_ANT_SEC    =   TIME_ANT_SEC - TIME_OUT_SEC   

         CLOSE(50)


         IF ( TIME_ANT_SEC .GT. 0 ) THEN
              WRITE(*,*) "  "
              WRITE(*,*) " *** PROCESSO INTERROMPIDO *** "
              WRITE(*,*) "     primeira data do arquivo ambiente eh mais recente que a data do blowout "   
              WRITE(*,*) "  " 
              WRITE(*,*) "     Reavalie data do blowout ou data do primeiro arquivo ambiente   "
              WRITE(*,*) "  "
              STOP
         ENDIF
         

   !  ---------------------------------------------------------------------------

      OPEN( 50 , FILE= FILE_AMBIENTE , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )    

 500  READ (50,*,IOSTAT = IIFF) LINE_POS ; IF( IIFF .NE. 0 ) GO TO 600

      TIME_POS_SEC   = TIME_ANT_SEC + (LINE_POS%TIME_S - LINE_ANT%TIME_S)

      IF (TIME_POS_SEC .GT.  TIME ) GO TO 600

      LINE_ANT     = LINE_POS
      TIME_ANT_SEC = TIME_POS_SEC

      GO TO 500

  600 CONTINUE


        CLOSE (50)

 !----------------------------
  END SUBROUTINE FILES_BETWEEN
 !----------------------------


!=========================================================
   SUBROUTINE JULDAY ( iyyyymmdd, jday)

!  -------------------------------------------------------
!  calcular o dia juliano de uma data
!  entra com ano, mes, dia, e sai o dia juliano
!  in  :  iyyyymmdd
!  out :  jday
!  -------------------------------------------------------

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
   INTEGER(KIND=8) , INTENT(IN) :: iyyyymmdd

! ------------------------------------------------------- parametros in
   INTEGER(KIND=8) , INTENT(OUT) :: jday

! ------------------------------------------------------- parametros locais
   INTEGER(KIND=8) :: iyyy , mm , id		! ano, mes e dia
   INTEGER(KIND=8) :: jy , jm , ja , igreg

!-----------------------------------------------------------------------------------
      PARAMETER (igreg=15+31*(10+12*1582))

      iyyy = iyyyymmdd/10000
      mm  = iyyyymmdd/100 - iyyy*100
      id = iyyyymmdd - iyyy*10000 - mm*100

      IF ( iyyy .EQ. 0 ) PAUSE 'there is no year zero.'
      IF ( iyyy .LT. 0) iyyy = iyyy + 1

      IF ( mm .GT. 2 ) THEN
           jy = iyyy

           jm = mm + 1
      ELSE
          jy = iyyy - 1

          jm = mm + 13
      END IF

      jday = INT(365.25D0*jy) + INT(30.6001D0*jm) + id + 1720995
      
      IF ( (id + 31*(mm + 12*iyyy)) .GE. igreg) THEN
            ja = INT(0.01D0*jy)
        
	    jday = jday + 2 - ja + INT(0.25D0*ja)
      END IF
 

 !----------------------
  END SUBROUTINE JULDAY
 !----------------------


!========================================================= 

    SUBROUTINE IJ_CENTRAL( LONOUT , LATOUT , I_CENTRAL , J_CENTRAL ) 

! descobrir os indices centrais do quadri-cubo 

    IMPLICIT NONE

! ------------------------------------------------------- parametros in
    DOUBLE PRECISION, INTENT(IN) :: LONOUT  , LATOUT 


! ------------------------------------------------------- parametros out
    INTEGER(KIND=8) , INTENT(OUT) :: I_CENTRAL , J_CENTRAL ! indices centrais do cubo ambiente

! ------------------------------------------------------- parametros locais
    DOUBLE PRECISION, DIMENSION (NX_OCEAN)   ::  X_COO 
    DOUBLE PRECISION, DIMENSION (NY_OCEAN)   ::  Y_COO 

    DOUBLE PRECISION         :: DIST                    ! na subroutina
    INTEGER(KIND=8)          :: IIFF , I , J            ! na subroutina

! ------------------------------------------------------- 


    OPEN ( 50, FILE= COORD_Y , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  Y_COO(I), I = 1 , NY_OCEAN )
    CLOSE( 50)

    OPEN ( 50, FILE= COORD_X , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  X_COO(J), J = 1 , NX_OCEAN )
    CLOSE( 50 )


    I_CENTRAL = 1
    DIST = ABS ( LONOUT -  X_COO(1) ) 
  
    DO I = 1 , NX_OCEAN
       IF (  ABS ( LONOUT -  X_COO(I) )  .LT.    DIST  ) THEN
           DIST =ABS ( LONOUT -  X_COO(I) )
           I_CENTRAL = I
       ENDIF

    ENDDO

    J_CENTRAL = 1
    DIST = ABS ( LATOUT -  Y_COO(1) )   
    DO J = 1 , NY_OCEAN
       IF (  ABS ( LATOUT -  Y_COO(J) )  .LT.    DIST  ) THEN
           DIST =ABS (LATOUT -  Y_COO(J) )
           J_CENTRAL = J
       ENDIF

    ENDDO


 !-------------------------
  END SUBROUTINE IJ_CENTRAL
 !-------------------------

!========================================================= 


   SUBROUTINE CORTA_MERCATOR ( I_CENTRAL , J_CENTRAL , LINE_FILE_INFO )

   IMPLICIT NONE


! ------------------------------------------------------- parametros in
   INTEGER(KIND=8) , INTENT(IN) :: I_CENTRAL , J_CENTRAL ! indices centrais do cubo ambiente

! ------------------------------------------------------- parametros out
   TYPE (LINE_FILE)  :: LINE_FILE_INFO      ! informacoes do arquivo ambiente anterior

! ------------------------------------------------------- parametros locais
   REAL, DIMENSION ( NX_OCEAN , NY_OCEAN , NZ_OCEAN )   ::  U_OCEAN , V_OCEAN , S_OCEAN , T_OCEAN
   REAL, DIMENSION ( NX_OCEAN , NY_OCEAN            )   ::  SSH_OCEAN

   REAL, DIMENSION ( -1:1 , -1:1, NZ_OCEAN  )           ::  U_RED , V_RED , S_RED, T_RED 
   REAL, DIMENSION ( -1:0 , -1:0, NZ_OCEAN  )           ::  W_RED
   REAL, DIMENSION ( -1:1 , -1:1 )                      ::  SSH_RED

   INTEGER(KIND=8)    :: IREC , I , J , K , IIFF
   INTEGER(KIND=8)    :: IX_1 , JY_1 
   DOUBLE PRECISION       :: X_MEIO , Y_MEIO , DX , DY  

   DOUBLE PRECISION, DIMENSION ( NX_OCEAN)   ::  X_COO 
   DOUBLE PRECISION, DIMENSION ( NY_OCEAN)   ::  Y_COO 
   DOUBLE PRECISION, DIMENSION ( NZ_OCEAN)   ::  Z_COO

   DOUBLE PRECISION       :: U_11 , U_21 , U_12 , U_22
   DOUBLE PRECISION       :: V_11 , V_21 , V_12 , V_22 
   DOUBLE PRECISION       :: DIV_U , DIV_V

   DOUBLE PRECISION       :: W_A  , W_B  , HALF_A , HALF_B

   DOUBLE PRECISION, DIMENSION ( NZ_OCEAN )   :: DIV_HOR

! ------------------------------------------------------- 

    OPEN ( 50, FILE= COORD_X , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  X_COO(J), J = 1 , NX_OCEAN )
    CLOSE( 50 )

    OPEN ( 50, FILE= COORD_Y , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  Y_COO(I), I = 1 , NY_OCEAN )
    CLOSE( 50)

    OPEN ( 50, FILE= COORD_Z , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  Z_COO(K), K = 1 , NZ_OCEAN )
    CLOSE( 50 )

 !-------------------------
!      WRITE(*,*) LINE_FILE_INFO  
! 2016012000 1641600 20 "/home/valdir/sub_mercator/d-mercator/outjunto_201601.dat"

      IREC = LINE_FILE_INFO%REC_N 

      OPEN(22, FILE = LINE_FILE_INFO%FILE , STATUS = 'UNKNOWN', FORM = 'UNFORMATTED'  ,           & 
               ACCESS = 'DIRECT',                                                                 & 
               RECL = (4* NX_OCEAN * NY_OCEAN * NZ_OCEAN + NX_OCEAN  * NY_OCEAN )* NBIT      )  

               READ( 22 , REC = IREC )  V_OCEAN , U_OCEAN , S_OCEAN , SSH_OCEAN, T_OCEAN

      CLOSE(22)

      DO J = -1 , 1 
      DO I = -1 , 1  
         DO K = 1 , NZ_OCEAN
            U_RED(I,J,K) = U_OCEAN(I_CENTRAL +I , J_CENTRAL +J, K) 
            V_RED(I,J,K) = V_OCEAN(I_CENTRAL +I , J_CENTRAL +J, K) 
            S_RED(I,J,K) = S_OCEAN(I_CENTRAL +I , J_CENTRAL +J, K) 
            T_RED(I,J,K) = T_OCEAN(I_CENTRAL +I , J_CENTRAL +J, K)

            IF ( U_RED(I,J,K) .LE. UNDEF ) U_RED(I,J,K)  = 0.D0  
            IF ( V_RED(I,J,K) .LE. UNDEF ) V_RED(I,J,K)  = 0.D0
         

         ENDDO
         SSH_RED(I,J) = SSH_OCEAN(I_CENTRAL +I , J_CENTRAL +J)
      ENDDO
      ENDDO

! -------------------------------------------------- velocidade vertical

      DO J = -1 , 0 !  sao 4 cubos, cada cubo eh identificado pelo indice inferior esquer
      DO I = -1 , 0

         IX_1 = I_CENTRAL + I

         JY_1 = J_CENTRAL + J
   

         X_MEIO = ( X_COO(IX_1) + X_COO(IX_1+1) )*0.5D0
         Y_MEIO = ( Y_COO(JY_1) + Y_COO(JY_1+1) )*0.5D0

         DX    = ( X_COO(IX_1+1) - X_COO(IX_1) ) * (PI/180.D0)* R_TERRA*COS(Y_MEIO*PI/180.D0) 
         DY    = ( Y_COO(JY_1+1) - Y_COO(JY_1) ) * (PI/180.D0)* R_TERRA 

         DO K = NZ_OCEAN , 1 , -1

            U_11 = U_RED(I   ,J   ,K)
            U_21 = U_RED(I+1 ,J   ,K)
            U_12 = U_RED(I   ,J+1 ,K)
            U_22 = U_RED(I+1 ,J+1 ,K)

            V_11 = V_RED(I   ,J   ,K)
            V_21 = V_RED(I+1 ,J   ,K)  
            V_12 = V_RED(I   ,J+1 ,K)
            V_22 = V_RED(I+1 ,J+1 ,K)




! --- calculo de du/dx e dv/dx no ponto medio entre IX_1 - IX_2 (JY_1 - IY_2)
            DIV_U = ( (U_22 + U_21)*0.5D0 - (U_12 + U_11)*0.5D0 ) /  DX
            DIV_V = ( (V_22 + V_12)*0.5D0 - (V_21 + V_11)*0.5D0 ) /  DY 

            DIV_HOR(K) =  DIV_U + DIV_V

         ENDDO  

         W_A = 0.D0
         HALF_A = ( Z_COO(NZ_OCEAN) - Z_COO(NZ_OCEAN-1) )*0.5D0   ! metade da distancia entre dois niveis



! ------------------------ 
         DO K = NZ_OCEAN , 2 , -1

            HALF_B = ( Z_COO(K) - Z_COO(K-1) )*0.5D0
            W_B = W_A - (HALF_B + HALF_A)*DIV_HOR(K)		 ! usando a divergencia


            W_RED(I,J,K) = (W_B + W_A)*0.5D0 

            W_A    = W_B
            HALF_A = HALF_B


         ENDDO 

         K = 1
         W_RED(I,J,K) =  W_A -  HALF_A *DIV_HOR(K)

      ENDDO
      ENDDO 

! ----------------------


      CALL SYSTEM("rm -f fort.50")
      OPEN(50, FILE = "fort.50" , STATUS = 'UNKNOWN', FORM = 'UNFORMATTED'  ,           & 
               ACCESS = 'DIRECT',                                                       & 
               RECL = (4 * 3*3*NZ_OCEAN + 3*3 + 2*2*NZ_OCEAN)* NBIT      )  
      WRITE(50,REC=1) V_RED , U_RED , S_RED , SSH_RED, T_RED , W_RED
      CLOSE(50)


  END SUBROUTINE CORTA_MERCATOR
 !-------------------------


!========================================================= 

  SUBROUTINE PESOS

! calcula os pesos das posicoes ao redor do ponto x, y , z

   IMPLICIT NONE

! ------------------------------------------------------- 

   DOUBLE PRECISION , DIMENSION(0:100 , 0:100 , 1:4) :: PESO_XY
   DOUBLE PRECISION , DIMENSION(4)                   :: PESO   , XY
   DOUBLE PRECISION                                  :: XYMIN  , XYMAX , PESO_DEN , SUM_P

   INTEGER(KIND=8) :: I , I_INI ,  I_FIM 
   INTEGER(KIND=8) :: J , J_INI ,  J_FIM 


 ! -------------------------------------------------------------------  pesos
      I_INI = 0
      I_FIM = 100
      
      J_INI = 0
      J_FIM = 100

    
      DO J = J_INI ,  J_FIM 
      DO I = I_INI ,  I_FIM 

          XY(1) = ( (I-I_INI)**2.D0 + (J-J_INI)**2.D0 )**0.5D0
          XY(2) = ( (I-I_FIM)**2.D0 + (J-J_INI)**2.D0 )**0.5D0
          XY(3) = ( (I-I_INI)**2.D0 + (J-J_FIM)**2.D0 )**0.5D0   
          XY(4) = ( (I-I_FIM)**2.D0 + (J-J_FIM)**2.D0 )**0.5D0

          XYMAX = MAXVAL(XY(:))
          XYMIN = MINVAL(XY(:))
 
          PESO_DEN = ( XYMIN + XYMAX )
          PESO(:) = ( PESO_DEN - XY(:))/PESO_DEN
 
          IF(PESO(1) .EQ. 1.D0 ) THEN ; PESO(2) = 0.D0 ; PESO(3) = 0.D0; PESO(4) = 0.D0 ; ENDIF
          IF(PESO(2) .EQ. 1.D0 ) THEN ; PESO(1) = 0.D0 ; PESO(3) = 0.D0; PESO(4) = 0.D0 ; ENDIF
          IF(PESO(3) .EQ. 1.D0 ) THEN ; PESO(1) = 0.D0 ; PESO(2) = 0.D0; PESO(4) = 0.D0 ; ENDIF  
          IF(PESO(4) .EQ. 1.D0 ) THEN ; PESO(1) = 0.D0 ; PESO(2) = 0.D0; PESO(3) = 0.D0 ; ENDIF

          SUM_P = SUM(PESO(:))
          PESO_XY(I,J,:) = PESO(:)/SUM_P

      ENDDO
      ENDDO

      CALL SYSTEM ("rm -f pesos.dat")
      OPEN(22, FILE= 'pesos.dat', STATUS= 'UNKNOWN', FORM= 'UNFORMATTED' , & 
         ACCESS= 'DIRECT', RECL = 2* 101*101*4 * NBIT ) 

         WRITE(22 , REC = 1) PESO_XY 
      CLOSE(22)
  

!----------------------
  END SUBROUTINE PESOS
!----------------------


!========================================================= 

  SUBROUTINE PROP_AMB_HIDTER ( LONPOS , LATPOS ,PROFPOS      		,    & 
                               I_CENTRAL    , J_CENTRAL         		,    & 
                               TIME_ANT_SEC , TIME_POS_SEC , TIME     	,    & 
		         U , V , W , T , SAL , IF_LOG    	)

! ---------------------------------------------------------------------------- 
!            
! - entra com 
! 
!     - a posicao do blow up da parcela (LAT_OUT, LON_OUT, PROF_OUT), 
! 
!     - a distancia da parcela (X , Y, Z) relativa a posicao do blowup, 
! 
!     - o valor dos indices horizontais mais proximo do blowup (I_CENTRAL ,
!       J_CENTRAL), i.e. indices da matriz da coordenada mercator  matriz coordenada  
! 
!     - TIME_ANT_SEC , TIME_POS_SEC , TIME : tempos pra fazer a interpolacao linear
!       i.e., tempo do mercator anterior, tempo do mercator posterior, e tempo 
!       atual
!    
! - sai com
!   U , V , W , T , SAL   interpolados para a posicao e tempo da parcela 
!   IF_LOG =  
! ---------------------------------------------------------------------------- 

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN)  :: LONPOS , LATPOS , PROFPOS    
   INTEGER(KIND=8)  , INTENT(IN)  :: I_CENTRAL , J_CENTRAL 		! indices centrais do cubo ambiente
   INTEGER(KIND=8)  , INTENT(IN)  :: TIME_ANT_SEC , TIME_POS_SEC   , TIME 

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(OUT) :: U       , V       , W   , T ,  SAL    
   INTEGER(KIND=8)  , INTENT(OUT) :: IF_LOG 

! ------------------------------------------------------- parametros locais
   INTEGER           ::  I , J , K , I1 , J1 , K1 , K2 , IT , I2 , J2 
   INTEGER           ::  IIFF , I_101 , J_101

   DOUBLE PRECISION, DIMENSION ( NX_OCEAN)   ::  X_COO 
   DOUBLE PRECISION, DIMENSION ( NY_OCEAN)   ::  Y_COO 
   DOUBLE PRECISION, DIMENSION ( NZ_OCEAN)   ::  Z_COO

   REAL, DIMENSION ( -1:1 , -1:1, NZ_OCEAN  )           ::  U_RED , V_RED , S_RED, T_RED 
   REAL, DIMENSION ( -1:0 , -1:0, NZ_OCEAN  )           ::  W_RED
   REAL, DIMENSION ( -1:1 , -1:1 )                      ::  SSH_RED

   DOUBLE PRECISION, DIMENSION (5,2)     ::  VAR_AMB
   DOUBLE PRECISION, DIMENSION ( 4 )     ::  VAR_XY
   DOUBLE PRECISION, DIMENSION ( 2 )     ::  VAR_Z
   DOUBLE PRECISION, DIMENSION ( 2 )     ::  PESO_Z , PESO_T

! --------------------------------------------------------------------------------------

! --- calcular os pesos no tempo

    PESO_T(2) =  FLOAT( TIME - TIME_ANT_SEC)/FLOAT(TIME_POS_SEC - TIME_ANT_SEC )
    PESO_T(1) = 1.D0 - PESO_T(2) 


! --- 


    OPEN ( 50, FILE= COORD_X , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  X_COO(J), J = 1 , NX_OCEAN )
    CLOSE( 50 )

    OPEN ( 50, FILE= COORD_Y , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  Y_COO(I), I = 1 , NY_OCEAN )
    CLOSE( 50)

    OPEN ( 50, FILE= COORD_Z , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )   
    READ ( 50, * , IOSTAT = IIFF ) (  Z_COO(J), J = 1 , NZ_OCEAN )
    CLOSE( 50 )

    I1 = 0
    J1 = 0

    IF(X_COO(I_CENTRAL) .GT. LONPOS ) I1 = -1
    IF(Y_COO(J_CENTRAL) .GT. LATPOS ) J1 = -1


  ! ----------------------------posicao da parcela transformada para uma grade 101X101
    I_101 =  NINT (100. * (LONPOS -X_COO(I_CENTRAL+I1)) / (X_COO(I_CENTRAL+I1+1) - X_COO(I_CENTRAL+I1) )) ! indices do peso
    J_101 =  NINT (100. * (LATPOS -Y_COO(J_CENTRAL+J1)) / (Y_COO(J_CENTRAL+J1+1) - Y_COO(J_CENTRAL+J1) )) 

    IF_LOG = 1

    IF( I_101 .LT.   0 )  GO TO 1000
    IF( I_101 .GT. 100 )  GO TO 1000
    IF( J_101 .LT.   0 )  GO TO 1000
    IF( J_101 .GT. 100 )  GO TO 1000

    IF_LOG = 0

  ! ------------------------- K1 : indice da altura acima da posicao da parcela 
    DO K = 1 , NZ_OCEAN

       IF (PROFPOS .LT.  Z_COO(K) ) GO TO 100

    ENDDO 


100 CONTINUE


    K1 = K-1
    K2 = K1+1
    IF ( K1 .EQ. 0) THEN
	K1 = 1
        PESO_Z(1) = 0.5D0
        PESO_Z(2) = 0.5D0
    ELSE
    	PESO_Z(1) = (Z_COO(K2)- PROFPOS )/ (Z_COO(K2) -Z_COO(K1))
    	PESO_Z(2) = (PROFPOS - Z_COO(K1)  )/ (Z_COO(K2) -Z_COO(K1))
    END IF



  ! ------------------------- 

    DO IT = 1 , 2 ! tempo anterior e tempo posterior do arquivo mercator

       IF ( IT .EQ. 1 ) THEN
           OPEN(50, FILE = "file_ant" , STATUS = 'UNKNOWN', FORM = 'UNFORMATTED'  ,     & 
                    ACCESS = 'DIRECT',                                                  & 
                    RECL = (4 * 3*3*NZ_OCEAN + 3*3 + 2*2*NZ_OCEAN)* NBIT      )  
       ELSE
           OPEN(50, FILE = "file_pos" , STATUS = 'UNKNOWN', FORM = 'UNFORMATTED'  ,     & 
                    ACCESS = 'DIRECT',                                                  & 
                    RECL = (4 * 3*3*NZ_OCEAN + 3*3 + 2*2*NZ_OCEAN)* NBIT      )  

       ENDIF

       READ(50,REC=1) V_RED , U_RED , S_RED , SSH_RED, T_RED , W_RED
       CLOSE(50)

  ! ------------------------- interpola nos dois niveis, e depois pra posicao Z

       CALL CORTA( V_RED(:,:,:) ,  I1 ,  J1 , K1, VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(1) ) 

       CALL CORTA( V_RED(:,:,:) ,  I1 ,  J1 , K2 , VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(2) ) 

       CALL INTERPOLA_2( VAR_Z(1) , VAR_Z(2) , PESO_Z(1) , PESO_Z(2) ,  VAR_AMB(1,IT) ) ! variavel V

     ! ----
       CALL CORTA( U_RED(:,:,:) ,  I1 ,  J1 , K1, VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(1) ) 

       CALL CORTA( U_RED(:,:,:) ,  I1 ,  J1 , K2 , VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(2) ) 

       CALL INTERPOLA_2( VAR_Z(1) , VAR_Z(2) , PESO_Z(1) , PESO_Z(2) ,  VAR_AMB(2,IT) ) ! variavel U

     ! ----
       CALL CORTA( S_RED(:,:,:) ,  I1 ,  J1 , K1, VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(1) ) 

       CALL CORTA( S_RED(:,:,:) ,  I1 ,  J1 , K2 , VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(2) ) 

       CALL INTERPOLA_2( VAR_Z(1) , VAR_Z(2) , PESO_Z(1) , PESO_Z(2) ,  VAR_AMB(3,IT) ) ! variavel S

     ! ----
       CALL CORTA( T_RED(:,:,:) ,  I1 ,  J1 , K1, VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(1) ) 

       CALL CORTA( T_RED(:,:,:) ,  I1 ,  J1 , K2 , VAR_XY(:) )
       CALL INTERPOLA( VAR_XY(:) ,I_101 , J_101 ,  VAR_Z(2) ) 

       CALL INTERPOLA_2( VAR_Z(1) , VAR_Z(2) , PESO_Z(1) , PESO_Z(2) ,  VAR_AMB(4,IT) ) ! variavel T


     ! ----

       I2 = 0 ; IF ( I1 .EQ. -1 ) I2 = -1

       J2 = 0 ; IF ( J1 .EQ. -1 ) J2 = -1

       VAR_Z(1) =   W_RED(I2,J2,K1) 
       VAR_Z(2) =   W_RED(I2,J2,K2) 

       CALL INTERPOLA_2( VAR_Z(1) , VAR_Z(2) , PESO_Z(1) , PESO_Z(2) ,  VAR_AMB(5,IT) ) ! variavel W



    ENDDO

    CALL INTERPOLA_2( VAR_AMB(1,1) , VAR_AMB(1,2) , PESO_T(1) , PESO_T(2) ,  V ) 
    CALL INTERPOLA_2( VAR_AMB(2,1) , VAR_AMB(2,2) , PESO_T(1) , PESO_T(2) ,  U ) 
    CALL INTERPOLA_2( VAR_AMB(3,1) , VAR_AMB(3,2) , PESO_T(1) , PESO_T(2) ,  SAL ) 
    CALL INTERPOLA_2( VAR_AMB(4,1) , VAR_AMB(4,2) , PESO_T(1) , PESO_T(2) ,  T   ) 
    CALL INTERPOLA_2( VAR_AMB(5,1) , VAR_AMB(5,2) , PESO_T(1) , PESO_T(2) ,  W   ) 


1000 CONTINUE


!----------------------
  END SUBROUTINE PROP_AMB_HIDTER 
!----------------------

!========================================================= 


 SUBROUTINE CORTA (VAR_IN ,  I_BELOW , J_BELOW , K1 ,VAR_XY ) 

! ---------------------------------------------------------------------------- 
!                                                              valdir 20180121 
! 
! - volta com uma matriz de 4 valores
!   - VAR_IN : matriz 3-D de quatro cubos (em cada z sao 3X3 valores)
!   - I_BELOW , J_BELOW , K1: indices da matriz VAR_IN do canto inferior esquerdo
!     e nivelK1
!  
!   - VAR_XY : variavel out, de 4 pontos tomados da matriz VAR_IN, no nivel K1
!  
! - estrategia
!     entra com uma matriz grande, indica um ponto destra matriz, e volta com apenas 
!     4 valores, no nivel K1, e nos pontos adjacentes :
!  
!                 3                        
!     J_BELOW+1  +-------------------------+ 4
!                |                         |  
!     J_BELOW    +-1-----------------------+ 2
!            I_BELOW               I_BELOW + 1
! 
!       I_BELOW      , J_BELOW     , K1 ---> VAR_XY(1)
!       I_BELOW + 1  , J_BELOW     , K1 ---> VAR_XY(2)
!       I_BELOW      , J_BELOW + 1 , K1 ---> VAR_XY(3)
!       I_BELOW + 1  , J_BELOW + 1 , K1 ---> VAR_XY(4)
! 
! ---------------------------------------------------------------------------- 

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
   REAL , INTENT(IN) , DIMENSION ( -1:1 , -1:1, NZ_OCEAN  ) ::  VAR_IN 
   INTEGER , INTENT(IN)    ::  I_BELOW , J_BELOW   , K1

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(OUT) , DIMENSION ( 4 ) ::  VAR_XY

! ------------------------------------------------------- 

      VAR_XY(1) = VAR_IN( I_BELOW    ,  J_BELOW     , K1) 
      VAR_XY(2) = VAR_IN( I_BELOW +1 ,  J_BELOW     , K1) 
      VAR_XY(3) = VAR_IN( I_BELOW    ,  J_BELOW  +1 , K1) 
      VAR_XY(4) = VAR_IN( I_BELOW +1 ,  J_BELOW  +1 , K1)


!----------------------
  END SUBROUTINE CORTA
!----------------------
!========================================================= 


   SUBROUTINE INTERPOLA(  VAR_IN , I_101 , J_101 ,  VAR ) 

! ---------------------------------------------------------------------------- 
!                                                              valdir 20180121 
! 
! - interpola a variavel
!   --------------------
!   - VAR_IN : entra com 4 variaveis, uma em cada cando de um quadrilatero
!   - I_101 , J_101 : indices de uma matriz peso que sera lida aqui
!   - VAR : variavel interpolada
! 
! - estrategia
!     entra com 4 variaveis, e o indice de uma matriz peso, calculada previamente
!     que  eh acessada de um arquivo. A matriz peso tem 101 posicoes, e em cada 
!     posicao ha 4 pesos, um peso relativo a cada vertice do retangulo.

!     A cada valor da variavel in que for indefinido ( < ou =  -999.9) o peso
!     correspondente eh zerado, e os pesos que nao forem indefinidos sao norma-
!     lizados para que a soma dos pesos seja 1.
! ---------------------------------------------------------------------------- 
!
!  
   IMPLICIT NONE

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN) , DIMENSION ( 4 ) ::  VAR_IN
   INTEGER, INTENT(IN)                             :: I_101,J_101

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(OUT) ::  VAR

! ------------------------------------------------------- parametros internos
   DOUBLE PRECISION , DIMENSION(0:100 , 0:100 , 1:4) :: PESO_XY
   DOUBLE PRECISION  :: SOMA_PESO
   INTEGER :: I

! ------------------------------------------------------- 

! ------------------------------------------------------- 

   OPEN(22, FILE= 'pesos.dat', STATUS= 'UNKNOWN', FORM= 'UNFORMATTED' , & 
         ACCESS= 'DIRECT', RECL = 2* 101*101*4 * NBIT ) 
         READ(22 , REC = 1) PESO_XY 
   CLOSE(22)

   DO I = 1 , 4
      IF ( VAR_IN(I) .LE. UNDEF ) PESO_XY(I_101,J_101,I )= 0.D0
   ENDDO


   SOMA_PESO = SUM(PESO_XY(I_101,J_101,:))

   IF( SOMA_PESO .LE. 0. ) THEN
       VAR  = UNDEF 
   ELSE
       VAR  = VAR_IN(1)*PESO_XY(I_101,J_101,1) +  &
              VAR_IN(2)*PESO_XY(I_101,J_101,2) +  &
              VAR_IN(3)*PESO_XY(I_101,J_101,3) +  & 
              VAR_IN(4)*PESO_XY(I_101,J_101,4)

       VAR = VAR / SOMA_PESO

   ENDIF

!----------------------
  END SUBROUTINE INTERPOLA
!----------------------
!====================================================================

   SUBROUTINE INTERPOLA_2( VAR_1 , VAR_2 , PESO_1 , PESO_2 ,  VAR ) 
! 
!  ------------------------------------------------------------------
!  in  : duas variaveis, dois pesos
!  out : valor interpolado 
!
!  usada para interpolar no tempo e na direcao z 
! 
!  estrategia:
!    se valor da variavel for menor que UNDEF, a variavel nao eh 
!    utilizada e os pesos sao recalculados peso nao eh utilizado
!   

!  ------------------------------------------------------------------

   IMPLICIT NONE

! ------------------------------------------------------- parametros in
   DOUBLE PRECISION , INTENT(IN)  :: VAR_1 , VAR_2 , PESO_1 , PESO_2

! ------------------------------------------------------- parametros out
   DOUBLE PRECISION , INTENT(OUT) ::  VAR

! ------------------------------------------------------- parametros internos
   DOUBLE PRECISION  :: SOMA_PESO , P_1 , P_2

! ------------------------------------------------------- 

   P_1 = PESO_1
   P_2 = PESO_2

   IF( VAR_1 .LE. UNDEF ) P_1  = 0.D0
   IF( VAR_2 .LE. UNDEF) P_2  = 0.D0
   
   SOMA_PESO =  P_1 +  P_2 

   IF( SOMA_PESO .LE. 0. ) THEN
       VAR = UNDEF
   ELSE
       VAR  = ( VAR_1*P_1  + VAR_2*P_2) /SOMA_PESO
   ENDIF


!----------------------
  END SUBROUTINE INTERPOLA_2

!========================================================= 
!========================================================= 

END MODULE module_ambiente



