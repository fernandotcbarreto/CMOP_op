Module dissolve_part

use processes
use vertical_dispersion

double precision, dimension(:), allocatable:: massdis

contains




SUBROUTINE DISSOLVE_OIL( MH2O, VIS_DIN_A, RO_A, RO_PL, DIAM_PL, Wp_OIL, VIS_DIN_PL, TS_OA, DT_STEP, NUMDROPS)


  IMPLICIT NONE
   DOUBLE PRECISION , INTENT(IN)    :: MH2O						! massa de agua em vc
   DOUBLE PRECISION :: VCRI , VMOL_EBN , COEF_WATER
   DOUBLE PRECISION :: KM_OIL , DOA							! coeficiente de transferencia de massa e difusividade do oleo na agua
   DOUBLE PRECISION , INTENT(IN)    :: VIS_DIN_A					! viscosidade dinamica do ambiente
   DOUBLE PRECISION :: CORRECAO								! correcao da solubilidade devido a salinidade de vc
   DOUBLE PRECISION :: CTE
   INTEGER :: I , J , ID_OG
   DOUBLE PRECISION , INTENT(IN)    :: RO_A, DIAM_PL, Wp_OIL						! densidade do ambiente
   DOUBLE PRECISION , INTENT(INOUT) :: RO_PL						! densidade da PL
   DOUBLE PRECISION :: RE , NSC , NSH , N_PEC						! numeros de Reynolds, Schimidt , Sherwood e Peclet
   DOUBLE PRECISION :: C1 , C2 , C3
   DOUBLE PRECISION , INTENT(IN)    :: VIS_DIN_PL
   DOUBLE PRECISION , INTENT(IN)    :: TS_OA						! tensao superficial do ambiente
  DOUBLE PRECISION , INTENT(IN) :: DT_STEP
   DOUBLE PRECISION :: MDISS_COMP,C_MDISS, FRAC_MOL			! concentracao do comp no ambiente	

   double precision :: NUMDROPS		     

!- calcular a difusividade do oleo na agua -> DOA (m2/s)
	 ! --- calculo de VMOL_EBN pelo Tyn and Calus Method (Reid et al 1987)
	  VCRI = VCRI_O * 1.D6 * PM_OIL_out *1.D-3		! em cm3/mol
	  VMOL_EBN = 0.285D0 * VCRI**1.048		! em cm3/mol


	  ! correlacao de Wilke,Chang(1955)
          COEF_WATER = 2.6D0					! 2.6 Wilke,Chang (1955)
	  DOA = 7.4D-8*T_A*(COEF_WATER*PM_H2O)**(0.5D0)/((VIS_DIN_A*1.D3)*(VMOL_EBN)**0.6D0)	! em cm2/s
	  DOA = DOA * 1.D-4									! em m2/s

!   print*, T_A, COEF_WATER, PM_H2O, VIS_DIN_A, VMOL_EBN
!   print*, 7.4D-8

          CORRECAO = 1.D0
	  CTE = 1.D0
	  ID_OG = 1							! identificar se eh oleo (1) ou gas (2)		
          IF(SAL_A .NE. 0.D0) CALL CORRECT_SAL (ID_OG , PROF_REF , PROFPOS , T_A , SAL_A , RO_A , &
 CTE , CTE , CTE , CTE,  RO_PL ,PM_OIL_out , CORRECAO)



   
!       print*, 'VIS_DIN_A',  WP_OIL, DIAM_PL, RO_A,VIS_DIN_A, DOA

!      print*, WP_OIL * DIAM_PL * RO_A 
!      print*, (WP_OIL * DIAM_PL * RO_A )/VIS_DIN_A

      RE = WP_OIL * DIAM_PL * RO_A /  VIS_DIN_A        		! numero de Reynolds

      NSC = VIS_DIN_A / ( DOA * RO_A)					! numero de Schimidt

      N_PEC = DIAM_PL * WP_OIL/DOA					! numero de peclet

      C1 = 50.D0 + 2.D0*(PI)**(-0.5D0)*N_PEC**(0.5D0)
      C2 = 5.26D-2*RE**(1.D0/3.D0 + 6.59D-2*RE**(1.D0/4.D0))*NSC**(1.D0/3.D0)*(WP_OIL*VIS_DIN_A/TS_OA)&
 **(1.D0/3.D0)*(1.D0/(1.D0 + (VIS_DIN_PL/VIS_DIN_A)**(1.1D0)))
      C3 = 2.D0 + 0.363D0*RE**(0.484D0)*NSC**(0.339D0)*(DIAM_PL*GRAVITY**(1.D0/3.D0)/DOA**(2.D0/3.D0))**(0.072D0)

      NSH = (C3 + C1*C2)/(1.D0 + C2)		


      KM_OIL = NSH * DOA / DIAM_PL					! coeficiente de trasnferencia de massa do oleo (m/s)

    MDISS_COMP=0




do comps=1, NCOMP_OIL

      C_MDISS   =  MDISS_COMP* RO_A * 1.D-3 / PM_COMP_OIL(comps)         		! concentracao de MDISS_COMP em kg/m**3

      FRAC_MOL = masscomp(jj,ii-1, comps) * PM_OIL_out / (massa(jj,ii-1) * PM_COMP_OIL(comps) )     	! fracao molar do componente na gota  

      
      massdis(comps) = KM_OIL*4.D0*PI*(DIAM_PL/2.D0)**(2.D0)*(FRAC_MOL*CORRECAO*SOLU_COMP_OIL(comps) - C_MDISS)*DT_STEP  	! massa de componente 

 !      massdis(comps) =  massdis(comps) * NUMDROPS

       massdis(comps) =  massdis(comps)   !back to particle reference

!      print*, 'aaaaaaaa', DIAM_PL


     MASSCOMP(jj,ii, comps)= MASSCOMP(jj,ii-1, comps) - massdis(comps)
     
     if ( MASSCOMP(jj,ii, comps) .le. 0 ) then
         MASSCOMP(jj,ii, comps) = 0
     endif


enddo  

   
 !     print*,'uuu', km_oil, nsh, doa, diam_pl
 !     print*, c1,c2,c3
 !     print*, WP_OIL
 !     print*, re, nsc, n_pec

 ! print*, 'massdiss', massdis(1), masscomp(jj,ii-1, 1), masscomp(jj,ii,1)

 
!   print*, km_oil, diam_pl, frac_mol
 
end subroutine  DISSOLVE_OIL







SUBROUTINE DISSOLVE_OIL_FASE2( MH2O, VIS_DIN_A, RO_A, RO_PL, DIAM_PL, Wp_OIL, VIS_DIN_PL, TS_OA, DT_STEP, NUMDROPS)


  IMPLICIT NONE
   DOUBLE PRECISION , INTENT(IN)    :: MH2O						! massa de agua em vc
   DOUBLE PRECISION :: VCRI , VMOL_EBN , COEF_WATER
   DOUBLE PRECISION :: KM_OIL , DOA							! coeficiente de transferencia de massa e difusividade do oleo na agua
   DOUBLE PRECISION , INTENT(IN)    :: VIS_DIN_A					! viscosidade dinamica do ambiente
   DOUBLE PRECISION :: CORRECAO								! correcao da solubilidade devido a salinidade de vc
   DOUBLE PRECISION :: CTE
   INTEGER :: I , J , ID_OG
   DOUBLE PRECISION , INTENT(IN)    :: RO_A, DIAM_PL, Wp_OIL						! densidade do ambiente
   DOUBLE PRECISION , INTENT(INOUT) :: RO_PL						! densidade da PL
   DOUBLE PRECISION :: RE , NSC , NSH , N_PEC						! numeros de Reynolds, Schimidt , Sherwood e Peclet
   DOUBLE PRECISION :: C1 , C2 , C3
   DOUBLE PRECISION , INTENT(IN)    :: VIS_DIN_PL
   DOUBLE PRECISION , INTENT(IN)    :: TS_OA						! tensao superficial do ambiente
  DOUBLE PRECISION , INTENT(IN) :: DT_STEP
   DOUBLE PRECISION :: MDISS_COMP,C_MDISS, FRAC_MOL			! concentracao do comp no ambiente	

   double precision :: NUMDROPS		     

!- calcular a difusividade do oleo na agua -> DOA (m2/s)
	 ! --- calculo de VMOL_EBN pelo Tyn and Calus Method (Reid et al 1987)
	  VCRI = VCRI_O * 1.D6 * PM_OIL_out *1.D-3		! em cm3/mol
	  VMOL_EBN = 0.285D0 * VCRI**1.048		! em cm3/mol


	  ! correlacao de Wilke,Chang(1955)
          COEF_WATER = 2.6D0					! 2.6 Wilke,Chang (1955)
	  DOA = 7.4D-8*T_A*(COEF_WATER*PM_H2O)**(0.5D0)/((VIS_DIN_A*1.D3)*(VMOL_EBN)**0.6D0)	! em cm2/s
	  DOA = DOA * 1.D-4									! em m2/s

!   print*, T_A, COEF_WATER, PM_H2O, VIS_DIN_A, VMOL_EBN
!   print*, 7.4D-8

          CORRECAO = 1.D0
	  CTE = 1.D0
	  ID_OG = 1							! identificar se eh oleo (1) ou gas (2)		
          IF(SAL_A .NE. 0.D0) CALL CORRECT_SAL (ID_OG , PROF_REF , PROFPOS , T_A , SAL_A , RO_A , &
 CTE , CTE , CTE , CTE,  RO_PL ,PM_OIL_out , CORRECAO)



   
!       print*, 'VIS_DIN_A',  WP_OIL, DIAM_PL, RO_A,VIS_DIN_A, DOA

!      print*, WP_OIL * DIAM_PL * RO_A 
!      print*, (WP_OIL * DIAM_PL * RO_A )/VIS_DIN_A

      RE = WP_OIL * DIAM_PL * RO_A /  VIS_DIN_A        		! numero de Reynolds

      NSC = VIS_DIN_A / ( DOA * RO_A)					! numero de Schimidt

      N_PEC = DIAM_PL * WP_OIL/DOA					! numero de peclet

      C1 = 50.D0 + 2.D0*(PI)**(-0.5D0)*N_PEC**(0.5D0)
      C2 = 5.26D-2*RE**(1.D0/3.D0 + 6.59D-2*RE**(1.D0/4.D0))*NSC**(1.D0/3.D0)*(WP_OIL*VIS_DIN_A/TS_OA)&
 **(1.D0/3.D0)*(1.D0/(1.D0 + (VIS_DIN_PL/VIS_DIN_A)**(1.1D0)))
      C3 = 2.D0 + 0.363D0*RE**(0.484D0)*NSC**(0.339D0)*(DIAM_PL*GRAVITY**(1.D0/3.D0)/DOA**(2.D0/3.D0))**(0.072D0)

      NSH = (C3 + C1*C2)/(1.D0 + C2)		


      KM_OIL = NSH * DOA / DIAM_PL					! coeficiente de trasnferencia de massa do oleo (m/s)

    MDISS_COMP=0




do comps=1, NCOMP_OIL

      C_MDISS   =  MDISS_COMP* RO_A * 1.D-3 / PM_COMP_OIL(comps)         		! concentracao de MDISS_COMP em kg/m**3

      FRAC_MOL = masscomp(jj,ii-1, comps) * PM_OIL_out / (massa(jj,ii-1) * PM_COMP_OIL(comps) )     	! fracao molar do componente na gota  

      
      massdis(comps) = KM_OIL*4.D0*PI*(DIAM_PL/2.D0)**(2.D0)*(FRAC_MOL*CORRECAO*SOLU_COMP_OIL(comps) - C_MDISS)*DT_STEP  	! massa de componente 

      massdis(comps) =  massdis(comps) * NUMDROPS


   !  print*, 'aaaaaaaa1',  VIS_DIN_PL,VIS_DIN_PL, VIS_DIN_A, VIS_DIN_PL/VIS_DIN_A


     MASSCOMP(jj,ii, comps)= MASSCOMP(jj,ii-1, comps) - massdis(comps)
     
     if ( MASSCOMP(jj,ii, comps) .le. 0 ) then
         MASSCOMP(jj,ii, comps) = 0
     endif


enddo  

   
 !     print*,'uuu', km_oil, nsh, doa, diam_pl
 !     print*, c1,c2,c3
 !     print*, WP_OIL
 !     print*, re, nsc, n_pec

 ! print*, 'massdiss', massdis(1), masscomp(jj,ii-1, 1), masscomp(jj,ii,1)

 
!   print*, km_oil, diam_pl, frac_mol
 
end subroutine  DISSOLVE_OIL_FASE2




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





end module
