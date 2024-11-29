program lagrange
USE ISO_FORTRAN_ENV, ONLY : REAL128
use white_space
use oil_fractions
use random_variables
use evaporation
use processes
use dissolve_part
use vertical_dispersion
use delft2oil
use era5_2_oil
use wind_variables
use dissolved_fase_mod
use coupling
use linear_interpolation
use datetime_module
!use iso_fortran_env, only: real64

  type (datetime) :: timenum


  double precision:: al, kl, spm, b(2), c(100,100), areat, diamt, ui(1,1), vi(1,1), ti(1,1), si(1,1),mwf, lo,wvp, rrr
  double precision:: scft, kf, rn1, rn2, voldis, odir, ovel, counttimeh_r0, probsl, rn10, ductwd
  integer i,a,j,cont, tsl, itsl, freq, ts_lim, cont_ts, outer_l, zlayer,  avel, reverse, depvi
  integer :: ts, numpalm, turn_off_mask, year11, month11, day11, hour11, minute11, secon11, deg_turn
  character:: g(100)
  double precision:: dt,vf, vft, x_model(100,100), y_model(100,100), u_model(100,100), x_mask(100,100), v_model(100,100)  !vf eh vel func
  double precision:: y_mask(100,100), mask(100,100), maskt(1,1)
  double precision, dimension(:, :), allocatable:: x, area, y, u,v, checkb
  double precision, dimension(:,:), allocatable :: diam, temps, salts
  double precision, dimension(:), allocatable:: dt_h_spr
  double precision :: time_lim, time_lim_wind, outp_h, watcontant
  character*20::  fmt1, x1, fmt, fmt2, fmt3, x2, fmt5, x5, fmt6, x3, fmt4, fmt7, fmt8, fmt9

!!!!!! RESHAPE
  double precision, dimension(:), allocatable:: res_ar
!!!!!

  double precision:: pvp, c21, k, k2, az, as

  real :: start_time, stop_time

  integer indexlon, indexlat

  NCOMP_OIL=25

  open(12,file='x_part.txt', status='UNKNOWN')
  open(13,file='y_part.txt', status='UNKNOWN')
  open(14,file='massa.txt', status='UNKNOWN')
  open(15,file='diam.txt', status='UNKNOWN')
  open(16,file='height.txt', status='UNKNOWN')
  open(199,file='vol.txt', status='UNKNOWN')
!  open(198,file='temp.txt', status='UNKNOWN')
!  open(197,file='salt.txt', status='UNKNOWN')
  open(198,file='u_wind.txt', status='UNKNOWN')
  open(197,file='v_wind.txt', status='UNKNOWN')

  open(11111,file='coastvalues.txt', status='UNKNOWN')

  open(141,file='massaf2.txt', status='UNKNOWN', FORM= 'FORMATTED ')

  OPEN(17 , FILE = "visc_e_CONT.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )  

  OPEN(18 , FILE = "time_hours2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' ) 

  OPEN(19 , FILE = "spmt.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
 
  OPEN(23 , FILE = "porc_evap_cont.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(231 , FILE = "mass_evap2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(232 , FILE = "mass_disp2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(233 , FILE = "mass_sedi.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(268 , FILE = "mass_deg.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(54 , FILE = "water_frac.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(221 , FILE = "lat_part_op.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(222 , FILE = "lon_part_op.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(223 , FILE = "lat_partf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(224 , FILE = "lon_partf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(227 , FILE = "lat_partf3.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(228 , FILE = "lon_partf3.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
 
  OPEN(225 , FILE = "zf1_192.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(2254 , FILE = "beached.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  
  OPEN(331 , FILE = "diam_drop.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  
  OPEN(332 , FILE = "spmtf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(226 , FILE = "zf3.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(400 , FILE = "prob_func_assu.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(402 , FILE = "num_part_assutxt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(500 , FILE = "probabilistic_marinha_9_315_3_10.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(6789 , FILE = "concentration_assu.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(6667 , FILE = "time_prob_assu.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(9967 , FILE = "lat_model_mod.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(9968 , FILE = "lon_model_mod.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
  OPEN(9969 , FILE = "depth_mod.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )


!  OPEN(20 , FILE = "u_model.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

!  OPEN(30 , FILE = "v_model.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )



!  OPEN(21 , FILE = "x_model.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

!  OPEN(22 , FILE = "y_model.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

 
!  OPEN(25 , FILE = "x_mask.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

!  OPEN(26 , FILE = "y_mask.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

!  OPEN(27 , FILE = "mask.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

 ! call init_delft

  zlayer = 1

  !coupling_ind=1

  !right_random = 0

  three_dim = 0

  PROBABILISTIC  = 0   !!!TURN ON in probabilistic analysis

  probabilistic_2= 0 

!  wind_theoretical = 1
   probsl=0.2

  linear_interp  = 1
 
  turn_off_mask = 0

  !api=36.5523D0  !!!!!!!!!!!!API Johansen Deep blow
  !API=36.7229556  !diesel
  !api=30.5410933707076
  !api=45.6367254222577
 ! api=45.6747897   !oscar jordbaer 2010 13
  !api=33.02352      !oscar mandalay batelle 20


!   api= 30.5410933707076
   
!  api= 45.6367254222577 !!!!!!!!!!!!API Johansen Deep blow

!   api  = 19.60d0
   
  scf=2.7  !handbook of oil spill scienc and technology

!  widfc=0.035
!  widfc=0.0

!  windspx=15800  !m/h 15800 is 15.8 km/h
!  windspy=0
!  windsp=(windspx**2 + windspy**2)**0.5
!  windspms=windsp/3600


!  !print*, windsp
  !kem=99
  !print*, trim(path)
  !stop
  open(6877,file='input.dat',  status='old')
  read(6877,*) year11, month11, day11, hour11, minute11, secon11
  read(6877,*) PATH
  read(6877,*) ts
  read(6877,*) numpal 
  read(6877,*) outp_h  
  read(6877,*) lon_ref
  read(6877,*) lat_ref
  read(6877,*) dt
  read(6877,*) voldis
  read(6877,*) DISSOLVE 
  read(6877,*) EMULSI 
  read(6877,*) EVAP_TURN 
  read(6877,*) ENTRAIN 
  read(6877,*) THEORETICAL 
  read(6877,*) DISSOLVED_FASE 
  read(6877,*) uwd 
  read(6877,*) vwd 
  read(6877,*) avel
  read(6877,*) odir 
  read(6877,*) ovel
  read(6877,*) wind_theoretical
  read(6877,*) reverse
  read(6877,*) coupling_ind
  read(6877,*) time_finish
  read(6877,*) right_random
  read(6877,*) path_ini_coup
  read(6877,*) path_frac
  read(6877,*) depvi
  read(6877,*) widfc
  read(6877,*) CDIF_HOR
  read(6877,*) apist
  read(6877,*) ductwd
  read(6877,*) watcontant
!  read(6877,*) deg_turn
  
 !print*, trim(apist)
 if (trim(apist) .eq. 'Diesel') then
   API=36.7229556  !diesel
 endif
 if (trim(apist) .eq. 'Jordbaer') then
   api=45.6747897   !oscar jordbaer 2010 13
 endif
 if (trim(apist) .eq. 'Mandalay') then
   api=33.02352      !oscar mandalay batelle 20
 endif  
 if (trim(apist) .eq. 'Linerle') then
   api=27.55361      !oscar Linerle
 endif  
 if (trim(apist) .eq. 'IF-30-bunker') then
   api=25.84299      !oscar Linerle
 endif  ! print*, api
 if (trim(apist) .eq. 'HAVIS_2012') then
   api=37.40500      !oscar Linerle
 endif  ! print*, api
 if (trim(apist) .eq. 'FOINAVEN_(IKU)') then
   api=25.88186      !oscar Linerle
 endif  ! print*, api    
 if (trim(apist) .eq. 'KUWAIT_2002') then
   api=36.09498      !oscar Linerle   
 endif  ! print*, api
 if (trim(apist) .eq. 'Emergency') then
   api=29.97929      !oscar Linerle   
 endif  ! print*, api  
 if (trim(apist) .eq. 'Oseberg') then
   api=33.05403      !oscar Linerle   
 endif  ! print*, api   
 if (trim(apist) .eq. 'Camorim') then
   api=30.43027      !oscar Linerle   
 endif  ! print*, api    
 if (trim(apist) .eq. 'Balder') then
   api=28.16920      !oscar Linerle   
 endif  ! print*, api   
 if (trim(apist) .eq. 'GASOIL') then
   api=37.98144      !oscar Linerle   
 endif  ! print*, api   
 if (trim(apist) .eq. 'VIGDIS') then
   api=33.9184      !oscar Linerle   
 endif  ! print*, api
 if (trim(apist) .eq. 'EMBLA_2000') then
   api=36.0713      !oscar Linerle   
 endif  ! print*, api 
 if (trim(apist) .eq. 'Njord_2002') then
   api=38.3566      !oscar Linerle   
 endif  ! print*, api  
 if (trim(apist) .eq. 'If_180_SHELL') then
   api=20.9615      !oscar Linerle   
 endif  ! print*, api   
 if (trim(apist) .eq. 'Arab_light_Battelle') then
   api=35.6941      !oscar Linerle   
 endif  ! print*, api    
 if (trim(apist) .eq. 'Gyda_2000') then
   api=33.3001      !oscar Linerle   
 endif  ! print*, api     !   stop
 if (trim(apist) .eq. 'Hago_4ss_LA_IKL') then
   api=28.6312      !oscar Linerle   
 endif  ! print*, api     !   stop 
 if (trim(apist) .eq. 'russian_crude') then
   api=31.6583      !oscar Linerle   
 endif  ! print*, api     !   stop  
  !print*, widfc, CDIF_HOR
  !print*, numparcels_dis
  !stop
!  print*, year11, month11, day11, hour11, minute11, secon11
!  print*, coupling_ind
 ! timenum = datetime(2020, 2, 28, 23, 59, 25)
 ! timenum = datetime(year11, month11, day11, hour11, minute11, secon11)
 ! print *, date2num(timenum)  - date2num(datetime(1970, 1, 1)) ! 734869.25000000000
 ! stop
  !print*, path_frac
  !stop
  if (coupling_ind.eq.1) then
	time_finish=0
  endif
    
    
  
  ucomp2=0
  vcomp2=0 
  !vel2=1
  !dir2=220
  !call vel_conv
  if (avel .eq. 1) then
    vel2=ovel
    dir2=odir
    call vel_conv   
  endif 
!  print*, avel, dir2, vel2, ucomp2,vcomp2 
  !print*, numpal, wind_theoretical, voldis, lon_ref,lat_ref
!  stop
 
  if (wind_theoretical.eq.1) then
   windsp = ((uwd**2 + vwd**2)**0.5)*3600   !wind velocity in m/h
   windspms = (uwd**2 + vwd**2)**0.5
  endif

!  print*, numpal, uwd, vwd, theoretical, windsp
!  stop
  
!  numpal=5!number of particles
  numparcels=numpal
  mwf= 933.45 !g/mol
!  dt=900! TS always in second
!  dt=3600! TS always in second

  dt_random = 1

!  ts = 3600 !number of hours you wanna run fernando noronh
!  ts = 24 !number of hours you wanna run

 ! time_finish = 10  ! number of hours that lasts the continuous release 10000000

  ts= (ts*60*60)/dt 

  tsl= dt/60. !minutes interval for realesing particle 240000000
  
  itsl=tsl

  tsl = (tsl*60) / dt

  freq=ts/tsl

  ts_lim=500
!  ts_lim=100

  if (time_finish .eq. 0) then
		freq=1
		!print*, freq
  else 
		freq=(time_finish*60*60)/dt
!		print*, freq
  endif 
!  stop
!   print*,numpal, ts - freq*tsl


!   print*, ts, tsl, itsl
!   stop

   if (coupling_ind .eq. 1) then
      go to 2534
   endif


   if (PROBABILISTIC.eq.0) then
          
 !    if ( (ts - freq*tsl) .eq.0) then
    if ((ts - freq * tsl) .eq. 0 .or. time_finish .eq. 0) then
	 
        !print*, '00000'
        allocate(x(int(numpal*freq ),int(ts_lim)), y(int(numpal*freq),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal*freq 
        write (x3,fmt) numparcels_dis
     else 
        !print*, 'not000'
        allocate(x(int(numpal*freq + numpal ),int(ts_lim)), y(int(numpal*freq + numpal ),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal*freq + numpal 
        write (x3,fmt) numparcels_dis
     endif
   else 

        allocate(x(int(numpal),int(ts_lim)), y(int(numpal),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal 
        write (x3,fmt) numparcels_dis
   endif  


  2534 continue



!print*, numpal, freq, tsl, ts
!stop
!  zf1 = 0
  contind=2

  call init_coupl

  if (coupling_ind.eq.1) then

    NCLAS_OIL = 8			! numero de classes de gota de oleo (!DEEPBLOW )
    call init_coupl
    open(1564, file=trim(path_ini_coup), status='old')
    open(1112, file=trim(path_frac), status='old')
    numb_lines_c=0
    DO
     READ(1564,*,iostat=io)
     IF (io/=0) EXIT
     numb_lines_c = numb_lines_c + 1
    END DO
    CLOSE (1564)
    allocate (time_res(numb_lines_c), lat_ref_res(numb_lines_c), lon_ref_res(numb_lines_c), vol_res(numb_lines_c), &
     res_par_in (numb_lines_c), zfc(numb_lines_c))

!    !print*, numb_lines_c

    open(1564, file=trim(path_ini_coup), status='old')

    DO i = 1, numb_lines_c
     READ(1564,*) res_par_in(i), time_res(i), lon_ref_res(i), lat_ref_res(i), zfc(i),  vol_res(i)
    END DO
    CLOSE (1564)    

!    print*, numpal, numb_lines_c
!	stop

    allocate (FRAC_MASS_OUT(NCOMP_OIL))
    allocate(x(int(numpal*numb_lines_c ),int(ts_lim)), y(int(numpal*numb_lines_c),int(ts_lim)))
    fmt= '(I10)'
    write (x1,fmt) numpal*numb_lines_c
    write (x3,fmt) numparcels_dis

    lon_ref = lon_ref_res(1)
    lat_ref = lat_ref_res(1)

!!print*, lon_ref, lat_ref

!    lon_ref = -90.78335556   ! golfo
!    lat_ref = 26.44927444    ! golfo

!    lon_ref = 4.83333    ! deep_spill
!    lat_ref = 65    ! deep_spill


  ALLOCATE ( V_COMP(NCOMP_OIL)    ) 
  ALLOCATE ( MOL_COMP(NCOMP_OIL)  )
  ALLOCATE ( FRAC_MASS_OUT_REF(NCOMP_OIL) )
  ALLOCATE ( PM_COMP_OIL  (NCOMP_OIL) )
  ALLOCATE ( SOLU_COMP_OIL(NCOMP_OIL) )
  ALLOCATE ( RO_COMP_OIL_15(NCOMP_OIL))
  ALLOCATE ( RO_COMP_OIL  (NCOMP_OIL) )
  ALLOCATE ( TEB_COMP_OIL (NCOMP_OIL) )     
  ALLOCATE ( CP_COMP_OIL  (NCOMP_OIL) )

  ALLOCATE ( SGCOMPS  (NCOMP_OIL) )
  ALLOCATE ( TGCOMPS  (NCOMP_OIL) )
  ALLOCATE ( PCRI_O_COMPS  (NCOMP_OIL) )
  ALLOCATE ( PC1  (NCOMP_OIL) )
  ALLOCATE ( PC2  (NCOMP_OIL) )
  ALLOCATE ( PC  (NCOMP_OIL) )
 
  ALLOCATE ( TCRI_O_COMPS  (NCOMP_OIL) )
  ALLOCATE ( TB_COMPS  (NCOMP_OIL) )

  ALLOCATE ( RED_TEMP (NCOMP_OIL) )
      ALLOCATE ( DIAM_OUT_OIL(NCLAS_OIL)  )
      ALLOCATE ( VOLU_OUT_OIL(NCLAS_OIL)  ) 

!   !print*, lat_ref, lon_ref, time_res

  endif


!!print*, shape(x)

 ! !print*, fmt1
 
!!  if ((ts-freq*tsl).eq.0) then
!      allocate(x(int(numpal*freq)+numpal,int(ts)), y(int(numpal*freq)+numpal,int(ts)))
!       fmt= '(I10)'
!!      write (x1,fmt) numpal*freq+pal
!      !print*, trim(x1)     
!        do i=1,dt, dt_random
!            !print*, i


!        enddo
!      fmt1='( '//trim(x1)//'f10.3)'    !"( 1000f10.3 )"
!      fmt1=trim(fmt1)
!      call StripSpaces(fmt1)  
!!      fmt2='( '//trim(x1)//'f10.6)'    !"( 1000f10.3 )"
!!      fmt2=trim(fmt2)
!      call StripSpaces(fmt2)
!  else
 
!    allocate(x(int(numpal*freq),int(ts)), y(int(numpal*freq),int(ts)))
  
!     fmt= '(I10)'

!    write (x1,fmt) numpal*freq
!    write (x2,fmt) numparcels

!   !print*, trim(x1)     
    fmt1='( '//trim(x1)//'f10.3)'    !"( 1000f10.3 )"
    fmt1=trim(fmt1)
    call StripSpaces(fmt1) 

    fmt2='( '//trim(x1)//'f16.6)'    !"( 1000f10.3 )"
 !   fmt2='( '//trim(x1)//'f10.3)'    !"( 1000f10.3 )"
    fmt2=trim(fmt2)
    call StripSpaces(fmt2)

    fmt3='( '//trim(x2)//'f13.6)'    !"( 1000f10.3 )"   !!!!!!!!!!!!!!!CHANGE IN THE CONTINUOS
    fmt3=trim(fmt3)
    call StripSpaces(fmt3)

    fmt4='( '//trim(x3)//'f13.6)'    !"( 1000f10.3 )"   !!!!!!!!!!!!!!!CHANGE IN THE CONTINUOS
    fmt4=trim(fmt4)
    call StripSpaces(fmt4)

 ! endif
!  allocate(x(int(numpal*freq),int(ts)), y(int(numpal*freq),int(ts)))

  b=shape(x)

 
  allocate( diam(int(b(1)),int(b(2))), massa(int(b(1)),int(b(2))), beached(int(b(1)),int(b(2))))

  allocate( vol(int(b(1)),int(b(2))), area(int(b(1)),int(b(2))), u(int(b(1)),int(b(2))), v(int(b(1)),int(b(2))) )
   
  
  allocate( temps(int(b(1)),int(b(2))), salts(int(b(1)),int(b(2))), checkb(int(b(1)),int(b(2))))
  

  allocate (massa2(int(b(1)),int(b(2))), massa1(int(b(1)),int(b(2))), height(int(b(1)),int(b(2))))

  allocate (massa3(int(b(1)),int(b(2))), massa4(int(b(1)),int(b(2))))


  allocate (massa5(int(b(1)),int(b(2))), massa6(int(b(1)),int(b(2))), massa7(int(b(1)),int(b(2))), massa8(int(b(1)),int(b(2))))

  allocate (massa9(int(b(1)),int(b(2))), massa10(int(b(1)),int(b(2))), massa11(int(b(1)),int(b(2))), massa12(int(b(1)),int(b(2))))

  allocate (massa13(int(b(1)),int(b(2))), massa14(int(b(1)),int(b(2))), massa15(int(b(1)),int(b(2))), massa16(int(b(1)),int(b(2))))

  allocate (massa17(int(b(1)),int(b(2))), massa18(int(b(1)),int(b(2))), massa19(int(b(1)),int(b(2))), massa20(int(b(1)),int(b(2))))

  allocate (massa21(int(b(1)),int(b(2))), massa22(int(b(1)),int(b(2))), massa23(int(b(1)),int(b(2))), massa24(int(b(1)),int(b(2))))

  allocate (massa25(int(b(1)),int(b(2))))
  allocate (massae(int(b(1)),int(b(2))), vole(int(b(1)),int(b(2))), watcont(int(b(1)),int(b(2))))

  allocate (visc_e(int(b(1)),int(b(2))), rho_e(int(b(1)),int(b(2))))

  allocate (mas_evap(int(b(1)),int(b(2))) , porc_evap(int(b(1)),int(b(2))), porc_evap_vol(int(b(1)),int(b(2))))
  
  allocate (vol_evap(int(b(1)),int(b(2))), vol_diss(int(b(1)),int(b(2))), porc_diss_vol(int(b(1)),int(b(2))))

  allocate (mass_diss(int(b(1)),int(b(2))), mass_sedi(int(b(1)),int(b(2))), mass_degr(int(b(1)),int(b(2))))

  allocate (area2(int(b(1)),int(b(2))), area1(int(b(1)),int(b(2))))
  

  allocate (area3(int(b(1)),int(b(2))), area4(int(b(1)),int(b(2))))


  allocate (area5(int(b(1)),int(b(2))), area6(int(b(1)),int(b(2))), area7(int(b(1)),int(b(2))), area8(int(b(1)),int(b(2))))

  allocate (area9(int(b(1)),int(b(2))), area10(int(b(1)),int(b(2))), area11(int(b(1)),int(b(2))), area12(int(b(1)),int(b(2))))

  allocate (area13(int(b(1)),int(b(2))), area14(int(b(1)),int(b(2))), area15(int(b(1)),int(b(2))), area16(int(b(1)),int(b(2))))

  allocate (area17(int(b(1)),int(b(2))), area18(int(b(1)),int(b(2))), area19(int(b(1)),int(b(2))), area20(int(b(1)),int(b(2))))

  allocate (area21(int(b(1)),int(b(2))), area22(int(b(1)),int(b(2))), area23(int(b(1)),int(b(2))), area24(int(b(1)),int(b(2))))

  allocate (area25(int(b(1)),int(b(2))))



  allocate (diam2(int(b(1)),int(b(2))), diam1(int(b(1)),int(b(2))))
  

  allocate (diam3(int(b(1)),int(b(2))), diam4(int(b(1)),int(b(2))))


  allocate (diam5(int(b(1)),int(b(2))), diam6(int(b(1)),int(b(2))), diam7(int(b(1)),int(b(2))), diam8(int(b(1)),int(b(2))))

  allocate (diam9(int(b(1)),int(b(2))), diam10(int(b(1)),int(b(2))), diam11(int(b(1)),int(b(2))), diam12(int(b(1)),int(b(2))))

  allocate (diam13(int(b(1)),int(b(2))), diam14(int(b(1)),int(b(2))), diam15(int(b(1)),int(b(2))), diam16(int(b(1)),int(b(2))))

  allocate (diam17(int(b(1)),int(b(2))), diam18(int(b(1)),int(b(2))), diam19(int(b(1)),int(b(2))), diam20(int(b(1)),int(b(2))))

  allocate (diam21(int(b(1)),int(b(2))), diam22(int(b(1)),int(b(2))), diam23(int(b(1)),int(b(2))), diam24(int(b(1)),int(b(2))))

  allocate (diam25(int(b(1)),int(b(2))))



allocate(MASSCOMP(int(b(1)),int(b(2)),NCOMP_OIL ))
allocate(MASSCOMPREF(NCOMP_OIL))


allocate( XAM(int(b(1)) ), XWM(int(b(1))), dt_h_spr(int(b(1))) )

allocate(evapmass(int(b(1)),NCOMP_OIL), porep(NCOMP_OIL), massdis(NCOMP_OIL))


allocate( FRAC_MASS_OUT_PART(int(b(1)),NCOMP_OIL )   )


allocate( watf ( int(b(1)),int(b(2)) ) )

allocate (dropdiam( int(b(1)),int(b(2)) ), numdrops( int(b(1)),int(b(2)) ), voldrops( int(b(1)),int(b(2)) ))

allocate(lat_part(int(b(1)),int(b(2))), lon_part(int(b(1)),int(b(2))),  zf1(int(b(1)),int(b(2))))

lat_part(:,:)=-999
lon_part(:,:)=-999

allocate(vec_tdelf(int(ts_lim)), vec_era(int(ts_lim)))

allocate(areaem(int(b(1)),int(b(2))), diamem(int(b(1)),int(b(2))), areadrop(int(b(1)),int(b(2))))


! if (ENTRAIN.EQ.1) THEN 
   allocate(xf2(int(numparcels),int(b(2))), yf2(int(numparcels),int(b(2))) , zf2(int(numparcels),int(b(2)))    )

   allocate(lon_partf2(int(numparcels),int(b(2))),lat_partf2(int(numparcels),int(b(2)))  )


   allocate(areaf2(int(numparcels),int(b(2))),volf2(int(numparcels),int(b(2))) , dropdiamf2(int(numparcels),int(b(2)))    )


   allocate( numdropsf2(int(numparcels),int(b(2))),voldropsf2(int(numparcels),int(b(2))) ,diamf2(int(numparcels),int(b(2)))    )


   allocate(MASSCOMPf2(int(numparcels),int(b(2)),NCOMP_OIL ), massaf2(int(numparcels),int(b(2))),spmtf2(int(numparcels),int(b(2))) )

   allocate(massadropsF2(int(numparcels),int(b(2))) , masscompdropF2(int(numparcels),int(b(2)),NCOMP_OIL )      )

! ENDIF


! if (dissolved_fase .eq. 1) then

   allocate(xf3(int(numparcels_dis),int(b(2))), yf3(int(numparcels_dis),int(b(2))) , zf3(int(numparcels_dis),int(b(2)))    )

   allocate(lon_partf3(int(numparcels_dis),int(b(2))),lat_partf3(int(numparcels_dis),int(b(2))), &
 massaf3(int(numparcels_dis),int(b(2)))  )
 
! endif



!allocate(tsevol(int(b(2)) *int(b(1)) ))
allocate(tsevol(int(b(2)) *int(b(1)) ), tsevol2(int(b(2)) *int(b(1)) ) )


allocate(lat_part_sum(int(b(1))),  lon_part_sum(int(b(1))), coord_sum_f2(int(b(1))) )





! call init_oil(porc1, porc2, scp1, scp2,mwf1, mwf2)

  call rand_var(vapp, r, temp)

  call celstokel(temp)
 
  call evap_vars

!!print*, r

!sto


!!print*, 'porep', porep(20)


!PRINT*, SHAPE(MASSCOMP)

!print*, b(1)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

 !print*, temp

   TEMP_OUT=temp    !temp at discharge



  if (coupling_ind .ne. 1) then
    zf1(:,:)=depvi
	if (depvi .lt. 0) then
	checkb(:,:)=19
	endif
   !print*, b(1), 'eee'
   !stop
!  spm=830  !kg/m3

  !   vol(1:NUMPAL,:)=voldis/numpal  !volume of each particle

     vol(1:NUMPAL,:)=voldis/b(1)  !volume of each particle

   ! print*, voldis/b(1), voldis, b(1), NUMPAL
!	stop
!   TEMP_OUT  = (273.15D0 + TEMP_OUT)   !The first calculations are performemed in Kelvin


!   CALL components( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
!	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT )   !----------------- 


!stop
! API=36.7229556  !diesel


     CALL components( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT )

 !stop
! api=88899
! print*, 'first', api, temp_out, spmt
! call calc_api(api)
! stop

     do i=1, numpal
       FRAC_MASS_OUT_PART(i,:)=FRAC_MASS_OUT
     enddo


   endif





   call emulsify_parameter

   call entrainment_param
   
   watf(:,:)=watcontant
!   print*, watf
!   stop




!!!!!!INITIALIZE DELFT INFO
 open(145,file=trim(path)//'time_delft.txt',  status='old')
 read(145,*) numtime, numlat, numlon, numz
! !print*, numtime, numlat, numlon, numz
 allocate(time_vec(numtime), lon_model(numlat, numlon)   ,   lat_model(numlat, numlon)) 
 allocate(u_model1(numlat, numlon, numz), u_model2 (numlat, numlon, numz), v_model1(numlat, numlon, numz),&
  v_model2(numlat, numlon, numz), t_model1(numlat, numlon, numz), t_model2 (numlat, numlon, numz), &
  s_model1(numlat, numlon, numz), s_model2 (numlat, numlon, numz)) 
  allocate(kz_model1(numlat, numlon, numz), kz_model2(numlat, numlon, numz))
! allocate(depth(numlat, numlon), lat_model_sum(numlat, numlon), lon_model_sum(numlat, numlon), &
 allocate(lat_model_sum(numlat, numlon), lon_model_sum(numlat, numlon), &
 coord_sum(numlat, numlon), levelsd(int(numz)))
! allocate(area_grid(numlat, numlon))
 allocate(deplevel(int(numz))) 
 allocate(w_model1(numlat, numlon, numz), w_model2(numlat, numlon, numz))
 
 do i=1,numtime
   read(145,*) time_vec(i) 
 enddo

 do i=1,numz
   read(145,*) levelsd(i) 
 enddo
  
 read(145,*)
 read(145,*) numlatm, numlonm
 allocate(depth(numlatm, numlonm), lat_modelm(numlatm, numlonm), lon_modelm(numlatm, numlonm))
 allocate(coastvalue(numlatm, numlonm), coastvalueb(numlatm, numlonm))
 allocate(lat_model_summ(numlatm, numlonm), lon_model_summ(numlatm, numlonm), &
 coord_summ(numlatm, numlonm)) 
! print*, numlatm, numlonm
! stop
 
 
  timenum = datetime(year11, month11, day11, hour11, minute11, secon11)
  counttimeh_r = (date2num(timenum)  - date2num(datetime(1970, 1, 1)))* (24*60)
 ! print*, counttimeh_r
!  stop
  if (counttimeh_r .lt. time_vec(1)) then
     print*, 'ERROR:INIT TIME BEFORE FORCING TIME', counttimeh_r, time_vec(1)
	 go to 8856
   endif

 time_lim = time_vec(numtime)*60 - itsl*60 - dt


counttimeh_r0=counttimeh_r
!print*, counttimeh_r
!!print*, time_vec(numtime) 
!!print*, time_vec
!  !print*, 'teste1',  time_lim, time_vec(numtime), numtime

 
 close(145)
 
 
  open(145,file=trim(path)//'lat_delft.txt',  status='old')
  open(146,file=trim(path)//'lon_delft.txt',  status='old')
  open(14555,file=trim(path)//'lat_dmask.txt',  status='old')
  open(14666,file=trim(path)//'lon_dmask.txt',  status='old')  
  open(111,file=trim(path)//'depth.txt',  status='old')
  open(1111,file=trim(path)//'coastvalue.txt',  status='old')

!  open(3766,file=trim(path)//'area_grid.txt',  status='old')
!  open(7535,file=trim(path)//'winter.dry',  status='old')
  !print*, shape(depth), numlon,numlat
 !!print*, numlat,numlon
 do i=1,numlat
     read (145,*) (lat_model(i,j), j=1,numlon)
     read (146,*) (lon_model(i,j), j=1,numlon) 
!     read (111,*) (depth(i,j), j=1,numlon) 
!     read (3766,*) (area_grid(i,j), j=1,numlon) 
 enddo
 do i=1,numlatm
     read (14555,*) (lat_modelm(i,j), j=1,numlonm)
     read (14666,*) (lon_modelm(i,j), j=1,numlonm)
     read (111,*) (depth(i,j), j=1,numlonm) 
     read (1111,*) (coastvalue(i,j), j=1,numlonm) 
 enddo
 close(145)
 close(146)
 close(111)
 close(14555)
 close(14666)
 
 coastvalueb=coastvalue
!print*, maxval(coastvalueb)
!stop
 !print*, lat_model
!  do i=1,110
!    read(7535,*) indexlon,indexlat
!      depth(indexlat, indexlon) = 0
!!      lat_model(indexlat, indexlon) = 0
!!      lon_model(indexlat, indexlon) = 0
! !   !print*, indexlon, indexlat, lat_model(indexlat, indexlon)
!  enddo

!!!!!!!!
!do i=1,numlat
!   write(9967, '(208f16.6)') (lat_model(i,:))
!   write(9968, '(208f16.6)') (lon_model(i,:))
!   write(9969, '(208f16.6)') (depth(i,:))
!!print*, 'hfhfhfh', -( 10 -  4 ) * 65 * 596 /180.D0
!!print*, 'hfhfhfh', -( 10 -  4 ) * (65/180.D0) * 596 
!!print*, pi
!enddo
!!print*, shape(lat_model), fmt2
! close(9967)
! close(9968)
!!!!


allocate(partcontmap(numlat, numlon  ) )

allocate(probmap(numlat, numlon  ) )

allocate(partcontmap_2(numlat, numlon  ) )

allocate(probmap_2(numlat, numlon  ) )

allocate(mass_dist(numlat, numlon  ) )

allocate(conc(numlat, numlon  ) )

allocate(time_sum(numlat, numlon  ) )

allocate(prob_time_sum(numlat, numlon  ) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INITIALIZE ERA 5 DATA
  if (wind_theoretical.ne.1) then

    open(687,file=trim(path)//'time_era_5.txt',  status='old')
    read(687,*) num_time_era, numlat_era, numlon_era
    allocate (time_vec_era(num_time_era), u10_era1(numlat_era, numlon_era), v10_era1(numlat_era, numlon_era) )
    allocate (u10_era2(numlat_era, numlon_era), v10_era2(numlat_era, numlon_era))
    allocate (lat_era(numlat_era, numlon_era), lon_era(numlat_era, numlon_era), lat_era_sum(numlat_era, numlon_era))
    allocate (lon_era_sum(numlat_era, numlon_era), coord_sum_era(numlat_era, numlon_era))
  
    do i=1,num_time_era
      read(687,*) time_vec_era(i) 
    enddo

 !   read(687,*) counttimeh_era
    close(687)

    open(687,file=trim(path)//'lat_era_5.txt',  status='old')
    open(787,file=trim(path)//'lon_era_5.txt',  status='old')
 
    do i=1,numlat_era
       read (687,*) (lat_era(i,j), j=1,numlon_era)
       read (787,*) (lon_era(i,j), j=1,numlon_era) 
    enddo

    close(687)
    close(787)

    time_lim_wind = time_vec_era(num_time_era)*60 - itsl*60 - dt
  endif


!  !print*, time_lim_wind, time_lim
!!print*, 'nickyyy'
! !print*, time_vec_era
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fmt= '(I10)'
    write (x5,fmt) numlon

  fmt5='( '//trim(x5)//'f15.5)'    !"( 1000f10.3 )"
  fmt5=trim(fmt5)
  call StripSpaces(fmt5)

  fmt6='( '//trim(x5)//'i7)'    !"( 1000f10.3 )"
  fmt6=trim(fmt6)
  call StripSpaces(fmt6)
  
      fmt= '(I10)'
    write (x5,fmt) numlonm

  fmt9='( '//trim(x5)//'f15.5)'    !"( 1000f10.3 )"
  fmt9=trim(fmt9)
  call StripSpaces(fmt9)

! !print*,  depth(50, :)

!!print*, lat_model(20,:)
  sal_a=salin                    ! salinity at discharge
  
 ! print*, 'sal1', salin
!  stop

  !print*, temp_out,sal_a

  call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   height_value = 0.0000005

   

  if (coupling_ind .ne. 1) then

    massref = massa(1,1)
    volreff = massa(1,1)/spmt
    MASSCOMPREF = MASSCOMP(1,1,:)

!print*, 'second', massref, volreff,  MASSCOMPREF


    aux=0
    do comps=1, NCOMP_OIL
       aux=aux + MASSCOMPREF(comps)/RO_COMP_OIL(comps)
    enddo

    spmtref=massref/aux

    viscst = (VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000

!!!!!! definition of droplet info
     dropdiamax =  cmax * ((viscst) ** (0.34) ) *  ( turbed**(-0.4)  ) !!! reed et al (1995)

     dropdiamin =  cmin * ((viscst) ** (0.34) ) *  ( turbed**(-0.4)  )

     dropdiam(1:NUMPAL,:) = ( ((dropdiamax - dropdiamin)/2.)   + dropdiamin ) * 0.000001

 !  do i = 1, int(b(1))
 !    CALL random_number(RN4)                      
   
 !    dropdiam(i,:) = ( ( (dropdiamax - dropdiamin) * RN4 ) + dropdiamin ) * 0.000001                !!! drop size distrubution reed et al(1995)

       areadrop(:,:) = pi*((dropdiam(:,:)/2)**2)
  
       voldrops(:,:) = ((dropdiam(:,:)/2.)**3.) *  pi * (4./3.)

       numdrops(:,:) =  (massa(:,:)/spmt) / voldrops (:,:)

!       area(:,:)=pi*((diam(:,:)/2)**2)

       time_ini_spread = ((1.45/1.14)**4.) * ( (  vol(NUMPAL,1) / ((VIS_DIN_A/RO_A) * gravity * ((RO_A &
- RO_OIL_OUT )/RO_A ) ) )**(1./3.)  )


       time_ini_spread = time_ini_spread/60

       area(1:NUMPAL,:) = 2.1*(pi)*  (   ( (( vol(1:NUMPAL,:)**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * &
((time_ini_spread*60)**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT))**(1./2.)  )**(1./3.)  ) 
 

        dt_h_spr(:) = time_ini_spread

!       areaem(:,:)=pi*((diamem(:,:)/2)**2)

       areaem(:,:) = area(1:NUMPAL,:)
 
 !      diam(:,:)=2*(((massa(:,:)/spmt)/(pi*height(:,:)))**0.5)  !meter

       diam(1:NUMPAL,:) = ((area(1:NUMPAL,:) / PI ) ** (1./2.)) * 2.


!       diamem(:,:)=2*(((massa(:,:)/spmt)/(pi*height(:,:)))**0.5)  !meter


       diamem(:,:) = diam(:,:)

       height(1:NUMPAL,:) = (vol(1:NUMPAL,:)) / area(1:NUMPAL,:)  !delft3d


       visc_e(1:NUMPAL,:) = VIS_DIN_OIL_OUT*1000   !in cP

       rho_e(1:NUMPAL,:) = RO_OIL_OUT

 !   !print*, areadrop(i,1) *numdrops (i,1), area(i,1), vol(1,1), massa(1,1)/spmt
  

 ! height(:,:)=0.01 !m

      deltad = ( (dropdiamax - dropdiamin)/10. ) * 0.000001

      call evap_mass


  endif


 ! diam(:,:)=2*(((massa(1,1)/spmt)/(pi*height(1,1)))**0.5)  !meter

!!!!!!!!!!!!!!!!!!!!!!!!1

!  diam(:,:)=2*(((vol(1,1)*3.)/(4.*pi))**(1./3.))   !meter
 ! area(:,:)=pi*((diam(1,1)/2)**2)

 ! areat=area(1,1)*numpal
 ! diamt=2*((areat/pi)**0.5)





!  !print*, diamt
!  !print*, diam(1,1)
!  !print*, diamt/diam(1,1)
  

!  !print*, diam
  
!  c(:,:)=(/((/(i, i=1,10)), j=1,100)/)  !metro
!  c=((/ ((i,i=1,10), j=1,100) /), (/ 10, 100 /))


 ! u(:,:)=0.05 !1m/s



  vf=0.1  ! root-mean-square value for the stochastic turbulent fluctuations of the horizontal velocity



!  X0=20
!  Y0=20


!  do i=1,numpal
!    CALL random_number(RN1)
!    CALL random_number(RN2)
!    x(i,1)=X0*RN1
!    y(i,1)=Y0*RN2
!  enddo

  X0=0
  Y0=0




  ZPOS=0
  PROFOUT=0
  opt=2
  circle_radius = 0   !radius of duct, following deltares

  if (coupling_ind .ne. 1) then

!  lon_ref=-40.278
!  lat_ref=-20.3168    !vitoria port


!  lon_ref=-40.985
!  lat_ref=-21.800   !acu port


!  lon_ref=-40.2422
!  lat_ref=-20.2986     !tubarao port

!  lon_ref=-31.5
!  lat_ref=-9     !Nordeste spill


!  lon_ref=-38.9516
!  lat_ref=-2.4761     !Nordeste spill north

!  lon_ref=-46
!  lat_ref=-27.8377     !Nordeste spill north

!   lon_ref=-42
!   lat_ref=-23.7

!  if (probabilistic.eq.1) then
       circle_radius = ((1.45**2.)/1.14) * ((   ( (voldis**5.)*gravity*((RO_A - RO_OIL_OUT)/RO_A)   &               
                   ) /((VIS_DIN_A/RO_A))**2.  )**(1./12.))   
!  endif
!  do i=1,numpal
!    CALL random_number(RN1)
!    CALL random_number(RN2)    !!OLD SQUARE SLICK
!    x(i,1)=X0*RN1
!	print*, circle_radius
!    y(i,1)=Y0*RN2
 !   circle_radius=0
!  enddo

  do i=1,numpal
    CALL random_number(RN1)
    CALL random_number(RN2)
    x(i,1)=(circle_radius*RN1) * COS(2*PI*RN2) +   X0             !!NEW CIRCULAR SLICK
    y(i,1)=(circle_radius*RN1) * SIN(2*PI*RN2) +   Y0

    call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(i,1), LAT_part(i,1), PROFOUT,  x(i,1),  y(i,1), ZPOS, OPT)

 !   !print*, x(i,1), y(i,1), lon_part(i,1), lat_part(i,1), opt
    
  enddo 
 

 endif
 ! !print*, x(10,1)
 ! !print*, (x_model(1,:)-x(10,1))
 ! !print*, minloc(abs(x_model(1,:)-x(10,1)))
 ! !print*, x_model(1,minloc(abs(x_model(1,:)-x(10,1))))

!!!!!!!!!!!!!!!!!!!!!!Variable height

  allocate(x_height(num_sp*2-1,num_sp*2-1), y_height(num_sp*2-1,num_sp*2-1 ), lat_height(num_sp*2-1,num_sp*2-1 ), &
lon_height(num_sp*2-1,num_sp*2-1 ), height_map(num_sp*2-1,num_sp*2-1 ))
  allocate(lat_model_sum_height(num_sp*2-1,num_sp*2-1), lon_model_sum_height(num_sp*2-1,num_sp*2-1), &
coord_sum_height(num_sp*2-1,num_sp*2-1))

 do i=num_sp+1, num_sp*2-1
    x_height(:,i) = x_height(i-1,i-1) + dx_h
    y_height(i,:) = y_height(i-1,i-1) + dy_h
 enddo
 do i=num_sp-1,1, -1
    x_height(:,i) = x_height(i+1,i+1) - dx_h
    y_height(i,:) = y_height(i+1,i+1) - dy_h
 enddo

 do i=1, num_sp*2-1
    do j=1, num_sp*2-1
      call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_height(i,j), LAT_height(i,j), PROFOUT,  x_height(i,j), &
  y_height(i,j), ZPOS, OPT)
    enddo
 enddo

  open(765,file='x_map.txt', status='UNKNOWN')
  open(766,file='y_map.txt', status='UNKNOWN')
  open(767,file='lon_map.txt', status='UNKNOWN')
  open(768,file='lat_map.txt', status='UNKNOWN')


    fmt= '(I10)'
    write (x5,fmt) num_sp*2-1

  fmt7='( '//trim(x5)//'f15.9)'    !"( 1000f10.3 )"
  fmt7=trim(fmt7)
  call StripSpaces(fmt7)

  fmt8='( '//trim(x5)//'f15.5)'    !"( 1000f10.3 )"
  fmt8=trim(fmt8)
  call StripSpaces(fmt8)


 do i=1, num_sp*2-1
   write(765, fmt8) x_height(i,:)
   write(766, fmt8) y_height(i,:)
   write(767, fmt8) lon_height(i,:)
   write(768, fmt8) lat_height(i,:)
 enddo

!!print*, lat_height(:,239)

!!!!!!!!!!!!!!!!!!!!!!Variable height

! lon_ref = 4.850
! lat_ref = 65.0000
! OPT=1
! LONOUT = 4.870
!  LATOUT = 65.0000
! call INI_AMB( lon_ref, lat_ref , PROF_REF, lonout, latout, PROFOUT,  X0,  Y0, ZPOS, OPT)
! !print*,  lon_ref, lat_ref, PROF_REF, lonout, latout, PROFOUT,  X0,  Y0, ZPOS, OPT
 ! !print*, 'pp'
 ! !print*, y(10,1)
 ! !print*, abs(y_model(:,1) - y(10,1))
 !   !print*, y_model(minloc(abs(y_model(:,1) - y(10,1))), 1)
!  !print*, y_model(1,:)


 do i=1, int(b(1))
    if (i.eq.int(b(1))) then
    write(12,'(i10)') i
    write(13,'(i10)') i
    write(14,'(i10)') i
    write(15,'(i10)') i     
    write(199,'(i10)') i     
    write(198,'(i10)') i     
    write(197,'(i10)') i     
    write(16,'(i10)') i 
    write(17,'(i10)') i 
    write(18,'(i10)') i
    write(19,'(i10)') i
    write(23,'(i10)') i
    write(231,'(i10)') i
    write(232,'(i10)') i
    write(233,'(i10)') i
    write(268,'(i10)') i
    write(221,'(i10)') i
    write(222,'(i10)') i
    write(54,'(i10)') i
	write(2254,'(i10)') i
    write(225,'(i10)') i

     else 
    write(12,'(i10)', advance='no') i
    write(13,'(i10)', advance='no') i
    write(14,'(i10)', advance='no') i
    write(15,'(i10)', advance='no') i
    write(199,'(i10)', advance='no') i
    write(198,'(i10)', advance='no') i
    write(197,'(i10)', advance='no') i
    write(16,'(i10)', advance='no') i
    write(17,'(i10)', advance='no') i
    write(18,'(i10)', advance='no') i
    write(19,'(i10)', advance='no') i
    write(23,'(i10)', advance='no') i
    write(231,'(i10)', advance='no') i
    write(232,'(i10)', advance='no') i
    write(233,'(i10)', advance='no') i
    write(268,'(i10)', advance='no') i
    write(221,'(i10)', advance='no') i
    write(222,'(i10)', advance='no') i
    write(54,'(i10)', advance='no') i
	write(2254,'(i10)', advance='no') i
    write(225,'(i10)', advance='no') i
    endif 
 enddo

 if(dissolved_fase .eq. 1) then
  do i = 1, int(numparcels_dis)

    if (i.eq.int(numparcels_dis)) then

       write(226,'(i10)') i
       write(227,'(i10)') i
       write(228,'(i10)') i

    else 
       write(226,'(i10)', advance='no') i
       write(227,'(i10)', advance='no') i
       write(228,'(i10)', advance='no') i

    endif

  enddo

 endif


 IF (ENTRAIN.EQ.1) THEN
  do i=1, int(b(1))
    if (i.eq.int(b(1))) then

       write(223,'(i10)') i
       write(224,'(i10)') i
       write(141,'(i10)') i
       write(331,'(i10)') i
       write(332,'(i10)') i
       
    else 
       write(223,'(i10)', advance='no') i
       write(224,'(i10)', advance='no') i
       write(141,'(i10)', advance='no') i
       write(331,'(i10)', advance='no') i
       write(332,'(i10)', advance='no') i
      
    endif

  enddo
 ENDIF
 
  write(12,fmt2) x(:,1)
  write(13,fmt2) y(:,1)
  write(14,fmt2) massa(:,1)
  write(231,fmt2) mas_evap(:,1)
  write(232,fmt2) mass_diss(:,1)
  write(233,fmt2) mass_sedi(:,1)
  write(268,fmt2) mass_degr(:,1)
  write(15,fmt2) diam(:,1)
  write(199,fmt2) vol(:,1)
  write(198,fmt2) temps(:,1)
  write(197,fmt2) salts(:,1)
  write(16,fmt2) height(:,1)
  write(17,fmt2) visc_e(:,1)
  write(221,fmt2) lat_part(:,1)
  write(222,fmt2) lon_part(:,1)
  write(54,fmt2) watf(:,1)
  write(18,fmt2) counttimeh_r/60. - counttimeh_r0/60.
  write(2254,fmt2) beached(:,1)
  write(225,fmt2) zf1(:,1)
  
!print*, lat_partf2(:,1)

     if (dissolved_fase .eq. 1) then
       write(227,fmt4) lat_partf3(:,1)
       write(228,fmt4) lon_partf3(:,1)
       write(226,fmt4) zf3(:,1)
     endif
	 
 IF (ENTRAIN.EQ.1) THEN
   write(223,fmt2) lat_partf2(:,1)
   write(224,fmt2) lon_partf2(:,1) 
   write(331,fmt2) dropdiam(:,1)
   write(332,fmt2) spmtf2(:,1)
   write(141,fmt2) massaf2(:,1)
 ENDIF


!print*, 'tnh',  counttimeh_r, counttimeh_era
!print*, shape(lat_part)
counttimeh_era=counttimeh_r

!print*, counttimeh_r, counttimeh_era
!print*, reverse
!stop
 if (reverse .ne. 1) then
  counttimeh = counttimeh_r + dt/60.
  counttimeh_era = counttimeh_era +dt/60
  print*, 'reverse deactivate'
  
 else 
  counttimeh = counttimeh_r - dt/60.
  counttimeh_era = counttimeh_era - dt/60
    print*, 'reverse activated'
endif


counttimeh_coup = dt/60.



counttime = dt


     lo=0

! numpalm=numpal

 if (coupling_ind .ne. 1) then
 
   NUMTOT = NUMPAL

 endif

 cont=1

tsevol2(:) = 0


! !print*, 'ts', numtot, ts
! !print*, shape(x)
  
 cont_ts = 1

 call cpu_time(start_time)

!dt_h=1

dt_h = (dt/(60.*60.))

dt_hf2 = limen

count_prob_2=1

do outer_l=1,1000000000

 !print*, 'EXECUTION FASE'
  do i=2,ts_lim
  print*, 'TS', cont_ts

   if (PROBABILISTIC.eq.1 .and. cont.eq.tsl) then 
     go to 468
   endif

     cont_ts = cont_ts + 1


     if (coupling_ind .eq. 1) then

       num_res_par = numtot   !I wrote before the loop to overcome the problem with droplets resurfacing

       do tc_in = inf_time, numb_lines_c
         if (counttimeh_coup .ge. time_res(tc_in)) then

            read (1112,*) dummy_val, (FRAC_MASS_OUT(frac_in), frac_in=1,NCOMP_OIL), xa, xw
! 		print*, xa, xw
!			print*, FRAC_MASS_OUT
!			print*, sum(FRAC_MASS_OUT)
!			stop

!        xa=0.4

!        xw=0.4  !fix

            numtot=numtot+numpal

!            vol(num_res_par+1:numtot, :)=vol_res(tc_in)/b(1)   ! tested and wrong
            vol(num_res_par+1:numtot, :)=vol_res(tc_in)/numpal
            
			zf1(num_res_par+1:numtot, :)=zfc(tc_in)
			
		    if (zfc(tc_in) .lt. 0) then
	              checkb(num_res_par+1:numtot,:)=19
	        endif
!  !print*, 'vol', vol(:,1), counttimeh
!print*, vol_res, b(1), numpal
!stop

            do coup_ct= num_res_par+1,numtot
              FRAC_MASS_OUT_PART(coup_ct,:)=FRAC_MASS_OUT
            enddo

!  !print*, 'frac', FRAC_MASS_OUT_PART(3,:)

            CALL COMPONENTS_COUPLING(VAZAO_OIL_OUT , DT , TEMP_OUT , &
 	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT, numtot, num_res_par)


            call evap_mass_coupling(numtot, num_res_par)

            call emulsify_parameter_coupling(numtot, num_res_par)


            call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

            OPT = 1
            call INI_AMB( LON_REF, LAT_REF, PROF_REF, lon_ref_res(tc_in), lat_ref_res(tc_in), PROFOUT,  X0,  Y0, ZPOS, OPT)

            OPT = 2

            circle_radius = ((1.45**2.)/1.14) * ((   ( ((vol_res(tc_in))**5.)*gravity*((RO_A - RO_OIL_OUT)/RO_A)   &               
                   ) /((VIS_DIN_A/RO_A))**2.  )**(1./12.))   

            !print*, circle_radius
			!stop
			!circle_radius=0
			
            do coup_ct= num_res_par+1,numtot
               CALL random_number(RN1)
               CALL random_number(RN2)
               x(coup_ct,i-1)=(circle_radius*RN1) * COS(2*PI*RN2) +   X0             !!NEW CIRCULAR SLICK
               y(coup_ct,i-1)=(circle_radius*RN1) * SIN(2*PI*RN2) +   Y0

              call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(coup_ct,i-1), LAT_part(coup_ct,i-1), PROFOUT, &
                   x(coup_ct,i-1),  y(coup_ct,i-1), ZPOS, OPT)
    
            enddo 


            viscst = (VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000

            dropdiamax =  cmax * ((viscst) ** (0.34) ) *  ( turbed**(-0.4)  ) !!! reed et al (1995)

            dropdiamin =  cmin * ((viscst) ** (0.34) ) *  ( turbed**(-0.4)  )

            dropdiam(num_res_par+1:numtot,:) = ( ((dropdiamax - dropdiamin)/2.)   + dropdiamin ) * 0.000001

            areadrop(num_res_par+1:numtot,:) = pi*((dropdiam(num_res_par+1:numtot,:)/2)**2)
         
            voldrops(num_res_par+1:numtot,:) = ((dropdiam(num_res_par+1:numtot,:)/2.)**3.) *  pi * (4./3.)

            numdrops(num_res_par+1:numtot,:) =  (massa(num_res_par+1:numtot,:)/RO_OIL_OUT) / voldrops (num_res_par+1:numtot,:)

!            diam(num_res_par+1:numtot,:)=2*(((massa(num_res_par+1:numtot,:)/RO_OIL_OUT)/(pi*height(num_res_par+1:numtot,:)))**0.5)  !meter

!            diamem(num_res_par+1:numtot,:)=2*(((massa(num_res_par+1:numtot,:)/RO_OIL_OUT)/(pi*height(num_res_par+1:numtot,:)))**0.5)  !meter

!            area(num_res_par+1:numtot,:)=pi*((diam(num_res_par+1:numtot,:)/2)**2)

!            areaem(num_res_par+1:numtot,:)=pi*((diamem(num_res_par+1:numtot,:)/2)**2)

!            height(num_res_par+1:numtot,:) = height_value  !delft3d

            time_ini_spread = ((1.45/1.14)**4.) * ( (  vol(numtot,1) / ((VIS_DIN_A/RO_A) * gravity * ((RO_A - &
 RO_OIL_OUT )/RO_A ) ) )**(1./3.)  )

            time_ini_spread = time_ini_spread/60

            area(num_res_par+1:numtot,:) = 2.1*(pi)*  (   ( ((vol(num_res_par+1:numtot, :)**2.)*gravity*((RO_A - &
 RO_OIL_OUT)/RO_A) * ((time_ini_spread*60)**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )

            dt_h_spr(num_res_par+1:numtot) = time_ini_spread

!            !print*, '333', vol(numtot, 1), massa(num_res_par+1:numtot,1)/RO_OIL_OUT

            areaem(num_res_par+1:numtot,:) = area(num_res_par+1:numtot,:)


            diam(num_res_par+1:numtot,:) =  ((area(num_res_par+1:numtot,:) / PI ) ** (1./2.)) * 2.

            diamem(num_res_par+1:numtot,:) = diam(num_res_par+1:numtot,:)


            height(num_res_par+1:numtot,:) = vol(num_res_par+1:numtot, :) / area(num_res_par+1:numtot,:) 

            visc_e(num_res_par+1:numtot,:) = VIS_DIN_OIL_OUT*1000   ! in cP
  
            rho_e(num_res_par+1:numtot,:) = RO_OIL_OUT


!             !print*, '3223333', diam(1,1)
           deltad = ( (dropdiamax - dropdiamin)/10. ) * 0.000001

!        !print*, 'area', area(:,1)
!           !print*, x(numtot,1), y(numtot,1), lon_part(numtot,1), lat_part(numtot,1)
!           !print*,  opt, LON_REF, LAT_REF
!           !print*, x0,y0, circle_radius

!           if (inf_time.eq.5) then

!               
!           endif
            num_res_par = numtot   
            inf_time = inf_time + 1
            start_index = 1

         endif

       enddo
     
       if (start_index .ne. 1) then          
          go to 1547
       endif
     endif

 !    !print*, 'TS', cont_ts

  
     do j=1,NUMTOT

       call random_seed()
       CALL random_number(RN1)
       CALL random_number(RN2)

!       !print*, x(j,i-1), y(j,i-1)
!       !print*, x_model(1,:)

!       !print*, y_model(:,1)
!       !print*, x_model(1,minloc(abs(x_model(1,:)-x(j,i-1))))


        if (( height(j,i-1) .gt. height_value)) then


          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )

              aux=0
              do comps=1, NCOMP_OIL
                aux=aux + masscomp(j,i-1,comps)/RO_COMP_OIL(comps)
              enddo

               spmt=massa(j,i-1)/aux

              do comps=1, NCOMP_OIL

                 MOL_COMP(comps)=masscomp(j,i-1, comps)* 1000.D0/PM_COMP_OIL(comps)

                 FRAC_MASS_OUT(comps)=masscomp(j,i-1, comps)/massa(j,i-1)

                 V_COMP(comps) = masscomp(j,i-1, comps)/ RO_COMP_OIL(comps)

              enddo


             call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


            MOL_TOT=sum(mol_comp)

            V_TOT=massa(j,i-1)/spmt
 
           VOL(j,i)=massa(j,i-1)/spmt

           moil=massa(j,i-1)

           FRAC_TOT=sum(FRAC_MASS_OUT)

          CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT,TS_A, TS_VC )


          if (EMULSI.eq.1) then

             areaem(j,i) = 2.1*(pi)*  (   ( ((  (massae(j,i-1)/rho_e(j,i-1))**2.)*gravity*((RO_A - &
 rho_e(j,i-1))/RO_A) * ((dt_h_spr(j)*60)**(3./2.)))/ (( (visc_e(j,i-1)/1000) / rho_e(j,i-1) )**(1./2.))  )**(1./3.)  )


            height(j,i) = (massae(j,i-1)/rho_e(j,i-1)) / areaem(j,i)  !delft3d

            area(j,i) = VOL(j,i-1) / height(j,i)
			

          else


             area(j,i) = 2.1*(pi)*  (   ( ((VOL(j,i-1)**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * ((dt_h_spr(j)* &
 60)**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )

             height(j,i) = VOL(j,i-1) / area(j,i)  !delft3d

          endif
  

  
           area(j,i-1) = area(j,i)

           diam(j,i) = ((area(j,i) / PI ) ** (1./2.)) * 2.

           diam(j,i-1) = diam(j,i)

           height(j,i-1) = height(j,i)

!     !print*, '1', height(j,i), j

          else
        
            area(j,i) = area(j,i-1)

            diam(j,i) = diam(j,i-1)

            height(j,i) = height(j,i-1)
 
        endif
!        !print*, y_model(minloc(abs(y_model(:,1) - y(j,i-1))), 1)

  tsevol2(contind) = i
  tsevol(contind) = i



  !print*, 'PARTICLE', J, 'TS', cont_ts, outer_l, i

  if (  (beached(j,i-1) .eq.1) .and. (dissolved_fase .eq. 0) ) then

    go to 77958
  endif
 
  IF (THEORETICAL.EQ.1) THEN
       ui=u_model(minloc(abs(y_model(:,1)-y(j,i-1))), minloc(abs(x_model(1,:)-x(j,i-1))))


       vi=v_model(minloc(abs(y_model(:,1)-y(j,i-1))), minloc(abs(x_model(1,:)-x(j,i-1))))

       ui=0.5
       vi=0.5
       kz=0.001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DELFT
  ELSE


    if (tsevol2(contind) .ne. tsevol2(contind-1)) then                   ! Condition to read this part only in the first particle

       tdelf1=minloc(abs(time_vec-counttimeh))

!print*, counttimeh, time_vec(tdelf1(1))

       if (tdelf1(1).eq.numtime) then
          tdelf2(1)=tdelf1(1)
       elseif (counttimeh.ge.time_vec(tdelf1(1)) ) then
          tdelf2(1)=tdelf1(1)+1
       elseif(counttimeh.lt.time_vec(tdelf1(1))) then
          tdelf2(1) = tdelf1(1)
          tdelf1(1) = tdelf1(1) - 1
       endif

       if (counttimeh.gt.time_vec(numtime)) then
         pre_end='end'
		               massa(:,i) =  massa(:,i-1)
              massae(:,i) = massae(:,i-1)
              masscomp(:,i,:) = masscomp(:,i-1,:)
              diam(:,i) = diam(:,i-1)
              diamem(:,i) = diamem(:,i-1)
              area(:,i) = area(:,i-1)
              areaem(:,i) = areaem(:,i-1)
              height(:,i) = height(:,i-1)
              vol(:,i)    = vol(:,i-1)  
	   vol_diss(:,i)=vol_diss(:,i-1)
	   porc_evap(:,i)=porc_evap(:,i-1)
	   mas_evap(:,i)=mas_evap(:,i-1)
	   mass_diss(:,i)=mass_diss(:,i-1)
	   rho_e(:,i)=rho_e(:,i-1)
	   visc_e(:,i)=visc_e(:,i-1)
	   watf(:,i)=watf(:,i-1)
	   zf1(:,i)=zf1(:,i-1)
	    lon_part(:,i)=lon_part(:,i-1)
		lat_part(:,i)=lat_part(:,i-1)
         go to 3000
       endif

       vec_tdelf(i) = tdelf1(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERA 5
      if (wind_theoretical.ne.1) then
         tera1 = minloc(abs(time_vec_era-counttimeh_era))

         if (tera1(1).eq.num_time_era) then
            tera2(1)=tera1(1)
         elseif (counttimeh_era.ge.time_vec_era(tera1(1)) ) then
            tera2(1)=tera1(1)+1
         elseif (counttimeh_era.lt.time_vec_era(tera1(1)) ) then
            tera2(1) = tera1(1)
            tera1(1) = tera1(1) - 1
         endif

         vec_era(i)= tera1(1) 

         write(bas_era,  '(I20)') tera1
         write(bas2_era, '(i20)') tera2
       
      endif
!      !print*, 'a', tera1, tera2, time_vec_era(tera1(1)), vec_era(i), counttimeh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       write(bas,  '(i20)') tdelf1
       write(bas2, '(i20)') tdelf2




!       !print*, 'y',trim(ADJUSTL(bas_era)), trim(ADJUSTL(bas2_era))
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! !print*, trim(path)//trim(ADJUSTL(bas))//'u.txt'
  !!print*, trim(path)//trim(ADJUSTL(bas2))//'u.txt'
! PRINT*, VEC_TDELF(I), VEC_TDELF(I-1), 'jgjgj'
  !print*, 'BAS', tdelf1, tdelf2, counttimeh
 ! print*, trim(path)//trim(ADJUSTL(bas))//'u.txt'
 ! print*, trim(path)//trim(ADJUSTL(bas2))//'u.txt'
 ! print*, vec_tdelf(i), vec_tdelf(i-1), counttimeh, time_vec(tdelf1(1))
 !print*, vec_tdelf(i), vec_tdelf(i-1)
 !stop

       if(vec_tdelf(i).ne. vec_tdelf(i-1)) then

      if(tdelf1(1).ne.tdelf2(1)) then

         open(145,file=trim(path)//trim(ADJUSTL(bas))//'u.txt',  status='old')
         open(146,file=trim(path)//trim(ADJUSTL(bas2))//'u.txt',  status='old')

         open(147,file=trim(path)//trim(ADJUSTL(bas))//'v.txt',  status='old')
         open(148,file=trim(path)//trim(ADJUSTL(bas2))//'v.txt',  status='old')
         
         open(151,file=trim(path)//trim(ADJUSTL(bas))//'temp.txt',  status='old')
         open(152,file=trim(path)//trim(ADJUSTL(bas2))//'temp.txt',  status='old')
         
         open(153,file=trim(path)//trim(ADJUSTL(bas))//'salt.txt',  status='old')
         open(154,file=trim(path)//trim(ADJUSTL(bas2))//'salt.txt',  status='old')

         if (three_dim.eq.1) then
           open(149,file=trim(path)//trim(ADJUSTL(bas))//'w.txt',  status='old')
           open(150,file=trim(path)//trim(ADJUSTL(bas2))//'w.txt',  status='old')
  
           open(160,file=trim(path)//trim(ADJUSTL(bas))//'kz.txt',  status='old')
           open(161,file=trim(path)//trim(ADJUSTL(bas2))//'kz.txt',  status='old')
         endif 
 

         do zl=1,numz
           do lati=1,numlat
             read (145,*) (u_model1(lati,lonj, zl), lonj=1,numlon)
             read (146,*) (u_model2(lati,lonj, zl), lonj=1,numlon) 
             read (147,*) (v_model1(lati,lonj, zl), lonj=1,numlon)
             read (148,*) (v_model2(lati,lonj, zl), lonj=1,numlon) 
             read (151,*) (t_model1(lati,lonj, zl), lonj=1,numlon)
             read (152,*) (t_model2(lati,lonj, zl), lonj=1,numlon) 
             read (153,*) (s_model1(lati,lonj, zl), lonj=1,numlon)
             read (154,*) (s_model2(lati,lonj, zl), lonj=1,numlon) 
             if (three_dim.eq.1) then
               read (149,*) (w_model1(lati,lonj, zl), lonj=1,numlon)
               read (150,*) (w_model2(lati,lonj, zl), lonj=1,numlon) 
               read (160,*) (kz_model1(lati,lonj, zl), lonj=1,numlon)
               read (161,*) (kz_model2(lati,lonj, zl), lonj=1,numlon)
             endif 
           enddo
         enddo

      else


         open(145,file=trim(path)//trim(ADJUSTL(bas))//'u.txt',  status='old')

         open(147,file=trim(path)//trim(ADJUSTL(bas))//'v.txt',  status='old')
         
         open(151,file=trim(path)//trim(ADJUSTL(bas))//'temp.txt',  status='old')

         open(153,file=trim(path)//trim(ADJUSTL(bas))//'salt.txt',  status='old')

         if (three_dim.eq.1) then

           open(149,file=trim(path)//trim(ADJUSTL(bas))//'w.txt',  status='old')
  
           open(160,file=trim(path)//trim(ADJUSTL(bas))//'kz.txt',  status='old')

         endif 
 

         do zl=1,numz
           do lati=1,numlat
             read (145,*) (u_model1(lati,lonj, zl), lonj=1,numlon)
             u_model2 = u_model1
             read (147,*) (v_model1(lati,lonj, zl), lonj=1,numlon)
             v_model2 = v_model1 
             read (151,*) (t_model1(lati,lonj, zl), lonj=1,numlon)
             t_model2=t_model1
             read (153,*) (s_model1(lati,lonj, zl), lonj=1,numlon)
             s_model2=s_model1
             if (three_dim.eq.1) then
               read (149,*) (w_model1(lati,lonj, zl), lonj=1,numlon)
               w_model2 =  w_model1
               read (160,*) (kz_model1(lati,lonj, zl), lonj=1,numlon)
               kz_model2 = kz_model1
             endif 
           enddo
         enddo


      endif
   

         close(145)
         close(146)
         close(147)
         close(148)
         close(149)
         close(150)
         close(151)
         close(152)
         close(153)
         close(154)
         close(160)
         close(161)
       endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERA 5

      if (wind_theoretical.ne.1) then

       if(vec_era(i).ne.vec_era(i-1)) then

       if(tera1(1).ne.tera2(1)) then

         open(145,file=trim(path)//trim(ADJUSTL(bas_era))//'u_wind.txt',  status='old')
         open(146,file=trim(path)//trim(ADJUSTL(bas2_era))//'u_wind.txt',  status='old')


         open(147,file=trim(path)//trim(ADJUSTL(bas_era))//'v_wind.txt',  status='old')
         open(148,file=trim(path)//trim(ADJUSTL(bas2_era))//'v_wind.txt',  status='old')

         do lati_era=1,numlat_era
             read (145,*) (u10_era1(lati_era,lonj_era), lonj_era=1,numlon_era)
             read (146,*) (u10_era2(lati_era,lonj_era), lonj_era=1,numlon_era) 
             read (147,*) (v10_era1(lati_era,lonj_era), lonj_era=1,numlon_era)
             read (148,*) (v10_era2(lati_era,lonj_era), lonj_era=1,numlon_era) 
         enddo

         close(145)
         close(146)
         close(147)
         close(148)

        else

         open(145,file=trim(path)//trim(ADJUSTL(bas_era))//'u_wind.txt',  status='old')

         open(147,file=trim(path)//trim(ADJUSTL(bas_era))//'v_wind.txt',  status='old')

         do lati_era=1,numlat_era
             read (145,*) (u10_era1(lati_era,lonj_era), lonj_era=1,numlon_era)
             u10_era2 = u10_era1
             read (147,*) (v10_era1(lati_era,lonj_era), lonj_era=1,numlon_era)
             v10_era2 =  v10_era1
         enddo

         close(145)
         close(146)
         close(147)
         close(148)



        endif


       endif
      endif

    endif
!  !print*, trim(path)//trim(ADJUSTL(bas_era))//'u_wind.txt'
!  !print*, trim(path)//trim(ADJUSTL(bas2_era))//'u_wind.txt'
!  !print*, trim(path)//trim(ADJUSTL(bas_era))//'v_wind.txt'
!  !print*, trim(path)//trim(ADJUSTL(bas2_era))//'v_wind.txt'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !print*, vec_tdelf(i), vec_tdelf(i-1), counttimeh


!!print*, 'c'

! lat_in = minloc(abs(lat_model - lat_part(j,i-1)))

! lon_in = minloc(abs(lon_model - lon_part(j,i-1)))

!!print*, lat_model, lat_part(j,i-1), lon_part(j,i-1), lon_in, lat_in
!!print*, lon_in, lat_in

       lat_model_sum=abs(lat_model-lat_part(j,i-1))
       lon_model_sum=abs(lon_model-lon_part(j,i-1))

       coord_sum=lat_model_sum + lon_model_sum

       lat_in =  minloc(coord_sum)
       lon_in=lat_in

       ! lat_model_summ=abs(lat_modelm-lat_part(j,i-1))
       ! lon_model_summ=abs(lon_modelm-lon_part(j,i-1))

       ! coord_summ=lat_model_summ + lon_model_summ

       ! lat_inm =  minloc(coord_summ)
       ! lon_inm=lat_inm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERA 5
     if (wind_theoretical.ne.1) then

       lat_era_sum=abs(lat_era-lat_part(j,i-1))
       lon_era_sum=abs(lon_era-lon_part(j,i-1))

       coord_sum_era = lat_era_sum + lon_era_sum

       lat_in_era =  minloc(coord_sum_era)
       lon_in_era = lat_in_era


       if(linear_interp.eq.1) then

          if ( (lat_in_era(1)+slic .gt. numlat_era) .or. (lat_in_era(1)- slic .lt. 1)  .or.  (lon_in_era(2)+slic &
 .gt. numlon_era) .or.  (lon_in_era(2)- slic .lt. 1)    ) then
             go to 85645
          endif

         slic_lat = lat_era(lat_in_era(1)-slic : lat_in_era(1)+slic, lon_in_era(2) - slic : lon_in_era(2) + slic )
         slic_lon = lon_era(lat_in_era(1)-slic : lat_in_era(1)+slic, lon_in_era(2) - slic : lon_in_era(2) + slic )

         call init_interpolation_wind(lat_in_era(1), lon_in_era(2))

         call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
         uwind1 = out_int(1)
         call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
         uwind2 = out_int(1)
         call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
         vwind1 = out_int(1)
         call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
         vwind2 = out_int(1)


!!print*, 'wind', vwind2

       else

         85645 continue

         uwind1 = u10_era1 (lat_in_era(1), lon_in_era(2))
         uwind2 = u10_era2 (lat_in_era(1), lon_in_era(2))
         vwind1 = v10_era1 (lat_in_era(1), lon_in_era(2))
         vwind2 = v10_era2 (lat_in_era(1), lon_in_era(2))

!!print*, 'windididi', vwind2
    
       endif

     endif


!!!!!!!!!!!!!!!!!!!!!!11
       if (zlayer .eq. 1) then

         deplevel = levelsd

! !print*, deplevel

       else
         deplevel = depth(lat_in(1), lon_in(2))  * levelsd

       endif 


       indexz = minloc(abs(deplevel-zf1(j,i-1)))


!print*, zf1

!numz=3
!deallocate(deplevel)
!allocate(deplevel(5))
!deplevel=(/0,-2,-3,-5,-6/)

!!print*, 'cachara', zf1, deplevel, numz, levelsd
!!print*, counttimeh_r, counttimeh

!PRINT*, indexz

!


!

  if(linear_interp.eq.1) then

     if ( (lat_in(1)+slic .gt. numlat) .or. (lat_in(1)- slic .lt. 1)  .or.  (lon_in(2)+slic .gt. numlon) .or. &
 (lon_in(2)- slic .lt. 1)    ) then
!        print*, '8888'
        go to 85641
     endif


     slic_lat = lat_model(lat_in(1)-slic : lat_in(1)+slic, lon_in(2) - slic : lon_in (2) + slic )
     slic_lon = lon_model(lat_in(1)-slic : lat_in(1)+slic, lon_in(2) - slic : lon_in (2) + slic )

     if (   (indexz(1).eq.1 .and. zf1(j,i-1).ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz &
	 .and. zf1(j,i-1).le. deplevel(indexz(1)))  .or. (zf1(j,i-1).eq.deplevel(indexz(1))) &
	 .or. (size(deplevel).eq.1) ) then

!!print*, 'abrac4', zf1, deplevel(indexz(1)), indexz
!!print*, deplevel, shape(u_model1)

    call init_interpolation(1, lat_in(1), lon_in(2), indexz(1), 0)


 !   call pwl_interp_2d ( size(lon_1d), size(lat_1d), lon_1d, lat_1d, slic_u1,1, inlon, inlat, out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        ui1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        ui2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        vi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        vi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        ti1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        ti2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        si1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        si2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        wi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        wi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        kz1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        kz2 = out_int(1)

!  !print*, 'luluaaaa', vi2

     else

        if (zf1(j,i-1).gt.deplevel(indexz(1))) then
           vert_index1 = indexz(1)
           vert_index2 = indexz(1) - 1
        else 

           vert_index1 = indexz(1)
           vert_index2 = indexz(1) + 1

        endif



       call init_interpolation(2, lat_in(1), lon_in(2), vert_index1, vert_index2)



        inter_depth(1) = deplevel(vert_index1)
        inter_depth(2) = deplevel(vert_index2)
        par_dep(1) = zf1(j,i-1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi2 = ui_vec(1)
                
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ti1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, t2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ti2 = ui_vec(1)
        
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        si1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, s2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        si2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d_2 , power, 1, lon_part(j,i-1), lat_part(j,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz2 = ui_vec(1)

!   !print*, 'luluzao', vi2
!   !print*, inter_depth
!   !print*, zf1
!   !print*, u2_1d
!   !print*, inter_depth
!   !print*, u2_1d_2 
  
     endif

!!print*, si1,si2, ti1, ti2

  else


       85641 continue

       ui1 = u_model1 (lat_in(1), lon_in(2), indexz(1))
       ui2 = u_model2 (lat_in(1), lon_in(2), indexz(1))
       vi1 = v_model1 (lat_in(1), lon_in(2), indexz(1))
       vi2 = v_model2 (lat_in(1), lon_in(2), indexz(1))
       ti1 = t_model1 (lat_in(1), lon_in(2), indexz(1))
       ti2 = t_model2 (lat_in(1), lon_in(2), indexz(1))
       si1 = s_model1 (lat_in(1), lon_in(2), indexz(1))
       si2 = s_model2 (lat_in(1), lon_in(2), indexz(1))       
       wi1 = w_model1 (lat_in(1), lon_in(2), indexz(1))
       wi2 = w_model2 (lat_in(1), lon_in(2), indexz(1))
       kz1 = kz_model1 (lat_in(1), lon_in(2), indexz(1))
       kz2 = kz_model2 (lat_in(1), lon_in(2), indexz(1))
	   

!  print*, vi2, depth(lat_in(1), lon_in(2))

  endif


       ! di = depth(lat_inm(1), lon_inm(2))
       ! coastad = coastvalue(lat_inm(1), lon_inm(2))

	 !  print*, di

       intertime(1) = time_vec(tdelf1(1))
       intertime(2) = time_vec(tdelf2(1))

       vel_interp(1) = ui1
       vel_interp(2) = ui2

       ts_vec(1) = counttimeh

       if (intertime(1).eq.intertime(2)) then                     !this option can only be satisfied if the end of the time series is reached
         intertime(1) = time_vec(numtime-1)
       endif


      ts_vec_w(1) = counttimeh_era
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!era 5
      if (wind_theoretical.ne.1) then
         intertime_era(1) = time_vec_era(tera1(1))
         intertime_era(2) = time_vec_era(tera2(1))

         if (intertime_era(1).eq.intertime_era(2)) then 
           intertime_era(1) = time_vec_era(num_time_era-1)
         endif


         vel_interp_era(1) = uwind1
         vel_interp_era(2) = uwind2

         call pwl_value_1d (size(intertime_era), intertime_era, vel_interp_era, 1, ts_vec_w , ui_vec )
         uwd= ui_vec(1)


       ! !print*, intertime_era
       ! !print*, vel_interp_era
!         !print*, ts_vec_w, counttimeh, counttimeh_era
       ! !print*, ui_vec


         vel_interp_era(1) = vwind1
         vel_interp_era(2) = vwind2

         call pwl_value_1d (size(intertime_era), intertime_era, vel_interp_era, 1, ts_vec_w , ui_vec )
         vwd= ui_vec(1)
    
         windsp = ((uwd**2 + vwd**2)**0.5)*3600   !wind velocity in m/h
         windspms = (uwd**2 + vwd**2)**0.5


     if (reverse .eq. 1) then
		vwd=-vwd
		uwd=-uwd
      endif

 ! print*, uwd, vwd
!!print*, intertime_era, ts_vec, counttimeh
       ! !print*, intertime_era
       ! !print*, vel_interp_era
       ! !print*, ts_vec
       ! !print*, ui_vec     
       !!print*, uwd, vwd
      endif

!       !print*, 'caos', uwd, vwd, ts_vec_w, counttimeh_era
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       ui= ui_vec(1)

!!print*, vec_tdelf(i), vec_tdelf(i-1)
! !print*, intertime, ts_vec
!!print*, vel_interp, ui_vec
 !!print*, tdelf1(1), tdelf2(1)

       vel_interp(1) = vi1
       vel_interp(2) = vi2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       vi=ui_vec(1)

       vel_interp(1) = ti1
       vel_interp(2) = ti2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       ti=ui_vec(1)
       
       vel_interp(1) = si1
       vel_interp(2) = si2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       si=ui_vec(1) 
       
!print*, si, ti
       
       vel_interp(1) = wi1
       vel_interp(2) = wi2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       wi=ui_vec(1)

       vel_interp(1) = kz1
       vel_interp(2) = kz2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       kz=ui_vec(1)
       
!	   temps(j,i)=ti(1,1)
!	   salts(j,i)=si(1,1)
 !ui=-1
 !vi=-1
      if (three_dim .ne. 1) then
         kz=0.001
      endif
! !print*, 'OWOWOWOW', ui, vi

!!print*, u_model2(1,:)
!
     if (reverse .eq. 1) then
		ui=-ui
		vi=-vi
      endif

  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111DELFT

!       !print*, ui, vi
!       !print*, ui
 !      !print*, minloc(abs(x_model(1,:)-x(j,i-1)))
 !      !print*, minloc(abs(y_model(:,1)-y(j,i-1)))
 !      !print*, x_model(1,:)-x(j,i-1)
       wanglex=wanglex*(pi/180)
       wangley=wangley*(pi/180)

       ow=(beta*exp(((-10.**(-8))*((windspx/3600)**3))/(10.*0.1)))    


!       uwd=0.03*(windspx/3600)*sin((wanglex-180+ow)*(3.14/180.))
 !      vwd=0.03*(windspx/3600)*cos((wangley-180+ow)*(3.14/180.))

 !      !print*, ow
!       !print*, beta
!       !print*, exp(((-0.38*(-10.**(-8)))*((windspx/3600)**3))/(10.*0.1))
!!       !print*, uwd
!       !print*, vwd
	   temps(j,i)=uwd
	   salts(j,i)=vwd
!       !print*, windspx/3600
      
!!print*,  massa1(j,i)/massa(1,1)
!!print*, sum(FRAC_MASS_OUT)
!!print*, evapmass(j,2)

77958 continue

 if (dissolved_fase .eq. 1) then   !! BEGIN FASE 3
 
! print*, "dissolved_fase", tsevol(contind-1), tsevol(contind)

  if (tsevol(contind) .ne. tsevol(contind-1)) then


    if(parcel_dis_cont.gt.0) then

     do m1_f3=1,parcel_dis_cont

       lat_model_sum=abs(lat_model-lat_partf3(m1_f3,i-1))
       lon_model_sum=abs(lon_model-lon_partf3(m1_f3,i-1))
       coord_sum=lat_model_sum + lon_model_sum
       lat_in_f3 =  minloc(coord_sum)
       lon_in_f3=lat_in_f3

       if (zlayer .eq. 1) then
         deplevel = levelsd
       else
         deplevel = depth(lat_in_f3(1), lon_in_f3(2))  * levelsd
       endif

       indexz = minloc(abs(deplevel-zf3(m1_f3,i-1)))
    !   dif3  = depth(lat_in_f3(1), lon_in_f3(2))

        ! if (dif3.le.0) then
         ! zf3(m1_f3,i) = zf3(m1_f3,i-1)
         ! xf3(m1_f3,i) = xf3(m1_f3,i-1)
         ! yf3(m1_f3,i) = yf3(m1_f3,i-1)
         ! lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         ! lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
 ! !       !print*, 'DI', di               
         ! cycle

        ! elseif(-zf3(m1_f3,i-1) .ge. dif3) then
         ! zf3(m1_f3,i) = zf3(m1_f3,i-1)
         ! xf3(m1_f3,i) = xf3(m1_f3,i-1)
         ! yf3(m1_f3,i) = yf3(m1_f3,i-1)
         ! lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         ! lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         ! cycle
       
        ! ! elseif ( (lat_model(lat_in_f3(1)-1, lon_in_f3(2)).eq.0) .or. (lat_model(lat_in_f3(1)+1, lon_in_f3(2)) .eq.0)  .or. &
                 ! ! (lat_model(lat_in_f3(1), lon_in_f3(2)+1).eq.0) .or. (lat_model(lat_in_f3(1), lon_in_f3(2)-1) .eq.0) ) then
         ! ! xf3(m1_f3,i)=xf3(m1_f3,i-1)
         ! ! yf3(m1_f3,i)=yf3(m1_f3,i-1)
         ! ! lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         ! ! lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         ! ! zf3(m1_f3,i) = zf3(m1_f3,i-1)
         ! ! cycle

        ! ! elseif ( (lon_model(lat_in_f3(1)-1, lon_in_f3(2)).eq.0) .or. (lon_model(lat_in_f3(1)+1, lon_in_f3(2)) .eq.0)  .or. &
                 ! ! (lon_model(lat_in_f3(1), lon_in_f3(2)+1).eq.0) .or. (lon_model(lat_in_f3(1), lon_in_f3(2)-1) .eq.0) ) then
         ! ! xf3(m1_f3,i)=xf3(m1_f3,i-1)
         ! ! yf3(m1_f3,i)=yf3(m1_f3,i-1)
         ! ! lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         ! ! lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         ! ! zf3(m1_f3,i) = zf3(m1_f3,i-1)
         ! ! cycle
        ! endif



  if(linear_interp.eq.1) then

     if ( (lat_in_f3(1)+slic .gt. numlat) .or. (lat_in_f3(1)- slic .lt. 1)  .or.  (lon_in_f3(2)+slic .gt. numlon) &
  .or.  (lon_in_f3(2)- slic .lt. 1)    ) then
        go to 85698
     endif


     slic_lat = lat_model(lat_in_f3(1)-slic : lat_in_f3(1)+slic, lon_in_f3(2) - slic : lon_in_f3(2) + slic )
     slic_lon = lon_model(lat_in_f3(1)-slic : lat_in_f3(1)+slic, lon_in_f3(2) - slic : lon_in_f3(2) + slic )

     if (   (indexz(1).eq.1 .and. zf3(m1_f3,i-1).ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz .and. &
 zf3(m1_f3,i-1).le. deplevel(indexz(1)))  .or. (zf3(m1_f3,i-1).eq.deplevel(indexz(1))) .or. (size(deplevel).eq.1)) then

!!print*, 'abrac4f3', zf3(m1_f3,i-1), deplevel(indexz(1)), indexz, levelsd, size(deplevel)


    call init_interpolation(1, lat_in_f3(1), lon_in_f3(2), indexz(1), 0)


 !   call pwl_interp_2d ( size(lon_1d), size(lat_1d), lon_1d, lat_1d, slic_u1,1, inlon, inlat, out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        ui1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        ui2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
        vi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        vi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        wi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        wi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        kz1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        kz2 = out_int(1)

!  !print*, 'luluaaaaf3', ui2

     else

        if (zf3(m1_f3,i-1).gt.deplevel(indexz(1))) then
           vert_index1 = indexz(1)
           vert_index2 = indexz(1) - 1
        else 

           vert_index1 = indexz(1)
           vert_index2 = indexz(1) + 1

        endif

       call init_interpolation(2, lat_in_f3(1), lon_in_f3(2), vert_index1, vert_index2)


        inter_depth(1) = deplevel(vert_index1)
        inter_depth(2) = deplevel(vert_index2)
        par_dep(1) = zf3(m1_f3,i-1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_partf3(m1_f3,i-1),&
 lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d_2 , power, 1, lon_partf3(m1_f3,i-1),&
 lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_partf3(m1_f3,i-1),&
 lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_partf3(m1_f3,i-1),&
 lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d_2 , power, 1, lon_partf3(m1_f3,i-1),&
 lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d_2 , power, 1, lon_partf3(m1_f3,i-1), &
lat_partf3(m1_f3,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz2 = ui_vec(1)
  
     endif

  else

       85698 continue

       ui1 = u_model1 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       ui2 = u_model2 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       vi1 = v_model1 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       vi2 = v_model2 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       wi1 = w_model1 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       wi2 = w_model2 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       kz1 = kz_model1 (lat_in_f3(1), lon_in_f3(2), indexz(1))
       kz2 = kz_model2 (lat_in_f3(1), lon_in_f3(2), indexz(1))

!  !print*, 'lulu1f3', ui2

  endif


       vel_interp(1) = ui1
       vel_interp(2) = ui2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       uif3= ui_vec(1)
       vel_interp(1) = vi1
       vel_interp(2) = vi2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       vif3=ui_vec(1)
       vel_interp(1) = wi1
       vel_interp(2) = wi2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       wif3=ui_vec(1)
       vel_interp(1) = kz1
       vel_interp(2) = kz2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       kzf3=ui_vec(1)

       if (three_dim .ne. 1) then
         kzf3=0.001
       endif



       if (wind_theoretical.ne.1) then

        lat_era_sum=abs(lat_era-lat_partf3(m1_f3,i-1))
        lon_era_sum=abs(lon_era-lon_partf3(m1_f3,i-1))
 
        coord_sum_era = lat_era_sum + lon_era_sum

        lat_in_era =  minloc(coord_sum_era)
        lon_in_era = lat_in_era



        if(linear_interp.eq.1) then

           if ( (lat_in_era(1)+slic .gt. numlat_era) .or. (lat_in_era(1)- slic .lt. 1)  .or.  &
 (lon_in_era(2)+slic .gt. numlon_era) .or.  (lon_in_era(2)- slic .lt. 1)    ) then
              !print*, '00000a'
              go to 85634
           endif

          slic_lat = lat_era(lat_in_era(1)-slic : lat_in_era(1)+slic, lon_in_era(2) - slic : lon_in_era(2) + slic )
          slic_lon = lon_era(lat_in_era(1)-slic : lat_in_era(1)+slic, lon_in_era(2) - slic : lon_in_era(2) + slic )

          call init_interpolation_wind(lat_in_era(1), lon_in_era(2))

          call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
          uwind1 = out_int(1)
          call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
          uwind2 = out_int(1)
          call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
          vwind1 = out_int(1)
          call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_partf3(m1_f3,i-1), &
 lat_partf3(m1_f3,i-1), out_int )
          vwind2 = out_int(1)


!!print*, 'windf3', vwind2

        else

          85634 continue

          uwind1 = u10_era1 (lat_in_era(1), lon_in_era(2))
          uwind2 = u10_era2 (lat_in_era(1), lon_in_era(2))
          vwind1 = v10_era1 (lat_in_era(1), lon_in_era(2))
          vwind2 = v10_era2 (lat_in_era(1), lon_in_era(2))

!!print*, 'windiif3', vwind2
    
        endif

     

        vel_interp_era(1) = uwind1
        vel_interp_era(2) = uwind2

        call pwl_value_1d (size(intertime_era), intertime_era, vel_interp_era, 1, ts_vec_w , ui_vec )
        uwdf3= ui_vec(1)

        vel_interp_era(1) = vwind1
        vel_interp_era(2) = vwind2
 
        call pwl_value_1d (size(intertime_era), intertime_era, vel_interp_era, 1, ts_vec_w , ui_vec )
        vwdf3= ui_vec(1)

       else

        uwdf3 = uwd

        vwdf3 = vwd
    
       endif



      if (zf3(m1_f3,i-1) .eq. 0) then

       if (right_random .eq. 0) then

         CALL random_number(RN1)
         CALL random_number(RN2)
         CALL random_number(RN3)

         randvert = -1. + 2.*RN3
       
         U_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * COS(2.D0*PI*RN2)    
         V_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * SIN(2.D0*PI*RN2)
         W_ALEA = randvert * ((2.D0*kzf3/dt)**(0.5D0))    !!  based on Reed et al., 1995   !! PAY ATTENTION KZF3

         xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt + U_ALEA*dt   +  ( widfc * (  (uwdf3 * cos(5*(pi/180))) &
 + (vwdf3 * sin(15*(pi/180))) )  )   *dt 

         yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt + V_ALEA*dt   +  ( widfc * (  (-uwdf3 * sin(5*(pi/180))) &
 + (vwdf3 * cos(15*(pi/180))) )  )   *dt

         zf3(m1_f3,i)=zf3(m1_f3,i-1) + wif3*dt +  W_ALEA*dt

         if ( zf3(m1_f3,i) .gt. 0) then
           zf3(m1_f3,i) = 0
         endif

         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(m1_f3,i), LAT_partf3(m1_f3,i), PROFOUT,  &
 xf3(m1_f3,i),  yf3(m1_f3,i), ZPOS, OPT)

       else 

        do step_random=1,dt, dt_random
          CALL random_number(RN1)
          CALL random_number(RN2)
          CALL random_number(RN3)


          randvert = -1. + 2.*RN3
       
           xrandom = xrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * COS(2.D0*PI*RN2) ) *  dt_random
           yrandom = yrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * SIN(2.D0*PI*RN2) ) *  dt_random   
           zrandom = zrandom + ( randvert * ((2.D0*kzf3/dt_random)**(0.5D0)) ) *dt_random  !!  based on Reed et al., 1995
        enddo


        xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt +  xrandom +  ( widfc * (  (uwdf3 * cos(5*(pi/180)))  + &
 (vwdf3 * sin(15*(pi/180))) )  )   *dt 

        yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt +  yrandom +  ( widfc * (  (-uwdf3 * sin(5*(pi/180)))  + &
 (vwdf3 * cos(15*(pi/180))) )  )   *dt

        zf3(m1_f3,i)=zf3(m1_f3,i-1) + wif3*dt +  zrandom

         if ( zf3(m1_f3,i) .gt. 0) then
           zf3(m1_f3,i) = 0
         endif

        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(m1_f3,i), LAT_partf3(m1_f3,i), PROFOUT,  &
xf3(m1_f3,i),  yf3(m1_f3,i), ZPOS, OPT)

        xrandom = 0 
        yrandom = 0
        zrandom = 0
 
       endif

             

      else  

       if (right_random .eq. 0) then


         CALL random_number(RN1)
         CALL random_number(RN2)
         CALL random_number(RN3)


         randvert = -1. + 2.*RN3
       
          U_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * COS(2.D0*PI*RN2)    
          V_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * SIN(2.D0*PI*RN2)
          W_ALEA = randvert * ((2.D0*kzf3/dt)**(0.5D0))    !!  based on Reed et al., 1995
 

         xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt +   U_ALEA*dt     !no wind
 
         yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt +   V_ALEA*dt     !no wind
 
  !     !print*, '444', spmt, spmtf2(m1,i-1)


         zf3(m1_f3,i)=zf3(m1_f3,i-1) + wif3*dt +  W_ALEA*dt

         if ( zf3(m1_f3,i) .gt. 0) then
           zf3(m1_f3,i) = 0
         endif


         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(m1_f3,i), LAT_partf3(m1_f3,i), PROFOUT,  &
 xf3(m1_f3,i),  yf3(m1_f3,i), ZPOS, OPT)

      

       else


         do step_random=1,dt, dt_random
          CALL random_number(RN1)
          CALL random_number(RN2)
          CALL random_number(RN3)


          randvert = -1. + 2.*RN3

       
          xrandom = xrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * COS(2.D0*PI*RN2) ) *  dt_random
          yrandom = yrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * SIN(2.D0*PI*RN2) ) *  dt_random   
          zrandom = zrandom + ( randvert * ((2.D0*kzf3/dt_random)**(0.5D0)) ) *dt_random  !!  based on Reed et al., 1995
         enddo


         xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt +   xrandom 

         yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt +   yrandom



         zf3(m1_f3,i)=zf3(m1_f3,i-1) + wif3*dt +  zrandom

         if ( zf3(m1_f3,i) .gt. 0) then
           zf3(m1_f3,i) = 0
         endif

         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(m1_f3,i), LAT_partf3(m1_f3,i), PROFOUT, &
 xf3(m1_f3,i),  yf3(m1_f3,i), ZPOS, OPT)


        xrandom = 0 
        yrandom = 0
        zrandom = 0


       endif



      endif


       lat_model_summ=abs(lat_modelm-lat_partf3(m1_f3,i))
       lon_model_summ=abs(lon_modelm-lon_partf3(m1_f3,i))

       coord_summ=lat_model_summ + lon_model_summ

       lat_inm =  minloc(coord_summ)
       lon_inm=lat_inm

       dif3 = depth(lat_inm(1), lon_inm(2))
	   
	   
	    if (dif3.le.0) then
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         xf3(m1_f3,i) = xf3(m1_f3,i-1)
         yf3(m1_f3,i) = yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
 !       !print*, 'DI', di               
         cycle

        elseif(-zf3(m1_f3,i-1) .ge. dif3) then
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         xf3(m1_f3,i) = xf3(m1_f3,i-1)
         yf3(m1_f3,i) = yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         cycle
       
        endif
     
     enddo
 
    endif

  endif


 endif   !!END FASE 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EVAPORATION SECTION
temp_out=ti(1,1) + 273.15
temp=temp_out
SAL_A=si(1,1)

!!print*, 'jfjf', ti, temp_out, sal_a, si, temp

ii=i
      
jj=j


aux=0


!print*,  '1ppp',  massa(j, i-1), zf1(j,i-1)
deg_turn=1
if (deg_turn .eq. 1) then
 if (mass_diss(j, i-1) .NE. 0) then
    ! Ação a ser executada se a e b não forem iguais

       tau = 12.0 * (3.0 ** ((20.0 - ti(1,1)) / 10.0))

    ! Calcular a fração biodegradada
       fraction_biodegraded = 1.0 - exp((-counttime/(60*60*24)) / tau)
	   
!      print*, counttime/(60*60*24), ti(1,1), mass_diss(j,i-1), mass_diss(j,i-1)*fraction_biodegraded
         
		mass_degr(j,i) =    mass_degr(j, i-1) + mass_diss(j,i-1)*fraction_biodegraded
        mass_diss(j,i) =    mass_diss(j, i-1)  - mass_diss(j,i-1)*fraction_biodegraded
		mass_diss(j, i-1)=mass_diss(j,i)
 end if		
endif 

if ( (massa(j,i-1).eq.0) ) then
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_sedi(j,i)=mass_sedi(j,i-1)
	   mass_degr(j,i)=mass_degr(j,i-1)
	   zf1(j,i)=zf1(j,i-1)
	   
   cycle
 endif 
 
 if (  (beached(j,i-1) .eq.1) ) then

!print*, 'rrrrrrrrrrrrrrrr', j
              massa(j,i) =  massa(j,i-1)
              massae(j,i) = massae(j,i-1)
              masscomp(j,i,:) = masscomp(j,i-1,:)
              diam(j,i) = diam(j,i-1)
              diamem(j,i) = diamem(j,i-1)
              area(j,i) = area(j,i-1)
              areaem(j,i) = areaem(j,i-1)
              height(j,i) = height(j,i-1)
              vol(j,i)    = vol(j,i-1)  
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_sedi(j,i)=mass_sedi(j,i-1)
	   mass_degr(j,i)=mass_degr(j,i-1)
	   rho_e(j,i)=rho_e(j,i-1)
	   visc_e(j,i)=visc_e(j,i-1)
	   watf(j,i)=watf(j,i-1)
	   zf1(j,i)=0
	    lon_part(j,i)=lon_part(j,i-1)
		lat_part(j,i)=lat_part(j,i-1)
		contind=contind+1		
       cycle
	  		stop
 
 endif 
 

	   
 if (  (zf1(j,i-1) .lt.0) ) then


       lat_model_summ=abs(lat_modelm-lat_part(j,i-1))
       lon_model_summ=abs(lon_modelm-lon_part(j,i-1))

       coord_summ=lat_model_summ + lon_model_summ

       lat_inm =  minloc(coord_summ)
       lon_inm=lat_inm

 !      di = depth(lat_inm(1), lon_inm(2))
!	   print*,'4', zf1(j,i-1), di
	   ! if (-zf1(j,i-1) .ge. di) then
	      ! zf1(j,i-1)=-di
	   ! endif

!print*, 'rrrrrrrrrrrrrrrr', j
              massa(j,i) =  massa(j,i-1)
              massae(j,i) = massae(j,i-1)
              masscomp(j,i,:) = masscomp(j,i-1,:)
              diam(j,i) = diam(j,i-1)
              diamem(j,i) = diamem(j,i-1)
               dropdiam(j,i) = dropdiam(j,i-1)			  
              area(j,i) = area(j,i-1)
              areaem(j,i) = areaem(j,i-1)
              height(j,i) = height(j,i-1)
              vol(j,i)    = vol(j,i-1)  
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   rho_e(j,i)=rho_e(j,i-1)
	   visc_e(j,i)=visc_e(j,i-1)
	   watf(j,i)=watf(j,i-1)

   go to 18767     !!! go to entrainment

 endif 
 



          CALL vapour_pressure(TEMP_OUT)

!print*, '5555555555555555', sum(masscomp(j,i-1,:))
!print*, "dia", counttime/(60*60*24)

 IF (EVAP_TURN.EQ.1) THEN
       do comps=1, NCOMP_OIL
	   
	!       print*, '11', masscomp(j,i-1, comps)

           if (masscomp(j,i-1, comps).eq.0) then
		     masscomp(j,i, comps)=0
              cycle
           endif
           
		   IF (windsp.le.0) THEN
		     windsp=0.1*3600
		     windsp=10
           endif
!		   print*, windsp

		   
           kf1=(0.0292*(windsp**0.78)*(diam(j,i-1)**(-0.11))*(scf**(-0.67))*(((PM_COMP_OIL(comps)+29.)/(PM_COMP_OIL(comps)))**0.5))  !unit grams, hours

           evap1=(((((kf1*PC(comps)*area(j,i-1))/(r*temp))*(PM_COMP_OIL(comps)*FRAC_MASS_OUT_PART(j,comps)))*0.001)/3600.)*dt   !0.001 is to convert from grams to kg    !betancout et al 2005+
 
            if (evap1.gt.evapmass(j, comps)) then
               evap1 = evapmass(j, comps)
            endif
  !         !print*, evap1, masscomp(j,i-1, comps), evapmass(j, comps)

  !           if (evap1.gt.masscomp(j,i-1, comps)) then
            !!print*, PC(comps)
  !               evap1 = masscomp(j,i-1, comps)
  !           endif
  !          print*, '2',evap1, evapmass(j, comps)
             if (counttime/(60*60*24).gt. 15) then
			  evap1=0
!			  print*, counttime/(60*60*24), aux
!			  stop
			 endif

            masscomp(j,i, comps)=masscomp(j,i-1, comps)-evap1

            if (masscomp(j,i, comps).le.0) then
               masscomp(j,i, comps)=0
            endif
        
   !        evapmass(j, comps)=evapmass(j, comps)-evap1   wrong

            evapmass(j, comps) = masscomp(j,i, comps)

   !         evapmass(j, comps) = ((maxwf-watf(j,i-1) )/maxwf) *  masscomp(j,i,comps)           !!!!EMULSIFICATION INTERFERENCE ON EVAPORATION

   !       !print*, watf(j,i-1), maxwf, evapmass(j, 10), masscomp(j,i,10), evap1
   
            aux=aux+evap1

       enddo

       mas_evap(j,i)  = mas_evap(j,i-1) + aux


!print*, '666666666666666666666666', temp_out
!print*, PC


!!print*, 'evapmass1', evapmass(23), masscomp(j,i, 23)
!!print*, '666', sum(masscomp(j,i, :))
!!!!!!!!!!!!!!!!!!!!EVAPORATION SECTION

!print*, 'o', massa(j,i)

massa(j,i)=sum(masscomp(j,i, :))

!print*, massa(j,i), massa(j,i-1)

!print*, '880', watf(j,i-1),  maxwf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


massa(j,i-1) =  massa(j,i)
                                            ! To overcome the problem of evaporation with dissolutino


masscomp(j, i-1, :) = masscomp(j, i, :)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111


porc_evap(j,i) = (mas_evap(j,i) / massref)*100


  CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


 
   aux=0
   do comps=1, NCOMP_OIL
     aux=aux + masscomp(j,i,comps)/RO_COMP_OIL(comps)
   enddo

   spmt=massa(j,i)/aux


     height(j,i)=massa(j,i)/(spmt*pi*((diam(j,i-1)/2)**2))

     height(j,i-1) = height(j,i)


do comps=1, NCOMP_OIL

   MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

   FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

   V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

   FRAC_MASS_OUT_PART(j,comps) = masscomp(j,i, comps)/massa(j,i)
enddo


call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


      MOL_TOT=sum(mol_comp)

       V_TOT=massa(j,i)/spmt
 
      VOL(j,i)=massa(j,i)/spmt

      moil=massa(j,i)

     FRAC_TOT=sum(FRAC_MASS_OUT)


    vol_evap(j,i)  = vol_evap(j,i-1) + ( VOL(J,I-1) -  VOL(j,i)  )
 
    porc_evap_vol(j,i) =  (  vol_evap(j,i) / volreff   ) * 100
    
!	print*, VOL(J,I-1),  VOL(j,i)
!   !print*, temp_out  


      VOL(j,i-1) = VOL(j,i)

   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT,TS_A, TS_VC )

!!!!!!!!!!!!!!!!!!!!recalculate droplets info


!!!!!! definition of droplet info
   dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

 !  dropdiam(j,i-1) =  dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)
   
   numdrops(j,i-1)= numdrops(j,i)

   if (massa(j,i).le.0) then
       massa(j,:)= 0
       lat_part(j,:)= 0
       lon_part(j,:)= 0
       dropdiam(j,:)= 0
       voldrops(j,:)= 0
       numdrops(j,:)= 0
       massae(j,:)  = 0
       area(j,:)    = 0
       areaem(j,:)  = 0
       diam(j,:)    = 0
       diamem(j,:)  = 0   
       height(j,:)  = 0
       vol(j,:)     = 0
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_sedi(j,i)=mass_sedi(j,i-1)	   
	   mass_degr(j,i)=mass_degr(j,i-1)	   
!       go to 18767     !!! go to entrainment
       cycle
   endif

!!!!!!!!!!!!!

!!print*,  'numdrop', numdrops(j,i), vol(j,i), voldrops (j,i)

else 


              massa(j,i) =  massa(j,i-1)
              massae(j,i) = massae(j,i-1)
              masscomp(j,i,:) = masscomp(j,i-1,:)
              diam(j,i) = diam(j,i-1)
              diamem(j,i) = diamem(j,i-1)
              area(j,i) = area(j,i-1)
              areaem(j,i) = areaem(j,i-1)
              height(j,i) = height(j,i-1)
              vol(j,i)    = vol(j,i-1)


          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


 
          aux=0
          do comps=1, NCOMP_OIL
            aux=aux + masscomp(j,i,comps)/RO_COMP_OIL(comps)
          enddo

          spmt=massa(j,i)/aux

       do comps=1, NCOMP_OIL

         MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

         FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

         V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

         FRAC_MASS_OUT_PART(j,comps) = masscomp(j,i, comps)/massa(j,i)
       enddo

       call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


       MOL_TOT=sum(mol_comp)

       V_TOT=massa(j,i)/spmt
 
       VOL(j,i)=massa(j,i)/spmt

       moil=massa(j,i)

       FRAC_TOT=sum(FRAC_MASS_OUT)

   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

   dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

  ! dropdiam(j,i-1) =  dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)
   
   numdrops(j,i-1)= numdrops(j,i)


!!print*, 'aaaaaaaaa', RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT

endif

!print*, '88', i, massa(j,i)


 	
	
	
 IF (EMULSI.EQ.1) THEN


   CALL EMULSIFY(DT, WINDSPMS, XAM(j), XWM(j) )            !!!!!!!!!!TURN ON EMULSIFICATION

!   watcont(j,i)=watcont(j-1,i) + watup + watout
   if (watf(j,i) .gt. maxwf) then
      watf(j,i) = maxwf
   endif
   if (watf(j,i) .lt. 0) then
      watf(j,i) =0
   endif   
   watcont(j,i)=massa(j,i-1)*watf(j,i)

!   !print*, massa(j,i), massae(j,i-1)


   massae(j,i)=massa(j,i)+watcont(j,i)

   massae(j,i-1) = massae(j,i)    
 !  call visc_emulsion(VIS_DIN_OIL_OUT)

   call visc_emulsion(VIS_DIN_OIL_OUT*1000)  !in cP

   call rho_emulsion(spmt)


   diamem(j,i)=2*(((massae(j,i)/rho_e(j,i))/(pi*height(j,i)))**0.5)  !meter

   diamem(j,i-1) = diamem(j,i)

   areaem(j,i)=pi*((diamem(j,i)/2)**2)

   areaem(j,i-1) = areaem(j,i)


 ! !print*, 'kgjd5', areaem(j,i), (massae(j,i)/rho_e(j,i))/height(j,i)

ENDIF 

!!print*, watcont(j,i), WINDSPMS,dt

!print*, 'API', (141.5D0 / (spmt / 1000.D0)) - 131.5D0
   !!!
!   abc=(1-0.658*(watf*0.65))*abc   !!!relation emulsification / evaporation
!   !print*, 'skhkjlsdkgjhldkfg', abc 
!!!!

!   !print*, watcont(j,i)
  

!   !print*,  VIS_DIN_OIL_OUT


     
!     !print*, FRAC_MASS_OUT(1)

 !    !print*, massa(j,i)



 DIAM_HYD=0



!  call VEL_ASCENSAO_BG ( diam(j,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )        !vertical velocity

call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

!spmt=1000

!print*, dropdiam(j,i-1) , DIAM_HYD , RO_A , spmt,  rho_e(j,i)

if (emulsi.eq.0) then

 call VEL_ASCENSAO_BG ( dropdiam(j,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )

else 

 call VEL_ASCENSAO_BG ( dropdiam(j,i) , DIAM_HYD , RO_A , rho_e(j,i) , VIS_DIN_A  , TS_VC  , Wvp  )

endif

!print*, spmt, rho_e(j,i), VIS_DIN_OIL_OUT, visc_e(j,i)

!stop
!print*, dropdiam(j,i-1)
!print*, Wvp, dropdiam(j,i-1) , DIAM_HYD , RO_A , rho_e(j,i) , VIS_DIN_A  , TS_VC

!print*,'glu', watcont(j,i), VIS_DIN_A, RO_A, spmt, diam(j,i-1), Wvp, VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
!                    dt, numdrops (j,i-1)
				
				
if  (DISSOLVE.EQ.1) THEN    !!DISSOLUTION

 !  call DISSOLVE_OIL(watcont(j,i), VIS_DIN_A, RO_A, spmt, dropdiam(j,i-1), Wvp, VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
  !                  dt, numdrops (j,i-1))

   call DISSOLVE_OIL(watcont(j,i), VIS_DIN_A, RO_A, spmt, diam(j,i-1), abs(Wvp), VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
                    dt, numdrops (j,i-1))


   massa(j,i) = sum(masscomp(j,i,:)) 
   if (massa(j,i).le.0) then
     massa(j,i) = 0
   endif


   mass_diss(j,i) =    mass_diss(j, i-1) + (massa(j,i-1) - massa(j,i)) 
   
   if (counttime/(60*60*24).gt. 15) then
	  mass_diss(j,i)=mass_diss(j, i-1)
	endif   
   mass_diss(j,i-1) = mass_diss(j,i)

!print*, 'massa', massa(j,i), Wvp
!stop

  if (DISSOLVED_FASE .EQ. 1) then

   if (parcel_dis_cont .lt. numparcels_dis) then

     parcel_dis_cont = parcel_dis_cont + 1

        xf3(parcel_dis_cont, i) = x(j,i-1)
        yf3(parcel_dis_cont, i) = y(j,i-1)
        zf3(parcel_dis_cont, i) = 0

!!print*, parcel_dis_cont, 'fase 1'


        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(parcel_dis_cont,i), LAT_partf3(parcel_dis_cont,i), PROFOUT, &
              xf3(parcel_dis_cont, i),  yf3(parcel_dis_cont, i), ZPOS, OPT)

        massaf3(parcel_dis_cont, :) = -(massa(j,i) - massa(j,i-1))    !! dissolved mass doest not change 
            
   
!      !print*, lat_part(j,i-1), LAT_partf3(parcel_dis_cont,i)
     
   endif

  endif



   massa(j,i-1) =  massa(j,i)                       !to overcome problem with entrainment

   masscomp(j, i-1, :) = masscomp(j, i, :)          !to overcome problem with entrainment


!!print*, 'masscomp', massa(j,i)

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


      aux=0
      do comps=1, NCOMP_OIL
       aux=aux + masscomp(j,i,comps)/RO_COMP_OIL(comps)
      enddo
 
    spmt=massa(j,i)/aux


  height(j,i)=massa(j,i)/(spmt*pi*((diam(j,i-1)/2)**2))    !consider initially that dissolution only affects height too

  height(j,i-1) = height(j,i)  !!update entrainment

!diam(j,i)=2*(((massa(j,i)/spmt)/(pi*height(j,i)))**0.5)  !meter
!!print*, spmt, sum(RO_COMP_OIL(:))

! !print*, 'dissolution_com_evap_conti', massa(j,i)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do comps=1, NCOMP_OIL

   MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

   FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

   V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

enddo


     call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

!!print*, 'visc', RO_A , VIS_DIN_A
!!print*, VIS_DIN_A/RO_A
!!print*, cos(pi)
   ! !print*, 'VIS_DIN_A', VIS_DIN_A


      MOL_TOT=sum(mol_comp)

       V_TOT=massa(j,i)/spmt
 
      VOL(j,i)=massa(j,i)/spmt


      moil=massa(j,i)

     FRAC_TOT=sum(FRAC_MASS_OUT)

    vol_diss(j,i)  = vol_diss(j,i-1) + ( VOL(J,I-1) -  VOL(j,i)  )
 
    porc_diss_vol(j,i) =  (  vol_diss(j,i) / volreff   ) * 100
	
	 VOL(j,i-1) = VOL(j,i)
     vol_diss(j,i-1)=vol_diss(j,i)
	 

!   !print*, temp_out  
!!print*, 1, v_tot
!!print*, 1, vol(j,i)
!!print*, 1, sum(v_comp)
! !print*, massa1(j,i)

   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

!!print*, 'aaaaaaaaa', RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!recalculate droplets info


!!!!!! definition of droplet info
   dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

 !  dropdiam(j,i-1) = dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)

   numdrops(j,i-1)= numdrops(j,i)



   if (massa(j,i).le.0) then
       massa(j,:)= 0
       lat_part(j,:)= 0
       lon_part(j,:)= 0
       dropdiam(j,:)= 0
       voldrops(j,:)= 0
       numdrops(j,:)= 0
       massae(j,:)  = 0
       area(j,:)    = 0
       areaem(j,:)  = 0
       diam(j,:)    = 0
       diamem(j,:)  = 0   
       height(j,:)  = 0
       vol(j,:)     = 0
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_sedi(j,i)=mass_sedi(j,i-1)	   
	   mass_degr(j,i)=mass_degr(j,i-1)	   
!       go to 18767      !!! go to entrainment
       cycle
   endif
!!!!!!!!!!!!!
!print*, dropdiamax, dropdiamin
!!print*, 'numdrop', j, numdrops(j,i), vol(j,i), voldrops (j,i), dropdiam(j,i), turbed, VIS_DIN_OIL_OUT, RO_OIL_OUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if 



 18767 continue



CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


      aux=0
      do comps=1, NCOMP_OIL
       aux=aux + masscomp(j,i-1,comps)/RO_COMP_OIL(comps)
      enddo
 
    spmt=massa(j,i-1)/aux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do comps=1, NCOMP_OIL

   MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

   FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

   V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

enddo

     call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

!!print*, 'visc', RO_A , VIS_DIN_A
!!print*, VIS_DIN_A/RO_A
!!print*, cos(pi)
   ! !print*, 'VIS_DIN_A', VIS_DIN_A


      MOL_TOT=sum(mol_comp)

       V_TOT=massa(j,i)/spmt
 
      VOL(j,i)=massa(j,i)/spmt

      VOL(j,i-1) = VOL(j,i)

      moil=massa(j,i)

     FRAC_TOT=sum(FRAC_MASS_OUT)



   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

!! SPMT  = RO_OIL_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!recalculate droplets info


!!!!!! definition of droplet info
 !  dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

 !  dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

  ! dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001
 ! print*, dropdiam(j,i), dropdiam(j,i-1), zf1(j,i-1)
 !  dropdiam(j,i-1) = dropdiam(j,i)
  !  print*, dropdiam(j,i), "1", RN1, checkb(j,i)

    if (checkb(j,i) .ne. 0) then
	  CALL random_number(RN1)	
!	 if ((ductwd .gt. 0) .and. (RN1 .gt. 0.5) .and. (checkb(j,i) .ne. 30)) then
	 if ((ductwd .gt. 0) .and. (RN1 .gt. 0.9) .and. (checkb(j,i) .ne. 30)) then
	   checkb(j,:) = 0
       dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
       dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
       dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001
	 else  
!	  dropdiam(j,i)=ductwd
	   dropdiam(j,:)=ductwd
	  checkb(j,:)=30
	  endif
    endif 

 !  print*, dropdiam(j,i), "fff", RN1, checkb(j,i), zf1(j,i)
!   stop

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

 !  voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)

!   numdrops(j,i-1)= numdrops(j,i)

 ! print*, numdrops(j,i), dropdiam(j,i)

! print*, 'third', dropdiam(j,i), zf1(j,i-1)
!print*, ro_oil_out,  rho_e(jj,ii)

!print*, '555',wvp, massa(j,i), ro_oil_out,  rho_e(jj,ii)

 IF (ENTRAIN.eq.1) then


 !  print*, checkb(j,i), dropdiam(j,i)

     ! if(j.eq.6) then
     !  !print*, 'caooooooooooooooooooooooooo', dropdiam(j,i-1), dropdiam(j,i), j
     !  STOP
     ! ENDIF
  ! !print*, 'ENTRAINMENT ON'
  
  ! !print*, 'jjj', temp_out, sal_a
  
 !  print*, VIS_DIN_OIL_OUT, RO_OIL_OUT, zf1(j,i-1)
 !  stop
   
   call  PROP_AMBIENTE(zf1(j,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


   
   viscst = (VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000  !!Absolute viscosity from current particle which will be the same of generated droplets(should change to emulsion)

   viscstem = ( (visc_e(jj,ii)/1000) / rho_e(jj,ii) ) * 1000000 
 !  viscst=100
 

    IF (EMULSI.EQ.1) THEN

 !   call vert_disp (viscstem, ro_A, windspms, deltad, dropdiam(j,i-1)/2, wvp, qd, zini, kz, seafrac)   !delvi
!print*,kz, wvp

      call vert_disp_li_2007 ( visc_e(jj,ii)/1000, TS_VC , ro_A,rho_e(jj,ii), windspms, qd, zini, seafrac, kz, wvp)  !lietal
      call size_distr_li_2007 ( visc_e(jj,ii)/1000, TS_VC , ro_A,rho_e(jj,ii), windspms, qd, zini, seafrac, kz, wvp)  !lietal
!     print *, 'Random samples 2:', rand_samples
	 

!	  stop

    else 


  !   call vert_disp (viscst/100, ro_A, windspms, deltad, dropdiam(j,i-1), wvp, qd, zini, kz, seafrac)

     call vert_disp_li_2007 ( VIS_DIN_OIL_OUT, TS_VC , ro_A,RO_OIL_OUT, windspms, qd, zini, seafrac, kz, wvp)  !lietal

      call size_distr_li_2007 ( VIS_DIN_OIL_OUT, TS_VC , ro_A,RO_OIL_OUT, windspms, qd, zini, seafrac, kz, wvp)  !lietal

	 ! if (rand_samples(1).eq.888898) then
        ! dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

        ! dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

        ! dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001	 
	 ! else
	   	 ! dropdiam(j,i) = rand_samples(1) 
	 ! endif
	 
    endif



 !  !print*, 10, tsevol(contind),  tsevol(contind-1)

   if(zf1(j,i-1).lt.0) then

          
      if (three_dim .ne. 1) then 
        wi=0  
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DISSOLUTION FASE 2

!     print*, 'ts f2', m1,temp_outf2, SAL_Af2
      m1=j

       call  PROP_AMBIENTE(zf1(j,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

!print*, '999', spmt, RO_OIL_OUT, rho_e(jj,ii)
	   
if (emulsi.eq.0) then

 call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )

else 

 call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , rho_e(j,i) , VIS_DIN_A  , TS_VC  , Wvp  )

endif	   
	   

  !     call VEL_ASCENSAO_BG ( dropdiam(M1,i-1) , DIAM_HYD , RO_A , spmt, VIS_DIN_A  , TS_VC  , Wvp  )

!print*, 'm1', m1, spmt

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


         do comps=1, NCOMP_OIL

           MOL_COMP(comps)=masscomp(m1,i-1, comps)* 1000.D0/PM_COMP_OIL(comps)

           FRAC_MASS_OUT(comps)=masscomp(m1,i-1, comps)/massa(m1,i-1)

           V_COMP(comps) = masscomp(m1,i-1, comps)/ RO_COMP_OIL(comps)

         enddo


         MOL_TOT=sum(mol_comp)

         V_TOT=massa(m1,i-1)/spmt
 
         vol(m1,i)=massa(m1,i-1)/spmt
		 
         VOL(J,I-1)=VOL(J,I)
		 
         moil=massa(m1,i-1)
  
         FRAC_TOT=sum(FRAC_MASS_OUT)

         CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)



        IF (DISSOLVE.EQ.1) THEN

  !		  call DISSOLVE_OIL_fase2(watcont(j,i), VIS_DIN_A, RO_A, spmt, dropdiam(M1,i-1), abs(Wvp), VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
  !                  dt, numdrops (j,i-1))
 
  ! 		  call DISSOLVE_OIL_fase2(watcont(j,i), VIS_DIN_A, RO_A, spmt, dropdiam(m1,i), abs(Wvp), VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
  !                  dt, numdrops (j,i-1))

   		  call DISSOLVE_OIL_fase2(watcont(j,i), VIS_DIN_A, RO_A, spmt, dropdiam(m1,i), abs(Wvp), VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
                    dt, numdrops (j,i)) 
       endif

!        print*,  1, massa(m1,i)
	!	print *, visc_e(jj,ii)

	   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

		! print *, visc_e(jj,ii)

! !		 stop
	       ! if (rand_samples(1).EQ.888898) then
             ! dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
             ! dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
            ! dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001	 
	       ! else
		   	! call size_distr_li_2007 ( visc_e(jj,ii)/1000, TS_VC , ro_A,rho_e(jj,ii), windspms, qd, zini, seafrac, kz, wvp)  !lietal
	   	    ! dropdiam(j,i) = rand_samples(1)
	       ! endif

		   
         massa(m1,i) = sum(masscomp(j,i,:)) 
		 
		 
		 mass_diss(j,i) =    mass_diss(j, i-1) + (massa(m1,i-1) - massa(m1,i))
		 
		 if (counttime/(60*60*24).gt. 15) then
	       mass_diss(j,i)=mass_diss(j, i-1)
	     endif 
       
	    mass_diss(j,i-1) = mass_diss(j,i)

  

 !        print*,  '2', massa(m1,i)
!		 stop

!         print*, '444444444444444444', sum(masscomp(j,i,:))
   if (massa(m1,i).le.0) then
       massa(j,:)= 0
       lat_part(j,:)= 0
       lon_part(j,:)= 0
       dropdiam(j,:)= 0
       voldrops(j,:)= 0
       numdrops(j,:)= 0
       massae(j,:)  = 0
       area(j,:)    = 0
       areaem(j,:)  = 0
       diam(j,:)    = 0
       diamem(j,:)  = 0   
       height(j,:)  = 0
       vol(j,:)     = 0
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_sedi(j,i)=mass_sedi(j,i-1)		   
	   mass_degr(j,i)=mass_degr(j,i-1)		   
 !      go to 18767      !!! go to entrainment
       cycle
   endif
 !        print*, 'forth', numdrops (j,i-1), massa(m1,i), VIS_DIN_OIL_OUT


       if (DISSOLVED_FASE .EQ. 1) then

        if (parcel_dis_cont .lt. numparcels_dis) then

         parcel_dis_cont = parcel_dis_cont + 1

         xf3(parcel_dis_cont, i) = x(m1,i-1)
         yf3(parcel_dis_cont, i) = y(m1,i-1)
         zf3(parcel_dis_cont, i) = zf1(m1,i-1)
 
!!print*, parcel_dis_cont, 'fase 2'


         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(parcel_dis_cont,i), LAT_partf3(parcel_dis_cont,i), PROFOUT, &
              xf3(parcel_dis_cont, i),  yf3(parcel_dis_cont, i), ZPOS, OPT)

         massaf3(parcel_dis_cont, :) = -(massa(m1,i) - massa(m1,i-1))   !!dissolved mass does not variate
            
   
!         !print*, lat_partf2(m1,i-1), LAT_partf3(parcel_dis_cont,i)
!         !print*, lon_partf2(m1,i-1), lon_partf3(parcel_dis_cont,i)
!         !print*, massaf2(m1,i), massaf2(m1, i-1),  massaf3(parcel_dis_cont, i) 
!         !print*, zf3(parcel_dis_cont, i)
!     !print*, 'lulu'
!stop     
       endif
 
      endif


        
         aux=0
         do comps=1, NCOMP_OIL
           aux=aux + masscomp(m1,i,comps)/RO_COMP_OIL(comps)
         enddo

         spmt=massa(m1,i)/aux
         vol(m1,i)=massa(m1,i)/spmt


    vol_diss(j,i)  = vol_diss(j,i-1) + ( VOL(J,I-1) -  VOL(j,i)  )
 
    porc_diss_vol(j,i) =  (  vol_diss(j,i) / volreff   ) * 100
	
	
 !   print*, 'fift2', m1, j, porc_diss_vol(j,i), numdrops (j,i-1), dropdiam(M1,i)
        


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF DISSOLUTION

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERACTION OF DROPLETS WITH HARD SURFACES 

   

     !!!!!!!!!!!!!!!!!!!!!!!!advection + random fase2
   if (right_random .eq. 0) then


        CALL random_number(RN1)
        CALL random_number(RN2)
        CALL random_number(RN3)


        randvert = -1. + 2.*RN3
       
         W_ALEA = randvert * ((2.D0*kz/dt)**(0.5D0))    !!  based on Reed et al., 1995
 
 
!        !print*, 'klk', temp_outf2 , SAL_Af2
              
        call  PROP_AMBIENTE(zf1(m1,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)   !!gotta update Z to zf1(m1,i-1)

!        call VEL_ASCENSAO_BG ( dropdiam(M1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )
        if (emulsi.eq.0) then

          call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )

        else 

         call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , rho_e(j,i) , VIS_DIN_A  , TS_VC  , Wvp  )


        endif	
  !     !print*, '444', spmt, spmtf2(m1,i-1)

        zf1(m1,i)=zf1(m1,i-1) + wi*dt +  W_ALEA*dt  + Wvp*dt
		!zf1(m1, i-1) = zf1(m1,i)
		go to 546

      

     else

     call  PROP_AMBIENTE(zf1(m1,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)   !!gotta update Z to zf1(m1,i-1)
!     call VEL_ASCENSAO_BG ( dropdiam(M1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )
        if (emulsi.eq.0) then

          call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )

        else 
 !         call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )
  !        dropdiam(m1,i)=10

         call VEL_ASCENSAO_BG ( dropdiam(m1,i) , DIAM_HYD , RO_A , rho_e(j,i) , VIS_DIN_A  , TS_VC  , Wvp  )

   !     print*, 'jdjdjdjdj', wvp, dropdiam(m1,i), rho_e(j,i)

        endif	
		
      do step_random=1,dt, dt_random
        CALL random_number(RN3)


        randvert = -1. + 2.*RN3
       

        zrandom = zrandom + ( randvert * ((2.D0*kz/dt_random)**(0.5D0)) ) *dt_random + Wvp*dt_random  !!  based on Reed et al., 1995
 !       zrandom = zrandom + ( randvert * ((2.D0*kz/dt_random)**(0.5D0)) ) *dt_random + Wvp*dt_random  !!  based on Reed et al., 1995

      enddo


!       print*, 'fift0', Wvp
!	  stop

        zf1(m1,i)=zf1(m1,i-1) + wi*dt +  zrandom 
        !zf1(m1, i-1) = zf1(m1,i)

     !  print*, wvp, zf1(m1,i), zf1(m1,i-1), dropdiam(m1,i) , randvert, randvert * ((2.D0*kz/dt_random)**(0.5D0)) , Wvp*dt_random
 !      print*, zf1(m1,i), zf1(m1,i-1), zrandom
 !print*, '444444', zf1(m1,i), zf1(m1,i-1), zrandom, di


       xrandom = 0 
       yrandom = 0
       zrandom = 0
       go to 546

     endif
     

!  493 continue

     
 ! print*, 'zf1', zf1(m1,i)

 ! if (zf1(m1,i) .ge. -0.1) then 
 !   zf1(m1,i) = 0
!	zf1(j, i-1) = zf1(m1,i)
!	go to 546
!  endif
  !     if ( dt_hf2 .ge.1) then

 !        !print*, 'WRITE DATA FASE2'porc_evap_cont.txt
  !       dt_hf2 = 0 
   !    endif

   endif
 
 
!  if(parcelcont.eq.1) then
!    !print*, '123'
!    go to 456

!  else

! !print*, qd, wvp * dt + (- zini)
!print*, dt_hf2, limen, limassem, massa(j,i)

!np.random.uniform(0, 1, len(self.elements.z))
!print*, 'qd', qd, zf1(j, i-1)

 !print*, VIS_DIN_OIL_OUT, RO_OIL_OUT, zf1(j,i-1)
 
 
 !print*, 'hhhhh', zf1(m1,i), zf1(m1,i-1), zrandom
   if ((qd.gt.0) .and. (zf1(j, i) .ge. 0)) then
          CALL random_number(RN1)
!		  oil_entrainment_probability = 1 - exp(-qd* 60)
		  oil_entrainment_probability = 1 - exp(-qd* dt)
!		  print*, RN1, exp(-qd* 60), dt
          if (RN1 .lt.  oil_entrainment_probability) then
!		    print*, 'rrr'
!			print*, RN1, oil_entrainment_probability
!            print*, "zini", zini
!			stop
            zf1(j, i) = - zini
	       if (rand_samples(1).EQ.888898) then
             dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
             dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
            dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001	 
	       else
!		    call size_distr_li_2007 ( visc_e(jj,ii)/1000, TS_VC , ro_A,rho_e(jj,ii), windspms, qd, zini, seafrac, kz, wvp)  !lietal
!            print*, rand_samples(1)
!			stop
	   	    dropdiam(j,i) = rand_samples(1)
	       endif
      !     print*, 	"44444", 	   rand_samples(1)
   !          dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
    !         dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )
       !      print*, 	"3333", 	   ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001	

 endif
		  
		   
!		  print*, oil_entrainment_probability, qd, RN1, TS_VC
		  
!		  stop
   endif
!  endif

 546 continue
if ( zf1(j, i) .ge. -0) then
   zf1(j,:)=0
   checkb(j,:)=0
endif

 !print*, 'qd', qd*dt*area(j,i), massa(j,i-1)
 !stop

 !      mass_degr(j,i-1) = mass_degr(j,i)
! 546 continue
 !  if (partcont .ge. 1) then

!   do m2=numpalm+1,numpalm+numpal

!    !print*, 'massa', massa(j,i), diam(j,i), area(j,
     !!!!PUT HERE PARTICLE UPDATED MASS!!!!!!!!!!!!!!!!! ALREAD PUT IN EVAPORATION




endif                                                           !!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT

    
  

!!print*, massa(j,i)

contind=contind+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 inicio interacao fundo
  ! if (turn_off_mask .eq. 1) then 

    ! go to 8867

  ! else  
     ! if(THEORETICAL.EQ.1) THEN


         ! maskt=mask(minloc(abs(y_mask(:,1)-y(j,i-1))), minloc(abs(x_mask(1,:)-x(j,i-1))))       !!!NEW MASK
      
         ! if (maskt(1,1).eq.1.) then
          ! x(j,i)=x(j,i-1)
          ! y(j,i)=y(j,i-1)                      
          ! go to 1000 
         ! endif

     ! ELSE     
                                                                                       ! !!!!!!!!!MASK WITH DELFT
! !print*, rn10, RN10, probsl
! !stop
! !print*, 'di', di
! ! print*, 'di', di, lon_model(lat_in(1), lon_in(2)),lat_model(lat_in(1), lon_in(2)-1)

       ! if (di.le.0) then
		! CALL random_number(rn10)
        ! OPT=1

! !        x(j,i)=x(j,i-1)
! !        y(j,i)=y(j,i-1)
! !        lon_part(j,i) = lon_part(j,i-1)
! !        lat_part(j,i) = lat_part(j,i-1)  
        ! if (rn10 .le.  probsl) then
         ! beached(j,:)=1
	     ! lon_part(j,i) = lon_part(j,i-1)
         ! lat_part(j,i) =  lat_part(j,i-1)		 
! !         lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
! !         lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
         ! call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
         ! OPT=2
		! ! print*, '111'
	! !	 print*, rn10, probsl
	! !	 stop
         ! go to 1000 
		! else
	     ! lon_part(j,i) = lon_part(j,i-1)
         ! lat_part(j,i) =  lat_part(j,i-1)
         ! call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
         ! OPT=2	
! !		 print*, '222'
! !		 print*, rn10, probsl, lon_part(j,i-1), lon_part(j,i-2)
! !		 stop
		 ! go to 8867
		! endif
       
! !       elseif ( (lat_model(lat_in(1)-1, lon_in(2)).eq.0) .or. (lat_model(lat_in(1)+1, lon_in(2)) .eq.0)  .or. &
! !             (lat_model(lat_in(1), lon_in(2)+1).eq.0) .or. (lat_model(lat_in(1), lon_in(2)-1) .eq.0) ) then
       ! ! elseif ( (lat_model(lat_in(1), lon_in(2)).eq.0) .or. (lat_model(lat_in(1), lon_in(2)) .eq.0)  .or. &
             ! ! (lat_model(lat_in(1), lon_in(2)).eq.0) .or. (lat_model(lat_in(1), lon_in(2)) .eq.0) ) then

        ! ! OPT=1

! ! !        x(j,i)=x(j,i-1)
! ! !        y(j,i)=y(j,i-1)
! ! !        lon_part(j,i) = lon_part(j,i-1)
! ! !        lat_part(j,i) = lat_part(j,i-1)  
        ! ! beached(j,:)=1
        ! ! lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
        ! ! lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
        ! ! call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
        ! ! OPT=2

! ! !        !print*, lon_part(j,i-1), lon_part(j,i)
! ! !        !print*, lat_part(j,i-1), lat_part(j,i)
! ! !        PRINT*, X(J,I-1), X(J,I)
! ! !        PRINT*, Y(J,I-1),Y(J,I)
! ! ! PRINT*, 'dryyyyy'
! ! !        STOP

        ! ! go to 1000 

! ! !       elseif ( (lon_model(lat_in(1)-1, lon_in(2)).eq.0) .or. (lon_model(lat_in(1)+1, lon_in(2)) .eq.0)  .or. &
! ! !             (lon_model(lat_in(1), lon_in(2)+1).eq.0) .or. (lon_model(lat_in(1), lon_in(2)-1) .eq.0) ) then
       ! ! elseif ( (lon_model(lat_in(1), lon_in(2)).eq.0) .or. (lon_model(lat_in(1), lon_in(2)) .eq.0)  .or. &
             ! ! (lon_model(lat_in(1), lon_in(2)).eq.0) .or. (lon_model(lat_in(1), lon_in(2)) .eq.0) ) then


        ! ! OPT=1

! ! !        x(j,i)=x(j,i-1)
! ! !        y(j,i)=y(j,i-1)
! ! !        lon_part(j,i) = lon_part(j,i-1)
! ! !        lat_part(j,i) = lat_part(j,i-1)  

        ! ! beached(j,:)=1
        ! ! lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
        ! ! lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
        ! ! call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
        ! ! OPT=2

! ! !        !print*, lon_part(j,i-1), lon_part(j,i)
! ! !        !print*, lat_part(j,i-1), lat_part(j,i)
! ! !        PRINT*, X(J,I-1), X(J,I)
! ! !        PRINT*, Y(J,I-1),Y(J,I)
! ! ! PRINT*, 'dryyyyy'
! ! !        STOP

        ! ! go to 1000 
       ! endif

     ! ENDIF

  ! endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 fim interacao fundo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
8867 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (zf1(m1,i) .ge. 0) then

   if (right_random .eq. 0) then

       CALL random_number(RN1)
       CALL random_number(RN2)
       CALL random_number(RN3)

       randvert = -1. + 2.*RN3
       
         U_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * COS(2.D0*PI*RN2)    
         V_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * SIN(2.D0*PI*RN2)
!         W_ALEA = randvert * ((2.D0*kz/dt)**(0.5D0))    !!  based on Reed et al., 1995


       x(j,i)=x(j,i-1) + ui(1,1)*dt + U_ALEA*dt   + (uwd * widfc)*dt 

       y(j,i)=y(j,i-1) + vi(1,1)*dt + V_ALEA*dt   + (vwd * widfc)*dt

       call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)


   else 

      do step_random=1,dt, dt_random
        CALL random_number(RN1)
        CALL random_number(RN2)
        CALL random_number(RN3)

        randvert = -1. + 2.*RN3
       
         xrandom = xrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * COS(2.D0*PI*RN2) ) *  dt_random
         yrandom = yrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * SIN(2.D0*PI*RN2) ) *  dt_random   
!         zrandom = zrandom + ( randvert * ((2.D0*kz/dt_random)**(0.5D0)) ) *dt_random  !!  based on Reed et al., 1995
      enddo

!print*, CDIF_HOR
!stop
!!print*, 'hh', uwd, vwd
!!print*, x(j,i-1) + ui(1,1)*dt +  xrandom, y(j,i-1) + vi(1,1)*dt +  yrandom
!!print*, uwd, (uwd * 0.03)*dt 
!!print*, x(j,i-1) + ui(1,1)*dt +  xrandom + (uwd * 0.03)*dt, y(j,i-1) + vi(1,1)*dt +  yrandom + (vwd * 0.03)*dt
!stop

!       x(j,i)=x(j,i-1) + 1*ui(1,1)*dt +  xrandom + ( widfc * (  (uwd * cos(-5*(pi/180)))  + (vwd * sin(-5*(pi/180))) ))*dt&
!	   + ucomp2*dt

!       y(j,i)=y(j,i-1) + 1*vi(1,1)*dt +  yrandom + ( widfc * (  (-uwd * sin(-5*(pi/180)))  + (vwd * cos(-5*(pi/180))) ))*dt&
!	   + vcomp2*dt

       x(j,i)=x(j,i-1) + 1*ui(1,1)*dt +  xrandom + ( widfc * (  (uwd * cos(0*(pi/180)))  + (vwd * sin(0*(pi/180))) ))*dt&
	   + ucomp2*dt

       y(j,i)=y(j,i-1) + 1*vi(1,1)*dt +  yrandom + ( widfc * (  (-uwd * sin(0*(pi/180)))  + (vwd * cos(0*(pi/180))) ))*dt&
	   + vcomp2*dt
       call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)

       xrandom = 0 
       yrandom = 0
       zrandom = 0

   endif
   
   
 else 
 
 
 
 
   if (right_random .eq. 0) then

       CALL random_number(RN1)
       CALL random_number(RN2)
       CALL random_number(RN3)

       randvert = -1. + 2.*RN3
       
         U_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * COS(2.D0*PI*RN2)    
         V_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * SIN(2.D0*PI*RN2)
!         W_ALEA = randvert * ((2.D0*kz/dt)**(0.5D0))    !!  based on Reed et al., 1995


       x(j,i)=x(j,i-1) + ui(1,1)*dt + U_ALEA*dt  

       y(j,i)=y(j,i-1) + vi(1,1)*dt + V_ALEA*dt  

       call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)


   else 

      do step_random=1,dt, dt_random
        CALL random_number(RN1)
        CALL random_number(RN2)
        CALL random_number(RN3)

        randvert = -1. + 2.*RN3
       
         xrandom = xrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * COS(2.D0*PI*RN2) ) *  dt_random
         yrandom = yrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * SIN(2.D0*PI*RN2) ) *  dt_random   
!         zrandom = zrandom + ( randvert * ((2.D0*kz/dt_random)**(0.5D0)) ) *dt_random  !!  based on Reed et al., 1995
      enddo

!!print*, 'hh', uwd, vwd
!!print*, x(j,i-1) + ui(1,1)*dt +  xrandom, y(j,i-1) + vi(1,1)*dt +  yrandom
!!print*, uwd, (uwd * 0.03)*dt 
!!print*, x(j,i-1) + ui(1,1)*dt +  xrandom + (uwd * 0.03)*dt, y(j,i-1) + vi(1,1)*dt +  yrandom + (vwd * 0.03)*dt
!stop

       x(j,i)=x(j,i-1) + 1*ui(1,1)*dt +  xrandom 

       y(j,i)=y(j,i-1) + 1*vi(1,1)*dt +  yrandom 

       call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)

       xrandom = 0 
       yrandom = 0
       zrandom = 0

   endif
   
   
 
 endif

! !print*, 'cool', uwd, vwd, CDIF_HOR
!print*, vol(j,i)
       lat_model_summ=abs(lat_modelm-lat_part(j,i))
       lon_model_summ=abs(lon_modelm-lon_part(j,i))

       coord_summ=lat_model_summ + lon_model_summ

       lat_inm =  minloc(coord_summ)
       lon_inm=lat_inm

       di = depth(lat_inm(1), lon_inm(2))
       coastad = coastvalue(lat_inm(1), lon_inm(2))
!print*, coastad
       if (di.le.0) then
		CALL random_number(rn10)
        OPT=1

!        x(j,i)=x(j,i-1)
!        y(j,i)=y(j,i-1)
!        lon_part(j,i) = lon_part(j,i-1)
!        lat_part(j,i) = lat_part(j,i-1)  
!        if (rn10 .le.  probsl) then                     !prob beach
!		print*, coastad
		if (coastad .gt.0) then
         beached(j,:)=1		 
		  if(coastad.le.vol(j,i)) then
		    coastvalue(lat_inm(1), lon_inm(2))=0
	!		print*, 'aa'
		  else 
		    coastvalue(lat_inm(1), lon_inm(2))= coastvalue(lat_inm(1), lon_inm(2)) - vol(j,i)
		  endif
    !   coastad = coastvalue(lat_inm(1), lon_inm(2))
		 
!	     lon_part(j,i) = lon_part(j,i-1)
!         lat_part(j,i) =  lat_part(j,i-1)		 
!         lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
!         lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
         OPT=2
		! print*, '111'
	!	 print*, rn10, probsl
	!	 stop
  !       go to 1000 
		else                  !prob beach
	     lon_part(j,i) = lon_part(j,i-1) 
         lat_part(j,i) =  lat_part(j,i-1)
         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
         OPT=2	
!		 print*, '222'
!		 print*, rn10, probsl, lon_part(j,i-1), lon_part(j,i-2)
!		 stop
	!	 go to 8867
		endif
       !
       endif
	   SEDIMENT=1
!	   print*, zf1(j,i), di
	   if ((-zf1(j,i) .ge. di) .and. (di .gt. 0)) then
	      zf1(j,i)=-di	
          if (SEDIMENT.EQ.1) then
	   mass_sedi(j,i)=massa(j,i)		   		  
       massa(j,:)= 0
       lat_part(j,:)= 0
       lon_part(j,:)= 0
       dropdiam(j,:)= 0
       voldrops(j,:)= 0
       numdrops(j,:)= 0
       massae(j,:)  = 0
       area(j,:)    = 0
       areaem(j,:)  = 0
       diam(j,:)    = 0
       diamem(j,:)  = 0   
       height(j,:)  = 0
       vol(j,:)     = 0
	   vol_diss(j,i)=vol_diss(j,i-1)
	   porc_evap(j,i)=porc_evap(j,i-1)
	   mas_evap(j,i)=mas_evap(j,i-1)
	   mass_diss(j,i)=mass_diss(j,i-1)
	   mass_degr(j,i)=mass_degr(j,i-1)
 !      go to 18767      !!! go to entrainment
       cycle		  
          endif		  
	   endif

  
1000   continue 

 !      !print*, cont, tsl

        lat_model_sum=abs(lat_model-lat_part(j,i))
        lon_model_sum=abs(lon_model-lon_part(j,i))

        coord_sum=lat_model_sum + lon_model_sum 

        lat_in =  minloc(coord_sum)
        lon_in=lat_in


        if (partcontmap_2(lat_in(1), lon_in(2)).eq.0) then
          partcontmap_2(lat_in(1), lon_in(2))=1
        endif
  

        time_sum(lat_in(1), lon_in(2)) = dt



       if (cont.eq.tsl-1) then

  !          !print*, 'entrou'
  
        lat_model_sum=abs(lat_model-lat_part(j,i))
        lon_model_sum=abs(lon_model-lon_part(j,i))

        coord_sum=lat_model_sum + lon_model_sum 

        lat_in =  minloc(coord_sum)
        lon_in=lat_in

        partcontmap(lat_in(1), lon_in(2)) = partcontmap(lat_in(1), lon_in(2)) + 1

        mass_dist(lat_in(1), lon_in(2)) = mass_dist(lat_in(1), lon_in(2)) + massa(j,i)

       endif



       !!!! HEIGHT LOOP
        lat_model_sum_height=abs(lat_height-lat_part(j,i))
        lon_model_sum_height=abs(lon_height-lon_part(j,i))
        coord_sum_height=lat_model_sum_height + lon_model_sum_height 

        lat_in_height =  minloc(coord_sum_height)
        lon_in_height=lat_in_height

        height_map(lat_in_height(1), lon_in_height(2))= height_map(lat_in_height(1), lon_in_height(2)) + vol(j,i)
!        !print*, height_map(lat_in_height(1), lon_in_height(2)), vol(j,i)
!        !print*, lat_in_height(1), lon_in_height(2), lat_part(j,i), lat_height(1,1), lat_height(239, 239)
!        !print*, lat_height(:,1)
!        !print*, lat_in_height(1), lon_in_height(2), lon_part(j,i), lon_height(1,1), lon_height(239, 239)
!        !print*, lon_height(1,:)

       !!!!!!!!!!!!!!


!!print*, i, CONT, TSL, NUMPALM, J
!PRINT*, X(:,I)
        726 continue


     dt_h_spr(j) = dt_h_spr(j) + (dt/60.)


     enddo                    !!! END PARTICLES LOOP

       
       prob_time_sum = prob_time_sum + time_sum
       time_sum = time_sum*0
     
       if (cont.eq.tsl-1) then
  
 
        probmap_2 = probmap_2 + partcontmap_2

        partcontmap_2= partcontmap_2*0

        contprob = contprob + 1

        conc = conc + mass_dist
       
        mass_dist = mass_dist*0
    
       endif


      height_map = height_map/(dx_h*dy_h)

!!print*, i


!     !print*, cont, tsl
!     !print*, tsl
!     !print*, cont
 468 continue

  !    print*, counttime_2, counttimeh, counttime_2_h, time_finish

  if (PROBABILISTIC.eq.0) then
     if (counttime_2_h.ge.time_finish) then
	!    print*, "rrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"
        go to 1547
     endif
  endif  


  if (cont.eq.tsl) then
       total_particles = total_particles + NUMTOT
   if (PROBABILISTIC.eq.0) then

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )

            aux=0
            do comps=1, NCOMP_OIL
              aux=aux +MASSCOMPREF(comps)/RO_COMP_OIL(comps)
            enddo

            spmt=sum(MASSCOMPREF)/aux

            do comps=1, NCOMP_OIL
                 MOL_COMP(comps)=MASSCOMPREF(comps)* 1000.D0/PM_COMP_OIL(comps)
                 FRAC_MASS_OUT(comps)=MASSCOMPREF(comps)/sum(MASSCOMPREF)
                 V_COMP(comps) = MASSCOMPREF(comps)/ RO_COMP_OIL(comps)
            enddo

            call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

            MOL_TOT=sum(mol_comp)

            V_TOT=sum(MASSCOMPREF)/spmt

            moil=sum(MASSCOMPREF)

            FRAC_TOT=sum(FRAC_MASS_OUT)


            CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

           time_ini_spread = ((1.45/1.14)**4.) * ( (  V_TOT / ((VIS_DIN_A/RO_A) * gravity * ((RO_A - RO_OIL_OUT &
 )/RO_A ) ) )**(1./3.)  )
         
           time_ini_spread = time_ini_spread/60


       do m=NUMTOT+1,NUMTOT+numpal
!           !print*, 1
!           CALL random_number(RN1)
!           CALL random_number(RN2)
!           x(m,i)=X0*RN1                              !!!checkout indices, they were wrong
!           y(m,i)=Y0*RN2
           CALL random_number(RN1)
           CALL random_number(RN2)
           x(m,i)=(circle_radius*RN1) * COS(2*PI*RN2) +   X0             !!NEW CIRCULAR SLICK
           y(m,i)=(circle_radius*RN1) * SIN(2*PI*RN2) +   Y0

           call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(m,i), LAT_part(m,i), PROFOUT,  x(m,i),  y(m,i), ZPOS, OPT)


           dt_h_spr(m) = time_ini_spread

           massa(m,i) = massref
           area(m,i) = 2.1*(pi)*  (   ( (( V_TOT**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * ((time_ini_spread*60) &
**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )
           diam(m,i) = ((area(m,i) / PI ) ** (1./2.)) * 2.
           height(m,i) = (V_tot) / area(m,i)
           massae(m,i) = massa(m,i)
           diamem(m,i) = diam(m,i)
           areaem(m,i) = area(m,i)
           vol(m,i) = volreff
           masscomp(m,i, :) = MASSCOMPREF
           FRAC_MASS_OUT_PART(m,:) = FRAC_MASS_OUT_REF
           visc_e(m,i) = VIS_DIN_OIL_OUT*1000   !in cP
           rho_e(m,i) = RO_OIL_OUT


!           !print*, x(m,1), y(m,1)
       enddo
       NUMTOT=NUMTOT+numpal   
       cont=0
       lo=lo+1
!       !print*, lo

!       !print*, freq
   else

     if (probabilistic_2.eq.1) then
         if (count_prob_2.eq.15) then
            end_prob='finish'
            goto 3888
         endif
     endif 


          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )
            aux=0
            do comps=1, NCOMP_OIL
              aux=aux +MASSCOMPREF(comps)/RO_COMP_OIL(comps)
            enddo

            spmt=sum(MASSCOMPREF)/aux

            do comps=1, NCOMP_OIL
                 MOL_COMP(comps)=MASSCOMPREF(comps)* 1000.D0/PM_COMP_OIL(comps)
                 FRAC_MASS_OUT(comps)=MASSCOMPREF(comps)/sum(MASSCOMPREF)
                 V_COMP(comps) = MASSCOMPREF(comps)/ RO_COMP_OIL(comps)
            enddo

            call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

            MOL_TOT=sum(mol_comp)

            V_TOT=sum(MASSCOMPREF)/spmt

            moil=sum(MASSCOMPREF)

            FRAC_TOT=sum(FRAC_MASS_OUT)


            CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)


           time_ini_spread = ((1.45/1.14)**4.) * ( (  V_TOT / ((VIS_DIN_A/RO_A) * gravity * ((RO_A - RO_OIL_OUT &
)/RO_A ) ) )**(1./3.)  )
         
           time_ini_spread = time_ini_spread/60

     do m=1,numpal
      CALL random_number(RN1)
      CALL random_number(RN2)
      x(m,i)=(circle_radius*RN1) * COS(2*PI*RN2) +   X0             !!NEW CIRCULAR SLICK
      y(m,i)=(circle_radius*RN1) * SIN(2*PI*RN2) +   Y0

      call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(m,i), LAT_part(m,i), PROFOUT,  x(m,i),  y(m,i), ZPOS, OPT)

 !   !print*, x(i,1), y(i,1), lon_part(i,1), lat_part(i,1), opt
           massa(m,i) = massref
           area(m,i) = 2.1*(pi)*  (   ( (( V_TOT**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * ((time_ini_spread*60) &
 **(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )
           diam(m,i) = ((area(m,i) / PI ) ** (1./2.)) * 2.
           height(m,i) = (V_tot) / area(m,i) 
           massae(m,i) = massa(m,i)
           diamem(m,i) = diam(m,i)
           areaem(m,i) = area(m,i)
           vol(m,i) = volreff
           masscomp(m,i, :) = MASSCOMPREF
           FRAC_MASS_OUT_PART(m,:) = FRAC_MASS_OUT_REF
           visc_e(m,i) = VIS_DIN_OIL_OUT*1000   !in cP
           rho_e(m,i) = RO_OIL_OUT

           lat_model_sum=abs(lat_model-lat_part(m,i))
           lon_model_sum=abs(lon_model-lon_part(m,i))
           coord_sum=lat_model_sum + lon_model_sum 
           lat_in =  minloc(coord_sum)
           lon_in=lat_in
           if (partcontmap_2(lat_in(1), lon_in(2)).eq.0) then
             partcontmap_2(lat_in(1), lon_in(2))=1
           endif
           time_sum(lat_in(1), lon_in(2)) = dt

     enddo 

      call evap_vars
      call evap_mass
      call emulsify_parameter

      NUMTOT = numpal
      parcelcont = 0 

      cont=0

      dt_h_spr(:) = time_ini_spread

     if (probabilistic_2.eq.1) then
 !      !print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       CALL random_number(RN5)
       count_prob_2 = count_prob_2+1
       counttime = time_lim*RN5 -dt   
       counttimeh_era = (time_lim_wind*RN5 -dt)/60.   
       counttimeh = counttime/60.

   
  !print*, rn5
     ENDIF

   endif

  endif

  cont=cont+1
 
 1547 continue
 !fmtr='f10.3'
! print*, "2", zf1(:,i)
!print*, counttimeh/(60)
  if (dt_h .ge. outp_h) then
!  if (dt_h .ge. 5) then
  ! if (cont.eq.tsl-1) then
     write(12,fmt2) x(:,i)
     write(13,fmt2) y(:,i)
     write(14,fmt2) massa(:,i)
     write(15,fmt2) diam(:,i)
     write(199,fmt2) vol(:,i)
     write(198,fmt2) temps(:,i)
     write(197,fmt2) salts(:,i) 
     write(16,fmt2) height(:,i)
     write(17,fmt2) visc_e(:,i)
     write(18,fmt1) counttimeh/(60) - counttimeh_r0/60.
     write(19,fmt2) spmt
     write(23,fmt2) porc_evap(:,i)
     write(231,fmt2) mas_evap(:,i)
     write(232,fmt2) mass_diss(:,i)
     write(233,fmt2) mass_sedi(:,i)
     write(268,fmt2) mass_degr(:,i)
     write(221,fmt2) lat_part(:,i)
     write(222,fmt2) lon_part(:,i)
     write(54,fmt2) watf(:,i) 
	 write(2254,fmt2) beached(:,i)
     write(225,fmt2) zf1(:,i)

!	 stop

    open(769,file='height_deep_55_0_time_1.txt', status='replace')
    do hh=1, num_sp*2-1
      write(769, fmt7) height_map(hh,:)
    enddo  
    close(769)

     if (dissolved_fase .eq. 1) then
       write(227,fmt4) lat_partf3(:,i)
       write(228,fmt4) lon_partf3(:,i)
       write(226,fmt4) zf3(:,i)
     endif

     if (ENTRAIN.EQ.1) THEN 
       write(223,fmt2) lat_partf2(:,i)
       write(224,fmt2) lon_partf2(:,i)
       write(331,fmt2) dropdiam(:,i)
       write(332,fmt2) spmtf2(:,i)       
       write(141,fmt2) massaf2(:,i)
     ENDIF
!     !print*, 'WRITE DATA', dt_h

!     if(counttime/(60*60) .ge. 3.84) then
!        stop
!     endif

     
     dt_h=0
   endif



!     !print*, x(:,i), y(:,i)
     counttime=counttime+dt


  if (reverse .ne. 1) then
   counttimeh = counttimeh + (dt/60.)
   counttimeh_era = counttimeh_era + (dt/60.) 
   !print*, 'reserve deactivate'
  
  else 
   counttimeh = counttimeh - (dt/60.)
   counttimeh_era = counttimeh_era - (dt/60.) 
   ! print*, 'reserve deactivate'
  endif

   !  counttimeh = counttimeh + (dt/60.)     !!!it's minute in fact

    ! counttimeh_era = counttimeh_era + (dt/60.)  

     counttimeh_coup = counttimeh_coup + (dt/60.)

     counttime_2 = counttime_2+dt

     counttime_2_h = counttime_2/(60*60)

     dt_h = dt_h + (dt/(60.*60.))

     if (dt_hf2 .ge. limen) then
       dt_hf2=0
     endif
     
     dt_hf2 = dt_hf2 + (dt/(60.*60.))
    
     height_map = height_map*0
!     !print*, 'count_era', counttimeh_era
 !    !print*, counttime/(60*60)
 !    !print*,  spmt

 !   !print*, i, CONT, TSL, NUMPALM, J
! !print*, massa(:,i)
! !print*, height(:,i)
! !print*, diam(:,i)



   if (i .eq. ts_lim) then
     lat_part(:,1) = lat_part(:,i)
     lon_part(:,1) = lon_part(:,i)
            x(:,1) = x(:,i)
            y(:,1) = y(:,i)
     vec_tdelf(1) = vec_tdelf(i)
   vec_era(1) = vec_era(i)
   massa(:,1) =  massa(:,i)
  massae(:,1) = massae(:,i)
    diam(:,1) =   diam(:,i)
  diamem(:,1) = diamem(:,i)
   rho_e(:,1) =  rho_e(:,i)
  visc_e(:,1) = visc_e(:,i)
    area(:,1) =   area(:,i)
  areaem(:,1) = areaem(:,i)
     vol(:,1) =    vol(:,i)
temps(:,1) =    temps(:,i)	 
salts(:,1) =    salts(:,i)	 
masscomp(:,1,:) = masscomp(:,i, :)
mas_evap(:,1) = mas_evap(:,i)
mass_diss(:,1) = mass_diss(:,i)
mass_sedi(:,1) = mass_sedi(:,i)
mass_degr(:,1) = mass_degr(:,i)
porc_evap(:,1) = porc_evap(:,i)
  height(:,1) =  height(:,i)
vol_evap(:,1) =  vol_evap(:,i)
porc_evap_vol(:,1) = porc_evap_vol(:,i)
dropdiam(:,1) = dropdiam(:,i)  
voldrops(:,1) = voldrops(:,i)
numdrops(:,1) = numdrops(:,i)
 watcont(:,1) =  watcont(:,i)
    watf(:,1) =     watf(:,i)
vol_diss(:,1) = vol_diss(:,i)
porc_diss_vol(:,1) =  porc_diss_vol(:,i)
  spmtf2(:,1) = spmtf2(:,i)
masscompf2(:,1,:) = masscompf2(:,i,:)
 massaf2(:,1) = massaf2(:,i)
   volf2(:,1) = volf2(:,i)
dropdiamf2(:,1) = dropdiamf2(:,i)
numdrops(:,1) = numdrops(:,i)
massadropsf2(:,1) = massadropsf2(:, i)
masscompdropf2(:,1,:) = masscompdropf2(:,i,:)
voldropsf2(:,1) = voldropsf2(:,i) 
     xf2(:,1) =  xf2(:,i)
     yf2(:,1) =  yf2(:,i)
     zf1(:,1) = zf1(:,i)
LON_partf2(:,1) = LON_partf2(:,i)
Lat_partf2(:,1) = Lat_partf2(:,i)
LON_partf3(:,1) = LON_partf3(:,i)
Lat_partf3(:,1) = Lat_partf3(:,i)
     xf3(:,1) =  xf3(:,i)
     yf3(:,1) =  yf3(:,i)
     zf3(:,1) = zf3(:,i)
    massaf3(:,1) = massaf3(:,i)
    contind=2

   !print*, 'ENTROU NO TS_LIM&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7'

 !  !print*, massa(:,i)
 !  !print*, shape(massa)
!!print*, shape(LON_partf2)

   endif

  if (cont_ts.eq. ts) then   !ts -1 causa it starts at 1
     print*, cont_ts, ts
     go to 78623
  endif

  enddo

enddo

 
3000 if (trim(pre_end).eq.'end') then
    print*, 'SIMULATION TIME GREATER THAN FORCING TIME'
  endif

3888 if (trim(end_prob).eq.'finish') then
    !print*, 'END OF PROBABILISTIC SIMULATION'
  endif

78623 continue

  probmap = (DBLE( partcontmap) / dble(total_particles)  )*100.

  probmap_2 = (DBLE(probmap_2)/DBLE(contprob))*100.

  conc=conc/contprob

!  !print*, '44444', sum(conc), massa(1,1)*500

!  conc= conc / ( 1 * area_grid )

  prob_time_sum = prob_time_sum/counttime_2

  do hh=1, numlat
   write(400,fmt5) probmap(hh,:)
  enddo

  do hh=1, numlat
   write(402,fmt6) partcontmap(hh,:)
  enddo

  do hh=1, numlat
   write(500,fmt5) probmap_2(hh,:)
  enddo


  do hh=1, numlat
   write(6789,fmt5) conc(hh,:)
  enddo


  do hh=1, numlat
   write(6667,fmt5) prob_time_sum(hh,:)
  enddo




                !print*, '77', maxval(probmap_2)

 call cpu_time(stop_time)


print *, "Simulation time:", &
      (stop_time - start_time)/60., "minutes"

!print*, 'total', total_particles
!print*, 'contprob', contprob
!print*, counttime

coastvalue=coastvalueb - coastvalue

  do hh=1, numlatm
   write(11111,fmt9) coastvalue(hh,:)
  enddo

!print last TS

     write(12,fmt2) x(:,i)
     write(13,fmt2) y(:,i)
     write(14,fmt2) massa(:,i)
     write(15,fmt2) diam(:,i)
     write(16,fmt2) height(:,i)
     write(17,fmt2) visc_e(:,i)
     write(18,fmt1) counttimeh/(60) - counttimeh_r0/60.
     write(19,fmt2) spmt
     write(23,fmt2) porc_evap(:,i)
     write(231,fmt2) mas_evap(:,i)
     write(232,fmt2) mass_diss(:,i)
     write(233,fmt2) mass_sedi(:,i)
     write(268,fmt2) mass_degr(:,i)
     write(221,fmt2) lat_part(:,i)
     write(222,fmt2) lon_part(:,i)
	 write(199,fmt2) vol(:,i)
     write(198,fmt2) temps(:,i)
     write(197,fmt2) salts(:,i) 
     write(54,fmt2) watf(:,i) 
	 write(2254,fmt2) beached(:,i)
	 write(225,fmt2) zf1(:,i)
!	 stop

    open(769,file='height_deep_55_0_time_1.txt', status='replace')
    do hh=1, num_sp*2-1
      write(769, fmt7) height_map(hh,:)
    enddo  
    close(769)

     if (dissolved_fase .eq. 1) then
       write(227,fmt4) lat_partf3(:,i)
       write(228,fmt4) lon_partf3(:,i)
       write(226,fmt4) zf3(:,i)
     endif

     if (ENTRAIN.EQ.1) THEN 
       write(223,fmt2) lat_partf2(:,i)
       write(224,fmt2) lon_partf2(:,i)
       write(331,fmt2) dropdiam(:,i)
       write(332,fmt2) spmtf2(:,i)       
       write(141,fmt2) massaf2(:,i)
     ENDIF
!     !print*, 'WRITE DATA', dt_h


!
  
8856  continue
PRINT*, 'END OF SIMULATION'
print*, widfc, CDIF_HOR
end program
  
