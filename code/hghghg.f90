program lagrange
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

  double precision:: al, kl, spm, b(2), c(100,100), areat, diamt, ui(1,1), vi(1,1), ti(1,1), si(1,1),mwf, lo,wvp, rrr
  double precision:: scft, kf, rn1, rn2, voldis 
  integer i,a,j,cont, tsl, itsl, freq, ts_lim, cont_ts, outer_l, zlayer
  integer :: ts, numpalm, turn_off_mask
  character:: g(100)
  double precision:: dt,vf, vft, x_model(100,100), y_model(100,100), u_model(100,100), x_mask(100,100), v_model(100,100)  !vf eh vel func
  double precision:: y_mask(100,100), mask(100,100), maskt(1,1)
  double precision, dimension(:, :), allocatable:: x, area, y, u,v
  double precision, dimension(:,:), allocatable :: diam
  double precision, dimension(:), allocatable:: dt_h_spr
  double precision :: time_lim, time_lim_wind
  character*20::  fmt1, x1, fmt, fmt2, fmt3, x2, fmt5, x5, fmt6, x3, fmt4, fmt7, fmt8

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

  open(141,file='massaf2.txt', status='UNKNOWN', FORM= 'FORMATTED ')

  OPEN(17 , FILE = "visc_e_CONT.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )  

  OPEN(18 , FILE = "time_hours.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' ) 

  OPEN(19 , FILE = "spmt.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
 
  OPEN(23 , FILE = "porc_evap_cont.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(54 , FILE = "water_frac.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(221 , FILE = "lat_part_inpe.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(222 , FILE = "lon_part_inpe.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(223 , FILE = "lat_partf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(224 , FILE = "lon_partf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(227 , FILE = "lat_partf3.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

  OPEN(228 , FILE = "lon_partf3.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )
 
  OPEN(225 , FILE = "zf2.txt" , STATUS= 'UNKNOWN', FORM= 'FORMATTED ' )

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



  DISSOLVE = 0    ! 0 TURN OFF DISSOLUTION
  EMULSI  = 0    !0 TURN OFF EMULSIFICATION
  EVAP_TURN = 0    !0 TURN OFF EVAPORATION
  ENTRAIN= 1
  THEORETICAL = 0   !IF THEORETICAL ACTIVATED DO NOT TURN ENTRAINMENT ON!!!!
  LOW_MEMORY  = 0     !IN THIS CASE THE PARCELS THAT THE REACH SURFACE ARE ADHERED TO NEAREST PARTICLE (USE IN SMALL PLACES)
  DISSOLVED_FASE=0

  zlayer = 1

  coupling_ind=0

  right_random = 1

  three_dim = 0

  PROBABILISTIC  = 0   !!!TURN ON in probabilistic analysis

  probabilistic_2= 0 

  wind_theoretical = 1

  linear_interp  = 1
 
  turn_off_mask = 0
!  api=36.5523D0  !!!!!!!!!!!!API Johansen Deep blow

  api= 45.6367254222577 !!!!!!!!!!!!API Johansen Deep blow


  scf=2.7  !handbook of oil spill scienc and technology

  widfc=0.030

!  windspx=15800  !m/h 15800 is 15.8 km/h
!  windspy=0
!  windsp=(windspx**2 + windspy**2)**0.5
!  windspms=windsp/3600

  if (wind_theoretical.eq.1) then
   uwd=2                              ! m/s
   vwd=3                                   ! m/s
   windsp = ((uwd**2 + vwd**2)**0.5)*3600   !wind velocity in m/h
   windspms = (uwd**2 + vwd**2)**0.5
  endif
!  print*, windsp


  numpal=100!number of particles
  mwf= 933.45 !g/mol
  dt=1800! TS always in second

  dt_random = 1

  ts = 24 !number of hours you wanna run

  time_finish = 240000  ! number of hours that lasts the continuous release

  ts= (ts*60*60)/dt 

  tsl= 240000 !minutes interval for realesing particle
  
  itsl=tsl

  tsl = (tsl*60) / dt

  freq=ts/tsl

  ts_lim=500

!   print*, freq, ts - freq*tsl

!   print*, ts
   if (coupling_ind .eq. 1) then
      go to 2534
   endif


   if (PROBABILISTIC.eq.0) then
          
     if ( (ts - freq*tsl) .eq.0) then
        print*, '00000'
        allocate(x(int(numpal*freq + numparcels),int(ts_lim)), y(int(numpal*freq + numparcels),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal*freq + numparcels
        write (x2,fmt) numparcels
        write (x3,fmt) numparcels_dis
     else 
        print*, 'not000'
        allocate(x(int(numpal*freq + numpal + numparcels),int(ts_lim)), y(int(numpal*freq + numpal + numparcels),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal*freq + numpal  + numparcels
        write (x2,fmt) numparcels
        write (x3,fmt) numparcels_dis
     endif
   else 

        allocate(x(int(numpal + numparcels),int(ts_lim)), y(int(numpal+ numparcels),int(ts_lim)))
        fmt= '(I10)'
        write (x1,fmt) numpal +  numparcels
        write (x2,fmt) numparcels
        write (x3,fmt) numparcels_dis
   endif  


  2534 continue



  zf1 = 0
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
     res_par_in (numb_lines_c))

!    print*, numb_lines_c

    open(1564, file=trim(path_ini_coup), status='old')

    DO i = 1, numb_lines_c
     READ(1564,*) res_par_in(i), time_res(i), lon_ref_res(i), lat_ref_res(i), dummy_val,  vol_res(i)
    END DO
    CLOSE (1564)    

!    print*, time_res

    allocate (FRAC_MASS_OUT(NCOMP_OIL))
    allocate(x(int(numpal*numb_lines_c + numparcels),int(ts_lim)), y(int(numpal*numb_lines_c+ numparcels),int(ts_lim)))
    fmt= '(I10)'
    write (x1,fmt) numpal*numb_lines_c +  numparcels
    write (x2,fmt) numparcels
    write (x3,fmt) numparcels_dis

    lon_ref = lon_ref_res(1)
    lat_ref = lat_ref_res(1)

!print*, lon_ref, lat_ref

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

!   print*, lat_ref, lon_ref, time_res

  endif


!print*, shape(x)

 ! print*, fmt1
 
!!  if ((ts-freq*tsl).eq.0) then
!      allocate(x(int(numpal*freq)+numpal,int(ts)), y(int(numpal*freq)+numpal,int(ts)))
!       fmt= '(I10)'
!!      write (x1,fmt) numpal*freq+pal
!      print*, trim(x1)     
!        do i=1,dt, dt_random
!            print*, i


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

!   print*, trim(x1)     
    fmt1='( '//trim(x1)//'f10.3)'    !"( 1000f10.3 )"
    fmt1=trim(fmt1)
    call StripSpaces(fmt1) 

    fmt2='( '//trim(x1)//'f16.6)'    !"( 1000f10.3 )"
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

!  print*, b


  allocate( diam(int(b(1)),int(b(2))), massa(int(b(1)),int(b(2))))

  allocate( vol(int(b(1)),int(b(2))), area(int(b(1)),int(b(2))), u(int(b(1)),int(b(2))), v(int(b(1)),int(b(2))) )

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

allocate(lat_part(int(b(1)),int(b(2))), lon_part(int(b(1)),int(b(2))) )

allocate(vec_tdelf(int(ts_lim)), vec_era(int(ts_lim)))

allocate(areaem(int(b(1)),int(b(2))), diamem(int(b(1)),int(b(2))), areadrop(int(b(1)),int(b(2))))


! if (ENTRAIN.EQ.1) THEN 
   allocate(xf2(int(numparcels),int(b(2))), yf2(int(numparcels),int(b(2))) , zf2(int(numparcels),int(b(2)))    )

   allocate(lon_partf2(int(numparcels),int(b(2))),lat_partf2(int(numparcels),int(b(2)))  )


   allocate(areaf2(int(numparcels),int(b(2))),volf2(int(numparcels),int(b(2))) , dropdiamf2(int(numparcels),int(b(2)))    )


   allocate( numdropsf2(int(numparcels),int(b(2))),voldropsf2(int(numparcels),int(b(2))) ,diamf2(int(numparcels),int(b(2)))    )


   allocate(MASSCOMPf2(int(numparcels),int(b(2)),NCOMP_OIL ), massaf2(int(numparcels),int(b(2))),spmtf2(int(numparcels),int(b(2))) )

   allocate(massadropsf2(int(numparcels),int(b(2))) , masscompdropf2(int(numparcels),int(b(2)),NCOMP_OIL )      )

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

!print*, r

!sto


!print*, 'porep', porep(20)


!PRINT*, SHAPE(MASSCOMP)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   TEMP_OUT=temp    



  if (coupling_ind .ne. 1) then

     voldis=100000 !m3

!  spm=830  !kg/m3

     vol(1:NUMPAL,:)=voldis/numpal  !volume of each particle

!   TEMP_OUT  = (273.15D0 + TEMP_OUT)   !The first calculations are performemed in Kelvin


!   CALL components( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
!	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT )   ----------------- 

!   call calc_api ( API  )            ---------------------------------------------------------------  generate oil with specific API


     CALL components( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT )


     do i=1, numpal
       FRAC_MASS_OUT_PART(i,:)=FRAC_MASS_OUT
     enddo


   endif



   call emulsify_parameter

   call entrainment_param

   call init_delft



!!!!!!INITIALIZE DELFT INFO
 open(145,file=trim(path)//'time_delft.txt',  status='old')
 read(145,*) numtime, numlat, numlon, numz
! print*, numtime, numlat, numlon, numz
 allocate(time_vec(numtime), lon_model(numlat, numlon)   ,   lat_model(numlat, numlon)) 
 allocate(u_model1(numlat, numlon, numz), u_model2 (numlat, numlon, numz), v_model1(numlat, numlon, numz),&
  v_model2(numlat, numlon, numz), t_model1(numlat, numlon, numz), t_model2 (numlat, numlon, numz), &
  s_model1(numlat, numlon, numz), s_model2 (numlat, numlon, numz)) 
  allocate(kz_model1(numlat, numlon, numz), kz_model2(numlat, numlon, numz))
 allocate(depth(numlat, numlon), lat_model_sum(numlat, numlon), lon_model_sum(numlat, numlon), &
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
  
  read(145,*) counttimeh_r

 time_lim = time_vec(numtime)*60 - itsl*60 - dt

!print*, counttimeh_r
!print*, time_vec(numtime) 
!print*, time_vec
!  print*, 'teste1',  time_lim, time_vec(numtime), numtime

 
 close(145)
 
  open(145,file=trim(path)//'lat_delft.txt',  status='old')
  open(146,file=trim(path)//'lon_delft.txt',  status='old')
  open(111,file=trim(path)//'depth.txt',  status='old')
!  open(3766,file=trim(path)//'area_grid.txt',  status='old')
!  open(7535,file=trim(path)//'winter.dry',  status='old')
  print*, shape(depth), numlon,numlat
 !print*, numlat,numlon
 do i=1,numlat
     read (145,*) (lat_model(i,j), j=1,numlon)
     read (146,*) (lon_model(i,j), j=1,numlon) 
     read (111,*) (depth(i,j), j=1,numlon) 
!     read (3766,*) (area_grid(i,j), j=1,numlon) 
 enddo
 close(145)
 close(146)

!  do i=1,110
!    read(7535,*) indexlon,indexlat
!      depth(indexlat, indexlon) = 0
!!      lat_model(indexlat, indexlon) = 0
!!      lon_model(indexlat, indexlon) = 0
! !   print*, indexlon, indexlat, lat_model(indexlat, indexlon)
!  enddo

!!!!!!!!
!do i=1,numlat
!   write(9967, '(208f16.6)') (lat_model(i,:))
!   write(9968, '(208f16.6)') (lon_model(i,:))
!   write(9969, '(208f16.6)') (depth(i,:))
!print*, 'hfhfhfh', -( 10 -  4 ) * 65 * 596 /180.D0
!print*, 'hfhfhfh', -( 10 -  4 ) * (65/180.D0) * 596 
!print*, pi
!enddo
!print*, shape(lat_model), fmt2
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

    read(687,*) counttimeh_era
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


!  print*, time_lim_wind, time_lim
!print*, 'nickyyy'
! print*, time_vec_era
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fmt= '(I10)'
    write (x5,fmt) numlon

  fmt5='( '//trim(x5)//'f15.5)'    !"( 1000f10.3 )"
  fmt5=trim(fmt5)
  call StripSpaces(fmt5)

  fmt6='( '//trim(x5)//'i7)'    !"( 1000f10.3 )"
  fmt6=trim(fmt6)
  call StripSpaces(fmt6)

! print*,  depth(50, :)

!print*, lat_model(20,:)
  sal_a=salin

  print*, temp_out,sal_a

  call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   height_value = 0.0000005

   

  if (coupling_ind .ne. 1) then

    massref = massa(1,1)
    volreff = massa(1,1)/spmt
    MASSCOMPREF = MASSCOMP(1,1,:)

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

 !   print*, areadrop(i,1) *numdrops (i,1), area(i,1), vol(1,1), massa(1,1)/spmt
  

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





!  print*, diamt
!  print*, diam(1,1)
!  print*, diamt/diam(1,1)
  

!  print*, diam
  
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
  circle_radius = 10   !radius of duct, following deltares

  if (coupling_ind .ne. 1) then

!  lon_ref=-40.278
!  lat_ref=-20.3168    !vitoria port


!  lon_ref=-40.985
!  lat_ref=-21.800   !acu port


!  lon_ref=-40.2422
!  lat_ref=-20.2986     !tubarao port

!  lon_ref=-31.5
!  lat_ref=-9     !Nordeste spill


  lon_ref=-38.9516
  lat_ref=-2.4761     !Nordeste spill north


!  if (probabilistic.eq.1) then
       circle_radius = ((1.45**2.)/1.14) * ((   ( (voldis**5.)*gravity*((RO_A - RO_OIL_OUT)/RO_A)   &               
                   ) /((VIS_DIN_A/RO_A))**2.  )**(1./12.))   
!  endif
!  do i=1,numpal
!    CALL random_number(RN1)
!    CALL random_number(RN2)    !!OLD SQUARE SLICK
!    x(i,1)=X0*RN1
!    y(i,1)=Y0*RN2
!  enddo

  do i=1,numpal
    CALL random_number(RN1)
    CALL random_number(RN2)
    x(i,1)=(circle_radius*RN1) * COS(2*PI*RN2) +   X0             !!NEW CIRCULAR SLICK
    y(i,1)=(circle_radius*RN1) * SIN(2*PI*RN2) +   Y0

    call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(i,1), LAT_part(i,1), PROFOUT,  x(i,1),  y(i,1), ZPOS, OPT)

 !   print*, x(i,1), y(i,1), lon_part(i,1), lat_part(i,1), opt
    
  enddo 
 

 endif
 ! print*, x(10,1)
 ! print*, (x_model(1,:)-x(10,1))
 ! print*, minloc(abs(x_model(1,:)-x(10,1)))
 ! print*, x_model(1,minloc(abs(x_model(1,:)-x(10,1))))

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

!print*, lat_height(:,239)

!!!!!!!!!!!!!!!!!!!!!!Variable height

! lon_ref = 4.850
! lat_ref = 65.0000
! OPT=1
! LONOUT = 4.870
!  LATOUT = 65.0000
! call INI_AMB( lon_ref, lat_ref , PROF_REF, lonout, latout, PROFOUT,  X0,  Y0, ZPOS, OPT)
! print*,  lon_ref, lat_ref, PROF_REF, lonout, latout, PROFOUT,  X0,  Y0, ZPOS, OPT
! print*, X0
!stop
 ! print*, 'pp'
 ! print*, y(10,1)
 ! print*, abs(y_model(:,1) - y(10,1))
 !   print*, y_model(minloc(abs(y_model(:,1) - y(10,1))), 1)
!  print*, y_model(1,:)


 do i=1, int(b(1))
    if (i.eq.int(b(1))) then
    write(12,'(i10)') i
    write(13,'(i10)') i
    write(14,'(i10)') i
    write(15,'(i10)') i     
    write(16,'(i10)') i 
    write(17,'(i10)') i 
    write(18,'(i10)') i
    write(19,'(i10)') i
    write(23,'(i10)') i
    write(221,'(i10)') i
    write(222,'(i10)') i
    write(54,'(i10)') i
     else 
    write(12,'(i10)', advance='no') i
    write(13,'(i10)', advance='no') i
    write(14,'(i10)', advance='no') i
    write(15,'(i10)', advance='no') i
    write(16,'(i10)', advance='no') i
    write(17,'(i10)', advance='no') i
    write(18,'(i10)', advance='no') i
    write(19,'(i10)', advance='no') i
    write(23,'(i10)', advance='no') i
    write(221,'(i10)', advance='no') i
    write(222,'(i10)', advance='no') i
    write(54,'(i10)', advance='no') i
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
  do i=1, int(numparcels)
    if (i.eq.int(numparcels)) then

       write(223,'(i10)') i
       write(224,'(i10)') i
       write(225,'(i10)') i
       write(141,'(i10)') i
    else 
       write(223,'(i10)', advance='no') i
       write(224,'(i10)', advance='no') i
       write(225,'(i10)', advance='no') i
       write(141,'(i10)', advance='no') i
    endif

  enddo
 ENDIF
 
  write(12,fmt2) x(:,1)
  write(13,fmt2) y(:,1)
  write(14,fmt2) massa(:,1)
  write(15,fmt2) diam(:,1)
  write(16,fmt2) height(:,1)
  write(221,fmt2) lat_part(:,1)
  write(222,fmt2) lon_part(:,1)
  write(54,fmt2) watf(:,1)
  write(18,fmt1) counttimeh_r/60.

 IF (ENTRAIN.EQ.1) THEN
   write(223,fmt3) lat_partf2(:,1)
   write(224,fmt3) lon_partf2(:,1)
   write(225,fmt3) zf2(:,1)
   write(141,fmt3) massaf2(:,1)
 ENDIF



counttimeh = counttimeh_r + dt/60.

counttimeh_coup = dt/60.

counttimeh_era = counttimeh_era +dt/60

!print*, 'tnh',  counttimeh

counttime = dt


     lo=0

! numpalm=numpal

 if (coupling_ind .ne. 1) then
 
   NUMTOT = NUMPAL

 endif

 cont=1

tsevol2(:) = 0


! print*, 'ts', numtot, ts
! print*, shape(x)
  
 cont_ts = 1

 call cpu_time(start_time)

dt_h=1

count_prob_2=1

do outer_l=1,1000000000

 print*, 'EXECUTION FASE'
  do i=2,ts_lim

   if (PROBABILISTIC.eq.1 .and. cont.eq.tsl) then 
     go to 468
   endif

     cont_ts = cont_ts + 1


     if (coupling_ind .eq. 1) then

       num_res_par = numtot   !I wrote before the loop to overcome the problem with droplets resurfacing

       do tc_in = inf_time, numb_lines_c
         if (counttimeh_coup .ge. time_res(tc_in)) then

            read (1112,*) dummy_val, (FRAC_MASS_OUT(frac_in), frac_in=1,NCOMP_OIL)


            numtot=numtot+numpal

            vol(num_res_par+1:numtot, :)=vol_res(tc_in)/numpal 

!  print*, 'vol', vol(:,1), counttimeh

            do coup_ct= num_res_par+1,numtot
              FRAC_MASS_OUT_PART(coup_ct,:)=FRAC_MASS_OUT
            enddo

!  print*, 'frac', FRAC_MASS_OUT_PART(3,:)

            CALL COMPONENTS_COUPLING(VAZAO_OIL_OUT , DT , TEMP_OUT , &
 	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT, numtot, num_res_par)


            call evap_mass_coupling(numtot, num_res_par)

            call emulsify_parameter_coupling(numtot, num_res_par)


            call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

            OPT = 1
            call INI_AMB( LON_REF, LAT_REF, PROF_REF, lon_ref_res(tc_in), lat_ref_res(tc_in), PROFOUT,  X0,  Y0, ZPOS, OPT)

            OPT = 2

            circle_radius = ((1.45**2.)/1.14) * ((   ( (vol_res(tc_in)**5.)*gravity*((RO_A - RO_OIL_OUT)/RO_A)   &               
                   ) /((VIS_DIN_A/RO_A))**2.  )**(1./12.))   


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

!            print*, '333', vol(numtot, 1), massa(num_res_par+1:numtot,1)/RO_OIL_OUT

            areaem(num_res_par+1:numtot,:) = area(num_res_par+1:numtot,:)


            diam(num_res_par+1:numtot,:) =  ((area(num_res_par+1:numtot,:) / PI ) ** (1./2.)) * 2.

            diamem(num_res_par+1:numtot,:) = diam(num_res_par+1:numtot,:)


            height(num_res_par+1:numtot,:) = vol(num_res_par+1:numtot, :) / area(num_res_par+1:numtot,:) 

            visc_e(num_res_par+1:numtot,:) = VIS_DIN_OIL_OUT*1000   ! in cP
  
            rho_e(num_res_par+1:numtot,:) = RO_OIL_OUT


!             print*, '3223333', diam(1,1)
           deltad = ( (dropdiamax - dropdiamin)/10. ) * 0.000001

!        print*, 'area', area(:,1)
!print*, '11111111111111111111111111111111111111111111111111111111111111111111111111111111'
!           print*, x(numtot,1), y(numtot,1), lon_part(numtot,1), lat_part(numtot,1)
!           print*,  opt, LON_REF, LAT_REF
!           print*, x0,y0, circle_radius

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

 !    print*, 'TS', cont_ts

  
     do j=1,NUMTOT

       call random_seed()
       CALL random_number(RN1)
       CALL random_number(RN2)

!       print*, x(j,i-1), y(j,i-1)
!       print*, x_model(1,:)

!       print*, y_model(:,1)
!       print*, x_model(1,minloc(abs(x_model(1,:)-x(j,i-1))))


        if ( height(j,i-1) .gt. height_value) then


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


!print*, '55', VOL(j,i), VOL(j,i-1)
!print*, '111111', vis_din_oil_out
          if (EMULSI.eq.1) then

             areaem(j,i) = 2.1*(pi)*  (   ( ((  (massae(j,i-1)/rho_e(j,i-1))**2.)*gravity*((RO_A - &
 rho_e(j,i-1))/RO_A) * ((dt_h_spr(j)*60)**(3./2.)))/ (( (visc_e(j,i-1)/1000) / rho_e(j,i-1) )**(1./2.))  )**(1./3.)  )


            height(j,i) = (massae(j,i-1)/rho_e(j,i-1)) / areaem(j,i)  !delft3d

            area(j,i) = VOL(j,i-1) / height(j,i)


!   print*, 'kgjd5', watf(j,i), areaem(j,i), (massae(j,i)/rho_e(j,i))/height(j,i)


          else

  
             area(j,i) = 2.1*(pi)*  (   ( ((VOL(j,i-1)**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * ((dt_h_spr(j)* &
 60)**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )

             height(j,i) = VOL(j,i-1) / area(j,i)  !delft3d

          endif
  

  
           area(j,i-1) = area(j,i)

           diam(j,i) = ((area(j,i) / PI ) ** (1./2.)) * 2.

           diam(j,i-1) = diam(j,i)

           height(j,i-1) = height(j,i)

!     print*, '1', height(j,i), j

          else
        
            area(j,i) = area(j,i-1)

            diam(j,i) = diam(j,i-1)

            height(j,i) = height(j,i-1)
 
        endif
!        print*, y_model(minloc(abs(y_model(:,1) - y(j,i-1))), 1)

  tsevol2(contind) = i
  tsevol(contind) = i



  print*, 'PARTICLE', J, 'TS', cont_ts, outer_l, i


  IF (THEORETICAL.EQ.1) THEN
       ui=u_model(minloc(abs(y_model(:,1)-y(j,i-1))), minloc(abs(x_model(1,:)-x(j,i-1))))


       vi=v_model(minloc(abs(y_model(:,1)-y(j,i-1))), minloc(abs(x_model(1,:)-x(j,i-1))))

       ui=0
       vi=0.05
       kz=0.001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DELFT
  ELSE


    if (tsevol2(contind) .ne. tsevol2(contind-1)) then                   ! Condition to read this part only in the first particle

       tdelf1=minloc(abs(time_vec-counttimeh))

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
!      print*, 'a', tera1, tera2, time_vec_era(tera1(1)), vec_era(i), counttimeh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       write(bas,  '(i20)') tdelf1
       write(bas2, '(i20)') tdelf2




!       print*, 'y',trim(ADJUSTL(bas_era)), trim(ADJUSTL(bas2_era))
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! print*, trim(path)//trim(ADJUSTL(bas))//'u.txt'
  !print*, trim(path)//trim(ADJUSTL(bas2))//'u.txt'
! PRINT*, VEC_TDELF(I), VEC_TDELF(I-1), 'jgjgj'
 !  print*, 'BAS', TRIM(ADJUSTL(BAS)), TRIM(ADJUSTL(BAS2))


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
!  print*, trim(path)//trim(ADJUSTL(bas_era))//'u_wind.txt'
!  print*, trim(path)//trim(ADJUSTL(bas2_era))//'u_wind.txt'
!  print*, trim(path)//trim(ADJUSTL(bas_era))//'v_wind.txt'
!  print*, trim(path)//trim(ADJUSTL(bas2_era))//'v_wind.txt'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print*, vec_tdelf(i), vec_tdelf(i-1), counttimeh


!print*, 'c'

! lat_in = minloc(abs(lat_model - lat_part(j,i-1)))

! lon_in = minloc(abs(lon_model - lon_part(j,i-1)))

!print*, lat_model, lat_part(j,i-1), lon_part(j,i-1), lon_in, lat_in
!print*, lon_in, lat_in

       lat_model_sum=abs(lat_model-lat_part(j,i-1))
       lon_model_sum=abs(lon_model-lon_part(j,i-1))

       coord_sum=lat_model_sum + lon_model_sum

       lat_in =  minloc(coord_sum)
       lon_in=lat_in


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


!print*, 'wind', vwind2

       else

         85645 continue

         uwind1 = u10_era1 (lat_in_era(1), lon_in_era(2))
         uwind2 = u10_era2 (lat_in_era(1), lon_in_era(2))
         vwind1 = v10_era1 (lat_in_era(1), lon_in_era(2))
         vwind2 = v10_era2 (lat_in_era(1), lon_in_era(2))

!print*, 'windididi', vwind2
    
       endif

     endif


!!!!!!!!!!!!!!!!!!!!!!11
       if (zlayer .eq. 1) then

         deplevel = levelsd

! print*, deplevel

       else
         deplevel = depth(lat_in(1), lon_in(2))  * levelsd

       endif 


       indexz = minloc(abs(deplevel-zf1))



!numz=3
!deallocate(deplevel)
!allocate(deplevel(5))
!zf1=-3.5
!deplevel=(/0,-2,-3,-5,-6/)
!indexz = minloc(abs(deplevel-zf1))

!print*, 'cachara', zf1, deplevel, numz, levelsd
!print*, counttimeh_r, counttimeh

!PRINT*, indexz

!
!lon_in(2) = 1
!if (   (indexz(1).eq.1 .and. zf1.ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz .and. zf1.le. deplevel(indexz(1)))  .or. (zf1.eq.deplevel(indexz(1))) ) then
!   print*, 'bababa',indexz, numz, deplevel(indexz(1)), zf1
!endif

!print*, 'sdfdfg', zf1, indexz, deplevel

!

  if(linear_interp.eq.1) then


!print*, 'aaa11444444444444444444'

     if ( (lat_in(1)+slic .gt. numlat) .or. (lat_in(1)- slic .lt. 1)  .or.  (lon_in(2)+slic .gt. numlon) .or. &
 (lon_in(2)- slic .lt. 1)    ) then
        go to 85641
     endif


     slic_lat = lat_model(lat_in(1)-slic : lat_in(1)+slic, lon_in(2) - slic : lon_in (2) + slic )
     slic_lon = lon_model(lat_in(1)-slic : lat_in(1)+slic, lon_in(2) - slic : lon_in (2) + slic )

     if (   (indexz(1).eq.1 .and. zf1.ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz .and. zf1.le. &
 deplevel(indexz(1)))  .or. (zf1.eq.deplevel(indexz(1))) .or. (size(deplevel).eq.1) ) then

!print*, 'abrac4', zf1, deplevel(indexz(1)), indexz
!print*, deplevel, shape(u_model1)

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

!  print*, 'luluaaaa', vi2

     else

        if (zf1.gt.deplevel(indexz(1))) then
           vert_index1 = indexz(1)
           vert_index2 = indexz(1) - 1
        else 

           vert_index1 = indexz(1)
           vert_index2 = indexz(1) + 1

        endif



       call init_interpolation(2, lat_in(1), lon_in(2), vert_index1, vert_index2)



        inter_depth(1) = deplevel(vert_index1)
        inter_depth(2) = deplevel(vert_index2)
        par_dep(1) = zf1

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

!   print*, 'luluzao', vi2
!   print*, inter_depth
!   print*, zf1
!   print*, u2_1d
!   print*, inter_depth
!   print*, u2_1d_2 
  
     endif

!print*, si1,si2, ti1, ti2

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

!  print*, 'lulu1', vi2

  endif


       di  = depth(lat_in(1), lon_in(2))

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


       ! print*, intertime_era
       ! print*, vel_interp_era
!         print*, ts_vec_w, counttimeh, counttimeh_era
       ! print*, ui_vec


         vel_interp_era(1) = vwind1
         vel_interp_era(2) = vwind2

         call pwl_value_1d (size(intertime_era), intertime_era, vel_interp_era, 1, ts_vec_w , ui_vec )
         vwd= ui_vec(1)
    
         windsp = ((uwd**2 + vwd**2)**0.5)*3600   !wind velocity in m/h
         windspms = (uwd**2 + vwd**2)**0.5


!print*, uwd, vwd, u10_era1, u10_era2
!print*, intertime_era, ts_vec, counttimeh
       ! print*, intertime_era
       ! print*, vel_interp_era
       ! print*, ts_vec
       ! print*, ui_vec     
       !print*, uwd, vwd
      endif

!       print*, 'caos', uwd, vwd, ts_vec_w, counttimeh_era
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
       ui= ui_vec(1)

!print*, vec_tdelf(i), vec_tdelf(i-1)
! print*, intertime, ts_vec
!print*, vel_interp, ui_vec
 !print*, tdelf1(1), tdelf2(1)

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
       
print*, si, ti
       
       vel_interp(1) = wi1
       vel_interp(2) = wi2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       wi=ui_vec(1)

       vel_interp(1) = kz1
       vel_interp(2) = kz2
       call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )

       kz=ui_vec(1)
 !ui=-1
 !vi=-1
      if (three_dim .ne. 1) then
         kz=0.001
      endif
! print*, 'OWOWOWOW', ui, vi
! print*, 'index', lat_in, lon_in
! print*, 'ooooooo', ui1,ui2,ui, intertime, ts_vec
 ! print*, 'UI', di

! print*, 'juju', lat_in, lat_part(j,i-1)
! print*, lon_in, lon_part(j,i-1)

! print*, ui1, ui2, vi1, vi2
!print*, lat_in(1), lon_in(2)

!print*, (time_vec(tdelf1(1)), time_vec(tdelf2(1))) , tdelf1(1)
! print*, lat_model(:,400)

!print*, 'dlldldldld'
!print*, intertime, vel_interp, 1, ts_vec , ui_vec
!print*, size(intertime)
!print*, u_model2(1,:)
!

  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111DELFT

!       print*, ui, vi
!       print*, ui
 !      print*, minloc(abs(x_model(1,:)-x(j,i-1)))
 !      print*, minloc(abs(y_model(:,1)-y(j,i-1)))
 !      print*, x_model(1,:)-x(j,i-1)
       wanglex=wanglex*(pi/180)
       wangley=wangley*(pi/180)

       ow=(beta*exp(((-10.**(-8))*((windspx/3600)**3))/(10.*0.1)))    


!       uwd=0.03*(windspx/3600)*sin((wanglex-180+ow)*(3.14/180.))
 !      vwd=0.03*(windspx/3600)*cos((wangley-180+ow)*(3.14/180.))

 !      print*, ow
!       print*, beta
!       print*, exp(((-0.38*(-10.**(-8)))*((windspx/3600)**3))/(10.*0.1))
!!       print*, uwd
!       print*, vwd
!       print*, windspx/3600
      
!print*, 111111111111111111, massa1(j,i)/massa(1,1)
!print*, sum(FRAC_MASS_OUT)
!print*, evapmass(j,2)
 if (dissolved_fase .eq. 1) then   !! BEGIN FASE 3

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
       dif3  = depth(lat_in_f3(1), lon_in_f3(2))

        if (dif3.le.0) then
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         xf3(m1_f3,i) = xf3(m1_f3,i-1)
         yf3(m1_f3,i) = yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
 !       print*, 'DI', di               
         cycle

        elseif(-zf3(m1_f3,i-1) .ge. dif3) then
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         xf3(m1_f3,i) = xf3(m1_f3,i-1)
         yf3(m1_f3,i) = yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         cycle
       
        elseif ( (lat_model(lat_in_f3(1)-1, lon_in_f3(2)).eq.0) .or. (lat_model(lat_in_f3(1)+1, lon_in_f3(2)) .eq.0)  .or. &
                 (lat_model(lat_in_f3(1), lon_in_f3(2)+1).eq.0) .or. (lat_model(lat_in_f3(1), lon_in_f3(2)-1) .eq.0) ) then
         xf3(m1_f3,i)=xf3(m1_f3,i-1)
         yf3(m1_f3,i)=yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         cycle

        elseif ( (lon_model(lat_in_f3(1)-1, lon_in_f3(2)).eq.0) .or. (lon_model(lat_in_f3(1)+1, lon_in_f3(2)) .eq.0)  .or. &
                 (lon_model(lat_in_f3(1), lon_in_f3(2)+1).eq.0) .or. (lon_model(lat_in_f3(1), lon_in_f3(2)-1) .eq.0) ) then
         xf3(m1_f3,i)=xf3(m1_f3,i-1)
         yf3(m1_f3,i)=yf3(m1_f3,i-1)
         lon_partf3(m1_f3,i) = lon_partf3(m1_f3,i-1)
         lat_partf3(m1_f3,i) = lat_partf3(m1_f3,i-1)  
         zf3(m1_f3,i) = zf3(m1_f3,i-1)
         cycle
        endif



  if(linear_interp.eq.1) then

     if ( (lat_in_f3(1)+slic .gt. numlat) .or. (lat_in_f3(1)- slic .lt. 1)  .or.  (lon_in_f3(2)+slic .gt. numlon) &
  .or.  (lon_in_f3(2)- slic .lt. 1)    ) then
        go to 85698
     endif


     slic_lat = lat_model(lat_in_f3(1)-slic : lat_in_f3(1)+slic, lon_in_f3(2) - slic : lon_in_f3(2) + slic )
     slic_lon = lon_model(lat_in_f3(1)-slic : lat_in_f3(1)+slic, lon_in_f3(2) - slic : lon_in_f3(2) + slic )

     if (   (indexz(1).eq.1 .and. zf3(m1_f3,i-1).ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz .and. &
 zf3(m1_f3,i-1).le. deplevel(indexz(1)))  .or. (zf3(m1_f3,i-1).eq.deplevel(indexz(1))) .or. (size(deplevel).eq.1)) then

!print*, 'abrac4f3', zf3(m1_f3,i-1), deplevel(indexz(1)), indexz, levelsd, size(deplevel)


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

!  print*, 'luluaaaaf3', ui2

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

!  print*, 'lulu1f3', ui2

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
              print*, '00000a'
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


!print*, 'windf3', vwind2

        else

          85634 continue

          uwind1 = u10_era1 (lat_in_era(1), lon_in_era(2))
          uwind2 = u10_era2 (lat_in_era(1), lon_in_era(2))
          vwind1 = v10_era1 (lat_in_era(1), lon_in_era(2))
          vwind2 = v10_era2 (lat_in_era(1), lon_in_era(2))

!print*, 'windiif3', vwind2
    
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

         xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt + U_ALEA*dt   +  ( widfc * (  (uwdf3 * cos(15*(pi/180))) &
 + (vwdf3 * sin(15*(pi/180))) )  )   *dt 

         yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt + V_ALEA*dt   +  ( widfc * (  (-uwdf3 * sin(15*(pi/180))) &
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


        xf3(m1_f3,i)=xf3(m1_f3,i-1) + uif3*dt +  xrandom +  ( widfc * (  (uwdf3 * cos(15*(pi/180)))  + &
 (vwdf3 * sin(15*(pi/180))) )  )   *dt 

        yf3(m1_f3,i)=yf3(m1_f3,i-1) + vif3*dt +  yrandom +  ( widfc * (  (-uwdf3 * sin(15*(pi/180)))  + &
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
 
  !     print*, '444', spmt, spmtf2(m1,i-1)


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

     
     enddo
 
    endif

  endif


 endif   !!END FASE 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EVAPORATION SECTION
temp_out=ti(1,1) + 273.15
temp=temp_out
SAL_A=si(1,1)

#print*, 'jfjf', ti, temp_out, sal_a, si, temp

ii=i
      
jj=j


aux=0

 if (massa(j,i-1).eq.0) then
   go to 18767     !!! go to entrainment
 endif 


          CALL vapour_pressure(TEMP_OUT)


 IF (EVAP_TURN.EQ.1) THEN
       do comps=1, NCOMP_OIL

           if (masscomp(j,i-1, comps).eq.0) then
              cycle
           endif

           kf1=(0.0292*(windsp**0.78)*(diam(j,i-1)**(-0.11))*(scf**(-0.67))*(((PM_COMP_OIL(comps)+29.)/(PM_COMP_OIL(comps)))**0.5))  !unit grams, hours

           evap1=(((((kf1*PC(comps)*area(j,i-1))/(r*temp))*(PM_COMP_OIL(comps)*FRAC_MASS_OUT_PART(j,comps)))*0.001)/3600.)*dt   !0.001 is to convert from grams to kg    !betancout et al 2005+
 
            if (evap1.gt.evapmass(j, comps)) then
               evap1 = evapmass(j, comps)
            endif
  !         print*, evap1, masscomp(j,i-1, comps), evapmass(j, comps)

  !           if (evap1.gt.masscomp(j,i-1, comps)) then
            !print*, PC(comps)
  !               evap1 = masscomp(j,i-1, comps)
  !           endif

            masscomp(j,i, comps)=masscomp(j,i-1, comps)-evap1

            if (masscomp(j,i, comps).le.0) then
               masscomp(j,i:, comps)=0
            endif
        
   !        evapmass(j, comps)=evapmass(j, comps)-evap1   wrong

            evapmass(j, comps) = masscomp(j,i, comps)

            evapmass(j, comps) = ((maxwf-watf(j,i-1) )/maxwf) *  masscomp(j,i,comps)           !!!!EMULSIFICATION INTERFERENCE ON EVAPORATION

   !       print*, watf(j,i-1), maxwf, evapmass(j, 10), masscomp(j,i,10), evap1
   
            aux=aux+evap1

       enddo

       mas_evap(j,i)  = mas_evap(j,i-1) + aux



!print*, 'evapmass1', evapmass(23), masscomp(j,i, 23)
!print*, '666', sum(masscomp(j,i, :))
!!!!!!!!!!!!!!!!!!!!EVAPORATION SECTION
!PRINT*, 'MASSA', windsp, vapp, massa(j,i), comps,area(j,i-1), kf1, temp, r, dt


massa(j,i)=sum(masscomp(j,i, :))


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
       go to 18767     !!! go to entrainment
   endif

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

!   print*, temp_out  
!print*, 1, v_tot
!print*, 1, vol(j,i)
!print*, 1, sum(v_comp)

! print*, massa1(j,i)

      VOL(j,i-1) = VOL(j,i)

   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT,TS_A, TS_VC )

!!!!!!!!!!!!!!!!!!!!recalculate droplets info


!!!!!! definition of droplet info
   dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

   dropdiam(j,i-1) =  dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)
   
   numdrops(j,i-1)= numdrops(j,i)

!!!!!!!!!!!!!

!print*,  'numdrop', numdrops(j,i), vol(j,i), voldrops (j,i)

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

   dropdiam(j,i-1) =  dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)
   
   numdrops(j,i-1)= numdrops(j,i)


!print*, 'aaaaaaaaa', RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT

endif


 IF (EMULSI.EQ.1) THEN


   CALL EMULSIFY(DT, WINDSPMS, XAM(j), XWM(j) )            !!!!!!!!!!TURN ON EMULSIFICATION

!   watcont(j,i)=watcont(j-1,i) + watup + watout

   watcont(j,i)=massa(j,i-1)*watf(j,i)

!   print*, massa(j,i), massae(j,i-1)


   massae(j,i)=massa(j,i)+watcont(j,i)

   massae(j,i-1) = massae(j,i)    
 !  call visc_emulsion(VIS_DIN_OIL_OUT)

   call visc_emulsion(VIS_DIN_OIL_OUT*1000)  !in cP

   call rho_emulsion(spmt)


   diamem(j,i)=2*(((massae(j,i)/rho_e(j,i))/(pi*height(j,i)))**0.5)  !meter

   diamem(j,i-1) = diamem(j,i)

   areaem(j,i)=pi*((diamem(j,i)/2)**2)

   areaem(j,i-1) = areaem(j,i)


 ! print*, 'kgjd5', areaem(j,i), (massae(j,i)/rho_e(j,i))/height(j,i)

ENDIF 

!print*, watcont(j,i), WINDSPMS,dt

!print*, 'API', (141.5D0 / (spmt / 1000.D0)) - 131.5D0

   !!!
!   abc=(1-0.658*(watf*0.65))*abc   !!!relation emulsification / evaporation
!   print*, 'skhkjlsdkgjhldkfg', abc 
!!!!

!   print*, watcont(j,i)
  

!   print*, 1111111, VIS_DIN_OIL_OUT


     
!     print*, FRAC_MASS_OUT(1)
!     print*, mol_tot
!      print*, vol(1,1)
!      print*, v_tot
!       print*, FRAC_MASS_OUT(3)
!     print*, 111,sum(v_comp)
 !    print*, spmt
 !    print*, massa(1,1)
 !    print*, massa(j,i)



 DIAM_HYD=0



!  call VEL_ASCENSAO_BG ( diam(j,i) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )        !vertical velocity

call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

call VEL_ASCENSAO_BG ( dropdiam(j,i-1) , DIAM_HYD , RO_A , spmt , VIS_DIN_A  , TS_VC  , Wvp  )




if  (DISSOLVE.EQ.1) THEN    !!DISSOLUTION

 !  call DISSOLVE_OIL(watcont(j,i), VIS_DIN_A, RO_A, spmt, dropdiam(j,i-1), Wvp, VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
  !                  dt, numdrops (j,i-1))

   call DISSOLVE_OIL(watcont(j,i), VIS_DIN_A, RO_A, spmt, diam(j,i-1), Wvp, VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
                    dt, numdrops (j,i-1))


   massa(j,i) = sum(masscomp(j,i,:)) 



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
       go to 18767      !!! go to entrainment
   endif

  if (DISSOLVED_FASE .EQ. 1) then

   if (parcel_dis_cont .lt. numparcels_dis) then

     parcel_dis_cont = parcel_dis_cont + 1

        xf3(parcel_dis_cont, i) = x(j,i-1)
        yf3(parcel_dis_cont, i) = y(j,i-1)
        zf3(parcel_dis_cont, i) = 0

!print*, parcel_dis_cont, 'fase 1'


        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(parcel_dis_cont,i), LAT_partf3(parcel_dis_cont,i), PROFOUT, &
              xf3(parcel_dis_cont, i),  yf3(parcel_dis_cont, i), ZPOS, OPT)

        massaf3(parcel_dis_cont, :) = -(massa(j,i) - massa(j,i-1))    !! dissolved mass doest not change 
            
   
!      print*, lat_part(j,i-1), LAT_partf3(parcel_dis_cont,i)
     
   endif

  endif



   massa(j,i-1) =  massa(j,i)                       !to overcome problem with entrainment

   masscomp(j, i-1, :) = masscomp(j, i, :)          !to overcome problem with entrainment


!print*, 'masscomp', massa(j,i)

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
!print*, spmt, sum(RO_COMP_OIL(:))

! print*, 'dissolution_com_evap_conti', massa(j,i)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do comps=1, NCOMP_OIL

   MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

   FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

   V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

enddo


     call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

!print*, 'visc', RO_A , VIS_DIN_A
!print*, VIS_DIN_A/RO_A
!print*, cos(pi)
   ! print*, 'VIS_DIN_A', VIS_DIN_A


      MOL_TOT=sum(mol_comp)

       V_TOT=massa(j,i)/spmt
 
      VOL(j,i)=massa(j,i)/spmt

      VOL(j,i-1) = VOL(j,i)

      moil=massa(j,i)

     FRAC_TOT=sum(FRAC_MASS_OUT)

    vol_diss(j,i)  = vol_diss(j,i-1) + ( VOL(J,I-1) -  VOL(j,i)  )
 
    porc_diss_vol(j,i) =  (  vol_diss(j,i) / volreff   ) * 100

!   print*, temp_out  
!print*, 1, v_tot
!print*, 1, vol(j,i)
!print*, 1, sum(v_comp)
! print*, massa1(j,i)

   CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)

!print*, 'aaaaaaaaa', RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!recalculate droplets info


!!!!!! definition of droplet info
   dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

   dropdiam(j,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

   dropdiam(j,i-1) = dropdiam(j,i)

   voldrops(j,i) = ((dropdiam(j,i)/2.)**3.) *  pi * (4./3.)

   voldrops(j,i-1) = voldrops(j,i)  

   numdrops(j,i) =  vol(j,i) / voldrops (j,i)

   numdrops(j,i-1)= numdrops(j,i)

!!!!!!!!!!!!!

!print*,  'numdrop', numdrops(j,i), vol(j,i), voldrops (j,i)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if 

 18767 continue



 IF (ENTRAIN .eq. 1) then


     ! if(j.eq.6) then
     !  print*, 'caooooooooooooooooooooooooo', dropdiam(j,i-1), dropdiam(j,i), j
     !  STOP
     ! ENDIF
  ! print*, 'ENTRAINMENT ON'
  
   print*, 'jjj', temp_out, sal_a

   call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)


  ! print*, '555', VIS_DIN_OIL_OUT, massa(j,i-1)
   
   viscst = (VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000  !!viscosity from current particle which will be the same of generated droplets(should change to emulsion)

   viscstem = ( (visc_e(jj,ii)/1000) / rho_e(jj,ii) ) * 1000000 
 !  viscst=100


  if (parcelcont.NE.numparcels) then
    IF (EMULSI.EQ.1) THEN

!    call vert_disp (viscstem, ro_A, windspms, deltad, dropdiam(j,i-1)/2, wvp, qd, zini, kz, seafrac)   !delvi

      call vert_disp_li_2007 ( visc_e(jj,ii)/1000, TS_VC , ro_A,rho_e(jj,ii), windspms, qd, zini, seafrac)  !lietal

!   print*, '444', VIS_DIN_OIL_OUT, visc_e(jj,ii)/1000, ((visc_e(jj,ii)/1000)/rho_e(jj,ii))*1.D6, rho_e(jj,ii)
!   print*,  'WATER', watf(jj,ii)
    else 


     call vert_disp (viscst, ro_A, windspms, deltad, dropdiam(j,i-1), wvp, qd, zini, kz, seafrac)


 !      call vert_disp_li_2007 (VIS_DIN_OIL_OUT, TS_VC , ro_A,RO_OIL_OUT, windspms, qd, zini,seafrac)

    endif

!   print*, RO_OIL_OUT, VIS_DIN_OIL_OUT 
!   print*, qd, dropdiam(j,i-1), moil, massref
  endif

  ! print*, viscst, ro_A, windspms, deltad, dropdiam(j,i-1), zini, wvp, kz
  ! stop

 !  print*, wvp * dt + (- zini)
  !print*, parcelcont



 if (tsevol(contind) .ne. tsevol(contind-1)) then

 !  print*, 10, tsevol(contind),  tsevol(contind-1)

   if(parcelcont.gt.0) then

    do m1=1,parcelcont
     
       if ((LON_partf2(M1,i-1).eq.0) .or. (Lat_partf2(M1,i-1).eq.0)) then
          cycle
       endif

      lat_model_sum=abs(lat_model-lat_partf2(m1,i-1))
      lon_model_sum=abs(lon_model-lon_partf2(m1,i-1))
      coord_sum=lat_model_sum + lon_model_sum
      lat_in_f2 =  minloc(coord_sum)
      lon_in_f2=lat_in_f2

       if (zlayer .eq. 1) then

         deplevel = levelsd
       else

        deplevel = depth(lat_in_f2(1), lon_in_f2(2))  * levelsd

       endif 


      indexz = minloc(abs(deplevel-zf2(m1,i-1)))
      dif2  = depth(lat_in_f2(1), lon_in_f2(2))
      
  if(linear_interp.eq.1) then

     print*,'fase 2'

     if ( (lat_in_f2(1)+slic .gt. numlat) .or. (lat_in_f2(1)- slic .lt. 1)  .or.  (lon_in_f2(2)+slic .gt. numlon) &
  .or.  (lon_in_f2(2)- slic .lt. 1)    ) then
        go to 45364
     endif


     slic_lat = lat_model(lat_in_f2(1)-slic : lat_in_f2(1)+slic, lon_in_f2(2) - slic : lon_in_f2(2) + slic )
     slic_lon = lon_model(lat_in_f2(1)-slic : lat_in_f2(1)+slic, lon_in_f2(2) - slic : lon_in_f2(2) + slic )

     if (   (indexz(1).eq.1 .and. zf2(m1,i-1).ge. deplevel(indexz(1))) .or.  (indexz(1).eq.numz .and. &
 zf2(m1,i-1).le. deplevel(indexz(1)))  .or. (zf2(m1,i-1).eq.deplevel(indexz(1))) .or. (size(deplevel).eq.1)) then

!print*, 'abrac4', zf2(m1,i-1), deplevel(indexz(1)), indexz, levelsd, size(deplevel)


    call init_interpolation(1, lat_in_f2(1), lon_in_f2(2), indexz(1), 0)


 !   call pwl_interp_2d ( size(lon_1d), size(lat_1d), lon_1d, lat_1d, slic_u1,1, inlon, inlat, out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        ui1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        ui2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_partf2(m1,i-1), &
 lat_partf2(m1,i-1), out_int )
        vi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        vi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        wi1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        wi2 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        kz1 = out_int(1)
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        kz2 = out_int(1)

!  print*, 'luluaaaa', ui2

     else

        if (zf2(m1,i-1).gt.deplevel(indexz(1))) then
           vert_index1 = indexz(1)
           vert_index2 = indexz(1) - 1
        else 

           vert_index1 = indexz(1)
           vert_index2 = indexz(1) + 1

        endif

       call init_interpolation(2, lat_in_f2(1), lon_in_f2(2), vert_index1, vert_index2)


        inter_depth(1) = deplevel(vert_index1)
        inter_depth(2) = deplevel(vert_index2)
        par_dep(1) = zf2(m1,i-1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d , power, 1, lon_partf2(m1,i-1), &
 lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u1_1d_2 , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, u2_1d_2 , power, 1, lon_partf2(m1,i-1), &
 lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        ui2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v1_1d_2 , power, 1, lon_partf2(m1,i-1), &
 lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d , power, 1, lon_partf2(m1,i-1), &
 lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, v2_1d_2 , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        vi2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d , power, 1, lon_partf2(m1,i-1),&
 lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w1_1d_2 , power, 1, lon_partf2(m1,i-1),&
 lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d , power, 1, lon_partf2(m1,i-1),&
 lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, w2_1d_2 , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        wi2 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d , power, 1, lon_partf2(m1,i-1),&
 lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz1_1d_2 , power, 1, lon_partf2(m1,i-1),&
 lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz1 = ui_vec(1)

        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int )
        call shepard_interp_2d ( size(lon_1d), lon_1d, lat_1d, kz2_1d_2 , power, 1, lon_partf2(m1,i-1), &
lat_partf2(m1,i-1), out_int_2 )
        vel_interp(1) = out_int(1)
        vel_interp(2) = out_int_2(1)       
        call pwl_value_1d (size(inter_depth), inter_depth, vel_interp, 1, par_dep , ui_vec )
        kz2 = ui_vec(1)
  
     endif


  else

       45364 continue


      ui1 = u_model1 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      ui2 = u_model2 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      vi1 = v_model1 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      vi2 = v_model2 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      wi1 = w_model1 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      wi2 = w_model2 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      kz1 = kz_model1 (lat_in_f2(1), lon_in_f2(2), indexz(1))
      kz2 = kz_model2 (lat_in_f2(1), lon_in_f2(2), indexz(1))

endif

      MM1 = M1

      vel_interp(1) = ui1
      vel_interp(2) = ui2
      call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
      uif2= ui_vec(1)
      vel_interp(1) = vi1
      vel_interp(2) = vi2
      call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
      vif2=ui_vec(1)
      vel_interp(1) = wi1
      vel_interp(2) = wi2
      call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
      wif2=ui_vec(1)
      vel_interp(1) = kz1
      vel_interp(2) = kz2
      call pwl_value_1d (size(intertime), intertime, vel_interp, 1, ts_vec , ui_vec )
      kzf2=ui_vec(1)

      if (three_dim .ne. 1) then 
        kzf2=0.0001
        wif2=0  
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DISSOLUTION FASE 2

       call  PROP_AMBIENTE(zf2(m1,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

       call VEL_ASCENSAO_BG ( dropdiamf2(M1,i-1) , DIAM_HYD , RO_A , spmtf2(m1,i-1), VIS_DIN_A  , TS_VC  , Wvp  )


          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )


         do comps=1, NCOMP_OIL

           MOL_COMP(comps)=masscompf2(m1,i-1, comps)* 1000.D0/PM_COMP_OIL(comps)

           FRAC_MASS_OUT(comps)=masscompf2(m1,i-1, comps)/massaf2(m1,i-1)

           V_COMP(comps) = masscompf2(m1,i-1, comps)/ RO_COMP_OIL(comps)

         enddo


         MOL_TOT=sum(mol_comp)

         V_TOT=massaf2(m1,i-1)/spmtf2(m1, i-1)
 
         volf2(m1,i)=massaf2(m1,i-1)/spmtf2(m1, i-1)

         moil=massaf2(m1,i-1)
  
         FRAC_TOT=sum(FRAC_MASS_OUT)

         CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)


         ! print*, '2', massaf2,  numdropsf2 (m1,i-1)


         call DISSOLVE_OIL_fase2(watcont(j,i), VIS_DIN_A, RO_A, spmtf2(m1,i-1),  dropdiamf2(M1,i-1), Wvp, VIS_DIN_OIL_OUT, TS_VC, &    !use last watcont and diam cause they ve been already modified, in thesis. 
                    dt, numdropsf2 (m1,i-1))


         massaf2 (m1,i) = sum(masscompf2(m1,i,:)) 

         if (massaf2(m1,i).le.0) then
             zf2(m1,:) = 0
             massaf2(m1,:) = 0
             LON_partf2(M1,:) = 0
             LAT_partf2(M1,:) = 0
             volf2(M1,:) = 0
             cycle
         endif



       if (DISSOLVED_FASE .EQ. 1) then

        if (parcel_dis_cont .lt. numparcels_dis) then

         parcel_dis_cont = parcel_dis_cont + 1

         xf3(parcel_dis_cont, i) = xf2(m1,i-1)
         yf3(parcel_dis_cont, i) = yf2(m1,i-1)
         zf3(parcel_dis_cont, i) = zf2(m1,i-1)
 
print*, parcel_dis_cont, 'fase 2'


         call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf3(parcel_dis_cont,i), LAT_partf3(parcel_dis_cont,i), PROFOUT, &
              xf3(parcel_dis_cont, i),  yf3(parcel_dis_cont, i), ZPOS, OPT)

         massaf3(parcel_dis_cont, :) = -(massaf2(m1,i) - massaf2(m1,i-1))   !!dissolved mass does not variate
            
   
!         print*, lat_partf2(m1,i-1), LAT_partf3(parcel_dis_cont,i)
!         print*, lon_partf2(m1,i-1), lon_partf3(parcel_dis_cont,i)
!         print*, massaf2(m1,i), massaf2(m1, i-1),  massaf3(parcel_dis_cont, i) 
!         print*, zf3(parcel_dis_cont, i)
!     print*, 'lulu'
!stop     
       endif
 
      endif



 !   PRINT*, 'MASSAF2', MASSAF2(M1,I), masscompf2(m1,i,1)

        
         aux=0
         do comps=1, NCOMP_OIL
           aux=aux + masscompf2(m1,i,comps)/RO_COMP_OIL(comps)
         enddo

         spmtf2(m1, i)=massaf2(m1,i)/aux
         volf2(m1,i)=massaf2(m1,i)/spmtf2(m1, i)


         massadropsf2(m1, i) = sum(masscompdropf2(m1,i,:))
   
         voldropsf2(m1,i) = massadropsf2(m1, i) / spmtf2(m1, i)

         dropdiamf2(m1,i) = 2 * (   (  ( 3* voldropsf2(m1,i) )  / (4. * pi)  )   ** (1./3.)       )


!       print*, '4',  massaf2(m1,i), massadropsf2(m1, i)


 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ENF OF DISSOLUTION

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERACTION OF DROPLETS WITH HARD SURFACES 

      if (turn_off_mask .eq. 1) then 
        go to 88689
      endif


        if (dif2.le.0) then
         zf2(m1,i) = zf2(m1,i-1)
         xf2(m1,i)=xf2(m1,i-1)
         yf2(m1,i)=yf2(m1,i-1)
         lon_partf2(m1,i) = lon_partf2(m1,i-1)
         lat_partf2(m1,i) = lat_partf2(m1,i-1)  
 !       print*, 'DI', di               
         go to 493 

        elseif(-zf2(m1,i-1) .ge. dif2) then
         zf2(m1,i) = zf2(m1,i-1)
         xf2(m1,i)=xf2(m1,i-1)
         yf2(m1,i)=yf2(m1,i-1)
         lon_partf2(m1,i) = lon_partf2(m1,i-1)
         lat_partf2(m1,i) = lat_partf2(m1,i-1)  
         go to 493 
       
        elseif ( (lat_model(lat_in_f2(1)-1, lon_in_f2(2)).eq.0) .or. (lat_model(lat_in_f2(1)+1, lon_in_f2(2)) .eq.0)  .or. &
                 (lat_model(lat_in_f2(1), lon_in_f2(2)+1).eq.0) .or. (lat_model(lat_in_f2(1), lon_in_f2(2)-1) .eq.0) ) then
         xf2(m1,i)=xf2(m1,i-1)
         yf2(m1,i)=yf2(m1,i-1)
         lon_partf2(m1,i) = lon_partf2(m1,i-1)
         lat_partf2(m1,i) = lat_partf2(m1,i-1)  
         zf2(m1,i) = zf2(m1,i-1)
         go to 493 

        elseif ( (lon_model(lat_in_f2(1)-1, lon_in_f2(2)).eq.0) .or. (lon_model(lat_in_f2(1)+1, lon_in_f2(2)) .eq.0)  .or. &
              (lon_model(lat_in_f2(1), lon_in_f2(2)+1).eq.0) .or. (lon_model(lat_in_f2(1), lon_in_f2(2)-1) .eq.0) ) then
         xf2(m1,i)=xf2(m1,i-1)
         yf2(m1,i)=yf2(m1,i-1)
         lon_partf2(m1,i) = lon_partf2(m1,i-1)
         lat_partf2(m1,i) = lat_partf2(m1,i-1)  
         zf2(m1,i) = zf2(m1,i-1)
         go to 493 
        endif

  88689 continue
 
     !!!!!!!!!!!!!!!!!!!!!!!!advection + random fase2
   if (right_random .eq. 0) then


        CALL random_number(RN1)
        CALL random_number(RN2)
        CALL random_number(RN3)


        randvert = -1. + 2.*RN3
       
         U_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * COS(2.D0*PI*RN2)    
         V_ALEA = RN1 * ((2.D0*CDIF_HOR/dt)**(0.5D0)) * SIN(2.D0*PI*RN2)
         W_ALEA = randvert * ((2.D0*kzf2/dt)**(0.5D0))    !!  based on Reed et al., 1995
 

        xf2(m1,i)=xf2(m1,i-1) + uif2*dt +   U_ALEA*dt  

        yf2(m1,i)=yf2(m1,i-1) + vif2*dt +  V_ALEA*dt 
 

        call  PROP_AMBIENTE(zf2(m1,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)   !!gotta update Z to zf2(m1,i-1)
        call VEL_ASCENSAO_BG ( dropdiamf2(M1,i) , DIAM_HYD , RO_A , spmtf2(m1, i) , VIS_DIN_A  , TS_VC  , Wvp  )

  !     print*, '444', spmt, spmtf2(m1,i-1)


        zf2(m1,i)=zf2(m1,i-1) + wif2*dt +  W_ALEA*dt  + Wvp*dt


        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf2(M1,i), LAT_partf2(M1,i), PROFOUT,  xf2(M1, i),  yf2(M1, i), ZPOS, OPT)

      

     else

     call  PROP_AMBIENTE(zf2(m1,i-1) , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)   !!gotta update Z to zf2(m1,i-1)
     call VEL_ASCENSAO_BG ( dropdiamf2(M1,i) , DIAM_HYD , RO_A , spmtf2(m1, i) , VIS_DIN_A  , TS_VC  , Wvp  )

      do step_random=1,dt, dt_random
        CALL random_number(RN1)
        CALL random_number(RN2)
        CALL random_number(RN3)


        randvert = -1. + 2.*RN3
       

        xrandom = xrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * COS(2.D0*PI*RN2) ) *  dt_random
        yrandom = yrandom + ( RN1 * ((2.D0*CDIF_HOR/dt_random)**(0.5D0)) * SIN(2.D0*PI*RN2) ) *  dt_random   
        zrandom = zrandom + ( randvert * ((2.D0*kzf2/dt_random)**(0.5D0)) ) *dt_random + Wvp*dt_random  !!  based on Reed et al., 1995

 
      enddo



        xf2(m1,i)=xf2(m1,i-1) + uif2*dt +   xrandom 

        yf2(m1,i)=yf2(m1,i-1) + vif2*dt +   yrandom
 

  !     print*, '444', spmt, spmtf2(m1,i-1)


        zf2(m1,i)=zf2(m1,i-1) + wif2*dt +  zrandom 


        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf2(M1,i), LAT_partf2(M1,i), PROFOUT,  xf2(M1, i),  yf2(M1, i), ZPOS, OPT)


       xrandom = 0 
       yrandom = 0
       zrandom = 0


     endif

  493 continue

 !       print*, 'cool'
  !      print*, lon_part(:,i), lon_part(:,i-2), LAT_part(:,i), LAT_part(:,i-2)
  !      print*, LON_partf2(M1,i), LAT_partf2(M1,i)
  !      print*, LON_partf2(M1,i-1), LAT_partf2(M1,i-1)  
  !      print*, shape(lat_part)    

  !     print*, 'minloc', pinf2, pinf2(1)

!       PRINT*, LAT_part(:,i-1) 
!       PRINT*, lon_part(:,i-1)
!       PRINT*, LAT_partf2(M1,i)
!       PRINT*, LON_partf2(M1,i)

!print*, 'kk', zf2(m1,i), zf2(m1,i-1), wif2, zrandom, wvp, kzf2

 !      print*, 'massa1', massa(pinf2(1), i),  zf2(m1,i), diam(pinf2(1), i), spmt
  !     print*, 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG', zf2(m1,i)

       if (zf2(m1,i) .ge. 0) then 

           IF (LOW_MEMORY.EQ.1) THEN
             zf2(m1,:) = 0
             LON_partf2(M1,:) = 0
             LAT_partf2(M1,:) = 0
             massaf2(m1,:) = 0

             lat_part_sum=abs(LAT_part(:,i-1)-LAT_partf2(M1,i))
             lon_part_sum=abs(lon_part(:,i-1)-LON_partf2(M1,i))
  
             coord_sum_f2=lat_part_sum + lon_part_sum  
   
             pinf2 = minloc(coord_sum_f2)
 !  print*, '5555'
             !!!!!!!!!!!!!!! MASS OF DROPLET WILL BE ADSORBED BY NEARBY PARTICLE
             do comps=1, NCOMP_OIL
                 masscomp(pinf2(1),i,comps)= masscomp(pinf2(1),i-1,comps) + masscompf2(m1,i,comps)
                 evapmass(pinf2(1), comps) = masscomp(pinf2(1),i-1,comps) + masscompf2(m1,i,comps)
             enddo
             massa(pinf2(1), i) = sum(masscomp(pinf2(1),i,:))

            massa(pinf2(1), i-1) = massa(pinf2(1), i)
            masscomp(pinf2(1),i-1,:) = masscomp(pinf2(1),i,:)

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )

            aux=0
            do comps=1, NCOMP_OIL
              aux=aux + masscomp(pinf2(1),i,comps)/RO_COMP_OIL(comps)
            enddo

            spmt=massa(pinf2(1),i)/aux

            diam(pinf2(1), i)=2*(((massa(pinf2(1), i)/spmt)/(pi*height(pinf2(1), i)))**0.5)  !meter      
  
            area(pinf2(1), i)=pi*((diam(pinf2(1), i)/2)**2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           do comps=1, NCOMP_OIL

              MOL_COMP(comps)=masscomp(pinf2(1),i,comps)* 1000.D0/PM_COMP_OIL(comps)

              FRAC_MASS_OUT(comps)=masscomp(pinf2(1),i,comps)/massa(pinf2(1),i)

              V_COMP(comps) = masscomp(pinf2(1),i,comps)/ RO_COMP_OIL(comps)

           enddo

           call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

           MOL_TOT=sum(mol_comp)

           V_TOT=massa(pinf2(1),i)/spmt
 
           VOL(pinf2(1),i)=massa(pinf2(1),i)/spmt

           moil=massa(pinf2(1),i)

           FRAC_TOT=sum(FRAC_MASS_OUT)


            CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!back to particle reference
  !           print*, 'massa2', massa(pinf2(1), i), zf2(m1,i), diam(pinf2(1), i), spmt
         ELSE 
 !            PRINT*, 'NUNTOT',numtot, m1, i, j
             NUMTOT  = NUMTOT + 1
             zf2(m1,:) = 0
             LON_partf2(M1,:) = 0
             LAT_partf2(M1,:) = 0
             massaf2(m1,:) = 0

             x(NUMTOT,i) = xf2(m1, i-1)
             y(NUMTOT,i) = yf2(m1, i-1)

             call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(NUMTOT,i), LAT_part(NUMTOT,i), PROFOUT, &
  x(NUMTOT,i),   y(NUMTOT,i), ZPOS, OPT)

             do comps=1, NCOMP_OIL
                 masscomp(NUMTOT,i,comps)= masscompf2(m1,i,comps)
                 evapmass(NUMTOT, comps) = masscompf2(m1,i,comps)
             enddo
             massa(NUMTOT, i) = sum(masscomp(NUMTOT,i,:))

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )

            aux=0
            do comps=1, NCOMP_OIL
              aux=aux + masscomp(NUMTOT,i,comps)/RO_COMP_OIL(comps)
            enddo

            spmt=massa(NUMTOT,i)/aux

            do comps=1, NCOMP_OIL
                 MOL_COMP(comps)=masscomp(NUMTOT,i, comps)* 1000.D0/PM_COMP_OIL(comps)
                 FRAC_MASS_OUT(comps)=masscomp(NUMTOT,i, comps)/massa(NUMTOT,i)
                 V_COMP(comps) = masscomp(NUMTOT,i, comps)/ RO_COMP_OIL(comps)
                 FRAC_MASS_OUT_PART(NUMTOT,comps) = masscomp(NUMTOT,i, comps)/massa(NUMTOT,i)
            enddo

            call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)

            MOL_TOT=sum(mol_comp)

            V_TOT=massa(NUMTOT,i)/spmt
 
            VOL(NUMTOT,i)=massa(NUMTOT,i)/spmt

            moil=massa(NUMTOT,i)

            FRAC_TOT=sum(FRAC_MASS_OUT)


            CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)


           visc_e(NUMTOT, i) = VIS_DIN_OIL_OUT*1000
  
           rho_e(NUMTOT, i) =  RO_OIL_OUT         

!print*, VOL(NUMTOT,i), VOL(1,i)
!            area(NUMTOT, i)=pi*((diam(NUMTOT, i)/2)**2)  

           time_ini_spread = ((1.45/1.14)**4.) * ( (  VOL(NUMTOT,i) / ((VIS_DIN_A/RO_A) * gravity * &
((RO_A - RO_OIL_OUT )/RO_A ) ) )**(1./3.)  )
         
           time_ini_spread = time_ini_spread/60

           dt_h_spr(NUMTOT) = time_ini_spread
  
            area(NUMTOT, i)=2.1*(pi)*  (   ( ((VOL(NUMTOT,i)**2.)*gravity*((RO_A - RO_OIL_OUT)/RO_A) * &
 ((time_ini_spread*60)**(3./2.)))/ ((VIS_DIN_OIL_OUT /RO_OIL_OUT)**(1./2.))  )**(1./3.)  )


            diam(NUMTOT, i) =  ((area(NUMTOT, i) / PI ) ** (1./2.)) * 2.
 !           diam(NUMTOT, i)=2*(((massa(NUMTOT, i)/spmt)/(pi*height(NUMTOT, i)))**0.5)  !meter   


            diamem(NUMTOT, i) = diam(NUMTOT, i)

            areaem(NUMTOT, i) = area(NUMTOT, i)

            height(NUMTOT, i) = vol(NUMTOT, i) / area(NUMTOT, i)

             masscomp(NUMTOT,i-1,:) = masscomp(NUMTOT,i,:)

             massa(NUMTOT, i-1) = massa(NUMTOT, i)

             diam(NUMTOT, i-1) =  diam(NUMTOT, i)

            dropdiamax =  cmax * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

            dropdiamin =  cmin * (((VIS_DIN_OIL_OUT /RO_OIL_OUT) * 1000000) ** (0.34) ) *  ( turbed**(-0.4)  )

            dropdiam(NUMTOT,i) = ((dropdiamax - dropdiamin)/2. + dropdiamin) * 0.000001

            voldrops(NUMTOT,i) = ((dropdiam(NUMTOT,i)/2.)**3.) *  pi * (4./3.)

            numdrops(NUMTOT,i) =  vol(NUMTOT,i) / voldrops (NUMTOT,i)

    
            teste=2

        ENDIF
            
       endif
 
  !     print*, '6656'
 !      print*, 'rtrt', zf2(m1,i), wif2*dt,  W_ALEA*dt, Wvp*dt, w_alea, wvp, dropdiamf2(M1,i-1)
 !      print*, i

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo

    
  !     if ( dt_hf2 .ge.1) then

 !        print*, 'WRITE DATA FASE2'
  !       dt_hf2 = 0 
   !    endif

   endif
 endif
 
!  if(parcelcont.eq.1) then
!    print*, '123'
!    go to 456

!  else
 if (massa(j,i-1) .eq. 0) then
   go to 726   
 endif 

! print*, qd, wvp * dt + (- zini)

   if (qd.gt.0) then
!       if ( wvp * dt + (- zini) .lt.0) then

          if (parcelcont.EQ.numparcels) then
              massa(j,i) =  massa(j,i-1)
              massae(j,i) = massae(j,i-1)
              masscomp(j,i,:) = masscomp(j,i-1,:)
              diam(j,i) = diam(j,i-1)
              diamem(j,i) = diamem(j,i-1)
              area(j,i) = area(j,i-1)
              areaem(j,i) = areaem(j,i-1)
              height(j,i) = height(j,i-1)
              vol(j,i)    = vol(j,i-1)
              goto 456
          endif


          parcelcont=parcelcont+1


          xf2(parcelcont, i) = x(j,i-1)
          yf2(parcelcont, i) = y(j,i-1)
          zf2(parcelcont, i) = - zini

 !         print*, 'iyyyyyyyyyyyyyyyyyyyyyyyyyy', zf2(parcelcont, i), i


          call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_partf2(parcelcont,i), LAT_partf2(parcelcont,i), &
 PROFOUT,  xf2(parcelcont, i),  yf2(parcelcont, i), ZPOS, OPT)


         massaf2 (parcelcont,i) =  qd*dt*area(j,i)

        if (massaf2 (parcelcont,i) .ge. massa(j,i-1) ) then 
            massaf2 (parcelcont,i) = massa(j,i-1)
            massa(j,i)= 0
        endif

         do comps=1, NCOMP_OIL
           masscompf2(parcelcont,i,comps)=massaf2 (parcelcont,i) * (  masscomp(j, i-1, comps) / massa(j,i-1) )
         enddo

          CALL PROPRIEDADES_COMP_OIL ( TEMP_OUT , PM_COMP_OIL(:)	 , &
		         TEB_COMP_OIL(:) , RO_COMP_OIL_15(:) , &
		         RO_COMP_OIL(:) , CP_COMP_OIL(:) )
            
         aux=0
         do comps=1, NCOMP_OIL
           aux=aux + masscompf2(parcelcont,i,comps)/RO_COMP_OIL(comps)
         enddo

         spmtf2(parcelcont, i)=massaf2(parcelcont,i)/aux


          
         volf2(parcelcont,i)=massaf2(parcelcont,i)/spmtf2(parcelcont, i)


         dropdiamf2(parcelcont,:) = dropdiam(j,i-1)

         voldropsf2(parcelcont,:) = ((dropdiamf2(parcelcont,i)/2.)**3.) *  pi * (4./3.)

         numdropsf2(parcelcont,:) =  volf2(parcelcont,i) / voldropsf2(parcelcont,i)
     
         massadropsf2(parcelcont, :) = voldropsf2(parcelcont,i)  * spmtf2(parcelcont, i)


         do comps=1, NCOMP_OIL
           masscompdropf2(parcelcont,i,comps)=massadropsf2 (parcelcont,i) * (  masscomp(j, i-1, comps) / massa(j,i-1) )
         enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  BACK TO PARTICLE REFERENCE

         if ( massa(j,i) .eq. 0) then
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
          go to 726   
         endif

         massa(j,i) = massa(j,i-1) - massaf2 (parcelcont,i) 

         do comps=1, NCOMP_OIL
           masscomp(j,i,comps)=massa(j,i) * (  masscomp(j, i-1, comps) / massa(j,i-1) ) !!mass of components also decrease, %of components did not change!!
           evapmass(j,comps) = masscomp(j,i,comps)
         enddo

 !        print*, 'yy', massa(j,i), sum(masscomp(j,i,:))

         spmt = spmtf2(parcelcont, i)

         height(j,i) = height(j,i-1)

         diam(j,i)=2*(((massa(j,i)/spmt)/(pi*height(j,i)))**0.5)  !meter      

         area(j,i)=pi*((diam(j,i)/2)**2)

!       print*, height(j,i-1), diam(j,i-1), diam(j,i)

!!!!! NEW properties of emulsion

         massae(j,i)=massa(j,i)+watcont(j,i) 

         watf(j,i) = watcont(j,i)/massa(j,i)
     
         call rho_emulsion(spmt)


         diamem(j,i)=2*(((massae(j,i)/rho_e(j,i))/(pi*height(j,i)))**0.5)  !meter


         areaem(j,i)=pi*((diamem(j,i)/2)**2)

  !!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do comps=1, NCOMP_OIL

        MOL_COMP(comps)=masscomp(j,i, comps)* 1000.D0/PM_COMP_OIL(comps)

        FRAC_MASS_OUT(comps)=masscomp(j,i, comps)/massa(j,i)

        V_COMP(comps) = masscomp(j,i, comps)/ RO_COMP_OIL(comps)

     enddo

     call  PROP_AMBIENTE(Z , PROF_REF , temp_out , SAL_A, RO_A , VIS_DIN_A , CP_A , TS_A)



      MOL_TOT=sum(mol_comp)

      V_TOT=massa(j,i)/spmt
 
      VOL(j,i)=massa(j,i)/spmt

      moil=massa(j,i)

      FRAC_TOT=sum(FRAC_MASS_OUT)


     CALL components_part2( API , VAZAO_OIL_OUT , DT , TEMP_OUT , &
	         	RO_OIL_OUT , PM_OIL_OUT , TS_OIL_OUT , VIS_DIN_OIL_OUT , TS_A, TS_VC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!back to particle reference
!    else 
!              massa(j,i) =  massa(j,i-1)
!              massae(j,i) = massae(j,i-1)
!              masscomp(j,i,:) = masscomp(j,i-1,:)
!              diam(j,i) = diam(j,i-1)
!              diamem(j,i) = diamem(j,i-1)
!              area(j,i) = area(j,i-1)
!              areaem(j,i) = areaem(j,i-1)
!              height(j,i) = height(j,i-1)
!              vol(j,i)    = vol(j,i-1)

!    endif


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

   endif
!  endif

456 continue
 !  if (partcont .ge. 1) then

!   do m2=numpalm+1,numpalm+numpal

!    print*, 'massa', massa(j,i), diam(j,i), area(j,
     !!!!PUT HERE PARTICLE UPDATED MASS!!!!!!!!!!!!!!!!! ALREAD PUT IN EVAPORATION




endif                                                           !!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT


!print*, massa(j,i)

contind=contind+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       if (x(j,i-1).ge.15)then
!        x(j,i)=15

!!        y(j,i)=y(j,i-1)                                    !very,very,very simplified sticking in solid surface
!        go to 1000
!       end if 
  if (turn_off_mask .eq. 1) then 

    go to 8867


  else  
     if(THEORETICAL.EQ.1) THEN


         maskt=mask(minloc(abs(y_mask(:,1)-y(j,i-1))), minloc(abs(x_mask(1,:)-x(j,i-1))))       !!!NEW MASK
      
         if (maskt(1,1).eq.1.) then
          x(j,i)=x(j,i-1)
          y(j,i)=y(j,i-1)                      
          go to 1000 
         endif

     ELSE     
                                                                                       !!!!!!!!!MASK WITH DELFT

!print*, 'di', di

       if (di.le.0) then

        OPT=1

!        x(j,i)=x(j,i-1)
!        y(j,i)=y(j,i-1)
!        lon_part(j,i) = lon_part(j,i-1)
!        lat_part(j,i) = lat_part(j,i-1)  

        lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
        lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
        OPT=2

        go to 1000 
       
       elseif ( (lat_model(lat_in(1)-1, lon_in(2)).eq.0) .or. (lat_model(lat_in(1)+1, lon_in(2)) .eq.0)  .or. &
             (lat_model(lat_in(1), lon_in(2)+1).eq.0) .or. (lat_model(lat_in(1), lon_in(2)-1) .eq.0) ) then

        OPT=1

!        x(j,i)=x(j,i-1)
!        y(j,i)=y(j,i-1)
!        lon_part(j,i) = lon_part(j,i-1)
!        lat_part(j,i) = lat_part(j,i-1)  

        lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
        lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
        OPT=2

!        print*, lon_part(j,i-1), lon_part(j,i)
!        print*, lat_part(j,i-1), lat_part(j,i)
!        PRINT*, X(J,I-1), X(J,I)
!        PRINT*, Y(J,I-1),Y(J,I)
! PRINT*, 'dryyyyy'
!        STOP

        go to 1000 

       elseif ( (lon_model(lat_in(1)-1, lon_in(2)).eq.0) .or. (lon_model(lat_in(1)+1, lon_in(2)) .eq.0)  .or. &
             (lon_model(lat_in(1), lon_in(2)+1).eq.0) .or. (lon_model(lat_in(1), lon_in(2)-1) .eq.0) ) then

        OPT=1

!        x(j,i)=x(j,i-1)
!        y(j,i)=y(j,i-1)
!        lon_part(j,i) = lon_part(j,i-1)
!        lat_part(j,i) = lat_part(j,i-1)  


        lon_part(j,i) = lon_model(lat_in(1), lon_in(2))
        lat_part(j,i) = lat_model(lat_in(1), lon_in(2))
        call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)
        OPT=2

!        print*, lon_part(j,i-1), lon_part(j,i)
!        print*, lat_part(j,i-1), lat_part(j,i)
!        PRINT*, X(J,I-1), X(J,I)
!        PRINT*, Y(J,I-1),Y(J,I)
! PRINT*, 'dryyyyy'
!        STOP

        go to 1000 
       endif

     ENDIF

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
8867 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

!print*, 'hh', uwd, vwd
!print*, x(j,i-1) + ui(1,1)*dt +  xrandom, y(j,i-1) + vi(1,1)*dt +  yrandom
!print*, uwd, (uwd * 0.03)*dt 
!print*, x(j,i-1) + ui(1,1)*dt +  xrandom + (uwd * 0.03)*dt, y(j,i-1) + vi(1,1)*dt +  yrandom + (vwd * 0.03)*dt
!stop

       x(j,i)=x(j,i-1) + 1*ui(1,1)*dt +  xrandom + ( widfc * (  (uwd * cos(0*(pi/180)))  + (vwd * sin(0*(pi/180))) )  )   *dt 

       y(j,i)=y(j,i-1) + 1*vi(1,1)*dt +  yrandom + ( widfc * (  (-uwd * sin(0*(pi/180)))  + (vwd * cos(0*(pi/180))) )  )   *dt 

       call INI_AMB( LON_REF, LAT_REF, PROF_REF, LON_part(j,i), LAT_part(j,i), PROFOUT,  x(j,i),  y(j,i), ZPOS, OPT)

       xrandom = 0 
       yrandom = 0
       zrandom = 0

   endif

! print*, 'cool', uwd, vwd, CDIF_HOR


1000   continue 

 !      print*, cont, tsl

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

  !          print*, 'entrou'
  
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
!        print*, height_map(lat_in_height(1), lon_in_height(2)), vol(j,i)
!        print*, lat_in_height(1), lon_in_height(2), lat_part(j,i), lat_height(1,1), lat_height(239, 239)
!        print*, lat_height(:,1)
!        print*, lat_in_height(1), lon_in_height(2), lon_part(j,i), lon_height(1,1), lon_height(239, 239)
!        print*, lon_height(1,:)

       !!!!!!!!!!!!!!


!print*, i, CONT, TSL, NUMPALM, J
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

!print*, i


!     print*, cont, tsl
!     print*, tsl
!     print*, cont
 468 continue

  if (PROBABILISTIC.eq.0) then
     if (counttime_2_h.ge.time_finish) then
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
!           print*, 1
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


!           print*, x(m,1), y(m,1)
       enddo
       NUMTOT=NUMTOT+numpal   
       cont=0
       lo=lo+1
!       print*, lo

!       print*, freq
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

 !   print*, x(i,1), y(i,1), lon_part(i,1), lat_part(i,1), opt
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
 !      print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       CALL random_number(RN5)
       count_prob_2 = count_prob_2+1
       counttime = time_lim*RN5 -dt   
       counttimeh_era = (time_lim_wind*RN5 -dt)/60.   
       counttimeh = counttime/60.


  print*, rn5
     ENDIF

   endif

  endif

  cont=cont+1

  if (dt_h .ge. 1) then
  ! if (cont.eq.tsl-1) then
     write(12,fmt2) x(:,i)
     write(13,fmt2) y(:,i)
     write(14,fmt2) massa(:,i)
     write(15,fmt2) diam(:,i)
     write(16,fmt2) height(:,i)
 !    write(17,fmt2) visc_e(:,i)
     write(18,fmt1) counttimeh/(60)
 !    write(19,fmt2) spmt
 !    write(23,fmt2) porc_evap(:,i)
     write(221,fmt2) lat_part(:,i)
     write(222,fmt2) lon_part(:,i)
     write(54,fmt2) watf(:,i) 

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
       write(223,fmt3) lat_partf2(:,i)
       write(224,fmt3) lon_partf2(:,i)
       write(225,fmt3) zf2(:,i)
       write(141,fmt3) massaf2(:,i)
     ENDIF
!     print*, 'WRITE DATA', dt_h

!     if(counttime/(60*60) .ge. 3.84) then
!        stop
!     endif

     
     dt_h=0
   endif


1547 continue
!     print*, x(:,i), y(:,i)
     counttime=counttime+dt

     counttimeh = counttimeh + (dt/60.)     !!!it's minute in fact

     counttimeh_era = counttimeh_era + (dt/60.)  

     counttimeh_coup = counttimeh_coup + (dt/60.)

     counttime_2 = counttime_2+dt

     counttime_2_h = counttime_2/(60*60)

     dt_h = dt_h + (dt/(60.*60.))

    
     height_map = height_map*0
!     print*, 'count_era', counttimeh_era
!      print*, di, lon_part(:,i)
 !    print*, counttime/(60*60)
 !    print*,  spmt

 !   print*, i, CONT, TSL, NUMPALM, J
! print*, massa(:,i)
! print*, height(:,i)
! print*, diam(:,i)



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
masscomp(:,1,:) = masscomp(:,i, :)
mas_evap(:,1) = mas_evap(:,i)
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
numdropsf2(:,1) = numdropsf2(:,i)
massadropsf2(:,1) = massadropsf2(:, i)
masscompdropf2(:,1,:) = masscompdropf2(:,i,:)
voldropsf2(:,1) = voldropsf2(:,i) 
     xf2(:,1) =  xf2(:,i)
     yf2(:,1) =  yf2(:,i)
     zf2(:,1) = zf2(:,i)
LON_partf2(:,1) = LON_partf2(:,i)
Lat_partf2(:,1) = Lat_partf2(:,i)
LON_partf3(:,1) = LON_partf3(:,i)
Lat_partf3(:,1) = Lat_partf3(:,i)
     xf3(:,1) =  xf3(:,i)
     yf3(:,1) =  yf3(:,i)
     zf3(:,1) = zf3(:,i)
    massaf3(:,1) = massaf3(:,i)
    contind=2

   print*, 'ENTROU NO TS_LIM&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7'

 !  print*, massa(:,i)
 !  print*, shape(massa)
!print*, shape(LON_partf2)

   endif

  if (cont_ts.eq. ts) then
     go to 78623
  endif

  enddo

enddo


3000 if (trim(pre_end).eq.'end') then
    print*, 'SIMULATION TIME GREATER THAN BOUNDARY CONDITION TIME - CONTINUOUS'
  endif

3888 if (trim(end_prob).eq.'finish') then
    print*, 'END OF PROBABILISTIC SIMULATION'
  endif

78623 continue

  probmap = (DBLE( partcontmap) / dble(total_particles)  )*100.

  probmap_2 = (DBLE(probmap_2)/DBLE(contprob))*100.

  conc=conc/contprob

!  print*, '44444', sum(conc), massa(1,1)*500

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




                print*, '77', maxval(probmap_2)

 call cpu_time(stop_time)


print *, "Simulation time:", &
      (stop_time - start_time)/60., "minutes"

print*, 'total', total_particles
print*, 'contprob', contprob
print*, counttime
  
PRINT*, 'END CONTINUOUS SIMULATION'
end program
  
