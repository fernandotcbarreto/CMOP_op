module evaporation

 use oil_fractions


implicit none

integer i2

double precision, dimension(:), allocatable ::  porep

double precision, dimension(:, :), allocatable :: evapmass

double precision:: kf1

contains


subroutine evap_vars


porep(1)=1
porep(2)=1
porep(3)=1
porep(4)=1
porep(5)=1
porep(6)=1
porep(7)=1
porep(8)=1
porep(9)=1
porep(10)=1
porep(11)=1
porep(12)=1
porep(13)=1
porep(14)=1
porep(15)=1
porep(16)=1
porep(17)=1
porep(18)=1   
porep(19)=1
porep(20)=1
porep(21)=1
porep(22)=1
porep(23)=1
porep(24)=1
porep(25)=1


end subroutine evap_vars



subroutine evap_mass


  do i2=1, NCOMP_OIL

    evapmass(:,i2)=MASSCOMPREF(i2)*porep(i2)

  enddo



end subroutine evap_mass


subroutine evap_mass_coupling(numtot, num_res_par)

  integer:: numtot, num_res_par


  do i2=1, NCOMP_OIL

    evapmass(num_res_par+1:numtot,i2)=masscomp(numtot,1, i2)*porep(i2)

  enddo



end subroutine evap_mass_coupling



end module evaporation
