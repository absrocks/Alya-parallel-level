program main

    ! use the Cantera module
    use cantera

    implicit none
    type(phase_t) :: gas
    integer nsp, nrxns
    double precision :: T, p, dnu, dlam

    integer  :: index_O2, index_N2, index_AR
    integer  :: ii
    double precision :: temp, rho, visc, cp, thermoCond
    double precision, allocatable, dimension(:) :: moleF_air, diff_k
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gas = importPhase('gri30.xml','gri30_mix')

    index_O2 = speciesIndex(gas, 'O2') - 1
    index_N2 = speciesIndex(gas, 'N2') - 1
    index_AR = speciesIndex(gas, 'AR') - 1

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nsp = nSpecies(gas)
    allocate(moleF_air(0:nsp-1)); moleF_air = 0.0
    allocate(diff_k   (0:nsp-1)); diff_k    = 0.0

    moleF_air(index_O2) = 0.2095
    moleF_air(index_N2) = 0.7809
    moleF_air(index_AR) = 0.0096

    temp = 298.0; p = 101325.0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call setState_TPX(gas, temp, p, moleF_air)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rho        = density            (gas)
    visc       = viscosity          (gas)
    cp         = cp_mass            (gas)
    thermoCond = thermalConductivity(gas)
    call getMixDiffCoeffs(gas, diff_k)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print'(A,F15.8)', "rho        =", rho
    print'(A,F15.8)', "visc       =", visc
    print'(A,F15.8)', "cp         =", cp
    print'(A,F15.8)', "thermoCond =", thermoCond



!!DMM    write(*,*)
!!DMM    write(*,*) '********   Fortran 90 Test Program   ********'
!!DMM    print*,''
!!DMM    !!gas = importPhase('sandiego20161214_mechCK.cti')
!!DMM    gas = importPhase('gri30.xml','gri30')
!!DMM
!!DMM    print*,''
!!DMM    print*,'Reading San Diego mechanism:'
!!DMM    nsp   = nSpecies(gas)    ! number of species
!!DMM    nrxns = nReactions(gas)  ! number of reactions
!!DMM   
!!DMM    print*,'# Species=',nsp
!!DMM    print*,'# Reactions=',nrxns
!!DMM    print*,''
!!DMM     
!!DMM    T = 300.0
!!DMM    p = 101325.0
!!DMM
!!DMM    call setState_TPX(gas, T, p, 'CH4:0.0445, O2:0.7328, N2:0.2227 ')
!!DMM    dnu  = viscosity(gas)
!!DMM    dlam = thermalConductivity(gas)
!!DMM
!!DMM    write(*,*) dnu, dlam


    !! call demo(gas, nsp, nrxns)

end program main

!--------------------------------------------------------

subroutine demo(gas, MAXSP, MAXRXNS)

  ! use the Cantera module
  use cantera

  implicit none

  ! declare the arguments
  type(phase_t), intent(inout) :: gas
  integer,       intent(in) :: MAXSP
  integer,       intent(in) :: MAXRXNS

  double precision q(MAXRXNS), qf(MAXRXNS), qr(MAXRXNS)
  double precision diff(MAXSP),diff_multi(MAXSP),y(MAXSP),x(MAXSP),W(MAXSP),D_t(MAXSP),wdot(MAXSP)

  character*80 eq
  character*20 name
  double precision :: dnu, dlam,W_mixt
  integer :: i, irxns, k

  !
  ! Define fuel composition in mass fraction
  !
  y = 0.0
  do k = 1, MAXSP
     call getSpeciesName(gas, k, name)
     if (name == 'CH4') then
       y(k) = 0.0445
     elseif (name == 'O2') then
       y(k) = 0.7328
     elseif (name == 'N2') then
       y(k) = 0.2227
     endif
  end do

  !
  ! Compute molecular weight mixture
  !
  call getMolecularWeights(gas,W)

  W_mixt = 0.0
  do k = 1, MAXSP
     W_mixt = W_mixt + y(k)/W(k)
  end do
  W_mixt = 1.0/W_mixt

  !
  ! Define fuel composition in gas phase
  !

  write(*,*) 'Unburnt conditions for CH4 at phi = 0.8:'
  write(*,*) '' 
  call setMassFractions(gas,y)

  y = 0.0
  print*,'Mass fractions:'
  do k = 1, MAXSP
     call getSpeciesName(gas, k, name)
     call getMassFractions(gas,y)
     write(*,*) name,'=', y(k)
  end do
  write(*,*) ''
  write(*,*)'Molecular weight mixture: W =',W_mixt

  write(*,*) ''
  write(*,10) temperature(gas), pressure(gas), density(gas), &
       enthalpy_mole(gas), entropy_mole(gas), cp_mass(gas)

  ! compute the equilibrium state holding the specific
  ! enthalpy and pressure constant
  call equilibrate(gas, 'HP')

  write(*,*) '----------------------------------------------------------- '
  write(*,*) 'Equilibrium conditions for CH4 at phi = 0.8:'
  write(*,*) ''

  y = 0.0
  write(*,*) 'Mass fractions equilibrium'
  do k = 1, MAXSP
     call getSpeciesName(gas, k, name)
     call getMassFractions(gas,y)
     write(*,*) name,'=', y(k)
  end do

  write(*,10) temperature(gas), pressure(gas), density(gas), &
       enthalpy_mole(gas), entropy_mole(gas), cp_mass(gas)

  print*,' ----------------------------------------------------------- '
  print*,' '
  print*,' '

  write(*,*) 'Transport properties at equilibrium'
!!  call setTemperature(gas,temperature(gas))
!!  call setDensity(gas,density(gas))
!!  call setMassFractions(gas,y)
!!  setState_TRY(gas,temperature(gas),density(gas),y)

  print*,'AFTER setState_TRY'

  write(*,10) temperature(gas), pressure(gas), density(gas), &
       enthalpy_mole(gas), entropy_mole(gas), cp_mass(gas)

  ! 
  ! Transport properties
  !
  dnu  = viscosity(gas)
  dlam = thermalConductivity(gas)

  write(*,30) dnu, dlam

!!  call getThermalDiffCoeffs(gas,D_t)
!!
!!  write(*,*) 'Thermal Coefficients'
!!  do k = 1, MAXSP
!!     call getSpeciesName(gas, k, name)
!!     write(*,41) name, D_t(k)
!!  end do
  write(*,*) ''

  call getMixDiffCoeffs(gas, diff)
  write(*,*) 'Species Diffusion Coefficients'
  do k = 1, MAXSP
     call getSpeciesName(gas, k, name)
     write(*,40) name, diff(k)
  end do

  print*,' ----------------------------------------------------------- '
  print*,'  '
  print*,' '
  write(*,*) 'Chemical kinetics information '
  print*,' '
  !     Reaction information
  irxns = nReactions(gas)

  ! forward and reverse rates of progress should be equal
  ! in equilibrium states
  call getFwdRatesOfProgress(gas, qf)
  call getRevRatesOfProgress(gas, qr)

  ! net rates of progress should be zero in equilibrium states
  call getNetRatesOfProgress(gas, q)

  write(*,*)''
  write(*,*)'Reaction rates: forward, backward and net'
  write(*,*)''
  ! for each reaction, print the equation and the rates of progress
  do i = 1,irxns
     call getReactionString(gas, i,eq)
     write(*,20) eq,qf(i),qr(i),q(i)
  end do
  write(*,*)' Note: net rates of progress should be zero at equilibrium.'

  write(*,*)''
  write(*,*)'Source terms of species:'

  ! for each reaction, print the source term for each species
  wdot = 0.0
  call getNetProductionRates(gas, wdot)
  do k = 1,MAXSP
     call getSpeciesName(gas, k, name)
     write(*,42) name, wdot(k)
  end do
  write(*,*)' Note: net rates of progress should be zero at equilibrium.'

  10 format(//'Temperature:   ',g14.5,' K'/ &
       'Pressure:      ',g14.5,' Pa'/ &
       'Density:       ',g14.5,' kg/m3'/ &
       'Molar Enthalpy:',g14.5,' J/kmol'/ &
       'Molar Entropy: ',g14.5,' J/kmol-K'/ &
       'Molar cp:      ',g14.5,' J/kmol-K'//)

  20   format(a27,3e14.5,' kmol/m3/s')

  30 format(//'Viscosity:             ',g14.5,'  Pa-s'/ &
       'Thermal conductivity:  ',g14.5,'  W/m/K'/)

  40   format(' ',a20,e14.5,' m2/s')
  41   format(' ',a20,e14.5,' m2/s')
  42   format(' ',a20,e14.5,' 1/s ')


end subroutine demo

