subroutine tem_reaphy()
  !------------------------------------------------------------------------
  !****f* Temper/tem_reaphy
  ! NAME 
  !    tem_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition for the
  !    temperature equation.
  ! USES
  ! USED BY
  !    tem_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_temper
  use def_domain
  use def_kermod
  use mod_ker_space_time_function
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_timei_tem = 0                                    ! Stationary flow
     kfl_advec_tem = 0                                    ! Convection is off
     kfl_joule_tem = 0                                    ! Joule effect is off
     kfl_radia_tem = 0                                    ! No radiation
     kfl_sourc_tem = 0                                    ! Sources are off
     kfl_cotur_tem = 0                                    ! Turbulence coupling
     kfl_tfles_tem = 0                                    ! Thickened flame model activation
     kfl_adiab_tem = 0                                    ! Adiabatic mixing problem activation (= 1, adiabatic mixing) 
     kfl_condu_tem = 1                                    ! Conductivity is on
     kfl_exint_tem = 0                                    ! No properties interpolation
     kfl_inter_tem = 0                                    ! No interpolation of arrays
     kfl_dynco_tem = 0                                    ! Dynamical coupling
     kfl_regim_tem = 0                                    ! Regime
     kfl_prope_tem = 0                                    ! properties update strategy CFI model
     kfl_parti_tem = 0                                    ! no particles in suspension
     kfl_flux_tem  = 0                                    ! Activate heat flux provided from fields
     turbu_tem     = 0.0_rp                               ! Turbulence parameters
     prtur_tem     = 1.0_rp                               ! Turbulent Prandtl number Prt = 0
     scond_tem     = 0.0_rp                               ! S conductivity
     react_tem     = 0.0_rp                               ! Reaction term
 
     !
     ! Reach the section
     !
     call ecoute('tem_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('tem_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('tem_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('tem_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='TEMPO') then                    ! Temporal evolution
                 if(words(2)=='ON   ') then
                    kfl_timei_tem = 1
                 else
                    kfl_timei_tem = 0
                 end if

              else if(words(1)=='REGIM') then              
                 if(words(2)=='INCOM') then
                    kfl_regim_tem=0                         ! Incompressible
                 else if(words(2)=='COMPR') then
                    kfl_regim_tem=1                         ! Compressible
                    if(exists('PRESS')) kfl_regim_tem=1
                    if(exists('DENSI')) kfl_regim_tem=2
                 else if(words(2)=='LOWMA') then
                    kfl_regim_tem = 3                       ! Low-Mach
                    if (words(3)=='ENTHA') then
                       kfl_regim_tem = 4                    ! Enthalpy equation
                       if(exists('ADIAB')) then             ! Input file --> REGIME:  LOWMACH, ENTHALPY, ADIABATIC: HMAX = 1000.0, HMIN = -300.0
                         cfi_hmax_tem = getrea('HMAX ',0.0_rp,'#Maximum enthalpy')
                         cfi_hmin_tem = getrea('HMIN ',0.0_rp,'#Minimum enthalpy')
                         kfl_adiab_tem = 1
                       end if
                    endif
                 end if

              else if(words(1)=='CONDU') then               ! Conductivity term
                 if(words(2)=='ON  ') then
                    kfl_condu_tem = 1 
                 else
                    kfl_condu_tem = 0
                 end if

              else if(words(1)=='DYNAM') then               ! Dynamical coupling
                 if(words(2)=='ON   ') kfl_dynco_tem=1

              else if(words(1)=='CONVE') then               ! Convective term
                 if(exists('ON   ')) then
                    kfl_advec_tem = 1
                    if(exists('FUNCT')) &
                         kfl_advec_tem = getint('FUNCT',0_ip,'#Velocity function')
                    if(words(3)=='VELOC') then
                       if(words(4)=='FUNCT') then
                          kfl_advec_tem = getint('FUNCT',0_ip,'#Velocity function')
                       else
                          kfl_advec_tem = 1                       
                       end if
                    else if(words(3)=='GRADI') then
                       kfl_advec_tem = -1                       
                    end if
                 end if

              else if(words(1).eq.'TURBU') then
                 if(words(2)=='LESMO') then
                    if(words(3)=='SMAGO') then
                       kfl_cotur_tem=-1
                       turbu_tem=getrea('PARAM',0.0_rp,'#Coefficient c')
                    elseif (words(3)=='WALE ') then
                       kfl_cotur_tem=-1
                       turbu_tem=getrea('PARAM',0.0_rp,'#Coefficient c')
                    end if
                 else if(words(2)=='RANSA') then
                    if(exists('MIXIN')) then
                       kfl_cotur_tem=-10
                       turbu_tem=getrea('PARAM',0.0_rp,'#Coefficient C L^2')
                    else if(exists('CONST')) then
                       kfl_cotur_tem=-11
                       turbu_tem=getrea('PARAM',0.0_rp,'#Constant turbulent viscosity')
                    else if(exists('XUCHE')) then
                       kfl_cotur_tem=-12
                    end if
                 else if(words(2)=='RANSD'.or.words(2)=='FROMT') then     ! We are using turmu
                    kfl_cotur_tem=1
                 end if

              else if( words(1) == 'SOURC' ) then               ! Source term

                 if( exists('SPACE') ) then
                    kfl_sourc_tem = space_time_function_number(getcha('SPACE','NONE ','#Space time function'))
                 else if( exists('FIELD') ) then 
                    kfl_sourc_tem = -getint('FIELD',1_ip,'#Element field number')
                 end if
                 if(exists('JOULE')) kfl_joule_tem = 1

              else if( words(1) == 'TFLES' ) then               ! Thickened flame model
                 if( exists('ON   ') ) kfl_tfles_tem = 1_ip

              else if( words(1) == 'HEATF' ) then               ! Source term
                    kfl_flux_tem = getint('FIELD',-1_ip,'#Boundary field number')

              else if(words(1)=='RADIA') then               ! Radiation
                 if(words(2)=='SURFA') kfl_radia_tem = 1

!!$              else if(words(1)=='PARTI' ) then  ! There are particles
!!$                 kfl_parti_tem=1
!!$                 ! We are going to request only one property to keep track with particles (temperature)
!!$                 idtem_tem=1_ip
!!$                 call lagdef(3_ip, idtem_tem)
!!$                 ! Now read properties for types
!!$                 call ecoute('tem_reaphy')
!!$                 if( words(1) == 'TYPE ' ) then
!!$                    itype = getint('TYPE ',1_ip,'#TYPE NUMBER OF LAGRANGIAN PARTICLE')
!!$                    if( itype < 1 .or. itype > mtyla ) call runend('TEMPER REAPHY: WRONG PARTICLE TYPE')
!!$                    if (parttyp(itype) % kfl_exist /= 1 ) call runend('TEMPER REAPHY: PARTICLE TYPE MUST ALSO BE DEFINED IN .DAT')
!!$                    parttyp(itype) % prova(idtem_tem) = ID_TEMPE  ! We assign the slot to temperature
!!$                    call ecoute('tem_reaphy')
!!$                    do while(words(1) /= 'ENDTY' )
!!$                       if( words(1) == 'TEMPE' ) then
!!$                          parttyp(itype) % prope(idtem_tem) = getrea('TEMPE',1.0_rp,'#PARTICLE DEFAULT TEMPERATURE')
!!$                       else if( words(1) == 'CALOR' ) then
!!$                          parttyp(itype) % calor = getrea('EMISI',1.0_rp,'#CALORIFIC CAPACITY')
!!$                       end if
!!$                       call ecoute('tem_reaphy')                 
!!$                    enddo
!!$                 else
!!$                    call runend('TEMPER REAPHY: IF PARTICLES REQUESTED NEED TO SPECIFY AT LEAST ONE TYPE')
!!$                 endif

              else if(words(1)=='REACT') then               
                 !
                 ! Reaction term (s)
                 !
                 react_tem = getrea('REACT',0.0_rp,'#Reaction parameter')
                 
              end if

              call ecoute('tem_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !
           ! Allocate memory
           !  

           call ecoute('tem_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='REACT') then               ! Reaction term (s)
                 react_tem = getrea('REACT',0.0_rp,'#Reaction parameter')

              else if(words(1)=='TURBU') then               ! Turbulent Prandtl number
                 prtur_tem = getrea('TURBU',1.0_rp,'#Turbulent Prandtl number')

              else if(words(1)=='UPDAT') then               ! properties update strategy CFI model
                 if(words(2)=='SYNCH') then
                   kfl_prope_tem = 1    ! syncronized update at the end of temper time step
                 else if(words(2)=='LOCAL') then 
                   kfl_prope_tem = 0    ! local update of properties
                 end if

              end if
              call ecoute('tem_reaphy')
           end do
        end if
     end do

  end if

end subroutine tem_reaphy
