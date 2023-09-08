subroutine rad_reaphy()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_reaphy
  ! NAME 
  !    rad_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition for the
  !    radiation heat transfer equation.
  ! USES
  ! USED BY
  !    rad_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_radiat
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ipara,ispec,itype,iprop
  real(rp)    :: remis_tmp,scatf_tmp

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     !
     ! Reach the section
     !
     call ecoute('rad_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('rad_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('rad_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           kfl_parti_rad=0  ! No particles in suspension
           call ecoute('rad_reaphy')
           do while(words(1)/='ENDPR')
!!$              if(words(1)=='PARTI' ) then  ! There are particles
!!$                 kfl_parti_rad=1
!!$                 ! We have to find out which variable carries temperature
!!$                 do iprop=1,mlapr
!!$                    if (parttyp(itype) % prova(iprop) == ID_TEMPE) idtem_rad=iprop
!!$                 enddo
!!$                 ! Now read properties for types
!!$                 call ecoute('rad_reaphy')
!!$                 if( words(1) == 'TYPE ' ) then
!!$                    itype = getint('TYPE ',1_ip,'#TYPE NUMBER OF LAGRANGIAN PARTICLE')
!!$                    if( itype < 1 .or. itype > mtyla ) call runend('RADIAT REAPHY: WRONG PARTICLE TYPE')
!!$                    if (parttyp(itype) % kfl_exist /= 1 ) call runend('RADIAT REAPHY: PARTICLE TYPE MUST ALSO BE DEFINED IN .DAT')
!!$                    call ecoute('rad_reaphy')
!!$                    do while(words(1) /= 'ENDTY' )
!!$                       if( words(1) == 'TEMPE' ) then
!!$                          parttyp(itype) % prope(idtem_rad) = getrea('TEMPE',1.0_rp,'#PARTICLE DEFAULT TEMPERATURE')
!!$                       else if( words(1) == 'EMISI' ) then
!!$                          parttyp(itype) % denpa = getrea('EMISI',1.0_rp,'#PARTICLE EMISIVITY')
!!$                       else if( words(1) == 'SCATT' ) then
!!$                          parttyp(itype) % denpa = getrea('SCATT',1.0_rp,'#PARTICLE SCATTERING FACTOR')
!!$                       end if
!!$                       call ecoute('rad_reaphy')                 
!!$                    enddo
!!$                 else
!!$                    call runend('RADIAT REAPHY: IF PARTICLES REQUESTED NEED TO SPECIFY AT LEAST ONE TYPE')
!!$                 endif
!!$              endif
              call ecoute('rad_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !
           ! Allocate memory
           !
           call rad_memphy(1_ip)

           kfl_atest_rad=0_ip
           do ispec=1,nspec_rad
              scatt_rad(ispec)=0.0_rp
              absor_rad(ispec)=0.0_rp
              aniso_rad(ispec)=0.0_rp
           enddo

           call ecoute('rad_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='SCATT') then               ! Scattering strength (molar)
                 scatt_rad=param(1:nspec_rad)
              else if(words(1)=='ABSOR') then               ! Molar Absorption coefficient
                 absor_rad=param(1:nspec_rad)
              else if(words(1)=='ANISO') then               ! Scattering anisotropy
                 aniso_rad=param(1:nspec_rad)
              else if(words(1)=='TESTP') then     !Test problem 
                 kfl_atest_rad = param(1)   ! Test number
                 !
                 ! Test parameters
                 !     TEST_PROBLEM 1 = Absorption coeff
                 !
                 expar_rad(1) = param(2)    
              end if
              call ecoute('rad_reaphy')
           end do
        end if

     end do

  end if

end subroutine rad_reaphy
