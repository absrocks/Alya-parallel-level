subroutine chm_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* partis/chm_output
  ! NAME 
  !    chm_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    chm_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_postpr
  use mod_ker_proper 
  use def_kermod

  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ipoin,iclas,iodes,ispec,dummi
  real(rp)                :: rutim,dummr(1),auxi
  character(5)            :: wopos(3)
  character(2)            :: wclas


  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim=cutim
  nullify(gesca)
  nullify(gevec)
  nullify(ger3p)

  select case (ivari)  

  case(1_ip)
     !
     ! Class concentration
     !
     wopos(2)=postp(1)%wopos(2,1)
     wopos(3)=postp(1)%wopos(3,1)
     
     if( kfl_model_chm == 2) then
        
        if( INOTMASTER ) call memgen(0_ip,npoin,0_ip)
        do iclas=1,nclas_chm
           wclas =  intost(iclas)
           if(iclas<10) then
              wopos(1)=postp(1)%wopos(1,1)(1:3)//'0'//trim(wclas)
           else
              wopos(1)=postp(1)%wopos(1,1)(1:3)//trim(wclas)
           end if
           !
           !  Convert to g/m3  
           !
           if( INOTMASTER ) then
              if(lawde_chm == 0) then
                 do ipoin = 1,npoin
                    gesca(ipoin) = conce(ipoin,iclas,1)*denma_chm*1.0e3_rp 
                 end do
              else
                 do ipoin = 1,npoin
                    gesca(ipoin) = conce(ipoin,iclas,1)*densi_chm(ipoin)*1.0e3_rp 
                 end do
              end if
           end if
           !
           call postpr(gesca,wopos,ittim,rutim)
        end do
        
        if( INOTMASTER ) call memgen(2_ip,npoin,0_ip)
        !
     else
        do ispec=1,min(8_ip,nspec_chm)
          if( INOTMASTER ) gesca => conce(1:,ispec,1) 
          if( kfl_model_chm == 4) then ! Combustion saves names of species
             wopos(1) = speci(ispec) % name
          else
             wclas =  intost(ispec)
             if(ispec<10) then
                wopos(1)=postp(1)%wopos(1,1)(1:3)//'0'//trim(wclas)
             else
                wopos(1)=postp(1)%wopos(1,1)(1:3)//trim(wclas)
             end if
          endif
          call postpr(gesca,wopos,ittim,rutim)

       end do
       if( INOTMASTER ) nullify(gesca)
    end if

     return

  case(2_ip)
     !
     ! Defects in clusters 
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = 0.0_rp
           do iodes = 1,nodes_chm
              ispec = iodes + nclas_chm
              gesca(ipoin) = gesca(ipoin) + real(ispec)*conce(ipoin,ispec,1)
           end do
        end do
     end if
         
  case(3_ip)
     !
     ! Total defects 
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = 0.0_rp
           do iodes = 1,nodes_chm
              ispec = iodes + nclas_chm
              gesca(ipoin) = gesca(ipoin) + real(ispec)*conce(ipoin,ispec,1)
           end do
           do iclas = 1,nclas_chm
              gesca(ipoin) = gesca(ipoin) + conce(ipoin,iclas,1)
           end do
        end do
     end if

  case(4_ip)
     !
     ! PROJE_CHM: Projection
     !
     if( INOTMASTER ) gesca => proje_chm(:,1)
         
  case(5_ip)
     !
     ! VELOC_CHM: Velocity
     !
     if( INOTMASTER ) gevec => veloc_chm
         
  case(6_ip)
     !
     ! COSGS_CHM: SGS
     !
     if( kfl_sgsti_chm == 0 ) return
     if( INOTMASTER ) then
           iclas_chm = 1
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        call runend('NOT PROGRAMMED FOR THIS MODEL')
        gesca => rhsid
     end if

  case(7_ip)
     !
     ! Limiter
     !
     if( kfl_limit_chm /= 0 ) then
        if( INOTMASTER ) then
           iclas_chm = 1
           do ipoin = 1,npoin
              rhsid(ipoin) = 0.0_rp
           end do
           call runend('NOT PROGRAMMED FOR THIS MODEL')
           call rhsmod(1_ip,rhsid)
           do ipoin = 1,npoin
              rhsid(ipoin) = rhsid(ipoin) / vmass(ipoin)
           end do
        end if
        gesca => rhsid
     end if

  case(8_ip)
     !
     ! Projection
     !
     if( kfl_stabi_chm > 1 ) then
        if( INOTMASTER ) then
           iclas_chm = 1
           do ipoin = 1,npoin
              rhsid(ipoin) = 0.0_rp
           end do
           call runend('NOT PROGRAMMED FOR THIS MODEL')
           call rhsmod(1_ip,rhsid)
           do ipoin = 1,npoin
              rhsid(ipoin) = rhsid(ipoin) / vmass(ipoin)
           end do
        end if
        gesca => rhsid
     end if

  case(9_ip)
     ! 
     ! Total concentration. Converted to (g/m3)
     !
     if(INOTMASTER ) then 
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp
           do iclas=1,nclas_chm
              gesca(ipoin) = gesca(ipoin) + conce(ipoin,iclas,1)
           end do
        end do
        !
        !  Convert to g/m3  
        !
        if(lawde_chm == 0) then
             do ipoin = 1,npoin
                gesca(ipoin) =  gesca(ipoin)*denma_chm*1.0e3_rp 
             end do
         else
             do ipoin = 1,npoin
                gesca(ipoin) =  gesca(ipoin)*densi_chm(ipoin)*1.0e3_rp 
             end do
         end if
     end if

  case(10_ip)
    ! 
    !  Ground accumulation (kg/m2)
    !
    if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = accum_chm(ipoin)
       end do
    end if

  case(11_ip)
     ! 
     ! Source term 
     !
     if (kfl_model_chm == 4 ) then
        !
        ! Class source
        !
        wopos(2)=postp(1)%wopos(2,11)
        wopos(3)=postp(1)%wopos(3,11)
        
        do ispec=1,min(8_ip,nspec_chm)
           if( INOTMASTER ) gesca => massk(1:,ispec) 
           wopos(1) = 'S'//speci(ispec)%name(1:4)
           call postpr(gesca,wopos,ittim,rutim)
        end do
        if( INOTMASTER ) nullify(gesca)
        return
     else if (kfl_model_chm == 5 ) then
        !
        ! CFI model, source only for class 1
        !
        if(INOTMASTER ) gesca => massk(:,1)
     else
        if(INOTMASTER ) then 
           !
           ! Global source
           !
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin) = 0.0_rp
              do iclas=1,nclas_chm
                 gesca(ipoin) = gesca(ipoin) + tmrat_chm(iclas,ipoin)
              end do
           end do
        end if
     end if
     
  case(12_ip)
    !
    ! Convection velocity of m
     if(INOTMASTER ) then
        gesca => proad_chm
     end if

     
  case(13_ip)
    !
    ! Viscosity
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('VISCO','NPOIN',dummi,dummi,gesca)
     end if

     
  case(14_ip)
    !
    ! Specific heat
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('SPHEA','NPOIN',dummi,dummi,gesca)
     end if
     
  case(15_ip)
    !
    ! Heat conductivity
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('CONDU','NPOIN',dummi,dummi,gesca)
     end if

  case(16_ip)
     !
     ! Enthalpy
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        gesca=0.0_rp
        do ispec=1,nspec_chm
           do ipoin=1,npoin
              gesca(ipoin) = gesca(ipoin) + entha_chm(ipoin,ispec) * conce(ipoin,ispec,1)
           enddo
        enddo
     end if


  case(17_ip)
     !
     ! Enthalpy heat source term
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call smooth (div_enthalpy_transport, gesca)      
     end if


  case(18_ip)
     !
     ! Chemical heat only source term 
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call smooth (chemical_heat, gesca)      
     end if

  case(19_ip)
     ! 
     !Sum of concentration
     !
     if(INOTMASTER ) then 
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp
           do iclas=1,nclas_chm
              gesca(ipoin) = gesca(ipoin) + conce(ipoin,iclas,1)
           end do
        end do
     end if

  case(20_ip)
     ! 
     ! Equivalence ratio
     !
     if(INOTMASTER ) then

       equiv_chm = 0.0_rp
       call chm_outeqr()
       call rhsmod(1_ip,equiv_chm)
       do ipoin = 1,npoin
          equiv_chm(ipoin) = equiv_chm(ipoin) / vmass(ipoin)
       end do

       gesca => equiv_chm(:)
     end if

  case(21_ip)
    !
    ! Molecular weight
     if(INOTMASTER ) then
        gesca => wmean(:,1)
     end if

  case(22_ip)
    !
    ! Flame front sensor
    !
     if(INOTMASTER ) then
        gesca => flsen_chm(:)
     end if
 
  case(23_ip)
    !
    ! Flame speed
    !
    if(INOTMASTER ) then

       flspe_chm = 0.0_rp
       call chm_outeqr()
       call rhsmod(1_ip,flspe_chm)
       do ipoin = 1,npoin
          flspe_chm(ipoin) = flspe_chm(ipoin) / vmass(ipoin)
       end do
       gesca => flspe_chm(:)
    end if

  case(24_ip)
    !
    ! Flame thickness
    !
    if(INOTMASTER ) then

       flthi_chm = 0.0_rp
       call chm_outeqr()
       call rhsmod(1_ip,flthi_chm)

       do ipoin = 1,npoin
          flthi_chm(ipoin) = flthi_chm(ipoin) / vmass(ipoin)
       end do

       gesca => flthi_chm(:)
    end if

  case(25_ip)
    !
    ! Subgrid scale wrinkling factor
    !
    if(INOTMASTER ) then

       flsgs_chm = 0.0_rp
       call chm_outeqr()
       call rhsmod(1_ip,flsgs_chm)

       do ipoin = 1,npoin
          flsgs_chm(ipoin) = flsgs_chm(ipoin) / vmass(ipoin)
       end do

       gesca => flsgs_chm(:)
    end if
 
    case(26_ip)
    !
    ! Subgrid scale wrinkling factor
    !
    if(INOTMASTER ) then

       flfac_chm = 0.0_rp
       call chm_outeqr()
       call rhsmod(1_ip,flfac_chm)

       do ipoin = 1,npoin
          flfac_chm(ipoin) = flfac_chm(ipoin) / vmass(ipoin)
       end do

       gesca => flfac_chm(:)
    end if

  case(27_ip)
    !
    ! Temperature for low-mach CFI combustion model
    !
    if(INOTMASTER ) then
       if(associated(tempe)) then
          gesca => tempe(:,1)
       end if
    end if

  case(28_ip)
    !
    ! AVTEM 
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avtem_chm(ipoin) / auxi
             avtem_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(29_ip)
    !
    ! Enthalpy scalar CFI model
     if(INOTMASTER ) then
        gesca => encfi(:)
     end if

  case(30_ip)
    !
    ! AVCON 
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avcon_chm(ipoin) / auxi
             avcon_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(31_ip)
    !
    ! AVCO2
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avco2_chm(ipoin) / auxi
             avco2_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif  

  case(32_ip)
    !
    ! AVVAR
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avvar_chm(ipoin) / auxi
             avvar_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(33_ip)
    !
    ! AVIME
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avime_chm(ipoin) / auxi
             avime_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(34_ip)
    !
    ! AVCHM
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avchm_chm(ipoin) / auxi
             avime_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(35_ip)
    !
    ! AVMIX 
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avmix_chm(ipoin) / auxi
             avmix_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(36_ip)
    !
    ! AVMI2
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avmi2_chm(ipoin) / auxi
             avmi2_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(37_ip)
    ! 
    !  Species post-processing for CFI, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(1,ipoin)
         end do
      end if
    end if

  case(38_ip)
    ! 
    !  Species post-processing for CFI, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(2,ipoin)
         end do
      end if
    end if

  case(39_ip)
     !
     ! Radiative source term for heat equation, CFI model 
     !
    if (kfl_radia_chm > 0) then
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call smooth (radiative_heat, gesca)
     end if
    end if

  case(40_ip)
     !
     ! Scaled reaction progress for partially premixed conditions 
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = yscale_chm(ipoin)
       end do
     end if

  case(41_ip)
     !
     ! Variance of scaled reaction progress for partially premixed conditions 
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = cvar_chm(ipoin)
       end do
     end if

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1)%wopos(1,ivari))

end subroutine chm_outvar 
