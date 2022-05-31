subroutine chm_outvar(ivari,imesh)
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
  use mod_memory
  use mod_arrays,          only : arrays
  use mod_solver,          only : solver_lumped_mass_system 
  use mod_chm_entropy,     only : chm_entropy_postprocess 
  use mod_interp_tab,      only : fw_scale_cont_var 
  use def_chemic,          only : table_fw
  use mod_chm_finiteRate,  only : get_mixture_fraction
  use mod_outvar,          only : outvar 
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ipoin,iclas,ipostvar,iodes,dummi,ielem,pelty,pgaus,imixf
  real(rp)                :: rutim,auxi,y_c_eq,y_c_0,retva(100_ip)
  character(5)            :: wopos(3)
  character(3)            :: wclasYk
  character(2)            :: wclas 
  type(r1p),pointer       :: aux_r1p(:)
  type(r2p),pointer       :: aux_r2p(:)
  real(rp), pointer       :: auxvar(:,:)
  
  real(rp)                  :: control(5_ip)   ! input of table lookup function 
  real(rp)                  :: scale_control(5_ip)
  real(rp)                  :: lim_control(5_ip,2_ip)

  character(len=10)         :: valZ, aux_name

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
     call arrays(ivari,'POSTPROCESS',conce)
     return

  case(2_ip)
     ! 
     ! Source term 
     !
     if (kfl_model_chm == 1 ) then
        !
        ! Flamelet model, source only for class 1
        !
        if (kfl_lookg_chm > 0) then
           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,max(1_ip,nelem))
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
              aux_r1p(ielem) % a = mass_gp(ielem) % a(:,1,1)
           end do
           call memgen(zero,max(1_ip,npoin),zero)
           call smooth (aux_r1p, gesca)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
        else
            if(INOTMASTER ) gesca => massk(:,1)
        endif
     
     else if (kfl_model_chm == 3 ) then
        !
        ! Finite Rate Chemistry - Chemical Sources
        !
        wopos(2)=postp(1)%wopos(2,2)
        wopos(3)=postp(1)%wopos(3,2)

#ifdef CANTERA
        do iclas=1,nclas_chm
           if( INOTMASTER ) gesca => src_chm(1:npoin,iclas) 
           if (iclas < nspec_chm ) then
              call getSpeciesName(gas_chm,iclas,wopos(1))
              wopos(1) = 'S_'//wopos(1)
           else
              wopos(1) = 'Q_'
           end if
           call postpr(gesca,wopos,ittim,rutim)
        end do
#endif
        if( INOTMASTER ) nullify(gesca)
        return
     end if
     
  case(3_ip)
    !
    ! Viscosity
    !
    if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       call ker_proper('VISCO','NPOIN',dummi,dummi,gesca)
    end if

     
  case(4_ip)
    !
    ! Specific heat
    !
    if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       call ker_proper('SPHEA','NPOIN',dummi,dummi,gesca)
    end if
     
  case(5_ip)
    !
    ! Heat conductivity
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('CONDU','NPOIN',dummi,dummi,gesca)
     end if

  case(6_ip)
     !
     ! Enthalpy
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        gesca=0.0_rp
        do iclas=1,nspec_chm
           do ipoin=1,npoin
              gesca(ipoin) = gesca(ipoin) + entha_chm(ipoin,iclas) * conce(ipoin,iclas,1)
           enddo
        enddo
     end if


  case(7_ip)
     !
     ! Enthalpy heat source term
     !
     if(INOTMASTER ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p(ielem)%a,pgaus)
           aux_r1p(ielem) % a = div_enthalpy_transport(ielem) % a(:,1,1)
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if


  case(8_ip)
     !
     ! Chemical heat only source term 
     !
     if(INOTMASTER ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p(ielem)%a,pgaus)
           aux_r1p(ielem) % a = chemical_heat(ielem) % a(:,1,1)
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if

  case(9_ip)
     ! 
     ! Sum of concentration
     !
     if(INOTMASTER ) then 
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp
           do iclas=1,nspec_chm
              gesca(ipoin) = gesca(ipoin) + conce(ipoin,iclas,1)
           end do
        end do
     end if

  case(10_ip)
     !
     ! Molecular weight
     !
     if(INOTMASTER ) then
        if (kfl_lookg_chm > 0 .and. kfl_model_chm /= 4) then
            nullify(aux_r1p)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
            do ielem = 1,nelem
               pelty = ltype(ielem)
               pgaus = ngaus(pelty)
               call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p(ielem)%a,pgaus)
               aux_r1p(ielem) % a = wmean_gp(ielem) % a(:,1,1)
            end do
            call memgen(zero,npoin,zero)
            call smooth (aux_r1p, gesca)
            call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
        else
            gesca => wmean(:,1)
        endif
     end if

  case(11_ip)
    !
    ! Temperature 
    !
    if(INOTMASTER ) then
        if(associated(tempe)) then
            gesca => tempe(:,1)
        end if
    end if

  case(12_ip)
    !
    ! AV Yc or C
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avY_chm(ipoin) / auxi
             avY_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(13_ip)
    !
    ! AVYv
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avYv_chm(ipoin) / auxi
             avYv_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif  

  case(14_ip)
    !
    ! AV variance of Z
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZv_chm(ipoin) / auxi
             avZv_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(15_ip)
    !
    ! AVCHM
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avchm_chm(ipoin) / auxi
             avchm_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(16_ip)
    !
    ! AVZ 
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZ_chm(ipoin) / auxi
             avZ_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(17_ip)
    !
    ! AVZ2
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZ2_chm(ipoin) / auxi
             avZ2_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case(18_ip)
    ! 
    !  Species post-processing for Flamelet, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(1,ipoin)
         end do
      end if
    end if

  case(19_ip)
    ! 
    !  Species post-processing for Flamelet, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(2,ipoin)
         end do
      end if
    end if

  case(20_ip)
     !
     ! Radiative source term for heat equation, Flamelet model 
     !
    if (kfl_radia_chm > 0) then
     if(INOTMASTER ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p(ielem)%a,pgaus)
           aux_r1p(ielem) % a = radiative_heat(ielem) % a(:,1,1)
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if
    end if

  case(21_ip)
     !
     ! Maximum gradient of the mixture fraction from diffusion flamelets
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = zgradmax_chm(ipoin)
       end do
     end if
  case(22_ip)
     !
     ! Weighting factor between premixed and diffusion flamelets
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = phi_chm(ipoin)
       end do
     end if

  case(23_ip)
     !
     ! Scalar dissipation rate of the reaction progress variable (resolved)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xYr_chm(ipoin)
       end do
     end if

  case(24_ip)
     !
     ! Scalar dissipation rate of the mixture fraction (resolved)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xZr_chm(ipoin)
       end do
     end if

  case(25_ip)
     !
     ! Scalar dissipation rate of the reaction progress variable (subgrid)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xYs_chm(ipoin)
       end do
     end if

  case(26_ip)
     !
     ! Scalar dissipation rate of the mixture fraction (subgrid)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xZs_chm(ipoin)
       end do
     end if

  case(27_ip,28_ip,29_ip,30_ip)
      call runend('chm_outvar: average scalar dissipation rates are not coded')

  case(31_ip)
     !
     ! Average progress variable squared
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avY2_chm(ipoin) / auxi
              avY2_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(32_ip)
     !
     ! Average liquid volume fraction
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm

        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avL_chm(ipoin) / auxi
              avL_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(33_ip)
     !
     ! Average liquid volume fraction squared
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avL2_chm(ipoin) / auxi
              avL2_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(34_ip)
     !
     ! Average interface surface density Sigma
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avS_chm(ipoin) / auxi
              avS_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(35_ip)
     !
     ! Average interface surface density Sigma_0
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avS0_chm(ipoin) / auxi
              avS0_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(36_ip)
     !
     ! Average Sauter mean diameter
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avd32_chm(ipoin) / auxi
              avd32_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(37_ip)
     !
     ! Average density
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avden_chm(ipoin) / auxi
              avden_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case(38_ip)
     !
     ! Interface surface density Sigma
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = Sigma_chm(ipoin)
       end do
     end if

  case(39_ip)
     !
     ! Interface surface density Sigma_0 or Sigma_min
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = Sigm0_chm(ipoin)
       end do
     end if

  case(40_ip)
     !
     ! Sauter mean diameter
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = d32_chm(ipoin)
       end do
     end if

  case(41_ip)
     !
     ! Gradient mass fractions
     !
     if( INOTMASTER )then 
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           gevec(1:ndime,ipoin) = grad_Yk(1,1:ndime,ipoin) 
        end do
     end if

  case(42_ip)
     !
     ! Enthalpy transport by diffusion 
     !
     if (kfl_model_chm /= 4) then
        if( INOTMASTER )then 
           nullify(aux_r2p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p(ielem)%a,pgaus,ndime)
              aux_r2p(ielem) % a = enthalpy_transport(ielem) % a(:,1:ndime,1)
           end do
           call memgen(zero,ndime,npoin)
           call smoot5(aux_r2p,gevec,ndime)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
        end if
     end if
  
  case(43_ip)
       !
       ! Elemental Mass Fraction - H
       !
       if( INOTMASTER )then 
          call memgen(zero,npoin,zero)
#ifdef CANTERA 
          do ipoin=1,npoin
             call cantera_elemh(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if
  
  case(44_ip)
       !
       ! Elemental Mass Fraction - O
       !
       if( INOTMASTER )then 
          call memgen(zero,npoin,zero)
#ifdef CANTERA
          do ipoin=1,npoin
             call cantera_elemo(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if
  
  case(45_ip)
       !
       ! Elemental Mass Fraction - C
       !
       if( INOTMASTER )then 
          call memgen(zero,npoin,zero)
#ifdef CANTERA
          do ipoin=1,npoin
             call cantera_elemc(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if


  case(46_ip)
     !
     ! Mass source from partis 
     !
     call memgen(zero,max(1_ip,npoin),zero)
     if (associated(mass_sink)) then
        do ipoin=1,npoin
           gesca(ipoin) = mass_sink(ipoin)
        enddo
        call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp 
        end do
     endif

  case(47_ip)
     !
     ! Entropy viscosity maximum among species 
     !
     call memgen(zero,max(npoin,1_ip),zero)
     call chm_entropy_postprocess(0_ip,gesca)

  case(48_ip)
     !
     ! Scaled progress variable
     !
     wopos(2)=postp(1)%wopos(2,48)
     wopos(3)=postp(1)%wopos(3,48)
     
     if( INOTMASTER ) then
        nullify(auxvar)
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, npoin, table_fw % main_table % ndim)
     endif    
     
     !
     ! Do lookup
     !
     do ipoin=1,npoin
         control = 0.0_rp
         do iclas = 1, table_fw % main_table % ndim
            select case (table_fw % main_table % coords(iclas) % name)
            case ('CMEAN','C    ')
                control(iclas) = conce(ipoin,1,1)
            case ('CVAR ')
                control(iclas) = conce(ipoin,2,1)
            case ('CHIST')
                control(iclas) = xZr_chm(ipoin) + xZs_chm(ipoin)
            case ('ZMEAN','Z    ')
                control(iclas) = conce(ipoin,3,1)
            case ('ZVAR ')
                control(iclas) = conce(ipoin,4,1)
            case ('IMEAN','I    ')
                control(iclas) = therm(ipoin,1)
            end select
         enddo

         call fw_scale_cont_var( control, scale_control, lim_control, table_fw)
         do iclas=1,table_fw % main_table % ndim
            auxvar(ipoin,iclas) = scale_control(iclas) 
         enddo

     end do
     
     do iclas=1,table_fw % main_table % ndim
        if( INOTMASTER ) then
            gesca => auxvar(:,iclas)
        endif
        
        wclas =  intost(iclas)
        if(iclas<10) then
           wopos(1)=postp(1)%wopos(1,48)(1:3)//'0'//trim(wclas)
        else
           wopos(1)=postp(1)%wopos(1,48)(1:3)//trim(wclas)
        end if
        call postpr(gesca,wopos,ittim,rutim)
     end do
     if( INOTMASTER ) nullify(gesca)
     
     if( INOTMASTER ) then
        call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
     endif    
     
     return

  case(49_ip)
     !
     ! Sum of Reactions (Finite Rate)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = sum(React_ind(ipoin,:))
       end do
     end if

  case(50_ip)
     !
     ! Mixture Fraction 
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          call get_mixture_fraction(conce(ipoin,:,1),gesca(ipoin))
       end do
     end if

  case(51_ip)
        !
        ! Finite Rate Chemistry - Instantaneous heat release
        !
        if( INOTMASTER ) gesca => hrr_chm

  case(52_ip)
        !
        ! Finite Rate Chemistry - Averaged heat release
        !
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           if ( INOTMASTER ) then
              do ipoin=1,npoin
                 unkno(ipoin)           = hrr_avg_chm(ipoin) / auxi
                 hrr_avg_chm(ipoin) = 0.0_rp
              end do
              gesca => unkno
           endif
        else
           return
        endif

  case(53_ip)
     !
     ! Average mass source from spray
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        do ipoin=1,npoin
           unkno(ipoin)     = avmsk_chm(ipoin) / auxi
           avmsk_chm(ipoin) = 0.0_rp
        end do
        gesca => unkno
     else
        return
     endif

  case(54_ip)
     !
     ! Post-processing lookup variables
     !
     wopos(2)=postp(1)%wopos(2,54)
     wopos(3)=postp(1)%wopos(3,54)
     
     if( INOTEMPTY ) then
        nullify(auxvar)
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar,  npoin, posttable_fw % main_table % nvar)
        !
        ! Lookup from postprocessing table
        !
        call chm_post_lookup(auxvar)
     endif    

    
     do ipostvar=1,posttable_fw % main_table % nvar
        if( INOTEMPTY ) then
            gesca => auxvar(:,ipostvar)
        endif
        wopos(1)=posttable_fw % main_table % varname(ipostvar)
        call postpr(gesca,wopos,ittim,rutim)
     end do
     if( INOTMASTER ) nullify(gesca)
     
     if( INOTEMPTY ) then
        call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
     endif    
     
     return

  case(55_ip)
    !
    ! AVPOT:  Average posttab properties
    !
    wopos(2)=postp(1)%wopos(2,55)
    wopos(3)=postp(1)%wopos(3,55)
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm

       do ipostvar=1,posttable_fw % main_table % nvar
          do ipoin=1,npoin
             unkno(ipoin)                  = avposttab_chm(ipoin,ipostvar) / auxi
             avposttab_chm(ipoin,ipostvar) = 0.0_rp
          end do
          if( INOTEMPTY ) then
              gesca => unkno
          endif
          if (posttable_fw % main_table % varname(ipostvar)(1:2) == 'YK') then
             wopos(1)='AV'//posttable_fw % main_table % varname(ipostvar)(3:5)
          else
             wopos(1)='AV'//posttable_fw % main_table % varname(ipostvar)(1:3)
          endif
          call postpr(gesca,wopos,ittim,rutim)
       end do
       if( INOTMASTER ) nullify(gesca)
    endif

    return

  case(56_ip)
      !
      ! Conditional mass fractions for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,56)
         wopos(3)=postp(1)%wopos(3,56)

#ifdef CANTERA
         do imixf = 1, nZ_CMC_chm
            do iclas = 1, nclas_chm
               if( INOTMASTER ) gesca => Yk_CMC_chm(imixf,1:npoin,iclas)
               call getSpeciesName(gas_chm,iclas,aux_name)
               !write(valZ,'(F7.5)') Z_CMC_chm(imixf)
               !wopos(1) = trim(aux_name) // '_Z' // valZ
               write(valZ,'(I2.2)') imixf
               !!!!wopos(1) = trim(aux_name) // '_' // valZ
               !!!!! ELIMINAR
               wopos(1) = trim(aux_name) // valZ

               call postpr(gesca,wopos,ittim,rutim)
            end do
         end do
#endif
         if( INOTMASTER ) nullify(gesca)
      end if


  case(57_ip)
      !
      ! Unconditional mass fractions for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,57)
         wopos(3)=postp(1)%wopos(3,57)

#ifdef CANTERA
         do iclas = 1, nclas_chm
            if( INOTMASTER ) gesca => Yk_int_CMC_chm(1:npoin,iclas)
            call getSpeciesName(gas_chm,iclas,wopos(1))
            call postpr(gesca,wopos,ittim,rutim)
         end do
#endif
         if( INOTMASTER ) nullify(gesca)
      end if


  case(58_ip)
      !
      ! Conditional enthalpy for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,58)
         wopos(3)=postp(1)%wopos(3,58)

         do imixf = 1, nZ_CMC_chm
            if( INOTMASTER ) then
                if (kfl_solve_enth_CMC_chm /= 0) then
                   gesca => enthalp_CMC_chm(imixf,1:npoin)
                end if
            end if
            !write(valZ,'(F7.5)') Z_CMC_chm(imixf)
            !wopos(1) = 'h_Z' // valZ
            write(valZ,'(I2.2)') imixf
            wopos(1) = 'h_' // valZ
            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if


  case(59_ip)
      !
      ! Unconditional enthalpy for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         if( INOTMASTER ) gesca => enthalp_int_CMC_chm(1:npoin)
      end if

  case(60_ip)
      !
      ! Conditional temperature for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,60)
         wopos(3)=postp(1)%wopos(3,60)

         do imixf = 1, nZ_CMC_chm
            if( INOTMASTER ) gesca => temp_CMC_chm(imixf,1:npoin)
            !write(valZ,'(F7.5)') Z_CMC_chm(imixf)
            !wopos(1) = 'T_Z' // valZ
            write(valZ,'(I2.2)') imixf
            wopos(1) = 'T_' // valZ

            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if

  case(61_ip)
      !
      ! Unconditional temperature for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         if( INOTMASTER ) gesca => temp_int_CMC_chm(1:npoin)
      end if


  case(62_ip)
      !
      ! Conditional mass fractions source terms for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,62)
         wopos(3)=postp(1)%wopos(3,62)

#ifdef CANTERA
         do imixf = 1, nZ_CMC_chm
            do iclas = 1, nclas_chm
               if( INOTMASTER ) gesca => src_Yk_CMC_chm(imixf,1:npoin,iclas)
               call getSpeciesName(gas_chm,iclas,aux_name)
               !write(valZ,'(F7.5)') Z_CMC_chm(imixf)
               !wopos(1) = 'w' // trim(aux_name) // '_Z' // valZ
               write(valZ,'(I2.2)') imixf
               wopos(1) = 'w' // trim(aux_name) // '_' // valZ

               call postpr(gesca,wopos,ittim,rutim)
            end do
         end do
#endif
         if( INOTMASTER ) nullify(gesca)
      end if

  case(63_ip)
      !
      ! Unconditional mass fractions source terms for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,63)
         wopos(3)=postp(1)%wopos(3,63)

#ifdef CANTERA
         do iclas = 1, nclas_chm
            if( INOTMASTER ) gesca => src_Yk_int_CMC_chm(1:npoin,iclas)
            call getSpeciesName(gas_chm,iclas,aux_name)
            wopos(1) = 'w' // trim(aux_name)
            call postpr(gesca,wopos,ittim,rutim)
         end do
#endif
         if( INOTMASTER ) nullify(gesca)
      end if

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1)%wopos(:,ivari),MESH_ID=imesh)

end subroutine chm_outvar 
