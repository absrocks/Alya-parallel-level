!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_averag.f90
!> @author  Daniel Mira
!> @date    13/06/204
!> @brief   Average variables: temperature
!> @details Average variables: temperature
!> @} 
!-----------------------------------------------------------------------
subroutine chm_averag()
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper
  use mod_postpr
  use def_kintyp,             only : ip,rp,r2p
  use mod_memory,             only : memory_alloca, memory_deallo
  use mod_solver,             only : solver_lumped_mass_system 
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip) :: ipoin,iclas_phi,dummi,ielem,igaus,pelty,pgaus,ipostvar
  real(rp)    :: zechm
  real(rp), pointer         :: prope_tmp(:)
  real(rp), pointer         :: auxvar(:,:)
  type(r1p),pointer         :: aux_r1p(:)

  zechm = epsilon(0.0_rp)

  if( cutim > avtim_chm ) then

     if( INOTMASTER ) then

        !
        ! AVY: average reaction progress Yc or C
        !
        if( output_postprocess_check_variable_postprocess(12_ip) ) then
          do ipoin=1,npoin
             avY_chm(ipoin) = avY_chm(ipoin)&
                                    + conce(ipoin,1,1) * dtime
          end do
        end if
        
        !
        ! AVYv: average variance of reaction progress Yc
        !
        if( output_postprocess_check_variable_postprocess(13_ip) ) then
          !
          ! Variance is transported directly
          !
          if ( kfl_varYc_chm == 1_ip ) then
             do ipoin=1,npoin
                avYv_chm(ipoin) = avYv_chm(ipoin)&
                                         + conce(ipoin,2,1) * dtime
             end do
          !
          ! Yc^2 is transported and variance must be computed
          !
          else if ( kfl_varYc_chm == 2_ip ) then
             do ipoin=1,npoin
                avYv_chm(ipoin) = avYv_chm(ipoin)&
                                + ( conce(ipoin,2,1) - conce(ipoin,1,1)*conce(ipoin,1,1) ) * dtime
             end do
          end if

        end if

        !
        ! AVZv: average variance of mixture fraction
        !
        if( output_postprocess_check_variable_postprocess(14_ip) ) then
          !
          ! Variance is transported directly
          !
          if ( kfl_varZ_chm == 1_ip ) then
             do ipoin=1,npoin
                avZv_chm(ipoin) = avZv_chm(ipoin)&
                                         + conce(ipoin,4,1) * dtime
             end do
          !
          ! Z^2 is transported and variance must be computed
          !
          else if ( kfl_varZ_chm == 2_ip ) then
             do ipoin=1,npoin
                avZv_chm(ipoin) = avZv_chm(ipoin)&
                                + ( conce(ipoin,4,1) - conce(ipoin,3,1)*conce(ipoin,3,1) ) * dtime
             end do
          end if

        end if

        !
        ! AVCHM: average chemical heat (CHEMICAL_HEAT)
        !
        if( output_postprocess_check_variable_postprocess(15_ip) ) then

           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_averag',aux_r1p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_averag',aux_r1p(ielem)%a,pgaus)
              aux_r1p(ielem) % a = chemical_heat(ielem) % a(:,1,1)
           end do

           nullify ( prope_tmp )
           allocate( prope_tmp(npoin) )
           call smooth (aux_r1p, prope_tmp)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_averag',aux_r1p)

           do ipoin=1,npoin
              avchm_chm(ipoin) = avchm_chm(ipoin)&
                                       + dtime*prope_tmp(ipoin)
           end do

           deallocate( prope_tmp )

        end if

        !
        ! AVZ: average mixture fraction Z
        !
        if( output_postprocess_check_variable_postprocess(16_ip) ) then
          do ipoin=1,npoin
             avZ_chm(ipoin) = avZ_chm(ipoin)&
                                    + conce(ipoin,3,1) * dtime
          end do
        end if

        !
        ! AVZ2: average squared of mixture fraction Z*Z
        !
        if( output_postprocess_check_variable_postprocess(17_ip) ) then
          do ipoin=1,npoin
             avZ2_chm(ipoin) = avZ2_chm(ipoin)&
                                      + conce(ipoin,3,1)*conce(ipoin,3,1) * dtime
          end do
        end if

        !
        ! AVY2: average squared of progress variable Yc*Yc or C*C
        !
        if( output_postprocess_check_variable_postprocess(31_ip) ) then
          do ipoin=1,npoin
             avY2_chm(ipoin) = avY2_chm(ipoin)&
                                      + conce(ipoin,1,1)*conce(ipoin,1,1) * dtime
          end do
        end if

        !
        ! AVL: average liquid volume fraction phi_L
        !
        if( output_postprocess_check_variable_postprocess(32_ip) ) then

          if ( kfl_spray_chm /= 0 ) then
             iclas_phi = nclas_chm - 1
             do ipoin=1,npoin
                   avL_chm(ipoin) = avL_chm(ipoin)&
                                      + conce(ipoin,iclas_phi,1) * dtime
             end do
          else
             call runend('CHEMIC CHM_AVERAG: Liquid volume fraction only valid with spray model')
          end if
        end if

        !
        ! AVL2: average liquid volume fraction squared phi_L*phi_L
        !
        if( output_postprocess_check_variable_postprocess(33_ip) ) then

          if ( kfl_spray_chm /= 0 ) then
             iclas_phi = nclas_chm - 1
             do ipoin=1,npoin
                   avL2_chm(ipoin) = avL2_chm(ipoin)&
                                            + conce(ipoin,iclas_phi,1)*conce(ipoin,iclas_phi,1) * dtime
             end do
          else
             call runend('CHEMIC CHM_AVERAG: Squared of liquid volume fraction only valid with spray model')
          end if
        end if

        !
        ! AVS: Average interface surface density Sigma
        !
        if( output_postprocess_check_variable_postprocess(34_ip) ) then

          do ipoin=1,npoin
             avS_chm(ipoin) = avS_chm(ipoin)&
                                      + dtime*Sigma_chm(ipoin)
          end do

        end if

        !
        ! AVS0: Average interface surface density Sigma_0
        !
        if( output_postprocess_check_variable_postprocess(35_ip) ) then

          do ipoin=1,npoin
             avS0_chm(ipoin) = avS0_chm(ipoin)&
                                      + dtime*Sigm0_chm(ipoin)
          end do

        end if

        !
        ! AVD32: Average Sauter mean diameter
        !
        if( output_postprocess_check_variable_postprocess(36_ip) ) then

          do ipoin=1,npoin
             avd32_chm(ipoin) = avd32_chm(ipoin)&
                                      + dtime*d32_chm(ipoin)
          end do

        end if

        !
        ! AVDEN: average density
        !
        if( output_postprocess_check_variable_postprocess(37_ip) ) then

           nullify ( prope_tmp )
           allocate( prope_tmp(npoin) )
           call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

           do ipoin=1,npoin
              avden_chm(ipoin) = avden_chm(ipoin) + prope_tmp(ipoin) * dtime
           end do

           deallocate( prope_tmp )

        end if

        !
        ! AVHRR: average heat release
        !
        if( output_postprocess_check_variable_postprocess(52_ip) ) then
           do ipoin=1,npoin
              hrr_avg_chm(ipoin) = hrr_avg_chm(ipoin) + hrr_chm(ipoin) * dtime
           end do
        end if

        !
        ! AVMSK: average mass source term from spray
        !
        if( output_postprocess_check_variable_postprocess(53_ip) ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp,npoin)
           if (associated(mass_sink)) then
              do ipoin=1,npoin
                 prope_tmp(ipoin) = mass_sink(ipoin)
              enddo
              call solver_lumped_mass_system(1_ip,prope_tmp,EXCHANGE=.false.)
           else
              do ipoin=1,npoin
                 prope_tmp(ipoin) = 0.0_rp 
              end do
           endif

           do ipoin=1,npoin
              avmsk_chm(ipoin) = avmsk_chm(ipoin) + prope_tmp(ipoin) * dtime
           end do
           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp)

        end if


        !
        ! AVPOT: average postprocessing table quantities
        !
        if( output_postprocess_check_variable_postprocess(55_ip) ) then
           if( INOTEMPTY ) then
              nullify(auxvar)
              call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar,  npoin, posttable_fw % main_table % nvar)
              !
              ! Lookup from postprocessing table
              !
              call chm_post_lookup(auxvar)
           endif    

           do ipostvar=1,posttable_fw % main_table % nvar
              do ipoin=1,npoin
                 avposttab_chm(ipoin,ipostvar) = avposttab_chm(ipoin,ipostvar)&
                                    + auxvar(ipoin,ipostvar) * dtime
              end do
           enddo
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
           endif    
        end if



     end if

  end if

end subroutine chm_averag
