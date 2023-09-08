!-----------------------------------------------------------------------
!> @addtogroup ExmediIonicurrents
!> @ingroup    Exmedi
!> @{
!> @file    mod_exm_ionicurrents.f90
!> @author  Mariano Vazquez
!> @date    28/11/2016
!> @brief   Module to compute ionicurrents  
!> @details Module to compute ionicurrents  
!-----------------------------------------------------------------------

module mod_exm_ionicurrents
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi
  use mod_messages, only : livinf
  use mod_timings,                     only : timings_assembly
  implicit none

  private
  real(rp) :: time1,time2

  integer(ip), parameter :: VECTOR = 1

  !
  ! Public stuff
  !
  public :: exm_nodalionicurrents


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    28/11/2016
  !> @brief   Compute nodal ionicurrents
  !> @details Compute nodal ionicurrents
  !-----------------------------------------------------------------------
  
  subroutine exm_nodalionicurrents
    implicit none

    integer(ip) :: ipoin,imate,imate_ipoin,kmodel_check, kmodel_ipoin,isubstep
    real(rp)    :: xioni,dioni,capmem,k1,k2,k3,k4,dv,dtimeEP,dtimeEPaux

    integer(ip), parameter :: EXPLICIT_SCHEME=1, IMPLICIT_SCHEME=2

    character(300)           :: messa
    
    time1     =  0.0_rp
    time2     =  0.0_rp

    kmodel_check= 0
    do imate= 1,nmate
       kmodel_check = max(kfl_cellmod(imate),kmodel_check)
    end do

    if (INOTSLAVE) then
       if (kmodel_check == 0)  then
            call livinf(165_ip,' ONLY DIFFUSIVE... ',0_ip)
       else
            call livinf(165_ip,' CELL... ',0_ip)
       endif
    end if

    dtimeEPaux= dtime !EP

    call cputim(time1)
    
    do isubstep= 1,msste_exm  ! ionic model substepping
       
       dtimeEP= dtimeEPaux / real(msste_exm)
       
       if (INOTSLAVE) then
          messa = '(' //trim(intost(isubstep))// ')'
          call livinf(165_ip,trim(messa),0_ip)
       end if
       
       if (INOTMASTER) then

          do ipoin=1,npoin              

             imate_ipoin  = nodemat(ipoin)
             kmodel_ipoin = kfl_cellmod(imate_ipoin)

             if (kmodel_ipoin == CELL_NOMOD_EXMEDI) then                       !no model
                xioni = 0.0_rp
                dioni = 0.0_rp

             else if (kmodel_ipoin == CELL_FITZHUGH_EXMEDI) then          !fitzhugh nagumo   

                call exm_mdfhn_ionicurrents(ipoin,xioni,dioni)

             else if (kmodel_ipoin == CELL_TT2006_EXMEDI) then  !TenTuscher 2006 Heterogeneo

                call exm_tt2006_ionicurrents(ipoin,xioni,dioni)  

             else if (kmodel_ipoin == CELL_OHARA_EXMEDI) then  !OHara-Rudy 2011
                call exm_oharaf_ionicurrents(ipoin,xioni,dioni)  !!! compute CURRENTS O'Hara-Rudy

             else if (kmodel_ipoin == CELL_FENTON_EXMEDI) then

                call runend ('EXM_IONICURRENTS: FENTON-KARMA MODELS TO BE PROGRAMMED')

             else if (kmodel_ipoin == CELL_SCATRIA_EXMEDI) then  !!Paci et al. 2013  hiPSC-Atrial like
               
                call exm_scatri_ionicurrents(ipoin,xioni,dioni)  !!! compute CURRENTS Atria-like
                
            else if (kmodel_ipoin == CELL_SCVENTRI_EXMEDI) then  !!Paci et al. 2013  hiPSC-Ventricular like
            
                call exm_scvent_ionicurrents(ipoin,xioni,dioni)  !!! compute CURRENTS Ventricular-like
            
            else

                call runend('EXM_IONICURRENTS: IONIC CURRENT FOR SELECTED MODEL NOT PROGRAMMED')

            end if


            ! (total cell ionic current * conversion_msecs_a_secs) / Cm
            !          ticel_exm(ipoin) =   xioni / xmccmmate_exm(kmodel_ipoin)
            ! (jacobian of the total cell ionic current * conversion_msecs_a_secs) / Cm
            jicel_exm(ipoin) =   dioni / xmccmmate_exm(imate_ipoin)

            !
            ! Update ionic currents, including applied ones
            !
            if (kfl_paced_exm == 0_ip) then
                capmem = 0.0_rp
            else
                capmem = appfi_exm(ipoin)
            end if

            if (kfl_appva_exm == 1_ip) capmem = 0.0_rp  ! when the stimuli is in voltage, it goes straightly to elmag

            if ((kmodel_ipoin > CELL_FENTON_EXMEDI) .and. (kmodel_ipoin < CELL_SCATRIA_EXMEDI)) then
                dtimeEP=dtime*1000.0_rp
            else
                dtimeEP = dtime
            end if

            dv = -(xioni + capmem) / xmccmmate_exm(imate_ipoin)          

            ! Forward Euler for the time strang splitting
            if (kfl_ticel_exm == 1_ip) then
               ticel_exm(ipoin) = dtimeEP * dv
            ! Runge-Kutta for the time strang splitting
            else if (kfl_ticel_exm == 2_ip) then
                   call runend("EXM_NODALIONICURRENTS: RUNGE-KUTTA TO BE PROGRAMMED. USE EULER.")
            end if

          end do

       end if

    end do

    call cputim(time2)
  
   
    !
    ! CPU times
    !
    cpu_ass_sol_exm(2) = time2 - time1
    call timings_assembly(cpu_ass_sol_exm(2), TYPE_OF_ASSEMBLY='NODES')

  end subroutine exm_nodalionicurrents

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

end module mod_exm_ionicurrents
!> @}
