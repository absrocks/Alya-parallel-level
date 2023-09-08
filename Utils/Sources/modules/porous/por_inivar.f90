!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_inivar.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Initialize some variables
!> @details Initialize some variables\n
!!          ITASK=1 ... When starting the run (from Turnon)\n
!!          ITASK=2 ... hhhpor_ old nsi First time step. This is needed as some variables
!!                      are not initialized before\n
!!          ITASK=3 ... hhhpor_ old nsiWhen starting a time step (from nsi_begste)\n
!> @} 
!------------------------------------------------------------------------
subroutine por_inivar(itask)

  use def_parame
  use def_master
  use def_domain
  use def_porous
  use def_solver
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'VELOC'
     postp(1) % wopos (1, 2) = 'PRESS'
     postp(1) % wopos (1, 3) = 'SATUR'
     postp(1) % wopos (1, 4) = 'IJKNO'
     postp(1) % wopos (1, 5) = 'POROS'
     postp(1) % wopos (1, 6) = 'WELLS'
     postp(1) % wopos (1, 7) = 'LINEL'
     postp(1) % wopos (1, 8) = 'GROUP'
     postp(1) % wopos (1, 9) = 'PERME'
 
     postp(1) % wopos (2, 1) = 'VECTO'
     postp(1) % wopos (2, 2) = 'SCALA'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'VECTO'
     postp(1) % wopos (2, 5) = 'SCALA'
     postp(1) % wopos (2, 6) = 'SCALA'
     postp(1) % wopos (2, 7) = 'SCALA'
     postp(1) % wopos (2, 8) = 'SCALA'
     postp(1) % wopos (2, 9) = 'VECTO'
     !
     ! Node set
     !
     postp(1) % wonse (1)    = 'WI   '  
     postp(1) % wonse (2)    = 'QT   '  
     postp(1) % wonse (3)    = 'QW   '  
     postp(1) % wonse (4)    = 'QO   '  
     postp(1) % wonse (5)    = 'PBH  '  
     postp(1) % wonse (6)    = 'DELTP'  
     postp(1) % wonse (7)    = 'SW   '  
     postp(1) % wonse (8)    = 'WELLN'  
     !
     ! Solver
     !     
     call soldef(-2_ip)                   ! Allocate memory : we will need 2 solvers
     solve(1) % kfl_solve = 1             ! Pressure:   Output flag
     solve(2) % kfl_solve = 1             ! Saturation: Output flag
     !
     ! Nullify pointers
     !
     nullify(nrows_por)
     nullify(nroww_por)
     nullify(iwell_por)
     nullify(ipwel_por)
     nullify(iheip_por)
     nullify(nheiw_por)
     nullify(kfl_wellc_por)

     nullify(tkrel_por)
     nullify(poro0_por)
     nullify(perme_por)
     nullify(satin_por)
     nullify(rwell_por)
     nullify(wvalu_por)
     nullify(dataw_por)
     nullify(tncod_por)
     nullify(kfl_fixno_por)
     nullify(bvess_por)
     nullify(nodpo_por)
     nullify(nodpe_por)
     nullify(winde_por)
     nullify(wmass_por)
     nullify(pboho_por)
 
  case(1_ip)
     !
     ! Dimensions
     !
     solve(1) % ndofn = 1
     solve(2) % ndofn = 1
     solve(1) % wprob = 'PRESSURE'
     solve(2) % wprob = 'SATURATION'
!    perhaps in the futre we may have a solver for veloc - we will need this if we smooth like Coutinho  

     do kprsa_por=1_ip,2_ip  ! 1 pres , 2 sat
        if(kfl_timei_por(kprsa_por)==1) then                      ! Number of velocity components
           if(kfl_tisch_por==1) then
              ncomp_por = max(ncomp_por,3)                              ! Trapezoidal rule
           else if(kfl_tisch_por==2) then
              ncomp_por = max(ncomp_por,2+kfl_tiacc_por(kprsa_por))               ! BDF scheme
           end if
        else
           ncomp_por = 2     
        end if
     end do
     !
     ! Solver fixity  
     !
!     if( INOTMASTER ) then        
!        solve(1) % bvess     => bvess_tem(:,:,1)
!        solve(1) % kfl_fixno => kfl_fixno_tem
!     end if

  case(2_ip)
     !
     ! Velocity subgrid scale
     !
     call runend('POR_INIVAR: velocity subgrid scales not ready')
!     do kprsa_por=1_ip,2_ip  ! 1 pres , 2 sat
!        if(associated(vesgs).and.kfl_advec_por(kprsa_por)==1) then
!           kfl_sgsve_por(kprsa_por)=1
!        else
!           kfl_sgsve_por(kprsa_por)=0
!        end if
!     end do

  end select

end subroutine por_inivar
