!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_inivar.f90
!> @author  Mariano Vazquez
!> @brief   Initialize data
!> @date    16/11/1966
!> @details Initialize data
!> @} 
!-----------------------------------------------------------------------
subroutine exm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/exm_inivar
  ! NAME 
  !    exm_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use def_solver
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iconc,iauxi,ivari

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !

     ! VARIABLE NAMES:
     
     postp(1) % wopos ( 1, 1) = 'FIBER'
     postp(1) % wopos ( 1, 2) = 'GRAFI'
     postp(1) % wopos ( 1, 4) = 'IOCON'   ! VCONC (ION CONCENTRATIONS), A MULTIDIMENSIONAL NODAL FIELD

     postp(1) % wopos ( 1,20) = 'INTRA'
     postp(1) % wopos ( 1,21) = 'EXTRA'
     postp(1) % wopos ( 1,22) = 'RECOV'
     postp(1) % wopos ( 1,23) = 'XFIBE'
     postp(1) % wopos ( 1,24) = 'YFIBE'
     postp(1) % wopos ( 1,25) = 'ZFIBE'

     postp(1) % wopos ( 1,26) = 'ISOCH'
     postp(1) % wopos ( 1,27) = 'CECO1'
     postp(1) % wopos ( 1,28) = 'CECO2'
     postp(1) % wopos ( 1,29) = 'CECO3'
!+MRV
     postp(1) % wopos ( 1,30) = 'FIBE2'     
!-MRV
     postp(1) % wopos ( 1,31) = 'CURRE'  ! VICEL (CURRENTS), A MULTIDIMENSIONAL NODAL FIELD
     postp(1) % wopos ( 1,32) = 'QNETO'  ! VICEL (CURRENTS), A MULTIDIMENSIONAL NODAL FIELD

     ! VARIABLE TYPES:

     postp(1) % wopos ( 2, 1) = 'VECTO'
     postp(1) % wopos ( 2, 2) = 'VECTO'
     postp(1) % wopos ( 2, 4) = 'MULTI'

     postp(1) % wopos ( 2,20) = 'SCALA'
     postp(1) % wopos ( 2,21) = 'SCALA'
     postp(1) % wopos ( 2,22) = 'SCALA'
     postp(1) % wopos ( 2,23) = 'SCALA'
     postp(1) % wopos ( 2,24) = 'SCALA'
     postp(1) % wopos ( 2,25) = 'SCALA'

     postp(1) % wopos ( 2,26) = 'SCALA'
     postp(1) % wopos ( 2,27) = 'VECTO'
     postp(1) % wopos ( 2,28) = 'VECTO'
     postp(1) % wopos ( 2,29) = 'VECTO'
!+MRV
     postp(1) % wopos ( 2,30) = 'VECTO'     
!-MRV

     postp(1) % wopos ( 2,31) = 'MULTI'
     postp(1) % wopos ( 2,32) = 'SCALA'

     
     do iconc= 1,nconc_exm
        ivari= 32 + iconc
        if (iconc.le.9) then
           postp(1) % wopos ( 1,ivari) = 'VCO0'//trim(intost(iconc))
        else
           postp(1) % wopos ( 1,ivari) = 'VCO'//trim(intost(iconc))
        end if
        postp(1) % wopos ( 2,ivari) = 'SCALA'
     end do
     
     do iauxi= 1,nauxi_exm
        ivari= 32 + nconc_exm + iauxi
        if (iconc.le.9) then
           postp(1) % wopos ( 1,ivari) = 'VAU0'//trim(intost(iconc))
        else
           postp(1) % wopos ( 1,ivari) = 'VAU'//trim(intost(iconc))
        end if
        postp(1) % wopos ( 2,ivari) = 'SCALA'
     end do     

     !
     ! Sets variables
     !

     postp(1) % woese (  1) = 'INTRA'
     postp(1) % woese (  2) = 'EXTRA'
     postp(1) % woese (  3) = 'RECOV' 
     postp(1) % woese (  4) = 'QNETO'         
     postp(1) % wonse (  1) = 'INTRA'
     postp(1) % wonse (  2) = 'EXTRA'
     postp(1) % wonse (  3) = 'RECOV'    
     postp(1) % wonse (  4) = 'QNETO'  
     postp(1) % wobse (  1) = 'INTRA'
     postp(1) % wobse (  2) = 'EXTRA'
     postp(1) % wobse (  3) = 'RECOV'    
     postp(1) % wobse (  4) = 'QNETO'   
     !
     ! Witness variables
     !
     postp(1) % wowit (1)     = 'INTRA'
     postp(1) % wowit (2)     = 'QNETO'
          !
     ! Pseudo-ecg
     !
     kcopeecg_exm = 0
     !
     ! Solvers
     !     
     call soldef(-1_ip)

     solve(1) % wprob     = 'ACTIVATION_POTENTIAL'
     solve(1) % kfl_solve = 1
     solve(1) % ndofn     = 1

     ! In the case of implicit, these values will be defined in nsa_reanut, when calling reasol
     ! These are default values, corresponding to explicit with global time step
     solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON            
     solve(1) % kfl_preco = SOL_CLOSE_MASS_MATRIX

     ! Time variables
     kfl_timei = 1      ! exmedi is always transient
     
     !
     ! Materials
     !
     nmate_exm =  nmate
     lmate_exm => lmate
     !lcell_exm => lcell

     epres = 0.0_rp

     sms_conversion_currents_exm = 1000.0_rp

     cpold_exm= 0.0_rp
     
  end select

end subroutine exm_inivar
