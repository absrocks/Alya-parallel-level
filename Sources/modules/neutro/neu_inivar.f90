!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_inivar.f90
!> @author  Guillaume Houzeaux
!> @brief   Initialize some variables
!> @details Initialize some variables\n
!>          ITASK=0 ... When starting the run (from Turnon)\n
!>          ITASK=1 ... After reading data\n
!> @} 
!------------------------------------------------------------------------
subroutine neu_inivar(itask)
  use def_parame
  use def_master
  use def_neutro
  use def_domain
  use def_solver
  use mod_ADR,   only : ADR_initialize_type
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)   
     !
     ! Postprocess Variable 
     ! name:   '.....' 
     ! type:   'SCALA'/'VECTO' 
     ! entity: 'NPOIN'/'NELEM'
     !
     postp(1) % wopos (1:3, 1) = (/ 'CURRE' , 'VECTO' , 'NPOIN' /)
     postp(1) % wopos (1:3, 2) = (/ 'FLUX ' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 3) = (/ 'RADIA' , 'SCALA' , 'NPOIN' /)
     !
     ! Element set variables 
     !
     postp(1) % woese ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % woese ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % woese ( 3)     = 'VARI3'   ! Variable 3
     !
     ! Boundary set variables
     !
     postp(1) % wobse ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wobse ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wobse ( 3)     = 'VARI3'   ! Variable 3
     !
     ! Node set variables
     !
     postp(1) % wonse ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wonse ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wonse ( 3)     = 'VARI3'   ! Variable 3
     !
     ! Witness variables 
     !
     postp(1) % wowit ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wowit ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wowit ( 3)     = 'VARI3'   ! Variable 3
     !
     ! Solvers
     !     
     call soldef(-1_ip)                           ! Allocate memory for NUM_SOLVERS solvers
     solve(1:1) % kfl_solve = 1                   ! Output flag
     solve(1) % ndofn       = 1                   ! dof
     solve(1) % wprob       = 'NEUTRON_RADIATION' ! Solver name
     !
     ! Nullify pointers
     !
     nullify(kfl_fixno_neu)
     nullify(kfl_fixbo_neu)
     nullify(bvess_neu)
     nullify(bvnat_neu)
     nullify(tncod_neu)    
     nullify(tbcod_neu)
     !
     ! Dimensions
     !
     nunkn_neu = num_energies_neu * num_directions_neu
     !
     ! ADR type
     !
     call ADR_initialize_type(ADR_neu)

  case(1_ip)   
     !
     ! Derived parameters
     !
     if( ADR_neu % kfl_time_integration /= 0 ) then
        ncomp_neu = 3
     else
        ncomp_neu = 2
     end if
     nprev_neu = min(3_ip,ncomp_neu)  ! Last time step or global iteration

  end select

end subroutine neu_inivar
