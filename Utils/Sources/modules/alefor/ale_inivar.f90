subroutine ale_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_inivar
  ! NAME 
  !    ale_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use def_solver
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'DISPM'
     postp(1) % wopos (1, 2) = 'VELOM'
     postp(1) % wopos (1, 3) = 'COALE'
     postp(1) % wopos (1, 4) = 'GROUP'
     postp(1) % wopos (1, 5) = 'BVESS'

     postp(1) % wopos (2, 1) = 'VECTO'
     postp(1) % wopos (2, 2) = 'VECTO'
     postp(1) % wopos (2, 3) = 'VECTO'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'VECTO'
     !
     ! Solver
     ! Displacement is solved separately
     !     
     call soldef(-1_ip)
     solve(1) % wprob     = 'DISPLACEMENT'    ! Equation name
     solve(1) % kfl_solve = 1                 ! Output flag
     solve(1) % ndofn     = ndime             ! Dimension
     !
     ! Nullify pointers
     !
     nullify(lnodb_ad)      
     nullify(ltypb_ad)       
     nullify(coord_ad)       
     nullify(kfl_funty_ale)  
     nullify(funpa_ale)      
     nullify(tncod_ale)      
     nullify(tgcod_ale)      
     nullify(tbcod_ale)      
     nullify(kfl_funno_ale)  
     nullify(kfl_funbo_ale)  
     nullify(kfl_fixbo_ale)  
     nullify(coord_ale)      
   
  case(1_ip)

  end select

end subroutine ale_inivar
