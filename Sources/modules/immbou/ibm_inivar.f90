subroutine ibm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_inivar
  ! NAME 
  !    ibm_inivar
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_solver

  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     postp(1) % wopos (1, 1) = 'DMESH'
     postp(1) % wopos (1, 2) = 'DISPM'
     postp(1) % wopos (1, 3) = 'LNDIB'
     postp(1) % wopos (1, 4) = 'WALLS'
     postp(1) % wopos (1, 5) = 'PREIB'
     postp(1) % wopos (1, 6) = 'VELIB'
     postp(1) % wopos (1, 7) = 'DISIB'
     postp(1) % wopos (1, 8) = 'DISTA'

     postp(1) % wopos (2, 1) = 'VECTO'
     postp(1) % wopos (2, 2) = 'VECTO'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'SCALA'
     postp(1) % wopos (2, 6) = 'VECTO'
     postp(1) % wopos (2, 7) = 'VECTO'
     postp(1) % wopos (2, 8) = 'SCALA'
     !
     ! Solver
     ! 
     call soldef(-1_ip) 
     solve(1) % ndofn     =  1
     solve(1) % kfl_solve =  1
     solve(1) % nrhss     =  ndime
     solve(1) % wprob     =  'MESH_DEFORMATION'
     !
     ! Others
     !
     kfl_stead_ibm        =  0 
     kfl_embed_ibm        =  0
     nwaib                =  0

  end select

end subroutine ibm_inivar
