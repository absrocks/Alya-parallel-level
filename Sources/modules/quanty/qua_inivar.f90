subroutine qua_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_inivar
  ! NAME 
  !    qua_inivar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use def_solver
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(0_ip)
     !
     ! Postprocess
     !
     wopos_qua(1,1) = 'EIGEN'
     wopos_qua(1,2) = 'DENSI'

     wopos_qua(2,1) = 'SCALA'
     wopos_qua(2,2) = 'SCALA'
     !
     ! Algebraic solver
     !
     call soldef(-1_ip)                           ! Initialize solver
     solve_qua => solve
     solve(1)%ndofn     = 1
     solve(1)%lun_cvgso = lun_conve_qua       ! Convergence unit
     solve(1)%kfl_solve = 1                   ! Output flag
     solve(1)%lun_solve = lun_solve_qua       ! Output unit
     solve(1)%wprob = 'QUANTY_EIGEN_SOLVER'
     !
     ! Eigen solver
     !
     allocate(eigen_qua(1),stat=istat)            ! Eigen Solver type
     eigen_sol => eigen_qua                       ! Eigen Solver pointer
     call eigdef(0_ip)                            ! Initialize eigen solver
     eigen_qua(1)%ndofn = 1
     eigen_sol(1)%lun_solei = lun_witne_qua

     ncomp_qua = 2                                ! Number of components
     kfl_symgr = 1

     ! estructuras de atomos y especies 
     allocate(atomo_qua(1),stat=istat)            ! atomo type


     atomo_qua(1)%wprob      = 'file NAME '
     atomo_qua(1)%natoms     = 1      ! nubmer of atoms
     atomo_qua(1)%nespec     = 1      ! nubmer of atoms


  case(1_ip)

  case(2_ip)

  case(5_ip)
     !
     ! Deflated Solver: groups composed by master
     !
     !solve_sol => solve_qua(1:)
     if( solve_sol(1)%kfl_algso == 2 ) then
        if( INOTMASTER ) solve_sol(1)%limpo => kfl_fixno_qua(1,1:)
        call cregro()
     end if

  end select

end subroutine qua_inivar
