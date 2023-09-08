subroutine qua_memall()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_memall
  ! NAME 
  !    qua_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    Shrodinger equation
  ! USES
  !    memchk
  !    mediso
  !    qua_memose
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_quanty
  use mod_memchk
  implicit none
  integer(ip) :: lodof,ielem,pelty,pgaus
  integer(4)  :: istat
  !
  ! Problem unknowns QUANTY initialization
  !

  if( INOTMASTER ) then
     !
     ! rhoon: density
     !

     allocate(rhoon(npoin,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHOON','qua_memall',rhoon)
     !
     ! Eigen solver
     !
     eigen_sol => eigen_qua
     call eigdef(4_ip)
     !
     ! Algebraic solver
     !
!     if( eigen_qua(1)%kfl_algso /= 2 ) then
        solve_sol => solve_qua
        call soldef(4_ip)  
             
!     else
!        solve_qua(1)%kfl_algso = 1
!        solve_qua(1)%kfl_symme = 0
!     end if

  else

     allocate(phion(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PHION','qua_memall',phion)
     allocate(rhoon(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RHOON','qua_memall',rhoon)

  end if

end subroutine qua_memall
      
