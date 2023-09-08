subroutine inieig()
  !-----------------------------------------------------------------------
  !****f* master/inieig
  ! NAME
  !    inieig
  ! DESCRIPTION
  !    This subroutine initializes the eigen solver arrays
  !    amatr ... Matrix stiffness
  !    cmass ... Consistent matrix mass
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_master 
  use def_solver
  implicit none  
  integer(ip)  :: izmat

  if( INOTMASTER ) then 
     
     do izmat = 1,eigen_sol(1)%nzmat 
        amatr(izmat) = 0.0_rp
     end do
     
  end if

end subroutine inieig
 
