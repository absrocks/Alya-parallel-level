subroutine ker_turnon(itask)
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_turnon
  ! NAME 
  !    ker_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the incompressible NS equations.
  !    - Write some info
  !    - Check errors
  !    - Allocate memory
  ! USES
  !    ker_openfi
  !    ker_reaphy
  !    ker_reabcs
  !    ker_reanut
  !    ker_reaous
  !    ker_outinf
  !    ker_outerr
  !    ker_memall
  ! USED BY
  !    Kermod
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_regularization, only: kfl_regularization, set_regularization_pars
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Just after reading the mesh (after readom.f90)
     !
     !-------------------------------------------------------------------
     !
     ! Initialize variables
     !
     call ker_inivar(0_ip)
     !
     ! Read the physical problem
     !
     call ker_readat()
     !
     ! Parall service
     !
     call ker_parall(1_ip)
     !
     ! Initialize variables
     !
     call ker_inivar(1_ip)  
     !
     ! Check errors
     !
     call ker_outerr()
     !
     ! set parameters if using regularization
     !
     if(kfl_regularization) call set_regularization_pars()

  case ( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Just after creating the domain (after domain.f90)
     !
     !-------------------------------------------------------------------
     !
     ! Allocate memory
     !
     call ker_memall()
     !
     ! Impose boundary conditions
     !
     call ker_inibcs()
     !
     ! Initialize variables
     !
     call ker_inivar(2_ip)
     !
     ! Parall service
     !
     call ker_parall(2_ip)

  end select

end subroutine ker_turnon
