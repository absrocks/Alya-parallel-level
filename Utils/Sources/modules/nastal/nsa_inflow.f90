subroutine nsa_inflow(ibset,ubulk)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_inflow
  ! NAME 
  !    nsa_inflow
  ! DESCRIPTION
  !    Inflow boundary condition that ensures a constant mass flow at the inlet
  ! USES
  !    nsa_bouset
  ! USED BY
  !    nsa_setvar
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use mod_iofile

  use def_nastal

  implicit none
  integer(ip)              :: dummi
  real(rp),    intent(out) :: ubulk 
  integer(ip), intent(in)  :: ibset

  !                                               
  ! Compute integrals over the inlet boundary: 
  !     rho*u*A = \int_{rho * u}dA
  !     rho*A   = \int_{rho }dA 
  !
  if( INOTMASTER ) then
    call nsa_bouset(lbsec(ibset),ibset)
  end if

  call posdef(22_ip,dummi) !! MPI_allreduce

  !
  ! Impose a conservation of the mass flow
  !
  if( INOTMASTER ) ubulk = inlet_nsa /vbset( 2,ibset)
 
end subroutine nsa_inflow
