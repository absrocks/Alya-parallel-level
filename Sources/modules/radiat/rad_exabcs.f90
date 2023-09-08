subroutine rad_exabcs(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_exabcs
  ! NAME 
  !    rad_exabcs
  ! DESCRIPTION
  !    This routine sets the boundary conditions for exact solution

  ! USED BY
  !    rad_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo
  real(rp)                :: dummr(3)
!!$
!!$  if(kfl_exacs_rad/=0) then
!!$
!!$     if(itask==1) then
!!$        !
!!$        ! Change fixity KFL_FIXNO_TEM: put Dirichlet bc on whole boundary
!!$        !
!!$        do ipoin=1,npoin
!!$           kfl_fixno_rad(1,ipoin)=0
!!$           ibopo=lpoty(ipoin)
!!$           if(ibopo/=0) kfl_fixno_rad(1,ipoin)=1
!!$        end do
!!$     end if
!!$     !
!!$     ! Prescribe rature on BVESS_TEM
!!$     !
!!$     do ipoin=1,npoin
!!$        if(kfl_fixno_rad(1,ipoin)==1)&
!!$             call rad_exacso(&
!!$             1_ip,coord(1,ipoin),dummr,dummr,&
!!$             dummr,dummr,dummr,dummr,bvess_rad(ipoin,1),dummr,dummr)
!!$     end do
!!$
!!$  end if

end subroutine rad_exabcs
