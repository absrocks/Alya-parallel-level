subroutine rad_outhfl()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_outhfl
  ! NAME 
  !    rad_outhfl
  ! DESCRIPTION
  !    This routine computes the heat flux
  ! USES
  ! USED BY
  !    rad_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_memchk
  use mod_postpr
  use mod_gradie
  implicit none
!!$  integer(ip)             :: ipoin,ibopo,idime
!!$  integer(4)              :: istat
!!$  real(rp), allocatable   :: gradt(:,:)
!!$  !
!!$  ! Allocate memory
!!$  ! 
!!$  allocate(gradt(ndime,npoin),stat=istat)
!!$  call memchk(zero,istat,mem_modul(1:2,modul),'GRADT','rad_outhfl',gradt)
!!$  !
!!$  ! Compute temp gradients
!!$  !
!!$  if(lawtc_rad(1)==0.and.nmate==1) then
!!$     call gradie(tempr(1:npoin,1),gradt)
!!$     gradt=gradt*tcond_rad(1,1)
!!$  else if(lawtc_rad(1)>0.and.nmate==1) then
!!$     call rad_heatfl(gradt)
!!$  end if
!!$  !
!!$  ! Compute heat flux
!!$  !
!!$  do ibopo=1,nbopo
!!$     gesca(ibopo)=0.0_rp
!!$  end do
!!$  do ipoin=1,npoin
!!$     ibopo=lpoty(ipoin)
!!$     if(ibopo>=1) then
!!$        do idime=1,ndime
!!$           gesca(ipoin)=gesca(ipoin)&
!!$                +gradt(idime,ipoin)*exnor(idime,1,ibopo)
!!$        end do
!!$     else
!!$        gesca(ipoin)=0.0_rp
!!$     end if
!!$  end do
!!$  !
!!$  ! Deallocate memory
!!$  !
!!$  call memchk(two,istat,mem_modul(1:2,modul),'GRADT','rad_outhfl',gradt)
!!$  deallocate(gradt,stat=istat)
!!$  if(istat/=0) call memerr(two,'GRADT','rad_outhfl',0_ip)

end subroutine rad_outhfl
