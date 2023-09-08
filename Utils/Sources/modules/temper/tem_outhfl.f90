subroutine tem_outhfl()
  !------------------------------------------------------------------------
  !****f* Temper/tem_outhfl
  ! NAME 
  !    tem_outhfl
  ! DESCRIPTION
  !    This routine computes the heat flux
  ! USES
  ! USED BY
  !    tem_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_memchk
  use mod_postpr
  use mod_gradie
  implicit none
  integer(ip)             :: ipoin,ibopo,idime
  integer(4)              :: istat
  real(rp), allocatable   :: gradt(:,:)
  !
  ! Allocate memory
  ! 
  allocate(gradt(ndime,npoin),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)
  !
  ! Compute temperature gradients
  !
 
  call tem_heatfl(gradt)

  !call grasca(tempe(1:npoin,1),gradt)
 
  !
  ! Compute heat flux
  !
  do ipoin=1,npoin
     ibopo=lpoty(ipoin)
     if(ibopo>=1) then
        do idime=1,ndime
           gesca(ipoin)=gesca(ipoin)&
                +gradt(idime,ipoin)*exnor(idime,1,ibopo)
        end do
     else
        gesca(ipoin)=0.0_rp
     end if
  end do
  !
  ! Deallocate memory
  !
  call memchk(two,istat,mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)
  deallocate(gradt,stat=istat)
  if(istat/=0) call memerr(two,'GRADT','tem_outhfl',0_ip)

end subroutine tem_outhfl
 
