subroutine arrzon(itask,ndim1,xarra)
  !
  ! Pass data from one zone to another through a contact zone
  !
  !
  ! Give force to solid nodes
  !
  ! 1.FORCF computed in nsi_solidz()
  !   x = assembled value
  !   y = local value (not complete)
  !
  ! o-------o o-------o
  ! |  NSI  | |  NSI  |   
  ! x-------y y-------x
  ! 0-------0 0-------0 
  ! |       | |       |   <= CONTACT
  ! 0-------0 0-------0
  ! 0-------0 0-------0
  ! |  SLD  | |  SLD  |  
  ! o-------o o-------o
  !
  ! => pararr('SLX')
  !
  ! o-------o o-------o
  ! |  NSI  | |  NSI  |   
  ! x-------x x-------x
  ! x-------x x-------x 
  ! |       | |       |   <= CONTACT
  ! 0-------0 0-------0
  ! 0-------0 0-------0
  ! |  SLD  | |  SLD  |  
  ! o-------o o-------o
  !
  ! => Copy only own boundary nodes from IPOIN (nastin) to JPOIN (solidz)
  !
  ! o-------o o-------o
  ! |  NSI  | |  NSI  |   
  ! x-------x x-------x
  ! x-------x x-------x 
  ! |       | |       |   <= CONTACT
  ! x-------x 0-------x
  ! 0-------0 0-------0
  ! |  SLD  | |  SLD  |  
  ! o-------o o-------o
  !
  ! => pararr('SLX')
  !
  ! o-------o o-------o
  ! |  NSI  | |  NSI  |   
  ! x-------x x-------x
  ! x-------x x-------x 
  ! |       | |       |   <= CONTACT
  ! x-------x x-------x
  ! x-------x x-------x
  ! |  SLD  | |  SLD  |  
  ! o-------o o-------o
  !
  ! The, solidz will only assemble its own-boundary contact nodes
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use def_inpout
  use mod_ker_proper 
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: ndim1
  real(rp),    intent(inout) :: xarra(ndim1,*)
  integer(ip)                :: isid1,isid2,incnt
  integer(ip)                :: ipoin,jpoin,idime

  if( INOTMASTER ) then

     if( itask == 1 ) then  ! Fluid to solid
        isid1 = 1   
        isid2 = 2
     else                   ! Solid to fluid
        isid1 = 2
        isid2 = 1
     end if

     call pararr('SLG',NPOIN_TYPE,ndime*npoin,xarra) 
     !
     ! If there is a contact element
     !
     call memgen(1_ip,npoin,0_ip)
     do incnt = 1,nncnt
        jpoin = lncnt(isid2,incnt)
        gisca(jpoin) = 1
     end do
     call parari('SLG',NPOIN_TYPE,npoin,gisca) 
     do incnt = 1,nncnt
        ipoin = lncnt(isid1,incnt) ! itask=1: Fluid node
        jpoin = lncnt(isid2,incnt) ! itask=1: Solid node
        do idime = 1,ndime 
           xarra(idime,jpoin) = xarra(idime,ipoin)
        end do
     end do
     call pararr('SLG',NPOIN_TYPE,ndime*npoin,xarra) 
     do ipoin = 1,npoin
        if( gisca(ipoin) /= 0 ) then
           do idime = 1,ndime
              xarra(idime,ipoin) = xarra(idime,ipoin) / real(gisca(ipoin),rp)
           end do
        end if
     end do
     call memgen(3_ip,npoin,0_ip)

  end if

end subroutine arrzon
 
