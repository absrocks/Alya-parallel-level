!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    nsa_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine rotates ahead or rotates back a nodal vector 
!> @details This routine rotates ahead or rotates back a nodal vector 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_rotunk(irose,istpo,ifipo,vetar)
!-----------------------------------------------------------------------
!****f* Nastal/nsm_rotunk
! NAME 
!    nsa_rotunk
! DESCRIPTION
!    This routine rotates ahead or rotates back a nodal vector 
!    using the appropiate rotation matrix. 
!    The sense of rotation is given by input integer irose:
!            irose = -1 --> rotate ahead, from GCB to LCB
!            irose =  1 --> rotate back, from LCB to GCB 
!    
!    Recall that:
!    GCB : Global Cartesian Basis
!    LCB : Local Curvilinear Basis
!    
! USES
!    mbvab0
! USED BY
!    
!***
!-----------------------------------------------------------------------
  use      def_domain
  use      def_master
  use      def_nastal
  implicit none
  integer(ip), intent(in) :: irose !< rotation sense
  integer(ip), intent(in) :: istpo !< starting point number
  integer(ip), intent(in) :: ifipo !< final point number
  real(rp)                :: vetar(3) !< rotated vector

  integer(ip) :: ipoin,ibopo,iroty
  real(rp)    :: venew(ndime),veold(ndime)

  if (kfl_local_nsa==0 .and. istpo < ifipo) return 
  
  if (INOTMASTER) then 

     do ipoin=istpo,ifipo
        
        ibopo=lpoty(ipoin)
        
        if(ibopo>0) then
           iroty=kfl_fixrs_nsa(ibopo)
           veold(1:ndime)=vetar(1:ndime)
           !
           ! Target Vector boundary conditions linked to exnor
           !           
           if(iroty==-1) then
              call nsa_rotvec(irose,veold,exnor(1,1,ibopo),venew,ndime)
              !
              ! Target Vector boundary conditions linked to skcos_nsa 
              ! (computed in nsa_autint)
              !           
           else if(iroty==-2) then
              call nsa_rotvec(irose,veold,skcos_nsa(1,1,ibopo),venew,ndime)
              !
              ! Target Vector boundary conditions linked to skcos
              !
           else if(iroty>=1) then
              call nsa_rotvec(irose,veold,skcos(1,1,iroty),venew,ndime)
           end if
           vetar(1:ndime)=venew(1:ndime)
        end if
        
     end do
  end if

end subroutine nsa_rotunk

