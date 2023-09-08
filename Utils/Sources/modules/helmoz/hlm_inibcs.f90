!-----------------------------------------------------------------------
!> @addtogroup HelmozInput
!> @ingroup    Helmoz
!> @{
!> @file    hlm_outerrr.f90
!> @author  Guillaume Houzeaux
!> @date    08/02/2016
!> @brief   Boundary conditions
!> @details Read boundary conditions
!> @} 
!-----------------------------------------------------------------------

subroutine hlm_inibcs()

  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_helmoz
  use mod_opebcs
  use def_kermod
  implicit none
  integer(ip) :: itotn

  if( INOTMASTER ) then
     !
     ! Allocate memory
     !
     call hlm_membcs(1_ip)

     if( kfl_edges_hlm == 1 ) then
        !
        ! Edge codes (requires boundary codes)
        !
        if( kfl_icodb > 0 ) then
           ifloc     =  0
           ifbop     =  0
           ifbes     =  0
           iffun     =  0
           kfl_fixno => kfl_fixno_hlm
           tncod     => tecod_hlm
           call reacod(IMPOSE_EDGE_CODES)
!do itotn = 1,meshe(ndivi)%nedge
!print*,kfl_fixno_hlm(1,itotn)
!end do
!stop
        end if
     else
        !
        ! Node codes
        !
        if( kfl_icodn > 0 ) then
           ifloc     =  0
           ifbop     =  0
           ifbes     =  0
           iffun     =  0
           kfl_fixno => kfl_fixno_hlm
           tncod     => tncod_hlm
           call reacod(IMPOSE_NODE_CODES)
        end if
     end if

     do itotn = 1,size(kfl_fixno_hlm,2)
        if( kfl_fixno_hlm(1,itotn) == -1 ) kfl_fixno_hlm(1,itotn) = 0
     end do

  end if

end subroutine hlm_inibcs
