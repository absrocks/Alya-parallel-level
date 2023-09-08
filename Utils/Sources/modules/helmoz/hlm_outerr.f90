!-----------------------------------------------------------------------
!> @addtogroup HelmozInput
!> @ingroup    Helmoz
!> @{
!> @file    hlm_outerrr.f90
!> @author  Guillaume Houzeaux
!> @date     08/02/2016
!> @brief   Error
!> @details Warning and erros
!> @} 
!-----------------------------------------------------------------------

subroutine hlm_outerr()

  use def_master
  use def_kermod
  use def_domain
  use def_helmoz
  use mod_outfor, only : outfor
  implicit none
  integer(ip)             :: ierro=0,iwarn=0,iboun,ipoin,idime,ibopo,ipara
  integer(ip)             :: iroty,dumm1,dumm2,ivara,ivar1,ivar2
  character(20)           :: messa
  character(200)          :: wmess
  !
  ! Set postprocess
  !
  if( kfl_edges_hlm /= 1 .and. kfl_edge_elements == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul) % lun_outpu,&
          'EDGE ELEMENTS MUST BE ACTIVATED IN KER.DAT FILE')
  end if

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,'MENSAJE VACIO???? ')

end subroutine hlm_outerr
