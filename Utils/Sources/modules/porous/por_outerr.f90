!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_outerr.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Checks if there are errros and warnings
!> @details Checks if there are errros and warnings
!> @} 
!------------------------------------------------------------------------
subroutine por_outerr()
  use def_master
  use def_domain
  use def_porous
  use def_kermod
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  !
  ! Properties
  !
!  if( kfl_prope == 0 ) then
!     ierro = ierro+1
!     call outfor(1_ip,momod(modul)%lun_outpu,&
!          'WHEN USING PORPER PROPERTIES SHOUD BE DECALRED IN KERMOD')
!  end if
  !
  ! Check the transient evolution
  !
!  if(kfl_timei/=0) then
!     if(kfl_timei_por == 0) then
!        iwarn=iwarn+1
!        call outfor(2_ip,momod(modul)%lun_outpu,&
!             'STEADY PORPERATURE IN A TRANSIENT CALCULATION')
!     end if
!  end if
  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------
  call errors(3_ip,ierro,iwarn,'POR_OUTERR')

end subroutine por_outerr
