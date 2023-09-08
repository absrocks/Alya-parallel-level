!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_outinf.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Information on equation solved
!> @details Information on equation solved
!> @} 
!------------------------------------------------------------------------
subroutine por_outinf()
  use def_master
  use def_domain
  use def_porous
  implicit none
!  character(60) :: equat
!  integer(ip)   :: ierhs,imate
!  character(2)  :: wcpcv

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then
     end if

  end if
  !
  ! Formats
  !
110 format(/,&
       & 10x,a,//,&
       & 10x,'Physical properties: ',i2,' material(s)',/,&
       & 10x,'------------------------------------------------------------ ')
111 format(&
       & 10x,a,10(e13.6,1x))
112 format(&
       & 10x,a,e13.6,a)

end subroutine por_outinf

