!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_openfi.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Gets ALL the file names for porous and opens files
!> @details Gets ALL the file names for porous and opens files
!> @} 
!------------------------------------------------------------------------
subroutine por_openfi(itask)
  use def_porous
  use def_master
  use def_domain
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if(kfl_paral<=0) then

     if(kfl_rstar==2) then
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)

     case (2_ip)

        if (kfl_naked==0) then

        else

        end if



     case(4_ip)


     end select

  end if

end subroutine por_openfi

