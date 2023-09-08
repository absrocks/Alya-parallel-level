!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_gtabl2.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Get Pbh from table
!> @details Get Pbh from table
!> @} 
!------------------------------------------------------------------------
subroutine por_gtabl2(ctime,Pbh,jwell)
!
  use def_kintyp, only       :  ip,rp
  use def_porous, only       :  wvalu_por,nroww_por
  implicit none 
  integer(ip), intent(in)    :: jwell        !> Number of well
  real(rp),    intent(in)    :: ctime        !> current time
  real(rp),    intent(out)   :: Pbh          !> Pbh

  integer(ip)                :: nrows,irows
  real(rp)                   :: delti,delpb

  nrows = nroww_por(jwell)
  if ( ctime <= wvalu_por(1,1,jwell) ) then
     Pbh = wvalu_por(2,1,jwell)


  else if ( ctime >= wvalu_por(1,nrows,jwell) ) then
     Pbh = wvalu_por(2,nrows,jwell)

  else
     rows : do irows = 2,nrows
        if ( ctime <= wvalu_por(1,irows,jwell) ) then
           delti = wvalu_por(1,irows,jwell) - wvalu_por(1,irows-1,jwell)
           if (delti == 0.0_rp) call runend ('POR_GTABL2: zero difference between 2 time values not valid') 
           delpb = wvalu_por(2,irows,jwell) - wvalu_por(2,irows-1,jwell)
           Pbh = wvalu_por(2,irows-1,jwell) + ( ctime - wvalu_por(1,irows-1,jwell) ) *delpb/delti
           exit rows
        end if
     end do rows
  end if

end subroutine por_gtabl2
