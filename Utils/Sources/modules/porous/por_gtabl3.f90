!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_gtabl3.f90
!> @date    26/08/2013
!> @author  Herbert Owen
!> @brief   Get Pbh (or flow rate) at all the nodes from table
!> @details Get Pbh (or flow rate) at all the nodes from table
!> @} 
!------------------------------------------------------------------------
subroutine por_gtabl3(ctime,Pbh,jwell,nheiw)
!
  use def_kintyp, only       :  ip,rp
  use def_porous, only       :  wvalu_por,nroww_por
  implicit none 
  integer(ip), intent(in)    :: jwell        !> Number of well
  integer(ip), intent(in)    :: nheiw        !> Number of nodes that form the well
  real(rp),    intent(in)    :: ctime        !> current time
  real(rp),    intent(out)   :: Pbh(nheiw)   !> Pbh for all the nodes that form the well

  integer(ip)                :: nrows,irows,iheig
  real(rp)                   :: delti,delpb

  nrows = nroww_por(jwell)
  if ( ctime <= wvalu_por(1,1,jwell) ) then
     do iheig = 1,nheiw
        Pbh(iheig) = wvalu_por(1+iheig,1,jwell)
     end do


  else if ( ctime >= wvalu_por(1,nrows,jwell) ) then
     do iheig = 1,nheiw
        Pbh(iheig) = wvalu_por(1+iheig,nrows,jwell)
     end do

  else
     rows : do irows = 2,nrows
        if ( ctime <= wvalu_por(1,irows,jwell) ) then
           delti = wvalu_por(1,irows,jwell) - wvalu_por(1,irows-1,jwell)
           if (delti == 0.0_rp) call runend ('POR_GTABL2: zero difference between 2 time values not valid') 
           do iheig = 1,nheiw
              delpb = wvalu_por(1+iheig,irows,jwell) - wvalu_por(1+iheig,irows-1,jwell)
              Pbh(iheig) = wvalu_por(1+iheig,irows-1,jwell) + ( ctime - wvalu_por(1,irows-1,jwell) ) *delpb/delti
           end do
           exit rows
        end if
     end do rows
  end if

end subroutine por_gtabl3
