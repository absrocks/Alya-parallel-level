!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_gtable.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Get kro & krw from table
!> @details Get kro & krw from table
!> @} 
!------------------------------------------------------------------------
subroutine por_gtable(gpswa,kro,krw,dkro,lmate,kfl_d_kro)
!
  use def_kintyp, only       :  ip,rp
  use def_porous, only       :  tkrel_por,nrows_por
  implicit none 
  integer(ip), intent(in)    :: lmate        !> material
  integer(ip), intent(in)    :: kfl_d_kro    !> calculate d kro/d Sw
  real(rp),    intent(in)    :: gpswa        !> gauss point watter saturation
  real(rp),    intent(out)   :: kro,krw,dkro !> kro,krw,dkro/dSw

  integer(ip)                :: nrows,irows
  real(rp)                   :: delsw,delkw,delko

  dkro = 0.0_rp
  nrows = nrows_por(lmate)
  if ( gpswa <= tkrel_por(1,1,lmate) ) then
     krw = tkrel_por(2,1,lmate)
     kro = tkrel_por(3,1,lmate)
     if(kfl_d_kro==1) dkro = ( tkrel_por(3,2,lmate) - tkrel_por(3,1,lmate) ) /  &
          ( tkrel_por(1,2,lmate) - tkrel_por(1,1,lmate) )   ! Extrapolation

  else if ( gpswa >= tkrel_por(1,nrows,lmate) ) then
     krw = tkrel_por(2,nrows,lmate)
     kro = tkrel_por(3,nrows,lmate)
     if(kfl_d_kro==1) dkro = ( tkrel_por(3,nrows,lmate) - tkrel_por(3,nrows-1,lmate) ) /  &
          ( tkrel_por(1,nrows,lmate) - tkrel_por(1,nrows-1,lmate) )   ! Extrapolation
  else
     rows : do irows = 2,nrows
        if ( gpswa <= tkrel_por(1,irows,lmate) ) then
           delsw = tkrel_por(1,irows,lmate) - tkrel_por(1,irows-1,lmate)
           if (delsw == 0.0_rp) call runend ('POR_GTABLE: zero difference between 2 water saturation values not valid') 
           delkw = tkrel_por(2,irows,lmate) - tkrel_por(2,irows-1,lmate)
           delko = tkrel_por(3,irows,lmate) - tkrel_por(3,irows-1,lmate)
           krw = tkrel_por(2,irows-1,lmate) + ( gpswa - tkrel_por(1,irows-1,lmate) ) *delkw/delsw
           kro = tkrel_por(3,irows-1,lmate) + ( gpswa - tkrel_por(1,irows-1,lmate) ) *delko/delsw
           if(kfl_d_kro==1) dkro = delko / ( tkrel_por(1,irows,lmate) - tkrel_por(1,irows-1,lmate) )
           exit rows
        end if
     end do rows
  end if

end subroutine por_gtable
