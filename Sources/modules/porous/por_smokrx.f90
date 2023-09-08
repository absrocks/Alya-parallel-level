!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_smokrx.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Smooth krw & kro
!> @details Smooth krw & kro
!> @} 
!------------------------------------------------------------------------
subroutine por_smokrx(smokr)
  use def_kintyp, only    :  ip,rp
  use def_domain, only    :  npoin,ltype,lnnod,ngaus,&
       &                     lnods,elmar,ndime,vmass,mnode,&
       &                     ntens,mgaus,lmate,coord,nelem
  use def_master, only    :  wasat,current_zone,INOTMASTER
  implicit none
  real(rp),   intent(out) :: smokr(2,npoin)            !> Smoothed values for krw and kro
  integer(ip)             :: pelty,pnode
  integer(ip)             :: pgaus,plapl
  integer(ip)             :: ielem,igaus,ipoin,kpoin
  real(rp)                :: elswa(mnode)              ! Gather 
  real(rp)                :: elrhs(2,mnode)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: gphes(ntens,mnode,mgaus)  ! dNk/dxidxj
  real(rp)                :: gpvol(mgaus)              ! |J|*w
  real(rp)                :: gpcar(ndime,mnode,mgaus)  ! dNk/dxj
  real(rp)                :: gpswa                     ! Water Saturation
  real(rp)                :: kro,krw,dummr,oovvm

  if( INOTMASTER ) then

     do ipoin = 1,npoin
        smokr(1,ipoin) = 0.0_rp
        smokr(2,ipoin) = 0.0_rp
     end do

     elements: do ielem = 1,nelem
        pelty = ltype(ielem)

        if( pelty > 0 ) then
           !
           ! Element dimensions
           !
           pnode = lnnod(ielem)
           pgaus = ngaus(pelty)
           plapl = 0
           !
           ! Gather operations
           !
           elswa(1:pnode)         = wasat(lnods(1:pnode,ielem),1)
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           elrhs(1,1:pnode)       = 0.0_rp
           elrhs(2,1:pnode)       = 0.0_rp
           !
           ! Cartesian derivatives and volume: GPCAR, PGVOL
           !
           call elmcar(&
                pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                gphes,ielem)

           do igaus = 1,pgaus
              !
              ! Water Saturation: GPSWA
              !
              gpswa = dot_product(elmar(pelty) % shape(1:pnode,igaus),elswa(1:pnode))
              !
              ! Obtain kro & krw from gpswa
              !
              call por_gtable(gpswa,kro,krw,dummr,lmate(ielem),0_ip) 

              elrhs(1,1:pnode) = elrhs(1,1:pnode) + gpvol(igaus) * krw * elmar(pelty) % shape(1:pnode,igaus)
              elrhs(2,1:pnode) = elrhs(2,1:pnode) + gpvol(igaus) * kro * elmar(pelty) % shape(1:pnode,igaus)

           end do

           call assrhs(2_ip,pnode,lnods(:,ielem),elrhs,smokr)
        end if

     end do elements
     !
     ! Periodicity and Parall service
     !
     call rhsmod(2_ip,smokr)
     !
     ! Update smoothed values
     !
     do ipoin = 1,npoin
        oovvm          = 1.0_rp / vmass(ipoin)
        smokr(1,ipoin) = smokr(1,ipoin) * oovvm
        smokr(2,ipoin) = smokr(2,ipoin) * oovvm
     end do

  end if

end subroutine por_smokrx
