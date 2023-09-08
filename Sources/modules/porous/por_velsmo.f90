!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_velsmo.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Obtains the smoothed velocity from the pressure
!> @details Obtains the smoothed velocity from the pressure
!!          Similar to nsi_assdia
!> @} 
!------------------------------------------------------------------------
subroutine por_velsmo()
  !
  use def_kintyp, only :  ip,rp
  use def_domain, only :  npoin,ltype,lnnod,ngaus,&
       &                  llapl,lorde,ltopo,lnods,elmar,lmate,&
       &                  ndime,vmass,mnode,ntens,mgaus,&
       &                  coord,nelem
  use def_porous, only :  prref_por,comwa_por,comoi_por,bwref_por,&
       &                  boref_por,muoil_por,muwat_por,perme_por,&
       &                  denoi_por,denwa_por,gravi_por,grnor_por,&
       &                  ncomp_por,denhy_por
  use def_master, only :  veloc,ID_POROUS,lzone,mem_modul,modul,&
       &                  INOTMASTER,wasat,press,current_zone
  use mod_memory, only :  memory_alloca,memory_deallo

  implicit none
  real(rp)          :: deltp,xx_w,Bw,xx_o,Bo,kro,krw
  real(rp)          :: lambd_o,lambd_w,perme(ndime),oovvm
  real(rp)          :: dummr(1),rho,lambd
  integer(ip)       :: pelty,pnode,kelem,inode,idime
  integer(ip)       :: pgaus,plapl,porde,ptopo
  integer(ip)       :: ielem,igaus,ipoin,kpoin
  real(rp)          :: elpre(mnode)            
  real(rp)          :: elswa(mnode)            
  real(rp)          :: elrhs(ndime,mnode)
  real(rp)          :: elcod(ndime,mnode)
  real(rp)          :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)          :: gpvol(mgaus)                          ! |J|*w
  real(rp)          :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)          :: gppre                                 ! Pressure
  real(rp)          :: gpgpr(3)                              ! Pressure Gradient
  real(rp)          :: gpswa                                 ! Water Saturation
  real(rp), pointer :: rhsid(:,:)

  if( INOTMASTER ) then

     nullify(rhsid)
     call memory_alloca(mem_modul(1:2,modul),'RHSID','por_velsmo',rhsid,ndime,npoin)    

     elements: do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = lnnod(ielem)
           pgaus = ngaus(pelty)
           plapl = 0
           porde = lorde(pelty)
           ptopo = ltopo(pelty)
           !
           ! Gather operations
           !
           elswa(1:pnode)         = wasat(lnods(1:pnode,ielem),1)
           elpre(1:pnode)         = press(lnods(1:pnode,ielem),1)
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
           elrhs(1:ndime,1:pnode) = 0.0_rp
           perme(1:ndime)         = perme_por(1:ndime,ielem)
           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
           !
           call elmcar(&
                pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                gphes,ielem)

           do igaus = 1,pgaus
              !
              ! Pressure: GPPRE & Water Saturation: GPSWA & Coordinates
              !
              gppre = 0.0_rp
              gpswa = 0.0_rp
              gpgpr = 0.0_rp
              do inode = 1,pnode
                 gppre          = gppre          + elmar(pelty) % shape(inode,igaus) * elpre(inode)
                 gpswa          = gpswa          + elmar(pelty) % shape(inode,igaus) * elswa(inode)
                 gpgpr(1:ndime) = gpgpr(1:ndime) + gpcar(1:ndime,inode,igaus)        * elpre(inode)
              end do
              !
              ! Preliminaries Bo, Bw, kro, krw, phi
              !
              deltp   = gppre - prref_por
              xx_o    = comoi_por * deltp
              xx_w    = comwa_por * deltp
              Bo      = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
              Bw      = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
              call por_gtable(gpswa,kro,krw,dummr,lmate(ielem),0_ip) ! obtain kro & krw
              lambd_o = kro / muoil_por
              lambd_w = krw / muwat_por

              do inode = 1,pnode
                 rho   = ( lambd_o * ( ( denoi_por / Bo ) - denhy_por ) ) + ( lambd_w * ( ( denwa_por / Bw ) - denhy_por ) ) 
                 lambd = lambd_o + lambd_w 
                 do idime = 1,ndime
                    elrhs(idime,inode) = elrhs (idime,inode) + &
                         gpvol(igaus) * perme(idime) * elmar(pelty) % shape(inode,igaus) &
                         * ( rho   * gravi_por(idime) * grnor_por   &
                         -   lambd * gpgpr(idime)  ) 
                 end do
              end do
           end do

           call assrhs(ndime,pnode,lnods(:,ielem),elrhs,rhsid)

        end if

     end do elements
     !
     ! Periodicity and Parall service
     !
     call rhsmod(ndime,rhsid)
     !
     ! Obtain velocity
     !
     do ipoin = 1,npoin
        oovvm = 1.0_rp / vmass(ipoin)
        veloc(1:ndime,ipoin,1) = rhsid(1:ndime,ipoin) * oovvm
     end do

     call memory_deallo(mem_modul(1:2,modul),'RHSID','por_velsmo',rhsid )

  end if

end subroutine por_velsmo
