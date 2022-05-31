subroutine nsi_elmset(iesec,ieset)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmset
  ! NAME 
  !    nsi_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are: 
  !    1. SETVO: set surface           = meas(W)=int_W
  !    2. SETVE: set mean vel. module  = int_W u^2 ]
  !    3. SETVR: set mean vort. module = int_W w^2 ]
  !                                      where w=dv/dx-du/dy
  !    SETVE and SETVR are normailzed further on in nsi_outset
  ! USES
  !    nsi_elmgat
  !    elmder
  ! USED BY
  !    nsi_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  use mod_ker_proper, only : ker_proper
  implicit none
  integer(ip), intent(in)  :: iesec,ieset
  real(rp),    pointer     :: setvo(:),setve(:),setvr(:)
  real(rp),    pointer     :: setki(:),divu2(:),setmo(:)
  integer(ip)              :: pnode,pgaus,pelty,nvabi,dummi
  integer(ip)              :: ielem,igaus,idime,inode,ipoin
  real(rp)                 :: gpcar(ndime,mnode,mgaus) 
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)                 :: elvel(ndime,mnode),elcod(ndime,mnode)
  real(rp)                 :: gpden(mgaus),divu
  real(rp)                 :: dvdx,dudy,gpvol,gpdet,gpvel(3),dummr
  !
  ! Initialization
  !
  nvabi =  postp(1) % nvaes+1
  setvo => veset( nvabi:nvabi ,ieset)  
  setvo =  0.0_rp  ! Set volume
  
  if( postp(1) % npp_setse(1) /= 0 ) setve => postp(1) % veset(1:1,ieset)
  if( postp(1) % npp_setse(2) /= 0 ) setvr => postp(1) % veset(2:2,ieset)
  if( postp(1) % npp_setse(3) /= 0 ) setki => postp(1) % veset(3:3,ieset)
  if( postp(1) % npp_setse(4) /= 0 ) divu2 => postp(1) % veset(4:4,ieset)
  if( postp(1) % npp_setse(5) /= 0 ) setmo => postp(1) % veset(5:7,ieset)

  if( postp(1) % npp_setse(1) /= 0 ) setve = 0.0_rp  ! Set mean velocity module
  if( postp(1) % npp_setse(2) /= 0 ) setvr = 0.0_rp  ! Set mean vorticity module
  if( postp(1) % npp_setse(3) /= 0 ) setki = 0.0_rp  ! Set kinetic turbulent energy
  if( postp(1) % npp_setse(4) /= 0 ) divu2 = 0.0_rp  ! rho * (div u) * u^2
  if( postp(1) % npp_setse(5) /= 0 ) setmo = 0.0_rp  ! rho * u 
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( leset(ielem) == iesec ) then
        ! 
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather operations
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           !
           ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
           !
           do igaus = 1,pgaus     
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
              gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
              setvo = setvo + gpvol
              gpvel = 0.0_rp
              do inode = 1,pnode
                 gpvel(1:ndime) = gpvel(1:ndime)          &
                      + elmar(pelty) % shape(inode,igaus) &
                      * elvel(1:ndime,inode)
              end do
              if( postp(1) % npp_setse(1) /= 0 ) then
                 !
                 ! Velocity average module: SETVE
                 !
                 do idime=1,ndime
                    setve=setve+gpvol*gpvel(idime)**2
                 end do
              end if
              if( postp(1) % npp_setse(2) /= 0 ) then
                 !
                 ! Vorticity average module: SETVR
                 !
                 dvdx = 0.0_rp
                 dudy = 0.0_rp
                 do inode = 1,pnode
                    dvdx = dvdx + gpcar(1,inode,igaus) * elvel(2,inode)
                    dudy = dudy + gpcar(2,inode,igaus) * elvel(1,inode)
                 end do
                 setvr = setvr + gpvol*(dvdx-dudy)**2
              end if
              if( postp(1) % npp_setse(3) /= 0 ) then
                 !
                 ! Kinetic energy: SETKI
                 !
                 dummr = 0.0_rp
                 do idime = 1,ndime
                    dummr = dummr + gpvel(idime)**2
                 end do
                 setki = setki + 0.5_rp * gpden(igaus) * gpvol * dummr
              end if
              if( postp(1) % npp_setse(4) /= 0 ) then
                 !
                 ! rho * (div u) * u^2
                 !
                 dummr = 0.0_rp
                 divu  = 0.0_rp
                 do inode = 1,pnode
                    divu = divu + dot_product(gpcar(:,inode,igaus),elvel(:,inode))
                 end do
                 do idime = 1,ndime
                    dummr = dummr + gpvel(idime)**2
                 end do
                 divu2 = divu2 + gpden(igaus) * gpvol * dummr * divu
              end if
              if( postp(1) % npp_setse(5) /= 0 ) then
                 !
                 ! Momentum
                 !
                 do idime = 1,ndime
                    setmo(idime) = setmo(idime)                        &
                                  + gpden(igaus) * gpvel(idime) * gpvol
                 end do
              end if
           end do
        end if
     end if

  end do elements

end subroutine nsi_elmset
