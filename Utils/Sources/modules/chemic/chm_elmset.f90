subroutine chm_elmset(iesec,ieset)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_elmset
  ! NAME 
  !    chm_elmset
  ! DESCRIPTION
  !    This routine computes variables on an element set W.
  !    The variable are: 
  !    1. SETVO: set surface           = meas(W)=int_W
  !    2. SETVE: set mean vel. module  = int_W u^2 ]
  !    3. SETVR: set mean vort. module = int_W w^2 ]
  !                                      where w=dv/dx-du/dy
  !    SETVE and SETVR are normailzed further on in chm_outset
  ! USES
  !    chm_elmgat
  !    elmder
  ! USED BY
  !    chm_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip), intent(in)  :: iesec,ieset
  real(rp),    pointer     :: setvo(:)
  real(rp),    pointer     :: setv1(:)
  real(rp),    pointer     :: setv2(:)
  real(rp),    pointer     :: setv3(:)
  real(rp),    pointer     :: setv4(:)
  integer(ip)              :: pnode,pgaus,pelty,nvabi,ispec,iodes
  integer(ip)              :: ielem,igaus,idime,inode,ipoin,iclas
  real(rp)                 :: gpcar(ndime,mnode,mgaus) 
  real(rp)                 :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: elcon(mnode,nspec_chm)
  real(rp)                 :: gpcon(mgaus,nspec_chm)
  real(rp)                 :: dvdx,dudy,gpvol,gpdet,dummr
  !
  ! Initialization
  !
  nvabi =  postp(1) % nvaes+1
  setvo => veset( nvabi:nvabi ,ieset)  
  setvo =  0.0_rp  ! Set volume

  if( postp(1) % npp_setse(1) /= 0 ) setv1 => postp(1) % veset(1:1,ieset)
  if( postp(1) % npp_setse(2) /= 0 ) setv2 => postp(1) % veset(2:2,ieset)
  if( postp(1) % npp_setse(3) /= 0 ) setv3 => postp(1) % veset(3:3,ieset)
  if( postp(1) % npp_setse(4) /= 0 ) setv4 => postp(1) % veset(4:4,ieset)

  if( postp(1) % npp_setse(1) /= 0 ) setv1 = 0.0_rp 
  if( postp(1) % npp_setse(2) /= 0 ) setv2 = 0.0_rp 
  if( postp(1) % npp_setse(3) /= 0 ) setv3 = 0.0_rp 
  if( postp(1) % npp_setse(4) /= 0 ) setv4 = 0.0_rp 
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( leset(ielem) == iesec ) then
        ! 
        ! Element properties and dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        !
        ! Gather operations
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do ispec = 1,nspec_chm
              elcon(inode,ispec) = conce(ipoin,ispec,1)
           end do
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        !
        ! Gauss point values
        !
        do igaus = 1,pgaus
           do ispec = 1,nspec_chm           
              gpcon(igaus,ispec) = 0.0_rp
           end do
        end do
        do ispec = 1,nspec_chm           
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gpcon(igaus,ispec) = gpcon(igaus,ispec) &
                      +elmar(pelty)%shape(inode,igaus)   &
                      *elcon(inode,ispec)
              end do
           end do
        end do
        !
        ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
        !
        do igaus = 1,pgaus     
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
           gpvol = elmar(pelty)%weigp(igaus)*gpdet                 ! |J|*wg
           setvo = setvo + gpvol

           if( postp(1) % npp_setse(1) /= 0 ) then
              !
              ! IPDE1: 1st class
              !
              setv1 = setv1 + gpvol * gpcon(igaus,1)
           end if

           if( postp(1) % npp_setse(2) /= 0 ) then
              !
              ! IODE1: 1st ODE
              !
              setv2 = setv2 + gpvol * gpcon(igaus,nclas_chm+1)
           end if

           if( postp(1) % npp_setse(3) /= 0 ) then
              !
              ! IODEN: remaining ODE's
              !
              dummr = 0.0_rp
              do iodes = 1,nodes_chm
                 ispec = iodes+nclas_chm
                 dummr = dummr + real(ispec) * gpcon(igaus,ispec)
              end do
              setv3 = setv3 + gpvol * dummr
           end if

           if( postp(1) % npp_setse(4) /= 0 ) then
              !
              ! Total mass
              !
              do iclas = 1,nclas_chm
                 setv4 = setv4 + gpvol * gpcon(igaus,iclas)
              end do

           end if

        end do
     end if

  end do elements

end subroutine chm_elmset

