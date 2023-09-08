subroutine chm_sizedi()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_sizedi
  ! NAME 
  !    chm_sizedi
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
  integer(ip) :: pnode,pgaus,pelty,ispec
  integer(ip) :: ielem,igaus,idime,inode,ipoin
  real(rp)    :: gpcar(ndime,mnode,mgaus) 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elcon(mnode,nspec_chm)
  real(rp)    :: gpcon(mgaus,nspec_chm)
  real(rp)    :: rsize(nspec_chm)
  real(rp)    :: gpvol,gpdet,xjacm(9),xjaci(9)

  if( kfl_sized_chm == 1 ) then
     !
     ! Initialization
     !
     do ispec = 1,nspec_chm
        rsize(ispec) = 0.0_rp
     end do

     if( INOTMASTER ) then
        !
        ! Loop over elements
        !
        elements: do ielem = 1,nelem
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

              do ispec = 1,nspec_chm
                 rsize(ispec) = rsize(ispec) + gpvol * gpcon(igaus,ispec)       
              end do

           end do

        end do elements

     end if

     call pararr('SUM',0_ip,nspec_chm,rsize)

     if( INOTSLAVE ) then
        do ispec = 1,nspec_chm
           write(lun_sized_chm,*) ispec,rsize(ispec)
        end do
     end if

  end if

end subroutine chm_sizedi

