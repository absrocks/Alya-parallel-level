subroutine chm_sizusr()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_sizusr
  ! NAME 
  !    chm_sizusr
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
  integer(ip) :: pnode,pgaus,pelty,ispec,ispev
  integer(ip) :: ielem,igaus,idime,inode,ipoin
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elcon(mnode,nspec_chm)
  real(rp)    :: gpcon(mgaus,nspec_chm)
  real(rp)    :: rsiz1(nspec_chm),rsiz2(nspec_chm)
  real(rp)    :: gpvol(mgaus),gpdet,xjacm(9)
  !
  ! Initialization
  !
  if( wprob_chm == 'COIN1' ) then
     do ispec = 1,nspec_chm
        rsiz1(ispec) = 0.0_rp
     end do
  else if( wprob_chm == 'INVN1' ) then
     do ispec = 1,nspec_chm
        rsiz1(ispec) = 0.0_rp
        rsiz2(ispec) = 0.0_rp
     end do
  else 
     return
  end if

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
           call jacdet(&
                ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                xjacm,gpdet)
           gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpdet 
        end do

        if( wprob_chm == 'COIN1' ) then
           do igaus = 1,pgaus
              do ispec = 2,nspec_chm
                 rsiz1(ispec) = rsiz1(ispec) + gpvol(igaus) * gpcon(igaus,ispec)
              end do
           end do
        else if( wprob_chm == 'INVN1' ) then
           do igaus = 1,pgaus
              ispev = 2 + nodes_chm/2
              do ispec = 3,2 + nodes_chm/2
                 ispev = ispev + 1
                 rsiz1(ispec) = rsiz1(ispec) + gpvol(igaus) * gpcon(igaus,ispec)
                 rsiz2(ispec) = rsiz2(ispec) + gpvol(igaus) * gpcon(igaus,ispev)
              end do
           end do
        end if


     end do elements

  end if
  
  if( wprob_chm == 'COIN1' ) then
     call pararr('SUM',0_ip,nspec_chm,rsiz1)
     if( INOTSLAVE ) then
        do ispec = 2,nspec_chm
           write(lun_resu2_chm,*) ispec,rsiz1(ispec)
        end do
     end if

  else if( wprob_chm == 'INVN1' ) then
     call pararr('SUM',0_ip,nspec_chm,rsiz1)
     call pararr('SUM',0_ip,nspec_chm,rsiz2)
     if( INOTSLAVE ) then
        do ispec = 3,2 + nodes_chm/2
           write(lun_resu2_chm,*) ispec-1,rsiz1(ispec)
           write(lun_resu3_chm,*) ispec-1,rsiz2(ispec)
        end do
     end if

  end if

end subroutine chm_sizusr

