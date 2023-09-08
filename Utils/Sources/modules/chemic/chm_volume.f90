subroutine chm_volume(energ_chm)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_volume
  ! NAME 
  !    chm_volume
  ! DESCRIPTION
  !    This routine computes the integral of the concentration
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
  real(rp),    intent(out) :: energ_chm(nspec_chm+1)
  integer(ip)              :: pnode,pgaus,pelty,ispec
  integer(ip)              :: ielem,igaus,idime,inode,ipoin
  real(rp)                 :: xjacm(ndime,ndime) 
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: elcon(mnode,nspec_chm)
  real(rp)                 :: gpcon(mgaus,nspec_chm)
  real(rp)                 :: gpvol,gpdet,xfact

  if( INOTMASTER ) then

     do ispec = 1,nspec_chm+1
        energ_chm(ispec) = 0.0_rp
     end do
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
        do ispec = 1,nspec_chm           
           do igaus = 1,pgaus
              gpcon(igaus,ispec) = 0.0_rp
           end do
        end do
        do ispec = 1,nspec_chm           
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gpcon(igaus,ispec) = gpcon(igaus,ispec)  &
                      + elmar(pelty) % shape(inode,igaus) &
                      * elcon(inode,ispec)
              end do
           end do
        end do
        !
        ! dV:=GPVOL=|J|*wg
        !
        do igaus = 1,pgaus    
           call jacdet(&
                ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                xjacm,gpdet)
           gpvol = elmar(pelty)%weigp(igaus)*gpdet 

           do ispec = 1,nspec_chm
              xfact                  =  gpvol * gpcon(igaus,ispec)
              energ_chm(ispec)       =  energ_chm(ispec)       + xfact
              energ_chm(nspec_chm+1) =  energ_chm(nspec_chm+1) + xfact
           end do

        end do

     end do elements

  end if

  call pararr('SUM',0_ip,nspec_chm+1,energ_chm)

end subroutine chm_volume
