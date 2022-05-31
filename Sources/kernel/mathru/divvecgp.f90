subroutine divvecgp(unkno,diunk)
  !------------------------------------------------------------------------
  !****f* Mathru/divvec
  ! NAME
  !    divvec
  ! DESCRIPTION
  !    This routine computes the divergence of a vector
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
    use def_kintyp
    use def_master, only             : INOTMASTER
    use def_domain, only             : ndime,npoin,nelem,nnode,mnode,&
         &                             lnods,ltype,coord,vmasc,elmar,&
         &                             ngaus,mgaus
    implicit none
    real(rp),   intent(in)          :: unkno(ndime,npoin)
    real(rp),   intent(out), target :: diunk(npoin)
    integer(ip)                     :: ipoin,idime,inode,ielem,igaus
    integer(ip)                     :: pnode,pelty,jnode,pgaus
    real(rp)                        :: detjm,gpvol,gpcar(ndime,mnode) 
    real(rp)                        :: eldiv
    real(rp)                        :: elunk(ndime,mnode),elcod(ndime,mnode)
    real(rp)                        :: xjaci(9),xjacm(9)
    real(rp)                        :: tragl(9),hleng(3),plapl
    real(rp)                        :: gp_divunk(mgaus)
    real(rp)                        :: xfact 

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       diunk     = 0.0_rp

       !
       ! Loop over elements
       !
       do ielem=1,nelem
          pelty = ltype(ielem) 
          pnode = nnode(pelty)
          pgaus = ngaus(pelty)

          !
          ! Gather vectors
          !
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime=1,ndime
                elcod(idime,inode) = coord(idime,ipoin)
                elunk(idime,inode) = unkno(idime,ipoin)
             end do
          end do

          !
          ! Loop over Gauss points
          !
          gp_divunk = 0.0_rp
          do igaus = 1,pgaus
             call elmder(&
                  pnode,ndime,elmar(pelty)%deric(1,1,igaus),&
                  elcod,gpcar,detjm,xjacm,xjaci)

             gpvol = elmar(pelty)%weigp(igaus)*detjm

             do inode=1,pnode
                do idime=1,ndime
                   gp_divunk(igaus) = gp_divunk(igaus) + gpcar(idime,inode) * elunk(idime,inode)
                end do
             end do

             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                xfact = elmar(pelty)%shape(inode,igaus) * gpvol

                diunk(ipoin) = diunk(ipoin) + xfact * gp_divunk(igaus)
             end do

          end do
       end do
       !
       ! Periodicity
       !
       call rhsmod(1_ip,diunk)
       !
       ! Solve diagonal system
       !
       do ipoin=1,npoin
          diunk(ipoin) = diunk(ipoin) / vmasc(ipoin)
       end do       
    end if

end subroutine divvecgp


