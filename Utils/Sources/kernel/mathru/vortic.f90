subroutine vortic(itask)
  !------------------------------------------------------------------------
  !***** mathru/vortic
  ! NAME 
  !    vortic
  ! DESCRIPTION
  !    This routine computes the vorticity and the invariants of the velocity 
  !    gradient.
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master, only    : veloc,vorti,INOTMASTER,kfl_paral
  use def_domain, only    : ndime,npoin,nelem,nnode,mnode,ntens,&
       &                    lnods,ltype,coord,vmasc,elmar,memor_dom
  use mod_memory, only    : memory_alloca
  use mod_memory, only    : memory_deallo
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,inode,ielem,jdime,jnode,jtask
  integer(ip)             :: pnode,pelty,ndofn
  integer(4)              :: istat,i,j
  real(rp)                :: detjm,gpvol,cartc(ndime,mnode),xauxi
  real(rp)                :: gvelo(ndime,ndime)
  real(rp)                :: elvel(ndime,mnode),elcod(ndime,mnode)
  real(rp)                :: xjaci(9),xjacm(9),dummr
  real(rp)                :: Q_crit,lambda2

  jtask = abs(itask)

  if( ndime == 1 ) return

  if( INOTMASTER ) then
     !
     ! Initialization
     !
     if( jtask == 1 ) then
        ndofn = ndime + 1
     else if( jtask == 2 ) then
        if( ndime == 3 ) then
           ndofn = ndime
        else
           ndofn = 1
        end if
     else if( jtask == 3 ) then
        call memory_deallo(memor_dom,'VORTI','vortic',vorti)
        return
     end if
     if( itask  < 0 .and. .not. associated(vorti) ) then
        call memory_alloca(memor_dom,'VORTI','vortic',vorti,ndofn,npoin)
     end if
     vorti = 0.0_rp
     !
     ! Loop over elements
     !
     elements: do ielem = 1,nelem
        pelty = ltype(ielem) 
        pnode = nnode(pelty)
        !
        ! Gather vectors
        !
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
        end do
        !
        ! Loop over Gauss points (which are nodes)
        ! 
        gauss_points: do inode=1,pnode
           ipoin=lnods(inode,ielem)
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod,cartc,detjm,xjacm,xjaci)
           gpvol=elmar(pelty)%weigc(inode)*detjm
           !
           ! Velocity gradients
           !
           do idime = 1,ndime
              do jdime = 1,ndime    
                 gvelo(idime,jdime) = 0.0_rp
                 do jnode = 1,pnode
                    gvelo(idime,jdime) = gvelo(idime,jdime) &
                         + cartc(idime,jnode)*elvel(jdime,jnode)
                 end do
              end do
           end do
           !
           ! Vorticity and 2nd. gradient velocity invariant
           !
           if( jtask == 1 ) then

              xauxi = gvelo(1,1)*gvelo(1,1) + gvelo(2,2)*gvelo(2,2) + 2.0_rp * gvelo(1,2) * gvelo(2,1)

              if (ndime == 3) then

                 vorti(1,ipoin) = vorti(1,ipoin) + gpvol*(gvelo(2,3) - gvelo(3,2))
                 vorti(2,ipoin) = vorti(2,ipoin) + gpvol*(gvelo(3,1) - gvelo(1,3))
                 vorti(3,ipoin) = vorti(3,ipoin) + gpvol*(gvelo(1,2) - gvelo(2,1))

                 xauxi = xauxi + gvelo(3,3) * gvelo(3,3) &
                      + 2.0_rp * gvelo(1,3) * gvelo(3,1) &
                      + 2.0_rp * gvelo(2,3) * gvelo(3,2)              

              else if (ndime == 2) then
                 vorti(1,ipoin) = vorti(1,ipoin) + gpvol*(gvelo(1,2) - gvelo(2,1))
              end if
              ! mariano way
              !vorti(ndime+1,ipoin) = vorti(ndime+1,ipoin) - gpvol * 0.5_rp * xauxi
              ! my way

              call veloc_grad_tensor(gvelo,gpvol,Q_crit,lambda2)

              !vorti(ndime+1,ipoin) = vorti(ndime+1,ipoin) + Q_crit
              vorti(ndime+1,ipoin) = vorti(ndime+1,ipoin) + lambda2
           else

              if (ndime == 3) then

                 vorti(1,ipoin) = vorti(1,ipoin) + gpvol*(gvelo(2,3) - gvelo(3,2))
                 vorti(2,ipoin) = vorti(2,ipoin) + gpvol*(gvelo(3,1) - gvelo(1,3))
                 vorti(3,ipoin) = vorti(3,ipoin) + gpvol*(gvelo(1,2) - gvelo(2,1))

              else if (ndime == 2) then

                 vorti(1,ipoin) = vorti(1,ipoin) + gpvol*(gvelo(1,2) - gvelo(2,1))

              end if

           end if
        end do gauss_points
     end do elements
     !
     ! Solve diagonal system
     !
     call rhsmod(ndofn,vorti) 
     do ipoin=1,npoin
        dummr = 1.0_rp/vmasc(ipoin)
        do idime= 1,ndofn
           vorti(idime,ipoin)=vorti(idime,ipoin)*dummr
        end do
     end do

  end if

end subroutine vortic
