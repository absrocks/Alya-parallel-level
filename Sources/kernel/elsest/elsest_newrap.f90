subroutine elsest_newrap(&
     coglo,coloc,ndime,pnode,elcod,&
     xjacm,xjaci,shapf,deriv)
  !------------------------------------------------------------------------
  !
  !    Calculate the inverse transformation (x,y,z)-->(s,t,r)
  !
  !    Iterate for f(s_i)=x_i: ds = J^{-1}.xpoin - J^{-1}.f(s_i)
  !                               = J^{-1}dx
  !                            ds = s_{i+1}-s_i  (deltas)
  !                            dx = xpoin-f(s_i) (deltax)
  !    where the s_i's are the coordinates in the local basis and
  !    xpoin(idime)'s the real ones.
  !
  !------------------------------------------------------------------------
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: coglo(3),elcod(ndime,pnode)
  real(rp),    intent(out) :: coloc(3)
  integer(ip)              :: inode,iiter,jdime,maxit,idime,ierro
  integer(ip)              :: ntens,jnode
  real(rp)                 :: shapf(pnode),deriv(ndime,pnode)
  real(rp)                 :: xjacm(ndime,ndime),xjaci(ndime,ndime)
  real(rp)                 :: deltx(3),delts(3),xnorm,detja
  real(rp)                 :: rnode,diame,coocg(3),dercg(ndime,pnode)
  real(rp)                 :: hessi(6,pnode),shacg(pnode)
  !
  ! Initial condition
  !
  coloc(1) = 0.0_rp
  coloc(2) = 0.0_rp
  coloc(3) = 0.0_rp
  ntens    = max(1_ip,3*ndime-3)
  call shader(coloc,ndime,pnode,ntens,shapf,deriv,hessi,ierro)
  !
  ! Element diameter
  !
  coocg(1:ndime) = 0.0_rp
  do inode = 1,pnode
     coocg(1:ndime) = coocg(1:ndime) + elcod(1:ndime,inode)
  end do
  rnode          = 1.0_rp / real(pnode)
  coocg(1:ndime) = rnode  * coocg(1:ndime)
  call shader(coocg,ndime,pnode,ntens,shacg,dercg,hessi,ierro)
  call mbmabt(xjacm,elcod,dercg,ndime,ndime,pnode)
  call invmtx(xjacm,xjaci,diame,ndime)
  if( diame <= 0.0_rp ) then
     diame = -1.0_rp
     do inode = 1,pnode
        do jnode = inode+1,pnode
           do idime = 1,ndime
              diame = max(diame,abs(elcod(idime,inode)-elcod(idime,jnode)))
           end do
        end do
     end do
  else
     diame = 1.0_rp / (diame**(2.0_rp/real(ndime)))
  end if
  !
  ! Initialize dx=coglo-f(s_1) with s_1=(0,0,0)
  !  
  do idime = 1,ndime
     deltx(idime) = 0.0_rp
     do inode = 1,pnode
        deltx(idime) = deltx(idime) + shapf(inode) * elcod(idime,inode)
     end do
  end do
  xnorm = 0.0_rp
  do idime = 1,ndime
     deltx(idime) = coglo(idime) - deltx(idime)
     xnorm        = xnorm + deltx(idime)*deltx(idime)
  end do
  xnorm =  xnorm * diame
  iiter =  0
  maxit = 10
  !
  ! Iterate for f(s_i)=x_i
  !
  do while( xnorm > 1e-8_rp .and. iiter <= maxit )
     iiter = iiter + 1
     !
     ! Compute J
     !
     do jdime = 1,ndime
        do idime = 1,ndime
           xjacm(idime,jdime) = 0.0_rp
           do inode = 1,pnode
              xjacm(idime,jdime) = xjacm(idime,jdime) + deriv(idime,inode) * elcod(jdime,inode)
           end do
        end do
     end do
     !
     ! Compute J^{-1}
     !
     call elsest_invmtx(xjacm,xjaci,detja,ndime)

     do idime = 1,ndime
        delts(idime) = 0.0_rp
        !
        ! Compute J^{-1}.dx
        !
        do jdime = 1,ndime
           delts(idime) = delts(idime) + deltx(jdime) * xjaci(jdime,idime)
        end do
     end do
     coloc(1:ndime) = coloc(1:ndime) + delts(1:ndime)
     if ( coloc(1) > 1d99 .or. coloc(2) > 1d99 .or. coloc(3) > 1d99 ) then
        iiter = maxit + 1
     else
        call shader(coloc,ndime,pnode,ntens,shapf,deriv,hessi,ierro)
     end if
     !
     ! Compute f_i
     !
     do idime = 1,ndime
        deltx(idime) = 0.0_rp
        do inode = 1,pnode
           deltx(idime) = deltx(idime) + shapf(inode) * elcod(idime,inode)
        end do
     end do
     !
     ! Compute dx=coglo-f
     !         xnorm=sum ds^2
     !
     xnorm = 0.0_rp
     do idime = 1,ndime
        deltx(idime) = coglo(idime) - deltx(idime)
        xnorm        = xnorm + delts(idime) * delts(idime)
     end do
     xnorm = xnorm * diame
  end do

  if( xnorm > 1e-8_rp ) coloc(1) = 2.0_rp

end subroutine elsest_newrap
