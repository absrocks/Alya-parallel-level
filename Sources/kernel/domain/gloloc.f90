subroutine gloloc(ndime,ntens,pnode,elcor,xpoin,posgp)
  !------------------------------------------------------------------------
  !****f* domain/gloloc
  ! NAME 
  !    gloloc
  ! DESCRIPTION
  !    Calculate the inverse transformation (x,y,z)-->(s,t,r)
  !    Iterate for f(s_i)=x_i: ds = J^{-1}.xpoin - J^{-1}.f(s_i)
  !                               = J^{-1}dx
  !                            ds = s_{i+1}-s_i  (deltas)
  !                            dx = xpoin-f(s_i) (deltax)
  !    where the s_i's are the elcorinates in the local basis and
  !    xpoin(idime)'s the real ones.    
  ! USES
  ! USED BY
  !    tem_elmbub
  !------------------------------------------------------------------------
  use def_kintyp
  implicit  none
  integer(ip), intent(in)  :: ndime,ntens,pnode
  real(rp),    intent(in)  :: xpoin(ndime),elcor(ndime,pnode)
  real(rp),    intent(out) :: posgp(ndime)
  integer(ip)              :: inode,iiter,jdime,maxit,idime,ierro
  real(rp)                 :: xnorm,detja
  real(rp)                 :: gpsha(pnode),gpder(ndime,pnode)
  real(rp)                 :: gphes(ntens,pnode)
  real(rp)                 :: xjacm(ndime,ndime),xjaci(ndime,ndime)
  real(rp)                 :: deltx(ndime),delts(ndime)
  !
  ! Initialization
  !
  do idime=1,ndime
     posgp(idime)=0.0_rp
  end do
  call shafun(&
       posgp(1:ndime),ndime,pnode,ntens,gpsha,gpder,gphes,ierro)
  do idime=1,ndime
     deltx(idime)=0.0_rp
     do inode=1,pnode
        deltx(idime)=deltx(idime)+gpsha(inode)*elcor(idime,inode)
     end do
  end do
  xnorm=0.0_rp
  !
  ! Initialize dx=xpoin-f(s_1) with s_1=(0,0,0)
  !
  do idime=1,ndime
     deltx(idime)=xpoin(idime)-deltx(idime)
     xnorm=xnorm+deltx(idime)*deltx(idime)
  end do
  iiter=0
  maxit=10
  !
  ! Iterate for f(s_i)=x_i
  !
  do while ( (xnorm>1e-8).and.(iiter<=maxit) )
     iiter=iiter+1
     !
     ! Compute J
     !
     do idime=1,ndime
        do jdime=1,ndime
           xjacm(idime,jdime)=0.0_rp
           do inode=1,pnode
              xjacm(idime,jdime)=xjacm(idime,jdime)&
                   +gpder(idime,inode)*elcor(jdime,inode)
           end do
        end do
     end do
     !
     ! Compute J^{-1}
     !
     call invmtx(xjacm,xjaci,detja,ndime)
     !
     ! Compute J^{-1}.dx
     !
     do idime=1,ndime
        delts(idime)=0.0_rp
        do jdime=1,ndime
           delts(idime)=delts(idime)+deltx(jdime)*xjaci(jdime,idime)
        end do
        posgp(idime)=posgp(idime)+delts(idime)
     end do
     if ((posgp(1)>1d99).or.(posgp(2)>1d99).or.(posgp(ndime)>1d99) ) then
        iiter=maxit+1
     else
        call shafun(&
             posgp(1:ndime),ndime,pnode,ntens,gpsha,gpder,gphes,ierro)
     end if
     !
     ! Compute f_i
     !
     do idime=1,ndime
        deltx(idime)=0.0_rp
        do inode=1,pnode
           deltx(idime)=deltx(idime)&
                +gpsha(inode)*elcor(idime,inode)
        end do
     end do
     xnorm=0.0_rp
     !
     ! Compute    dx=xpoin-f
     !            xnorm=sum dx^2
     do idime=1,ndime
        deltx(idime)=xpoin(idime)-deltx(idime)
        xnorm=xnorm+deltx(idime)*deltx(idime)
     end do
  end do

end subroutine gloloc
