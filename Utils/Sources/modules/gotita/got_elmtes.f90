subroutine got_elmtes(&
     pnode,pgaus,gpsha,gpcar,gppor,gpvdr,gpugu,gpvel,&
     gpdiv,gpst1,gpst2,tesmv,tesmb,tescv,tescb) 
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmtes
  ! NAME 
  !    got_elmtes
  ! DESCRIPTION
  !    Compute the test function at the Gauss points
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  mnode,ndime,ntens
  use def_gotita, only     :  kfl_staty_got,kfl_coupl_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gppor(pgaus),gpvdr(ndime,pgaus)
  real(rp),    intent(in)  :: gpugu(ndime,pgaus),gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpdiv(pgaus)
  real(rp),    intent(in)  :: gpst1(pgaus),gpst2(pgaus)
  real(rp),    intent(out) :: tesmv(pnode,pgaus)
  real(rp),    intent(out) :: tesmb(ndime,pnode,pgaus)
  real(rp),    intent(out) :: tescv(ndime,pnode,pgaus)
  real(rp),    intent(out) :: tescb(pnode,pgaus)
  integer(ip)              :: idime,inode,igaus
  real(rp)                 :: fact1,fact2,fact3

  if(kfl_staty_got==0) then
     !
     ! NO STABILIZATION
     !
     do igaus=1,pgaus
        do inode=1,pnode
           tesmv(inode,igaus)=gpsha(inode,igaus)
           tescb(inode,igaus)=gpsha(inode,igaus)
        end do
     end do

  else if(kfl_staty_got==1) then
     !
     ! SUPG
     !
     do igaus=1,pgaus
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime                                               ! (u.grad)v_h and
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)beta_h
           end do
           tesmv(inode,igaus)=gpsha(inode,igaus)+gpst1(igaus)*fact1 
           tescb(inode,igaus)=gpsha(inode,igaus)+gpst2(igaus)*fact1 
        end do
     end do

  else if(kfl_staty_got==2.or.kfl_staty_got==3) then
     !
     ! ASGS
     !
     ! TESMV= [1-tau1'*sig]*v_h      + tau1'*(u.grad)v_h
     ! TESCB= [1-tau2*div(u)]*beta_h + tau2*(u.grad)beta_h
     !
     do igaus=1,pgaus
        fact2=1.0_rp-gpst1(igaus)*gppor(igaus)                            ! 1-tau1'*sig
        fact3=1.0_rp-gpst2(igaus)*gpdiv(igaus)                            ! 1-tau2*div(u)
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime                                               ! (u.grad)v_h and
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)beta_h
           end do
           tesmv(inode,igaus)=gpsha(inode,igaus)*fact2+gpst1(igaus)*fact1 
           tescb(inode,igaus)=gpsha(inode,igaus)*fact3+gpst2(igaus)*fact1 
        end do
     end do
     if(kfl_coupl_got==1) then
        if(kfl_staty_got==3) then
           !
           ! TESMB=tau1'*grad(beta_h)
           !
           do igaus=1,pgaus
              do inode=1,pnode
                 do idime=1,ndime
                    tesmb(idime,inode,igaus)=gpst1(igaus)*gpcar(idime,inode,igaus)
                 end do
              end do
           end do
           !
           ! TESCV=tau2*[sig*(ua-u)-(u.grad)u].v_h
           ! TESCV=tau2*sig*ua.v_h
           !
           if(kfl_staty_got==3) then
              do igaus=1,pgaus
                 do inode=1,pnode
                    fact1=gpst2(igaus)*gpsha(inode,igaus)
                    do idime=1,ndime 
                       tescv(idime,inode,igaus)=fact1&
                            *(gppor(igaus)*(gpvel(idime,igaus)-gpvdr(idime,igaus))&
                            -gpugu(idime,igaus))
                       tescv(idime,inode,igaus)=fact1*gppor(igaus)*gpvel(idime,igaus)
                       tescv(idime,inode,igaus)=0.0_rp
                    end do
                 end do
              end do
           end if
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 do idime=1,ndime
                    tesmb(idime,inode,igaus)=0.0_rp
                    tescv(idime,inode,igaus)=0.0_rp
                 end do
              end do
           end do
        end if
     end if

  end if

end subroutine got_elmtes
