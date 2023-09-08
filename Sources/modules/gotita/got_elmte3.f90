subroutine got_elmte3(&
     pnode,pgaus,gpsha,gpcar,gpvdr,gpdiv,gpst2,tescb) 
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmte3
  ! NAME 
  !    got_elmte3
  ! DESCRIPTION
  !    Compute the test function at the Gauss point
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  mnode,ndime,ntens
  use def_gotita, only     :  kfl_staty_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpvdr(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(in)  :: gpst2(pgaus)
  real(rp),    intent(out) :: tescb(pnode,pgaus)
  integer(ip)              :: idime,inode,igaus
  real(rp)                 :: fact1,fact2

  if(kfl_staty_got==0) then
     !
     ! NO STABILIZATION: TESCB= beta_h
     !
     do igaus=1,pgaus
        do inode=1,pnode
           tescb(inode,igaus)=gpsha(inode,igaus)
        end do
     end do

  else if(kfl_staty_got==1) then
     !
     ! SUPG: TESCB= tau2*(u.grad)beta_h
     !
     do igaus=1,pgaus
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)beta_h
           end do
           tescb(inode,igaus)=gpsha(inode,igaus)+gpst2(igaus)*fact1 
        end do
     end do

  else if(kfl_staty_got>=2) then
     !
     ! ASGS: TESCB= [1-tau2*div(u)]*beta_h + tau2*(u.grad)beta_h
     !
     do igaus=1,pgaus
        fact2=1.0_rp-gpst2(igaus)*gpdiv(igaus)                            ! 1-tau2*div(u)
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)beta_h
           end do
           tescb(inode,igaus)=gpsha(inode,igaus)*fact2+gpst2(igaus)*fact1 
        end do
     end do

  end if

end subroutine got_elmte3
