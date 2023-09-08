subroutine got_elmte2(&
     pnode,pgaus,gpsha,gpcar,gppor,gpvdr,gpvel,gpst1,tesmv) 
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmtes
  ! NAME 
  !    got_elmtes
  ! DESCRIPTION
  !    Compute the test function at the Gauss point
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
  real(rp),    intent(in)  :: gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpst1(pgaus)
  real(rp),    intent(out) :: tesmv(pnode,pgaus)
  integer(ip)              :: idime,inode,igaus
  real(rp)                 :: fact1,fact2

  if(kfl_staty_got==0) then
     !
     ! NO STABILIZATION
     !
     do igaus=1,pgaus
        do inode=1,pnode
           tesmv(inode,igaus)=gpsha(inode,igaus)
        end do
     end do

  else if(kfl_staty_got==1) then
     !
     ! SUPG: TESMV= v_h + tau1*(u.grad)v_h
     !
     do igaus=1,pgaus
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)v_h
           end do
           tesmv(inode,igaus)=gpsha(inode,igaus)+gpst1(igaus)*fact1 
        end do
     end do

  else if(kfl_staty_got>=2) then
     !
     ! ASGS: TESMV= (1-tau1*sig)*v_h + tau1*(u.grad)v_h
     !
     do igaus=1,pgaus
        fact2=1.0_rp-gpst1(igaus)*gppor(igaus)                            ! 1-tau1*sig
        do inode=1,pnode
           fact1=0.0_rp
           do idime=1,ndime
              fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)     ! (u.grad)v_h
           end do
           tesmv(inode,igaus)=gpsha(inode,igaus)*fact2+gpst1(igaus)*fact1 
        end do
     end do

  end if

end subroutine got_elmte2
