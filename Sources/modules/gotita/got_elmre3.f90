subroutine got_elmre3(&
     pnode,pgaus,gpsha,gpcar,gpvdr,gpdiv,resic)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmre3
  ! NAME 
  !    got_elmre3
  ! DESCRIPTION
  !    Compute the residual of the momentum and continuity equations
  ! OUTPUT
  !    RESIM(k,i,j) ... k-component of momentum residual of node i
  !                     Gauss point j
  ! CONTINUITY EQUATION: 1/(theta*dt)*alpha_h + grad(alpha_h).u
  !                      + alpha_h*div(u) + grad(alpha).u_h
  !                      + alpha*div(u_h)
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_gotita, only     :  kfl_timei_got,dtinv_got,ndofn_got,&
       &                      kfl_forme_got,penal_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpvdr(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: resic(pnode,pgaus)
  integer(ip)              :: idime,inode,igaus
  real(rp)                 :: fact1,fact2
  !
  ! Divergence: alpha_h*div(u) + u.grad(alpha_h)
  !
  do igaus=1,pgaus
     do inode=1,pnode
        resic(inode,igaus)=gpdiv(igaus)*gpsha(inode,igaus)
        do idime=1,ndime
           resic(inode,igaus)=resic(inode,igaus)&
                +gpvdr(idime,igaus)*gpcar(idime,inode,igaus)
        end do
     end do
  end do
  !
  ! Time integration: alpha_h/(dt*theta)
  !
  if(kfl_timei_got==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           resic(inode,igaus)=resic(inode,igaus)&
                +dtinv_got*gpsha(inode,igaus)
        end do
     end do
  end if

end subroutine got_elmre3
