subroutine got_elmre2(&
     pnode,pgaus,ndofr,gpsha,gpcar,gppor,gpvdr,gpgvd,resim)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmre2
  ! NAME 
  !    got_elmre2
  ! DESCRIPTION
  !    Compute the residual of the momentum and continuity equations
  !    1. For Picard linearization:
  !       1/(theta*dt)*u_h + (u.grad)u_h + sig*u_h
  !    2. For Newton_Raphson linearization:
  !       1/(theta*dt)*u_h + (u.grad)u_h + (uh.grad)u + sig*u_h
  !                     
  ! OUTPUT
  !    RESIM(k,i,j) ... k-component of momentum residual of node i
  !                     Gauss point j
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_gotita, only     :  kfl_timei_got,dtinv_got,ndofn_got,&
       &                      kfl_forme_got,penal_got,kfl_coupl_got,&
       &                      kfl_linea_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,ndofr
  real(rp),    intent(in)  :: gppor(pgaus)
  real(rp),    intent(in)  :: gpvdr(ndime,pgaus)
  real(rp),    intent(in)  :: gpgvd(ndime,ndime,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: resim(ndofr,pnode,pgaus)
  integer(ip)              :: idime,jdime,inode,igaus,idofn
  real(rp)                 :: fact1
  !
  ! Initialization
  !
  if(kfl_linea_got==2) then
     do igaus=1,pgaus
        do inode=1,pnode
           do idofn=1,ndofr
              resim(idofn,inode,igaus)=0.0_rp
           end do
        end do
     end do
  end if
  !
  ! Advection and porosity: (u.grad)u_h + sig*u_h
  !
  do igaus=1,pgaus
     do inode=1,pnode
        fact1=0.0_rp
        do idime=1,ndime
           fact1=fact1+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)
        end do
        resim(1,inode,igaus)=fact1+gppor(igaus)*gpsha(inode,igaus)
     end do
  end do
  !
  ! Time integration: 1/(dt*theta)*u_h
  !
  if(kfl_timei_got==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           resim(1,inode,igaus)=resim(1,inode,igaus)&
                +dtinv_got*gpsha(inode,igaus)
        end do
     end do
  end if
  !
  ! Newton-Raphson: (uh.grad)u
  !
  if(kfl_linea_got==2) then
      do igaus=1,pgaus
        do inode=1,pnode
           do idime=2,ndime
              idofn=(idime-1)*ndime+idime
              resim(idofn,inode,igaus)=resim(1,inode,igaus)
           end do
        end do
     end do
     do igaus=1,pgaus
        do inode=1,pnode
           do idime=1,ndime
              idofn=(idime-1)*ndime
              do jdime=1,ndime
                 idofn=idofn+1
                 resim(idofn,inode,igaus)=resim(idofn,inode,igaus)&
                      +gpsha(inode,igaus)*gpgvd(jdime,idime,igaus)
              end do
           end do
        end do
     end do
  end if

end subroutine got_elmre2
