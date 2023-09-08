subroutine got_elmres(&
     pnode,pgaus,gpsha,gpcar,gppor,gpcdr,gpvdr,gpgcd,&
     gpdiv,gpvel,gpugu,elcdr,elvdr,resim,resic)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmres
  ! NAME 
  !    got_elmres
  ! DESCRIPTION
  !    Compute the residual of the momentum and continuity equations
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
  use def_gotita, only     :  dtinv_got,ndofn_got,kfl_forme_got,&
       &                      penal_got,kfl_coupl_got,kfl_timem_got,&
       &                      kfl_timec_got,almax_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gppor(pgaus),gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpcdr(pgaus),gpgcd(ndime,pgaus)
  real(rp),    intent(in)  :: gpvdr(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(in)  :: gpugu(ndime,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: elcdr(pnode),elvdr(ndime,pnode)
  real(rp),    intent(out) :: resim(ndofn_got(3),pnode,pgaus)
  real(rp),    intent(out) :: resic(ndofn_got(3),pnode,pgaus)
  integer(ip)              :: idime,inode,igaus
  real(rp)                 :: fact1,fact2
  !
  ! Initialization
  !
  do igaus=1,pgaus                
     do inode=1,pnode
        do idime=1,ndofn_got(3)
           resim(idime,inode,igaus)=0.0_rp
           resic(idime,inode,igaus)=0.0_rp
        end do
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! MOMENTUM EQUATIONS: alpha/(theta*dt)*u_h + (alpha*u.grad)u_h
  !                     + sig*alpha*u_h + (sig*(u-ua)+(u.grad)u)*alpha_h
  !
  !----------------------------------------------------------------------
  !
  ! Time integration
  !
  if(kfl_timem_got==1) then
     if(kfl_coupl_got==1) then
        !
        ! Coupled eqns: (alpha+penal)/(dt*theta)*u_h
        !
        do igaus=1,pgaus
           fact1=dtinv_got*(gpcdr(igaus)+penal_got)
           do inode=1,pnode
              resim(1,inode,igaus)=resim(1,inode,igaus)&
                   +fact1*gpsha(inode,igaus)
           end do
        end do
     else
        !
        ! Uncoupled eqns: alpha_max/(dt*theta)*u_h
        !
        do igaus=1,pgaus
           fact1=dtinv_got*almax_got
           do inode=1,pnode
              resim(1,inode,igaus)=resim(1,inode,igaus)&
                   +fact1*gpsha(inode,igaus)
           end do
        end do
     end if
  end if
  !
  ! Advection and porosity 
  !     
  if(kfl_coupl_got==1) then
     !
     ! Coupled eqns: (alpha+eps)*[ (u.grad)u_h + sig*u_h ]
     !     
     do igaus=1,pgaus
        fact1=gpcdr(igaus)+penal_got
        do inode=1,pnode
           fact2=0.0_rp
           do idime=1,ndime
              fact2=fact2+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)
           end do
           resim(1,inode,igaus)=resim(1,inode,igaus)&
                +fact1*(fact2+gppor(igaus)*gpsha(inode,igaus))
        end do
     end do
  else 
     !
     ! Uncoupled eqns: alpha_max*[ (u.grad)u_h + sig*u_h ]
     !
     do igaus=1,pgaus
        fact1=almax_got
        do inode=1,pnode
           fact2=0.0_rp
           do idime=1,ndime
              fact2=fact2+gpvdr(idime,igaus)*gpcar(idime,inode,igaus)
           end do
           resim(1,inode,igaus)=resim(1,inode,igaus)&
                +fact1*(fact2+gppor(igaus)*gpsha(inode,igaus))
        end do
     end do
  end if
  !
  ! Pressure-like and convection terms: [sig*(u-ua)+(u.grad)u]*alpha_h
  !     
  if(kfl_coupl_got==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           fact1=gppor(igaus)*gpsha(inode,igaus)
           do idime=1,ndime
              !resim(1+idime,inode,igaus)=gpsha(inode,igaus)&
              !     *(gppor(igaus)*(gpvdr(idime,igaus)-gpvel(idime,igaus))&
              !     +gpugu(idime,igaus))
              !resim(1+idime,inode,igaus)=-fact1*gpvel(idime,igaus)
           end do
        end do
     end do
  end if
  !
  !----------------------------------------------------------------------
  !
  ! CONTINUITY EQUATION: 1/(theta*dt)*alpha_h + grad(alpha_h).u
  !                      + alpha_h*div(u) + grad(alpha).u_h
  !                      + alpha*div(u_h)
  ! 
  !----------------------------------------------------------------------
  !
  ! Time integration: alpha_h/(dt*theta)
  !
  if(kfl_timec_got==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           resic(ndofn_got(3),inode,igaus)=resic(ndofn_got(3),inode,igaus)&
                +dtinv_got*gpsha(inode,igaus)
        end do
     end do
  end if
  !
  ! Divergence: alpha_h*div(u) + grad(alpha_h).u
  !
  do igaus=1,pgaus
     do inode=1,pnode
        resic(ndofn_got(3),inode,igaus)=resic(ndofn_got(3),inode,igaus)&  
             +gpdiv(igaus)*gpsha(inode,igaus)
        do idime=1,ndime
           resic(ndofn_got(3),inode,igaus)=resic(ndofn_got(3),inode,igaus)&
                +gpvdr(idime,igaus)*gpcar(idime,inode,igaus)
        end do
     end do
  end do
  !
  ! Divergence (only coupled eqns): alpha*div(u_h) + grad(alpha).u_h
  !
  if(kfl_coupl_got==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           do idime=1,ndime           
              !resic(idime,inode,igaus)=resic(idime,inode,igaus)&
              !     +gpcar(idime,inode,igaus)*gpcdr(igaus)&
              !     +gpgcd(idime,igaus)*gpsha(inode,igaus)
           end do
        end do
     end do
  end if
  !
  ! Conservative form: u_h*(continuity residual)
  !
  if(kfl_forme_got==2.and.kfl_coupl_got==1) then
     do igaus=1,pgaus
        fact1=0.0_rp
        do inode=1,pnode
           fact1=fact1+resic(ndofn_got(3),inode,igaus)*elcdr(inode)
           do idime=1,ndime
              fact1=fact1+resic(idime,inode,igaus)*elvdr(idime,inode)
           end do
        end do
        do inode=1,pnode
           resim(1,inode,igaus)=&
                resim(1,inode,igaus)+gpsha(inode,igaus)*fact1   
        end do
     end do
  end if

end subroutine got_elmres
