subroutine adrtes(&
     pnode,pgaus,plapl,kfl_advec,kfl_grdif,kfl_react,&
     kfl_taust,gpsha,gpcar,gplap,gprea,gpadv,gpgrd,gpdif,&
     gpcod,gpstp,gpstt,gptes)
  !-----------------------------------------------------------------------
  !****f* mathru/adrtes
  ! NAME
  !   adrtes
  ! DESCRIPTION
  !    Compute the residual GPRES of an ADR equation:
  !    Galerkin:          tau^{-1}*tau'*Ni
  !    Advection:         rho*a.grad(Ni)
  !    Diffusion:         grad(k).grad(Ni)-tau'*k*Lapl(Ni)
  !    Reaction           -tau'*s*Ni
  !    Axisymmetric flow: tau1'*k/r*dNi/dr
  ! OUTPUT 
  !    GPRES(PNODE,PGAUS) ... Residual(Ni) at Gauss points
  !    Compute the adjoint operator
  !    -[rho*a+grad(k)].grad(v) -k*Lapl(v) -k/r*dv/dr+ r*v
  ! OUTPUT 
  !    GPADG ... Adjoint operator at Gauss point
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,kfl_naxis
  implicit none
  integer(ip), intent(in)  :: kfl_advec,kfl_grdif,kfl_react,kfl_taust
  integer(ip), intent(in)  :: pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gplap(pnode,pgaus),gpgrd(ndime,pgaus)
  real(rp),    intent(in)  :: gpdif(pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(in)  :: gprea(pgaus)
  real(rp),    intent(in)  :: gpadv(pnode,pgaus)
  real(rp),    intent(in)  :: gpstp(pgaus),gpstt(pgaus)
  real(rp),    intent(out) :: gptes(pnode,pgaus)
  integer(ip)              :: inode,idime,igaus
  real(rp)                 :: fact1

  if(kfl_taust==0) then
     !
     ! Galerkin
     !
     do igaus=1,pgaus
        do inode=1,pnode
           gptes(inode,igaus)=gpsha(inode,igaus)
        end do
     end do
     
  else
     !
     ! Source term: -s*Ni
     !
     if(kfl_react==0) then
        do igaus=1,pgaus
           do inode=1,pnode
              gptes(inode,igaus)=0.0_rp
           end do
        end do
     else
        do igaus=1,pgaus
           do inode=1,pnode
              gptes(inode,igaus)=-gpsha(inode,igaus)*gprea(igaus)
           end do
        end do
     end if
     !
     ! Advection: rho*a.grad(Ni)
     !
     if(kfl_advec/=0) then
        do igaus=1,pgaus
           do inode=1,pnode
              gptes(inode,igaus)=gptes(inode,igaus)&
                   +gpadv(inode,igaus)
           end do
        end do
     end if
     !
     ! Diffusion: grad(k).grad(Ni)
     !
     if(kfl_grdif==1) then
        do igaus=1,pgaus
           do inode=1,pnode
              do idime=1,ndime
                 gptes(inode,igaus)=gptes(inode,igaus)&
                      +gpgrd(idime,igaus)*gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end if
     !
     ! Diffusion: k*Lapl(Ni)
     !
     if(plapl==1) then     
        do igaus=1,pgaus
           do inode=1,pnode
              gptes(inode,igaus)=gptes(inode,igaus)&
                   +gplap(inode,igaus)
           end do
        end do
     end if
     !
     ! Axisymmetric: k/r*dNi/dr
     !
     if(kfl_naxis==1) then     
        do igaus=1,pgaus
           fact1=gpdif(igaus)/gpcod(1,igaus)
           do inode=1,pnode
              gptes(inode,igaus)=gptes(inode,igaus)&
                   +fact1*gpcar(1,inode,igaus)
           end do
        end do
     end if
     !
     ! Multiply by tau' and add Galerkin term
     !
     do igaus=1,pgaus
        do inode=1,pnode
           gptes(inode,igaus)=&
                gpstt(igaus)*gpsha(inode,igaus)&
                +gpstp(igaus)*gptes(inode,igaus)
        end do
     end do
     
  end if

end subroutine adrtes
