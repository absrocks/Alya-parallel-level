subroutine adrres(&
     pnode,pgaus,plapl,kfl_timei,kfl_advec,kfl_diffu,&
     kfl_react,kfl_grdif,kfl_tisch,kfl_taust,dtinv,pabdf,&
     gpden,gpgrd,gprea,gpadv,gpdif,gpcod,gpsha,gpcar,&
     gplap,gpres)
  !-----------------------------------------------------------------------
  !****f* mathru/adrres
  ! NAME
  !    addres
  ! DESCRIPTION
  !    Compute the residual GPRES of an ADR equation:
  !    Time:              rho/(dt*theta)*Ni
  !    Advection:         rho*a.grad(Ni)
  !    Diffusion:         -k*lapl(Ni)-grad(Ni).grad(k)
  !    Reaction           s*Ni
  !    Axisymmetric flow: -k/r*dNi/dr
  ! OUTPUT 
  !    GPRES(PNODE,PGAUS) ... Residual(Ni) at Gauss points
  ! USES
  ! USED BY
  !    *_elmope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  mnode,ndime,kfl_naxis
  implicit none 
  integer(ip), intent(in)  :: pnode,pgaus,plapl
  integer(ip), intent(in)  :: kfl_timei,kfl_advec,kfl_diffu,kfl_react
  integer(ip), intent(in)  :: kfl_grdif,kfl_tisch,kfl_taust
  real(rp),    intent(in)  :: dtinv,pabdf(*)
  real(rp),    intent(in)  :: gpden(pgaus),gpgrd(ndime,pgaus)
  real(rp),    intent(in)  :: gprea(pgaus),gpadv(pnode,pgaus)
  real(rp),    intent(in)  :: gpdif(pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gplap(pnode,pgaus)
  real(rp),    intent(out) :: gpres(pnode,pgaus)
  integer(ip)              :: igaus,inode,idime
  real(rp)                 :: fact1,factt
  !
  ! Reaction: s*Ni
  !
  if(kfl_react==1) then
     do igaus=1,pgaus                
        do inode=1,pnode
           gpres(inode,igaus)=gprea(igaus)*gpsha(inode,igaus)
        end do
     end do
  else
     do igaus=1,pgaus                
        do inode=1,pnode
           gpres(inode,igaus)=0.0_rp
        end do
     end do
  end if
  !
  ! Advection: rho*a.grad(Ni)
  !
  if(kfl_advec==1) then
     do igaus=1,pgaus                
        do inode=1,pnode
           gpres(inode,igaus)=gpres(inode,igaus)&
                +gpden(igaus)*gpadv(inode,igaus)
        end do
     end do
  end if
  !
  ! Time integration: rho/(dt*theta)*Ni
  !
  if(kfl_timei==1) then
     if(kfl_tisch==1) then
        factt=dtinv                 ! Trapezoidal rule
     else
        factt=dtinv*pabdf(1)        ! BDF scheme
     end if
     do igaus=1,pgaus
        fact1=factt*gpden(igaus)
        do inode=1,pnode
           gpres(inode,igaus)=gpres(inode,igaus)&
                +fact1*gpsha(inode,igaus)
        end do
     end do
  end if
  !
  ! Axisymmetric flow: -k/r*dNi/dr
  !
  if(kfl_naxis==1) then
     do igaus=1,pgaus
        fact1=-gpdif(igaus)/gpcod(1,igaus)
        do inode=1,pnode
           gpres(inode,igaus)=gpres(inode,igaus)&
                +fact1*gpcar(1,inode,igaus)
        end do
     end do     
  end if
  !
  ! Diffusion term: Only stabilized formulation
  !
  if(kfl_taust/=0) then

     if(plapl==1) then
        !
        ! Laplacian: -k*lapl(Ni) 
        !
        do igaus=1,pgaus
           do inode=1,pnode
              gpres(inode,igaus)=gpres(inode,igaus)-gplap(inode,igaus)
           end do
        end do
     end if

     if(kfl_grdif==1) then
        !
        ! Conductivity gradient: -grad(Ni).grad(k)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              do idime=1,ndime
                 gpres(inode,igaus)=gpres(inode,igaus)&
                      -gpgrd(idime,igaus)*gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end if

  end if

end subroutine adrres
