subroutine adrmat(&
     pnode,pgaus,plapl,kfl_grdif,kfl_diffu,&
     kfl_taust,bemol,gpdif,gpgrd,gprhs,gpvol,gpsha,&
     gpcar,gplap,gpcod,gpres,gptes,gpadv,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* mathru/adrmat
  ! NAME
  !   elmmat
  ! DESCRIPTION
  !    Compute elemental matrix and RHS
  ! OUTPUT
  !    ELMAT ... LHS matrix
  !    ELRHS ... RHS vector
  ! USES
  ! USED BY
  !    elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,ndime,kfl_naxis
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,plapl
  integer(ip), intent(in)    :: kfl_grdif,kfl_diffu,kfl_taust
  real(rp),    intent(in)    :: bemol
  real(rp),    intent(in)    :: gpdif(pgaus),gpgrd(ndime,pgaus) 
  real(rp),    intent(in)    :: gprhs(pgaus),gpvol(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gplap(pnode,pgaus)
  real(rp),    intent(in)    :: gpcod(ndime,pgaus)
  real(rp),    intent(in)    :: gptes(pnode,pgaus),gpres(pnode,pgaus)
  real(rp),    intent(in)    :: gpadv(pnode,pgaus)
  real(rp),    intent(out)   :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,jnode,kdime,igaus
  real(rp)                   :: fact1,fact2,fact3
  !
  ! Initialization
  !
  do inode=1,pnode 
     elrhs(inode)=0.0_rp
     do jnode=1,pnode
        elmat(inode,jnode)=0.0_rp
     end do
  end do

  if(kfl_taust==0) then
     !
     ! Galerkin
     !
     do igaus=1,pgaus
        fact1=gprhs(igaus)*gpvol(igaus)
        do inode=1,pnode
           elrhs(inode)=elrhs(inode)+fact1*gpsha(inode,igaus)
           fact2=gpvol(igaus)*gpsha(inode,igaus)
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode)&
                   +fact2*gpres(jnode,igaus)
           end do
        end do
     end do
     do igaus=1,pgaus
        fact1=gpdif(igaus)*gpvol(igaus)
        do inode=1,pnode
           do jnode=1,pnode
              do kdime=1,ndime
                 elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                      *gpcar(kdime,inode,igaus)&
                      *gpcar(kdime,jnode,igaus)
              end do
           end do
        end do
     end do

  else
     !
     ! ASGS stabilization: GPTES(Ni)=(tau^-1*tau')*Nj-tau'*L*(Nj)
     !
     do igaus=1,pgaus
        fact1=gprhs(igaus)*gpvol(igaus)
        do inode=1,pnode
           elrhs(inode)=elrhs(inode)+fact1*gptes(inode,igaus)
           fact2=gpvol(igaus)*gptes(inode,igaus)
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode)&
                   +fact2*gpres(jnode,igaus)
           end do
        end do
     end do
  end if

  if(kfl_diffu==1) then

     do igaus=1,pgaus
        !
        ! Diffusion term:  ( k*grad(Nj) , grad(Ni) )
        !
        fact1=gpdif(igaus)*gpvol(igaus)
        do inode=1,pnode
           do jnode=1,pnode
              do kdime=1,ndime
                 elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                      *gpcar(kdime,inode,igaus)&
                      *gpcar(kdime,jnode,igaus)
              end do
           end do
        end do
     end do

     if(kfl_grdif==1) then
        !
        ! Diffusion term: ( grad(k).grad(Nj) , Ni )
        !
        do igaus=1,pgaus
           do inode=1,pnode
              fact1=gpvol(igaus)*gpsha(inode,igaus)
              do jnode=1,pnode
                 do kdime=1,ndime
                    elmat(inode,jnode)=elmat(inode,jnode)&
                         +fact1*gpcar(kdime,jnode,igaus)&
                         *gpgrd(kdime,igaus)
                 end do
              end do
           end do
        end do
     end if

     if(plapl==1) then
        !
        ! Diffusion term: ( k lapl(Nj) , Ni )
        !
        do igaus=1,pgaus
           do inode=1,pnode
              fact1=gpvol(igaus)*gpsha(inode,igaus)
              do jnode=1,pnode
                 elmat(inode,jnode)=elmat(inode,jnode)&
                      +fact1*gplap(jnode,igaus)
              end do
           end do
        end do
     end if

     if(kfl_naxis==1) then
        !
        ! Axisymmetric term: ( k/r*grad(Nj) , Ni )
        !
        do igaus=1,pgaus
           fact1=gpvol(igaus)*gpdif(igaus)/gpcod(1,igaus)
           do inode=1,pnode
              fact2=fact1*gpsha(inode,igaus)
              do jnode=1,pnode
                 elmat(inode,jnode)=elmat(inode,jnode)&
                      +fact2*gpcar(1,jnode,igaus)
              end do
           end do
        end do

     end if

  end if
  !
  ! Convective term: - b*( Nj , rho*a.grad(Ni) ) 
  !                  - b*( rho*a.grad(Nj) , Ni )
  !
  if(bemol>0.0_rp) then
     do igaus=1,pgaus
        fact1=-bemol*gpvol(igaus)
        do inode=1,pnode
           fact2=fact1*gpsha(inode,igaus)
           do jnode=1,pnode
              fact3=gpadv(jnode,igaus)*fact2
              elmat(inode,jnode)=elmat(inode,jnode)+fact3
              elmat(jnode,inode)=elmat(jnode,inode)+fact3
           end do
        end do
     end do
  end if

end subroutine adrmat
