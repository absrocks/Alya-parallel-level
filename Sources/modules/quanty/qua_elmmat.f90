subroutine qua_elmmat(&
     pnode,pgaus,gpdif,gposc,gpvol,gpsha,gpcar,&
     gpcod,gprhs,gpcou,gpaxi,elmat,gpelf,gpesp,elrhs,gpion,&
     gphar,gppxc)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmmat
  ! NAME
  !   qua_elmmat
  ! DESCRIPTION
  !    Compute elemental matrix and RHS
  ! OUTPUT
  !    ELMAT ... LHS matrix
  !    ELRHS ... RHS vector
  ! USES
  ! USED BY
  !    qua_elmope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,kfl_naxis,mnode
  use def_quanty
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus
  real(rp),    intent(in)    :: gpdif(pgaus)
  real(rp),    intent(in)    :: gposc(pgaus),gpvol(pgaus),gpelf(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus),gpesp(pgaus)
  real(rp),    intent(in)    :: gpcou(pgaus),gpaxi(pgaus)
  real(rp),    intent(in)    :: gpcod(ndime,pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus),gpion(pgaus),gphar(pgaus),gppxc(pgaus)
  real(rp),    intent(out)   :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,jnode,kdime,igaus
  real(rp)                   :: fact1,fact2
  !
  ! Initialization
  !
  do inode = 1,pnode 
     do jnode = 1,pnode
        elmat(inode,jnode) = 0.0_rp
     end do
     elrhs(inode) = 0.0_rp
  end do
  !
  ! Diffusion term:  ( gpdif * grad(Nj) , grad(Ni) )
  !

  ! termino independiente para Poisson si no es cero
  !               gprhs * Ni

  do igaus = 1,pgaus
     fact1 = gpdif(igaus) * gpvol(igaus)
     fact2 = gprhs(igaus) * gpvol(igaus)
     do inode = 1,pnode
        elrhs(inode) = elrhs(inode) + fact2 * gpsha(inode,igaus)
        do jnode = 1,pnode
           do kdime = 1,ndime
              elmat(inode,jnode) = elmat(inode,jnode) + fact1 * &
                   gpcar(kdime,inode,igaus) * gpcar(kdime,jnode,igaus)
           end do
        end do
     end do
  end do
  !
  ! termino axisymetrico
  !
  do igaus = 1,pgaus
     fact1 = gpvol(igaus) * gpaxi(igaus)
     do inode = 1,pnode
        do jnode = 1,pnode
           elmat(inode,jnode) = elmat(inode,jnode) + fact1 &
                * gpsha(inode,igaus) * gpsha(jnode,igaus)
        end do
     end do
  end do
  !
  ! termino espheric
  !
  do igaus=1,pgaus
     fact1=gpvol(igaus)*gpesp(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do
  !
  ! termino de coulomb  (gpcou()*Nj*Ni)
  !
  do igaus=1,pgaus
     fact1=gpvol(igaus)*gpcou(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do
  !
  ! termino de oscilador  (gpcou()*Nj*Ni)
  !
  do igaus=1,pgaus
     fact1=gpvol(igaus)*gposc(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do
  !
  ! termino de E field  (gpelf()*Nj*Ni)
  !
  do igaus=1,pgaus
     fact1=gpvol(igaus)*gpelf(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do



  ! termino de potencial ionico 
  do igaus=1,pgaus
     fact1= gpvol(igaus)*gpion(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do

  ! termino de potencial XC
  do igaus=1,pgaus
     fact1=gpvol(igaus)*gppxc(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do

  ! termino de Hartree

  do igaus=1,pgaus
     fact1=gpvol(igaus)*gphar(igaus)
     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+fact1&
                *gpsha(inode,igaus)*gpsha(jnode,igaus)
        end do
     end do
  end do




end subroutine qua_elmmat
