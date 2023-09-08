subroutine chm_levvel(&
     pnode,pgaus,elcon,gpDik,gpgDk,gpsha,gpcar,gplap,gpvel,gpgve,gpmol,gpgmo,gphmo,gpgac,gplac)

  !-----------------------------------------------------------------------
  !****f* partis/chm_elmprc
  ! NAME 
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute correction velocity (ensures conservation of total mass)  
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus,mnode
  use def_chemic, only      :  nspec_chm,dtinv_chm,pabdf_chm,kfl_corve_chm
  use def_master, only      :  speci
  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nspec_chm)  ! NOTICE THIS CAN BE PREVIOUS ITERATION VALUES
  real(rp),     intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus),gplap(mgaus,pnode)  ! 
  real(rp),     intent(in)  :: gpDik(pgaus,nspec_chm)    ! Difussion coefficients
  real(rp),     intent(in)  :: gpgDk(ndime,pgaus,nspec_chm) ! Gradient of difussion coeffs
  real(rp),     intent(inout) :: gpvel(ndime,pgaus) ! Corrected advection velocity
  real(rp),     intent(inout) :: gpgve(pgaus) ! div u
  real(rp),intent(in) :: gpmol(pgaus)
  real(rp),intent(in) :: gpgmo(ndime,pgaus)
  real(rp),intent(in) :: gphmo(pgaus)
  real(rp),intent(in)    :: gpgac(ndime,pgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp),intent(in)    :: gplac(pgaus,nspec_chm)                ! Laplacian(activity)/activity
  integer(ip)               :: igaus,ispec,inode,jspec,idime

  do ispec = 1,nspec_chm             
     do igaus = 1, pgaus
        do inode = 1,pnode
           gpgve(igaus) = gpgve(igaus) &
                + gpDik(igaus,ispec) * gplap(igaus,inode) * elcon(inode,ispec) * gpmol(igaus)  &
                / speci(ispec) % weigh &
                + gpDik(igaus,ispec) * elcon(inode,ispec) * gphmo(igaus) / speci(ispec) % weigh
           do idime = 1,ndime          
              gpvel(idime,igaus) = gpvel(idime,igaus)  &
                   + gpdik(igaus,ispec) * gpcar(idime,inode,igaus) * elcon(inode,ispec) * gpmol(igaus)  &
                   / speci(ispec) % weigh &
                   + gpdik(igaus,ispec) * elcon(inode,ispec) * gpgmo(idime,igaus)  / speci(ispec) % weigh
              gpgve(igaus) = gpgve(igaus) &
                   + gpgdk(idime,igaus,ispec) * gpcar(idime,inode,igaus) * elcon(inode,ispec) * gpmol(igaus)  &
                   / speci(ispec) % weigh &
                   + gpgdk(idime,igaus,ispec) * gpgmo(idime,igaus) * elcon(inode,ispec) / speci(ispec) % weigh &
                   + gpdik(igaus,ispec) * gpgmo(idime,igaus) * gpcar(idime,inode,igaus) * elcon(inode,ispec) &
                   / speci(ispec) % weigh
           end do
        enddo
     enddo
  enddo

end subroutine chm_levvel
