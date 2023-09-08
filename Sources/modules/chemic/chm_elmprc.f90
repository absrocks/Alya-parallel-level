subroutine chm_elmprc(&
  iclas,pgaus,gpcon,gpden,gpgde,gpDik,gpgDk,gpmas,gpvel,gpvec,gpgve,&
       gpadv,gpdif,gpgrd,gprea,gprhs,gpmol,gpgmo,gphmo,gpgac,gplac)

  !-----------------------------------------------------------------------
  !****f* partis/chm_elmprc
  ! NAME 
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute terms for each species ADR equation 
  ! USES
  ! USED BY
  !    chm_elmcom
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus
  use def_chemic, only      :  nspec_chm,kfl_react_chm,kfl_diffu_chm, kfl_activ_chm
  use def_master, only      :  speci
  implicit none
  integer(ip),  intent(in)  :: iclas,pgaus
  real(rp),     intent(in)  :: gpden(pgaus),gpgde(ndime,pgaus)
  real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm)
  real(rp),     intent(in)  :: gpDik(pgaus,nspec_chm)  ! Species diffusion coeff
  real(rp),     intent(in)  :: gpgDk(ndime,pgaus,nspec_chm) ! Grad of diffusion coeff
  real(rp),     intent(in)  :: gpmas(pgaus,nspec_chm)  ! Species source term (reaction rate)
  real(rp),     intent(in)  :: gpvel(ndime,pgaus),gpvec(ndime,pgaus) ! Velocity and correction velocity
  real(rp),     intent(in)  :: gpgve(pgaus) ! Gradient of correction velocity

  real(rp), intent(in) :: gpmol(pgaus)
  real(rp), intent(in) :: gpgmo(ndime,pgaus)
  real(rp), intent(in) :: gphmo(pgaus)
  real(rp), intent(in)    :: gpgac(ndime,pgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp), intent(in)    :: gplac(pgaus,nspec_chm)                ! Laplacian(activity)/activity

  real(rp),     intent(out) :: gpadv(ndime,pgaus)
  real(rp),     intent(out) :: gpdif(pgaus)
  real(rp),     intent(out) :: gpgrd(ndime,pgaus)
  real(rp),     intent(out) :: gprea(pgaus)
  real(rp),     intent(out) :: gprhs(pgaus)
  integer(ip)               :: igaus,idime
  
  real(rp) :: dif_temp

  do igaus = 1,pgaus
     !
     ! Mass source term
     !
     gprhs(igaus) = gpmas(igaus,iclas)
     !
     ! Diffusion coefficient
     !
     gpdif(igaus) = gpden(igaus) * gpDik(igaus,iclas)
     !
     ! Diffusion gradient
     !
     do idime = 1,ndime
        gpgrd(idime,igaus)  =  & !0.0_rp
            +  gpden(igaus) * gpgDk(idime,igaus,iclas) &
            +  gpDik(igaus,iclas) * gpgde(idime,igaus)
     enddo
     !
     ! Reaction
     !
     if (kfl_react_chm == 1) then
        gprea(igaus) = gpgve(igaus) * gpden(igaus)
        do idime=1,ndime
           gprea(igaus) = gprea(igaus) + gpvec(idime,igaus) * gpgde(idime,igaus)
        enddo
     else ! We assemble reaction term on the rhs
        gprea(igaus) = 0.0_rp
        gprhs(igaus) = gprhs(igaus) - gpcon(igaus,iclas) * gpgve(igaus)* gpden(igaus)
        do idime=1,ndime
           gprhs(igaus) = gprhs(igaus) - gpcon(igaus,iclas) * gpvec(idime,igaus) * gpgde(idime,igaus)
        enddo
     endif
     ! 
     ! Advection
     !
     do idime = 1,ndime
        gpadv(idime,igaus) = gpvel(idime,igaus) + gpvec(idime,igaus)
     enddo
  enddo

  !
  ! Correction because of complete diffusion term
  !
  if (kfl_diffu_chm==2) then
     dif_temp=gpdif(igaus) ! Before change is needed for one eq. below (*)
     do idime=1,ndime
        gprea(igaus) = gprea(igaus) + 2.0_rp * gpgmo(idime,igaus)*gpgrd(idime,igaus)*gpdif(igaus)/speci(iclas)%weigh
     enddo
     gprea(igaus) = gprea(igaus) - gpdif(igaus) * gphmo(igaus)/speci(iclas)%weigh
     do idime=1,ndime
        gpgrd(idime,igaus) =  (gpgrd(idime,igaus)*gpmol(igaus) + 2.0_rp*gpgmo(idime,igaus)*gpdif(igaus)) / speci(iclas)%weigh
     enddo
     gpdif(igaus)=gpdif(igaus)*gpmol(igaus)/speci(iclas)%weigh
     if (kfl_activ_chm == 1) then ! Activity coefficients present
        do idime=1,ndime
           gprea(igaus) = gprea(igaus) + gpgac(idime,igaus,iclas)*  &
                ( gpdif(igaus) * gpgac(idime,igaus,iclas) - gpgrd(idime,igaus) + dif_temp * gpgmo(idime,igaus)/speci(iclas)%weigh) ! (*) here
        enddo
        gprea(igaus) = gprea(igaus) - gpdif(igaus) * gplac(igaus,iclas)
        do idime=1,ndime
           gpgrd(idime,igaus) =  gpgrd(idime,igaus) + gpdif(igaus) * gpgac(idime,igaus,iclas)
        enddo
     endif
  endif


end subroutine chm_elmprc
