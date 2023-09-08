subroutine chm_elmvel(&
     pnode,pgaus,elcon,gpDik,gpgDk,gpsha,gpcar,gplap,gpvel,gpgve,gpmol,gpgmo,gphmo)

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
  use def_chemic, only      :  nspec_chm,dtinv_chm,pabdf_chm,kfl_corve_chm,kfl_diffu_chm
  use def_master, only      :  speci
  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nspec_chm)       ! NOTICE THIS CAN BE PREVIOUS ITERATION VALUES
  real(rp),     intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus),gplap(mgaus,pnode)  ! 
  real(rp),     intent(in)  :: gpDik(pgaus,nspec_chm)       ! Difussion coefficients
  real(rp),     intent(in)  :: gpgDk(ndime,pgaus,nspec_chm) ! Gradient of difussion coeffs
  real(rp),     intent(out) :: gpvel(ndime,pgaus)           ! Corrected advection velocity
  real(rp),     intent(out) :: gpgve(pgaus)                 ! div u
  real(rp),     intent(in) :: gpmol(pgaus)
  real(rp),     intent(in) :: gpgmo(ndime,pgaus)
  real(rp),optional,intent(in) :: gphmo(pgaus)
  integer(ip)               :: igaus,ispec,inode,jspec,idime
  real(rp)                  :: gpden(pgaus)


  ! Various initializations
  do igaus=1,pgaus
     gpgve(igaus) = 0.0_rp
     do idime = 1,ndime
        gpvel(idime,igaus) = 0.0_rp
     enddo
  enddo
  
  if (kfl_corve_chm==1) then !Compute correction velocity

     select case (kfl_diffu_chm)
        case(1)  ! Approximate mole fraction with constant density * mass fraction
           !
           ! Here we compute correction velocity to keep the total mass conservation ok
           !
           ! u_Total = u + sum_k Dk grad(Y_k) = u+uc
           ! div(uc) = sum_k grad(Dk).grad(Y_k) + sum_k Dk lapl(Y_k)
           !
           ! NOTICE that we might be using the concentration from the previous step, 
           ! which one is decided in assembly routine chm_elmcom by the kfl_gauss_chm flag
           !
           do ispec = 1,nspec_chm             
              do igaus = 1, pgaus
                 do inode = 1,pnode
                    gpgve(igaus) = gpgve(igaus) &
                         + gpDik(igaus,ispec) * gplap(igaus,inode) * elcon(inode,ispec) 
                    do idime = 1,ndime          
                       gpvel(idime,igaus) = gpvel(idime,igaus)  &
                            + gpdik(igaus,ispec) * gpcar(idime,inode,igaus) * elcon(inode,ispec)
                       gpgve(igaus) = gpgve(igaus) &
                            + gpgdk(idime,igaus,ispec) * gpcar(idime,inode,igaus) * elcon(inode,ispec)
                    end do
                 enddo
              enddo
           enddo
        case(2)
           !
           ! IF we use the advanced for of the diffusion term (mole fraction) then we need to add terms
           !
           if (kfl_diffu_chm==2) then
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
           endif
        end select

     endif

end subroutine chm_elmvel
