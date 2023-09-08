subroutine chm_elmpre_species(&
     pnode,pgaus,elcon,elden,elvel,elDik,elmas,elmut,elsen,gpsha,gpcar,gplap,&
     gpcon,gpvel,gpDik,gpgDk,gpmas,elmol,gpmol,gpgmo,gphmo,gpfar,&
     gpvol,gpdiv,gpspe,gpthi,sgsef,gpfac,gpsen,gptur) 
  
  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmpre_species
  ! NAME 
  !    chm_elmpre_species
  ! DESCRIPTION
  !    Compute quantities at Gauss points
  ! USES
  ! USED BY
  !    chm_elmpre_species
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nspec_chm,kfl_diffu_chm,kfl_tfles_chm,kfl_model_chm, &
                               diffu_chm,kfl_cotur_chm,nclas_chm,     &
                               ADR_chm
  use def_master, only      :  speci,kfl_paral
  use mod_ADR,    only      :  BDF
  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nspec_chm,*)
  real(rp),     intent(in)  :: elden(pnode)
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: elDik(pnode,nspec_chm)
  real(rp),     intent(in)  :: elmas(pnode,nspec_chm)
  real(rp),     intent(in)  :: elmut(pnode)
  real(rp),     intent(in)  :: elsen(pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gplap(mgaus,mnode)                  
  real(rp),     intent(in)  :: gpvol(mgaus)
  real(rp),     intent(out) :: gpcon(pgaus,nspec_chm,*)
  real(rp),     intent(out) :: gpvel(ndime,pgaus)
  real(rp),     intent(out) :: gpDik(pgaus,nspec_chm)        ! Difussion coefficients
  real(rp),     intent(out) :: gpgDk(ndime,pgaus,nspec_chm)  ! Gradient of difussion coeffs
  real(rp),     intent(out) :: gpmas(pgaus,nspec_chm)        ! Mass source term
  real(rp),     intent(in)  :: elmol(pnode)
  real(rp),     intent(out) :: gpmol(pgaus)
  real(rp),     intent(out) :: gpgmo(ndime,pgaus)
  real(rp),     intent(out) :: gphmo(pgaus)
  real(rp),     intent(out) :: gpfar(mgaus)
  real(rp),     intent(out) :: gpdiv(pgaus)                   ! Velocity divergence
  real(rp),     intent(out) :: gpspe(mgaus)                   ! Flame speed at gauss points
  real(rp),     intent(out) :: gpthi(mgaus)                   ! Flame thickness at gauss points
  real(rp),     intent(out) :: sgsef(mgaus)                   ! Efficiency function for the wrinkling of the flame front
  real(rp),     intent(out) :: gpfac(mgaus)                   ! Dynamic thickening factor F(c) for DTFLES 
  real(rp),     intent(out) :: gpsen(mgaus)                   ! Flame sensor for DTFLES
  real(rp),     intent(out) :: gptur(mgaus)                   ! Turbulence viscosity at gauss points
  real(rp)                  :: gpden(mgaus)                   ! Density at gauss points
  real(rp)                  :: aux,vol,vol2,rvol2,new
  integer(ip)               :: igaus,ispec,iclas,inode,idime,ifuel,ioxo2,ioxn2,icoun,itime

  !
  ! Initialization
  !
  gpDik = 0.0_rp
  gpmas = 0.0_rp
  gpgDk = 0.0_rp 
  gpfar = 0.0_rp
  gpspe = 0.0_rp
  gpthi = 0.0_rp
  gptur = 0.0_rp
  sgsef = 0.0_rp
  gpfac = 0.0_rp
  gpsen = 0.0_rp
  gpden = 0.0_rp
  gpvel = 0.0_rp
  !
  ! Concentration
  !
  do iclas = 1,nclas_chm
     do igaus = 1,pgaus
        gpcon(igaus,iclas,1) = 0.0_rp
        do inode = 1,pnode
           gpcon(igaus,iclas,1) = gpcon(igaus,iclas,1)&
                                + gpsha(inode,igaus) * elcon(inode,iclas,1)
        end do
     end do
  end do
  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then
     do iclas = 1,nclas_chm
        do igaus = 1,pgaus
           gpcon(igaus,iclas,2) = 0.0_rp
           do inode = 1,pnode
              gpcon(igaus,iclas,2) = gpcon(igaus,iclas,2)&
                                   + gpsha(inode,igaus) * elcon(inode,iclas,2)
           end do
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do igaus = 1,pgaus
                 gpcon(igaus,iclas,itime) = 0.0_rp
                 do inode = 1,pnode
                    gpcon(igaus,iclas,itime) = gpcon(igaus,iclas,itime)&
                                             + gpsha(inode,igaus) * elcon(inode,iclas,itime)
                 end do
              end do
           end do
        end do
     end if
  end if

  !
  ! Fluid velocity
  !
  do igaus = 1, pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gpvel(idime,igaus) = gpvel(idime,igaus) &
                +  gpsha(inode,igaus) * elvel(idime,inode) 
        end do
     enddo
  enddo
  !
  ! Fluid velocity divergence
  !
  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gpdiv(igaus) = gpdiv(igaus) + gpcar(idime,inode,igaus) * elvel(idime,inode)
        end do
     end do
  end do
  !
  ! Mass source terms
  !
  do ispec = 1,nspec_chm
     do igaus = 1,pgaus
        do inode = 1,pnode
           gpmas(igaus,ispec) = gpmas(igaus,ispec) + gpsha(inode,igaus) * elmas(inode,ispec)
        enddo
     enddo
  enddo
  !
  ! turbulent viscosity
  !
  if (kfl_cotur_chm /= 0) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           gptur(igaus) = gptur(igaus) + gpsha(inode,igaus) * elmut(inode)
        end do
     end do
  end if 
  !
  ! DTFLES combustion model
  !
  if (kfl_tfles_chm == 1_ip) then
     ifuel = 1  ! Fuel
     ioxo2 = 2  ! O2
     ioxn2 = 5  ! N2
     !
     ! Flame/turbulence interactions
     !
     do igaus = 1,pgaus
        gpfar(igaus) = gpcon(igaus,ifuel,1) / ( gpcon(igaus,ioxo2,1) + gpcon(igaus,ioxn2,1))
        do inode = 1,pnode
           gptur(igaus) = max(gptur(igaus),0.0_rp)
           gpsen(igaus) = gpsen(igaus) + gpsha(inode,igaus) * elsen(inode)
        end do
        call chm_flame_prop(gpfar(igaus),gpspe(igaus),gpthi(igaus))
        call chm_tfles_sgs(gpspe(igaus),gpthi(igaus),&
                           gptur(igaus),gpvol(igaus),sgsef(igaus),&
                           gpfac(igaus),gpsen(igaus))
     end do
     !
     ! Density
     !
     do igaus = 1,pgaus
        do inode = 1,pnode
           gpden(igaus) = gpden(igaus)  &
                        + gpsha(inode,igaus) * elden(inode)
        end do
     end do

  end if
  !
  ! Diffusion coefficients
  !
  if (kfl_tfles_chm >= 1) then
    do ispec = 1,nspec_chm          !!For DTFLES => D_eff = D_lam * E * F_dyn + D_tur * (1 - OMEGA); F_dyn = 1+(F_max-1)OMEGA
        do igaus = 1,pgaus      
           do inode = 1,pnode
              gpDik(igaus,ispec) = gpDik(igaus,ispec) & 
                                 + gpsha(inode,igaus) * elDik(inode,ispec) * sgsef(igaus) * gpfac(igaus)  &
                                 + gptur(igaus) / diffu_chm(1,1) * (1.0_rp - gpsen(igaus))   

              do idime = 1,ndime
                 gpgDk(idime,igaus,ispec) = gpgDk(idime,igaus,ispec)  &
                                          + gpcar(idime,inode,igaus) * elDik(inode,ispec) * sgsef(igaus) * gpfac(igaus) &
                                          + gpcar(idime,inode,igaus) * gptur(igaus) / diffu_chm(1,1)
              enddo
           end do
        end do
     enddo 
  else
     do ispec = 1,nspec_chm
        do igaus = 1,pgaus
           do inode = 1,pnode
              gpDik(igaus,ispec) = gpDik(igaus,ispec)&
                                 + gpsha(inode,igaus) * elDik(inode,ispec)
              do idime = 1,ndime
                 gpgDk(idime,igaus,ispec)= gpgDk(idime,igaus,ispec) + gpcar(idime,inode,igaus) * elDik(inode,ispec)
              enddo
           end do
        end do
     enddo
  endif
  !
  ! Molar mass
  !
  if (kfl_diffu_chm==2) then
     gpmol=0.0_rp
     gpgmo=0.0_rp
     gphmo=0.0_rp
     do igaus = 1,pgaus
        do inode = 1,pnode
           gpmol (igaus) = gpmol(igaus) + gpsha(inode,igaus) * elmol(inode)
           gphmo (igaus) = gphmo(igaus) + gplap(igaus,inode) * elmol(inode)
           do idime = 1,ndime
              gpgmo(idime,igaus)=gpgmo(idime,igaus) +  gpcar(idime,inode,igaus) * elmol(inode)
           enddo
        enddo
     enddo
  endif

end subroutine chm_elmpre_species
