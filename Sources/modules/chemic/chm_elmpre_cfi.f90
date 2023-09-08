subroutine chm_elmpre_cfi(&
     pnode,pgaus,elcon,elden,elvel,elDik,elmas,elkey,eleps,elrrt,gpsha,gpcar,gplap,&
     gpcon,gpvel,gpDik,gpgDk,gpmas,gpmol,gpgmo,gphmo,&
     gpvol,gpdiv,gpdis,gpprd,gprrt,gptur,gpden,gphco,gpsph) 
  
  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmpre_cfi
  ! NAME 
  !    chm_elmpre
  ! DESCRIPTION
  !    Compute quantities at Gauss points for CFI combustion model
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nspec_chm,kfl_diffu_chm,kfl_tfles_chm,kfl_model_chm, &
                               diffu_chm,kfl_cotur_chm,kfl_tucfi_chm,nclas_chm,     &
                               ADR_chm,kfl_uncfi_chm
  use def_master, only      :  speci,kfl_paral
  use mod_ADR,    only      :  BDF

  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nspec_chm,*)
  real(rp),     intent(in)  :: elden(pnode)
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: elDik(pnode,nspec_chm)
  real(rp),     intent(in)  :: elmas(pnode,nspec_chm)
  real(rp),     intent(in)  :: elkey(pnode)
  real(rp),     intent(in)  :: eleps(pnode)
  real(rp),     intent(in)  :: elrrt(pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gplap(mgaus,mnode)                  
  real(rp),     intent(in)  :: gpvol(mgaus)
  real(rp),     intent(in)  :: gpden(mgaus)                  ! Density at gauss points
  real(rp),     intent(in)  :: gphco(mgaus)                  ! Conductivity at gauss points
  real(rp),     intent(in)  :: gpsph(mgaus)                  ! Specific heat at gauss points
  real(rp),     intent(in)  :: gptur(mgaus)                  ! Turbulence viscosity at gauss points
  real(rp),     intent(out) :: gpcon(pgaus,nspec_chm,*)
  real(rp),     intent(out) :: gpvel(ndime,pgaus)
  real(rp),     intent(out) :: gpDik(pgaus,nspec_chm)        ! Difussion coefficients
  real(rp),     intent(out) :: gpgDk(ndime,pgaus,nspec_chm)  ! Gradient of difussion coeffs
  real(rp),     intent(out) :: gpmas(pgaus,nspec_chm)        ! Mass source term
  real(rp),     intent(out) :: gpmol(pgaus)
  real(rp),     intent(out) :: gpgmo(ndime,pgaus)
  real(rp),     intent(out) :: gphmo(pgaus)
  real(rp),     intent(out) :: gpdiv(pgaus)                  ! Velocity divergence
  real(rp),     intent(out) :: gpdis(mgaus)                  ! Dissipation rate for CFI model
  real(rp),     intent(out) :: gpprd(pgaus,nspec_chm)        ! Production term of the variance of c and f in the CFI model
  real(rp),     intent(out) :: gprrt(mgaus)                  ! {w_c * c} - {w_c}*{c} 
  real(rp)                  :: gpkey(mgaus)                  ! TKE from K-EPS or K-OMEGA models (RANS)
  real(rp)                  :: gpeps(mgaus)                  ! Dissipation from K-EPS or K-OMEGA models (RANS)
  real(rp)                  :: aux,vol,vol2,rvol2,new
  real(rp)                  :: prod_turb,prod_lam
  integer(ip)               :: igaus,ispec,iclas,inode,idime,ifuel,ioxo2,ioxn2,icoun,itime
  integer(ip)               :: ivari,nvari
  !
  ! Initialization
  !
  gpDik = 0.0_rp
  gpmas = 0.0_rp
  gpgDk = 0.0_rp 
  gprrt = 0.0_rp
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
  ! CFI-RANS modelling 
  !
  if (kfl_cotur_chm > 0_ip .and. kfl_tucfi_chm > 0_ip) then                
    !
    ! Transport of reaction rate fluctuations: {w_c * c} 
    !
    if (kfl_tucfi_chm /= 2_ip) then     ! Deactivate when VRPV = OFF & VMF = ON
      do igaus = 1,pgaus
         do inode = 1,pnode
            gprrt(igaus) = gprrt(igaus) + gpsha(inode,igaus) * elrrt(inode)
         end do
         gprrt(igaus) = 2.0_rp * (gprrt(igaus) - gpmas(igaus,1_ip) * gpcon(igaus,1_ip,1))
      end do
    end if
    !
    ! Dissipation term (RPV & MF): D_k = - C_x * rho * eps / k * Phi, rho is added in chm_elmprc_cfi 
    !
    do igaus = 1,pgaus
       gpkey(igaus) = 0.0_rp
       gpeps(igaus) = 0.0_rp
       do inode = 1,pnode
          gpkey(igaus) = gpkey(igaus) + gpsha(inode,igaus) * elkey(inode) 
          gpeps(igaus) = gpeps(igaus) + gpsha(inode,igaus) * eleps(inode) 
       end do
       if(gpkey(igaus) <= 0.0_rp) gpkey(igaus) = 1.0e-10_rp
       gpdis(igaus) = 2.0_rp * gpeps(igaus) / gpkey(igaus)     !! Reaction term (dissipation rate) 
    end do
    !
    ! Production term (RPV & MF): P_k = 2 * mu_t / Sc_t * |grad Phi|**2    
    ! 
    gpprd = 0.0_rp
    if     (kfl_tucfi_chm == 1_ip .or. kfl_tucfi_chm == 2_ip ) then    !! Either RPV or MF variance
       nvari = 1_ip
       ivari = kfl_tucfi_chm
    elseif (kfl_tucfi_chm == 3_ip ) then                                 !! RPV and MF variances activated
       nvari = 2_ip
    end if

    do icoun = 1,nvari 
       if (kfl_tucfi_chm == 3_ip ) ivari = icoun
       ispec = 2_ip * ivari
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                gpprd(igaus,ispec) = gpprd(igaus,ispec) + gpcar(idime,inode,igaus) * elcon(inode,ispec-1_ip,1)
             end do
          end do
          aux = gpprd(igaus,ispec)
          gpprd(igaus,ispec) = 2.0_rp * aux * aux * gptur(igaus) / diffu_chm(1,1) 
       end do
    end do

  !
  ! CFI-LES modelling 
  !
  elseif (kfl_cotur_chm < 0_ip .and. kfl_tucfi_chm > 0_ip) then
    !
    ! Transport of reaction rate fluctuations: {w_c * c} - {w_c}*{c} 
    ! - Domingo et al. (2005) Combust. Flame
    !
    if (kfl_tucfi_chm /= 2_ip) then     ! Deactivate when VRPV = OFF & VMF = ON
       do igaus = 1,pgaus
          do inode = 1,pnode
             gprrt(igaus) = gprrt(igaus) + gpsha(inode,igaus) * elrrt(inode)
          end do
          if (kfl_uncfi_chm == 1_ip) then
             gprrt(igaus) = 2.0_rp * gprrt(igaus)
          else
             gprrt(igaus) = 2.0_rp * (gprrt(igaus) - gpmas(igaus,1_ip) * gpcon(igaus,1_ip,1))
          endif
       end do
    end if
    !
    ! Dissipation term: 2 nu_t / (Delta^2*Sc_t) 
    ! - Domingo et al. (2005) Combust. Flame
    !
    do igaus = 1,pgaus
       vol   = gpvol(igaus)**(0.33333_rp)
       vol2  = vol*vol
       rvol2 = 1.0_rp / vol2
       gpdis(igaus) = 2.0_rp * gptur(igaus) * rvol2 / diffu_chm(1,1) 
    end do
    !
    ! Production term: 2 rho nu_t/Sc_t * (grad c)^2  
    ! - Domingo et al. (2005) Combust. Flame
    !
    gpprd = 0.0_rp
    if     (kfl_tucfi_chm == 1_ip .or. kfl_tucfi_chm == 2_ip ) then    !! Either RPV or MF variance
       nvari = 1_ip
       ivari = kfl_tucfi_chm
    elseif (kfl_tucfi_chm == 3_ip ) then                               !! RPV and MF variances activated
       nvari = 2_ip
    end if

    if (kfl_uncfi_chm == 1_ip) then
       !
       ! Unscale model: Y_c^2 & Z_v
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                gpprd(igaus,2_ip) = gpprd(igaus,ispec) + gpcar(idime,inode,igaus) * elcon(inode,1_ip,1)
                gpprd(igaus,4_ip) = gpprd(igaus,ispec) + gpcar(idime,inode,igaus) * elcon(inode,3_ip,1)
             end do
          end do
          aux               = gpprd(igaus,4_ip)
          gpprd(igaus,4_ip) = 2.0_rp * gptur(igaus) * gpden(igaus) * aux * aux / diffu_chm(1,1)

          aux               = gpprd(igaus,2_ip)
          vol               = gpvol(igaus)**(0.33333_rp)
          vol2              = vol*vol
          rvol2             = 1.0_rp / vol2

          prod_turb = gptur(igaus) * gpden(igaus) * gpcon(igaus,1_ip,1) * gpcon(igaus,1_ip,1) * rvol2 / diffu_chm(1,1)
          prod_lam  = gphco(igaus) / gpsph(igaus) * aux * aux

          gpprd(igaus,2_ip) = 2.0_rp * (  prod_turb - prod_lam )
 
       end do
    else
       !
       ! Scale model: c_v & Z_v
       !
       do icoun = 1,nvari
          if (kfl_tucfi_chm == 3_ip ) ivari = icoun 
          ispec = 2_ip * ivari
          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   gpprd(igaus,ispec) = gpprd(igaus,ispec) + gpcar(idime,inode,igaus) * elcon(inode,ispec-1_ip,1)
                end do
             end do
             aux = gpprd(igaus,ispec)
             gpprd(igaus,ispec) = 2.0_rp * gptur(igaus) * gpden(igaus) * aux * aux / diffu_chm(1,1) 
          end do
       end do
    end if

  else
    !
    ! variances are not used
    !
    gpdis = 0.0_rp
    gpprd = 0.0_rp
  end if
  !
  ! Diffusion to compute the gradient of diffusion coefficient gpgDk
  !
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
  
end subroutine chm_elmpre_cfi
