subroutine chm_elmpre_flamLet(&
     ielem,pnode,pgaus,elcon,elvel,elmas,elkey,eleps,elrrt,gpsha,gpcar,gplap,&
     gpcon,gpvel,gpmas,hleng,gpdiv,gpdis,gpprd,gprrt,gptur,gpden,gphco,&
     gpsph) 
  
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_elmpre_flamLet
  ! NAME 
  !    chm_elmpre
  ! DESCRIPTION
  !    Compute quantities at Gauss points for flamelet combustion model
  ! USES
  ! USED BY
  !    mod_chm_element_operations
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nclas_chm,&
                               diffu_chm,kfl_cotur_chm,nclas_chm,     &
                               ADR_chm,kfl_premix_chm,kfl_lookg_chm,mass_gp,rrt_gp
  use def_chemic, only      :  kfl_varZ_chm,kfl_varYc_chm
  use mod_ADR,    only      :  BDF

  implicit none
  integer(ip),  intent(in)  :: ielem,pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nclas_chm,*)
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: elmas(pnode,nclas_chm)
  real(rp),     intent(in)  :: elkey(pnode)
  real(rp),     intent(in)  :: eleps(pnode)
  real(rp),     intent(in)  :: elrrt(pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gplap(mgaus,mnode)                  
  real(rp),     intent(in)  :: hleng(3)
  real(rp),     intent(in)  :: gpden(mgaus)                  ! Density at gauss points
  real(rp),     intent(in)  :: gphco(mgaus)                  ! Conductivity at gauss points
  real(rp),     intent(in)  :: gpsph(mgaus)                  ! Specific heat at gauss points
  real(rp),     intent(in)  :: gptur(pgaus)                  ! Turbulence viscosity at gauss points
  real(rp),     intent(out) :: gpcon(pgaus,nclas_chm,*)
  real(rp),     intent(out) :: gpvel(ndime,pgaus)
  real(rp),     intent(out) :: gpmas(pgaus,nclas_chm)        ! Mass source term
  real(rp),     intent(out) :: gpdiv(pgaus)                  ! Velocity divergence
  real(rp),     intent(out) :: gpdis(pgaus,nclas_chm)        ! Dissipation rate for combustion model
  real(rp),     intent(out) :: gpprd(pgaus,nclas_chm)        ! Production term of the variance of c and z in the combustion model
  real(rp),     intent(out) :: gprrt(pgaus,nclas_chm)        ! {w_c * c} - {w_c}*{c} 
  real(rp)                  :: gpkey(mgaus)                  ! TKE from K-EPS or K-OMEGA models (RANS)
  real(rp)                  :: gpeps(mgaus)                  ! Dissipation from K-EPS or K-OMEGA models (RANS)
  real(rp)                  :: aux,delta,delta2,rdelta2,seci4
  real(rp)                  :: gpgve(ndime,ndime,pgaus)
  real(rp)                  :: sgs_Chi_Yc,res_Chi_Yc
  real(rp)                  :: sgs_Chi_z,res_Chi_z
  real(rp)                  :: rtau_turb                     ! Inverse of turbulent time scale
  real(rp)                  :: z_var,Yc_var
  real(rp)                  :: grad_Z(pgaus,ndime),mod_grad_Z(pgaus) 
  real(rp)                  :: grad_Yc(pgaus,ndime),mod_grad_Yc(pgaus)
  integer(ip)               :: igaus,iclas,inode,idime,jdime,icoun,itime
  integer(ip)               :: ivari,nvari
  
  !
  ! Initialization
  !
  gpmas = 0.0_rp
  gprrt = 0.0_rp
  gpvel = 0.0_rp
  gpdis = 0.0_rp
  gpprd = 0.0_rp
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
  gpmas = 0.0_rp 
  if (kfl_lookg_chm > 0) then
     gpmas = mass_gp(ielem) % a(:,:,1)
  else
      do iclas = 1,nclas_chm
        do igaus = 1,pgaus
           do inode = 1,pnode
              gpmas(igaus,iclas) = gpmas(igaus,iclas) + gpsha(inode,igaus) * elmas(inode,iclas)
           enddo
        enddo
     enddo
  end if
  
  !
  ! LES modelling 
  !
  if (kfl_cotur_chm < 0_ip ) then

    !
    ! Fluctuations of the reaction progress variable: Ycv, Y_c^2 or c_var 
    ! 
    if ( kfl_varYc_chm /= 0_ip ) then

       grad_Yc     = 0.0_rp
       mod_grad_Yc = 0.0_rp

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                grad_Yc(igaus,idime) = grad_Yc(igaus,idime) + gpcar(idime,inode,igaus) * elcon(inode,1,1)
             end do
          end do
       end do

       !
       ! Transport of reaction rate fluctuations: {w_c * c} or {w_Yc * Yc} 
       !
       if (kfl_lookg_chm > 0) then
          gprrt = rrt_gp(ielem) % a(:,:,1)
       else 
          do igaus = 1,pgaus
            do inode = 1,pnode
               gprrt(igaus,2) = gprrt(igaus,2) + gpsha(inode,igaus) * elrrt(inode)
            end do
          end do
       end if

       do igaus = 1,pgaus
          
          delta   = ( hleng(1) * hleng(2) * hleng(3) )** 0.3333333_rp
          delta2  = delta*delta
          rdelta2 = 1.0_rp / delta2

          !
          ! Compute |grad Yc|^2
          !
          mod_grad_Yc(igaus) = dot_product(grad_Yc(igaus,:),grad_Yc(igaus,:))

          !
          ! Compute strain rate square |S_ij|*|S_ij|
          !
          gpgve = 0.0_rp 
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                end do
             end do
          end do

          seci4 = 0.0_rp
          do idime = 1,ndime                     ! |S_ij|*|S_ij|
             do jdime = 1,ndime         
                seci4 = seci4 + 0.5_rp * gpgve(idime,jdime,igaus) * (gpgve(idime,jdime,igaus) + gpgve(jdime,idime,igaus))
             end do
          end do

          !
          ! Compute inverse of turbulent time scale: 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
          !
          rtau_turb = ( 3.24_rp * gptur(igaus) * seci4 * rdelta2 ) ** (0.33333_rp)

          !
          ! Non-premixed combustion based on Y_c^2 transport 
          ! - Domingo et al. (2008) Combust. Flame
          !
          if (kfl_premix_chm == 0_ip .and. kfl_varYc_chm == 2_ip ) then

             !
             ! Transport of reaction rate fluctuations: W_Y_c = 2 * {w_Y_c * Y_c}
             !
             gprrt(igaus,2) = 2.0_rp * gprrt(igaus,2)

             !
             ! Dissipation term: 2 * ( rho * D * |grad Y_c|^2  + C_x * rho * 1 / tau * Y_cv )
             !
             Yc_var = gpcon(igaus,2,1) - gpcon(igaus,1,1) * gpcon(igaus,1,1)
             res_Chi_Yc = gphco(igaus) / gpsph(igaus) * mod_grad_Yc(igaus)
             sgs_Chi_Yc = gpden(igaus) * rtau_turb * Yc_var

             gpdis(igaus,2) = 2.0_rp * ( res_Chi_Yc + sgs_Chi_Yc )

             !
             ! Production term: 2 mu_t / (Delta^2*Sc_t) * (Y_c)^2
             !
             gpprd(igaus,2) = 0.0_rp 

          !
          ! Premixed combustion or non-premixed combustion based on Y_c variance transport
          ! - Domingo et al. (2005) Combust. Flame (PREMIXED, c)
          ! - Galpin et al. (2008) Combust. Flame  (NON-PREMIXED, Yc)
          !
          else if ( kfl_premix_chm == 1_ip .or. ( kfl_premix_chm == 0_ip .and. kfl_varYc_chm == 1_ip ) ) then

             !
             ! Transport of reaction rate fluctuations: W_c  = 2 * ( {w_c * c} - {w_c}*{c} ) 
             !                                          W_Yc = 2 * ( {w_Yc * Yc} - {w_Yc}*{Yc} )
             !
             gprrt(igaus,2) = 2.0_rp * (gprrt(igaus,2) - gpmas(igaus,1) * gpcon(igaus,1,1))

             !
             ! Dissipation term: D_c  = 2 * C_x * rho * 1 / tau * c_var 
             !                   D_Yc = 2 * C_x * rho * 1 / tau * Y_cv
             !
             gpdis(igaus,2) = 2.0_rp * rtau_turb * gpden(igaus) * gpcon(igaus,2,1)

             !
             ! Production term: P_c  = 2 * rho * nu_t/Sc_t * |grad c|^2 
             !                  P_Yc = 2 * rho * nu_t/Sc_t * |grad Yc|^2 
             !
             gpprd(igaus,2) = 2.0_rp * gptur(igaus) * gpden(igaus) * mod_grad_Yc(igaus) / diffu_chm(1,1) 

          endif

       end do

    end if

    !
    ! Fluctuations of the mixture fraction: Z_var or Z^2
    ! 
    if ( kfl_varZ_chm /= 0_ip ) then

       grad_Z     = 0.0_rp
       mod_grad_Z = 0.0_rp

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                grad_Z(igaus,idime) = grad_Z(igaus,idime) + gpcar(idime,inode,igaus) * elcon(inode,3,1)
             end do
          end do
       end do

       do igaus = 1,pgaus

          delta   = ( hleng(1) * hleng(2) * hleng(3) )** 0.3333333_rp
          delta2  = delta*delta
          rdelta2 = 1.0_rp / delta2

          !
          ! Compute |grad Z|^2
          !
          mod_grad_Z(igaus) = dot_product(grad_Z(igaus,:),grad_Z(igaus,:))
  
          !
          ! Compute strain rate square |S_ij|*|S_ij|
          !
          gpgve = 0.0_rp 
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                end do
             end do
          end do

          seci4 = 0.0_rp
          do idime = 1,ndime                     ! |S_ij|*|S_ij|
             do jdime = 1,ndime         
                seci4 = seci4 + 0.5_rp * gpgve(idime,jdime,igaus) * (gpgve(idime,jdime,igaus) + gpgve(jdime,idime,igaus))
             end do
          end do

          !
          ! Compute inverse of turbulent time scale: 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
          !
          rtau_turb = ( 3.24_rp * gptur(igaus) * seci4 * rdelta2 ) ** (0.33333_rp) 

          !
          ! Transport Z variance equation           
          ! - Domingo et al. (2005) Combust. Flame 
          !
          if (kfl_varZ_chm == 1_ip) then
             !
             ! Production term: P_Zv = 2 * D_t * |grad Z|^2
             !       
             gpprd(igaus,4) = 2.0_rp * gptur(igaus) * gpden(igaus) * mod_grad_Z(igaus) / diffu_chm(1,1)

             !
             ! Dissipation term: D_Zv = 2 * rho * 1 / tau * Z_var
             !       
             gpdis(igaus,4) = 2.0_rp * rtau_turb * gpden(igaus) * gpcon(igaus,4,1) 

          !
          ! Transport Z^2 equation           
          ! - Knudsen et al. (2012) Phys. Fluids 
          !
          elseif (kfl_varZ_chm == 2_ip) then
             !
             ! Production term: P_Z*Z = 0 
             !
             gpprd(igaus,4) = 0.0_rp
             !
             ! Dissipation term: D_Z = 2 D |grad Z|^2 + X_z,sgs, X_z,sgs = 2 * rho * 1 / tau * Z_var
             !
             z_var = gpcon(igaus,4,1) - gpcon(igaus,3,1) * gpcon(igaus,3,1) 
             res_Chi_z = 2.0_rp * gphco(igaus) / gpsph(igaus) * mod_grad_Z(igaus)
             sgs_Chi_z = 2.0_rp * rtau_turb * gpden(igaus) * z_var

             gpdis(igaus,4) = res_Chi_z + sgs_Chi_z  
             
          endif

       end do

    end if


  else
    !
    ! variances are not used
    !
    gpdis = 0.0_rp
    gpprd = 0.0_rp

  end if
  
end subroutine chm_elmpre_flamLet
