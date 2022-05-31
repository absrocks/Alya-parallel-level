subroutine chm_calc_scalar_dissip(&
     pnode,pgaus,elcon,elvel,gpsha,gpcar,&
     hleng,gptur,gphco,gpsph,gpden,res_Chi_Yc,sgs_Chi_Yc,res_Chi_z,&
     sgs_Chi_z) 
  
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_calc_scalar_dissip
  ! NAME 
  !    chm_calc_scalar_dissip
  ! DESCRIPTION
  !    Compute scalar dissipation rate at Gauss points for flamelet combustion model
  ! USES
  ! USED BY
  !    chm_post_scalar_dissipation_rate
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nclas_chm, &
                               ADR_chm,kfl_premix_chm
  use def_chemic, only      :  kfl_varZ_chm,kfl_varYc_chm
  use mod_ADR,    only      :  BDF

  implicit none
  integer(ip),  intent(in)  :: pnode,pgaus
  real(rp),     intent(in)  :: elcon(pnode,nclas_chm,*)
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: hleng(3)
  real(rp),     intent(in)  :: gptur(pgaus)                  ! Turbulence viscosity at gauss points
  real(rp),     intent(in)  :: gphco(mgaus)                  ! Conductivity at gauss points
  real(rp),     intent(in)  :: gpsph(mgaus)                  ! Specific heat at gauss points
  real(rp),     intent(in)  :: gpden(mgaus)                  ! Density at gauss points
  real(rp),     intent(out) :: res_Chi_Yc(pgaus)
  real(rp),     intent(out) :: sgs_Chi_Yc(pgaus)
  real(rp),     intent(out) :: res_Chi_z(pgaus)
  real(rp),     intent(out) :: sgs_Chi_z(pgaus)
  real(rp)                  :: gpcon(pgaus,nclas_chm,4)
  real(rp)                  :: seci4,delta,delta2,rdelta2
  real(rp)                  :: gpgve(ndime,ndime,pgaus)
  real(rp)                  :: rtau_turb                     ! Inverse of turbulent time scale
  real(rp)                  :: z_var,Yc_var
  real(rp)                  :: grad_Z(pgaus,ndime),mod_grad_Z(pgaus)
  real(rp)                  :: grad_Yc(pgaus,ndime),mod_grad_Yc(pgaus)
  real(rp)                  :: diffu                         ! diffusivity assuming Le=1 
  integer(ip)               :: igaus,iclas,inode,idime,jdime,itime

  !
  ! Initialization
  !
  res_Chi_z  = 0.0_rp
  res_Chi_Yc = 0.0_rp
  sgs_Chi_z  = 0.0_rp
  sgs_Chi_Yc = 0.0_rp

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
  ! Scalar dissipation rates of the reaction progress variable Yc and
  !   mixture fraction Z: X_Yc & X_Z
  ! 
  grad_Yc     = 0.0_rp
  mod_grad_Yc = 0.0_rp

  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           grad_Yc(igaus,idime) = grad_Yc(igaus,idime) + gpcar(idime,inode,igaus) * elcon(inode,1,1)
        end do
     end do
  end do

  if (kfl_premix_chm == 0_ip) then
     grad_Z     = 0.0_rp
     mod_grad_Z = 0.0_rp
     do igaus = 1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              grad_Z(igaus,idime) = grad_Z(igaus,idime) + gpcar(idime,inode,igaus) * elcon(inode,3,1)
           end do
        end do
     end do
  end if

  do igaus = 1,pgaus
     
     !
     ! Compute scalar dissipation rate (resolved) X_Yc = 2 * D * |grad Yc|^2
     !

     diffu = gphco(igaus) /max((gpsph(igaus) * gpden(igaus)),1e-10_rp)

     mod_grad_Yc(igaus) = dot_product(grad_Yc(igaus,:),grad_Yc(igaus,:))
     res_Chi_Yc(igaus) = 2.0_rp * diffu * mod_grad_Yc(igaus)


     !
     ! Compute scalar dissipation rate (resolved) X_Z = 2 * D * |grad Z|^2
     !
     if (kfl_premix_chm == 0_ip) then 
        mod_grad_Z(igaus) = dot_product(grad_Z(igaus,:),grad_Z(igaus,:))
        res_Chi_z(igaus) = 2.0_rp * diffu * mod_grad_Z(igaus)
     end if

     !
     ! Compute subgrid scalar dissipation rates: X_Yc^sgs & X_Z^sgs
     !

     if (kfl_varYc_chm /= 0_ip .or. kfl_varZ_chm /= 0_ip ) then
        
        delta   = ( hleng(1) * hleng(2) * hleng(3) )** 0.3333333_rp
        delta2  = delta*delta
        rdelta2 = 1.0_rp / delta2


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
        ! Compute X_Yc (subgrid) from transport Y_c^2 equation 
        !

        if (kfl_premix_chm == 0_ip .and. kfl_varYc_chm == 2_ip ) then
          Yc_var = gpcon(igaus,2,1) - gpcon(igaus,1,1) * gpcon(igaus,1,1)
          sgs_Chi_Yc(igaus) = 2.0_rp * rtau_turb * Yc_var

        !
        ! Compute X_Yc (subgrid) from transport Y_c variance equation
        !
        else if ( kfl_premix_chm == 1_ip .or. ( kfl_premix_chm == 0_ip .and. kfl_varYc_chm == 1_ip ) ) then
          sgs_Chi_Yc(igaus) = 2.0_rp * rtau_turb * gpcon(igaus,2,1)

        end if

        !
        ! Compute X_Z (subgrid) from transport Z variance equation           
        !
        if (kfl_varZ_chm == 1_ip) then
          sgs_Chi_z(igaus) = 2.0_rp * rtau_turb * gpcon(igaus,4,1)

        !
        ! Compute X_Z (subgrid) from transport Z^2 equation           
        !
        elseif (kfl_varZ_chm == 2_ip) then
          z_var = gpcon(igaus,4,1) - gpcon(igaus,3,1) * gpcon(igaus,3,1) 
          sgs_Chi_z(igaus) = 2.0_rp * rtau_turb * z_var

        endif

     end if

  end do

end subroutine chm_calc_scalar_dissip
