subroutine chm_upcpmu(itask) 
  !-----------------------------------------------------------------------
  ! Update
  ! NAME 
  !    chm_upcpmu
  ! DESCRIPTION
  !    Update each of the species specific heat and viscosity
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,npoin,mnode,coord
  use def_chemic, only      :  nspec_chm,entha_chm,tempe_chm,lawte_chm,temma_chm,kfl_warni_chm, &
                               sponge_chm, visco_factor_chm, visco_axis, visco_range,kfl_field_chm
  use def_master, only      :  conce,wmean,inotmaster,speci,visck,condk,sphek

  implicit none
  integer(ip),intent(in)    :: itask   ! 1 = update local variables, 2 = update also global variables
  integer(ip)               :: ipoin,idime,icoef,ispec
  real(rp)                  :: T,c(0:7),T_min_cp,T_max_cp
  integer(ip)               :: i_max_cp(1) ! Range that reaches max temperature
  real(rp)                  :: axis_projection,visco_factor

  if (INOTMASTER) then
     !
     ! We also compute enthalpies
     !
     ! h_k = enthalpy_of_formation_k - int(T_0 -> T) C_p,k dT
     !
     ! Cp is computed from the Shomate equation
     ! Cp = A + B*t + C*t^2 + D*t^3 + E/t^2 + X*t^0.5
     !           with t = T/1000
     !
     ! h_k = H0f(298.15K) + A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H
     !
     ! and cpcoe = (A, B, C, D, E, X, F, H)
     !
     ! Temperature below and above the defined ranges end up in constant Cp and h
     !
     ! Viscosity follows a Sutherland Law
     !
     !                    T0+110
     ! mu/mu_0=(T/T0)^1.5 ------ 
     !                    T +110
     !
     do ispec= 1_ip,nspec_chm
        T_min_cp = speci(ispec)%trang(1)/1000.0_rp
        i_max_cp = maxloc( speci(ispec)%trang ) 
        T_max_cp = speci(ispec)%trang(i_max_cp(1)) /1000.0_rp
        do ipoin = 1, npoin
           select case (lawte_chm)
           case(0)
              T=temma_chm(1)/1000.0_rp
           case(-2)
              T=tempe_chm(ipoin)/1000.0_rp  ! This 1000 is just from the Shomate eq.
           end select

           if ( T < T_min_cp) then ! Temp below minimum, assume constant
              T=T_min_cp
              icoef=1
           else if ( T > T_max_cp) then ! Temp above max range, assume constant
              T=T_max_cp
              icoef=i_max_cp(1)-1
           else
              ! Search the temperature range for the Cp coefficients
              icoef=1
              do while (  icoef < 6 .and. 1000.0_rp *T >= speci(ispec)%trang(icoef+1) ) !!We only have up to 4 ranges=5 temperatures
                 icoef=icoef+1
              enddo
           endif
           c=speci(ispec)%cpcoe(1:8,icoef)
           sphek(ipoin,ispec) = c(0) + c(1)*T + c(2)*T**2 + c(3)*T**3 + c(4)/T**2 !+ c(5)*sqrt(T) 
           !hCp0 = c(0)*Th + c(1)*Th**2/2.0_rp + c(2)*Th**3/3.0_rp + c(3)*Th**4/4.0_rp - c(4)/Th + (2.0_rp/3.0_rp)*c(5)*Th**1.5
           entha_chm(ipoin,ispec) = 1000.0_rp * &  ! For some reason Shomate eq. gives enthalpy in kJ so we transform to J
                (speci(ispec)%entha(1) + &
                c(0)*T + c(1)*T**2/2.0_rp + c(2)*T**3/3.0_rp + c(3)*T**4/4.0_rp - c(4)/T + c(6) - c(7) ) ! + (2.0_rp/3.0_rp)*c(5)*T**1.5
           if (speci(ispec) % lawvi .eq. 1_ip) then ! Sutherland Law
              visck(ipoin,ispec) = speci(ispec)%visco(1) * (1000.0_rp*T/speci(ispec)%visco(2))**1.5_rp &
                   * (speci(ispec)%visco(2)+110.0_rp) / (1000.0_rp*T+110.0_rp)
           else if (speci(ispec) % lawvi .eq. 2_ip) then ! Linear law
              visck(ipoin,ispec) = speci(ispec)%visco(1) + (1000.0_rp*T) * speci(ispec)%visco(2)
           endif
           ! Check if we have sponge layer
           if (sponge_chm == 1_ip) then
              axis_projection = 0.0_rp
              do idime = 1,ndime
                 axis_projection = axis_projection + coord(idime,ipoin) * visco_axis(idime)
              enddo
              if (axis_projection <= visco_range(1) ) then
                 visco_factor=1.0_rp
              else if  (axis_projection >= visco_range(2) ) then
                 visco_factor=visco_factor_chm
              else
                 visco_factor= 1.0_rp + (visco_factor_chm-1.0_rp) * (axis_projection-visco_range(1))/(visco_range(2)-visco_range(1))
              endif
              visck(ipoin,ispec)=  visck(ipoin,ispec)*visco_factor
           endif
           !
           ! Conductivity = mu * Cp / Prandtl
           !
           condk(ipoin,ispec) = visck(ipoin,ispec) *sphek(ipoin,ispec) / speci(ispec)%prand 
        enddo
     enddo
     
endif

end subroutine chm_upcpmu
