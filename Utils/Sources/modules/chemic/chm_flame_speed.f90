subroutine chm_flame_speed() 
  !-----------------------------------------------------------------------
  ! ****f* partis/chm_flame_speed
  ! NAME 
  !    chm_flame_speed
  ! DESCRIPTION
  !    Compute the laminar flame speed and flame thickness for a given fuel
  !    based on local properties.
  !    The correlation is based on Gottgens et al (1992):
  !     -->  Gottgens J, Mauss F and Peters N. Analytic approximations of burning velocities
  !          and flame thickness of lean hydrogen, methane, ethylene, ethane, acetylene
  !          and propane flames. Twentyfourth Symp. (Int.) on Combustion, Pittsburgh,
  !          120135 (1992).
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_chemic
  use def_master 
  use mod_ker_proper 
  use def_kermod

  implicit none

  integer(ip)            :: ipoin,ispec,fuel_ispec,air_o2,air_n2,dummi
  real(rp)               :: adiab_temp,rst_air_fuel,unburnt_temp,unburnt_dens,a_const,   &
                            d_const,c_const,b_const,x_const,y_const,z_const,densl,temsl, &
                            presl,consl(nspec_chm),phisl,epsil,temp_0,gamma_0,dummr(1)
  real(rp)               :: vissl(npoin)

  epsil  = 0.0000000001_rp

  if (INOTMASTER) then

     !! METHANE PROPERTIES
        adiab_temp    = 2242.000_rp   ! 1 atm & 298.15 K
        rst_air_fuel  =   17.160_rp   ! For methane 
        unburnt_temp  =  298.000_rp   ! Unburnt temperature (fuel injection temperature)
        unburnt_dens  =    1.104_rp   ! Unburnt density (fuel injection density)      

        ! GOTTGENS PROPERTIES
        a_const =    22.176_rp ! [cm/s]
        d_const =  3.1557e8_rp ! [bar]
        c_const =   23873.0_rp ! [K]
        b_const =   6444.27_rp ! [K]
        x_const =  0.565175_rp
        y_const =    2.5158_rp
        z_const =     0.994_rp
     !!
  
     !! DMM To be removed
        fuel_ispec = 1
        air_o2     = 2
        air_n2     = 5
     !! DMM

     flspe_chm = 0.0_rp
     flthi_chm = 0.0_rp

     call ker_proper('VISCO','NPOIN',dummi,dummi,vissl,dummi)

     do ipoin = 1,npoin

        densl = densi(ipoin,1)
        temsl = tempe(ipoin,1)
        presl = press(ipoin,1) / 101325.0_rp ! Conversion to bar
        consl = conce(ipoin,:,1)

        !
        ! Equivalence ratio phisl
        !
   
        do ispec = 1,nspec_chm

!           if ( speci(ispec)%name == 'CH4  ' ) then ! Fuel
!              fuel_ispec = ispec
!              write(11111,*)'FUEL',ispec
!           else if ( speci(ispec)%name == 'O2   ' ) then ! Oxygen
!              air_o2     = ispec
!              write(11111,*)'O2  ',ispec
!           else if ( speci(ispec)%name == 'N2   ' ) then ! Nitrogen
!              air_n2     = ispec
!              write(11111,*)'N2  ',ispec
!           endif
        enddo
 
        phisl = rst_air_fuel * consl(fuel_ispec) / (consl(air_o2) + consl(air_n2) + epsil)
        if (phisl < 0.0_rp ) then
           phisl = 0.0_rp    ! Correction to avoid negative equivalence ratio
        endif

        !
        ! Local reference temperature To
        !
        temp_0  = - c_const / log(presl/d_const)  
       

        !
        ! Gamma function evaluated at To
        !
        if (temp_0 > 0.0_rp ) gamma_0 = 0.0000258_rp * (temp_0 / 298.0_rp)**0.7_rp  

        !
        ! Laminar flame speed in [m/s]
        !
         flspe_chm(ipoin) = 0.0_rp
!        flspe_chm(ipoin) = 0.01_rp * a_const * (phisl / (rst_air_fuel + phisl) )**x_const  &
!                           * exp(b_const / temp_0) * (unburnt_temp / temp_0)     &
!                           * ((adiab_temp - temp_0)/(adiab_temp - unburnt_temp))**y_const
        !
        ! Laminar flame thickness in [m]
        !
        if ( (flspe_chm(ipoin) - 0.0_rp) < 0.0001_rp ) then
           flthi_chm(ipoin) = 0.0_rp       
        else
           !! flthi_chm(ipoin) = 0.001_rp * z_const * gamma_0 * (adiab_temp - unburnt_temp) &       !! Proposed by Gottgens et al (1992)
           !!                   / (unburnt_dens * flspe_chm(ipoin) * (temp_0 - unburnt_temp))
           flthi_chm(ipoin) = 0.001_rp * vissl(ipoin) / (0.7_rp * unburnt_dens * flspe_chm(ipoin))  !! Proposed by H.M. Heravi et al (2007)
        endif
     
     enddo

  endif
 
end subroutine chm_flame_speed
