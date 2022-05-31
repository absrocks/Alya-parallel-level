!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_ctcalc.f90
!> @author  Matias
!> @brief   Calculates thrust coeff for actuator disk in terms of hub velocity
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_ctcalc(                     &
    xforc_material_nsi, pmate)

  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi
  use def_master, only     :  kfl_paral
  use def_nastin, only     :  velta_nsi, thrta_nsi, powta_nsi
  implicit none
  
  real(rp), intent(inout)  :: xforc_material_nsi(*)
  integer(ip), intent(in)  :: pmate
  integer(ip)              :: iters, index
  real(rp)                 :: thrco, toler, oldve, powco, uturb(3)
  real(rp)                 :: veinf, uhub1, u_ref(3), n_disk(3), afact, muhub
  real(rp)                 :: uinfi, uinfj, ct1, error, zefun, uhub(3)
  integer(ip)              :: method=2
  character(20)            :: messa
  character(100)           :: wmess
 
  iters = 0


  ! convergence tolerance between infinite velocity and Ct
  toler = 1.0e-7_rp
  !  disk normal 
  n_disk(1:3) = xforc_material_nsi( 3:5)
  ! choose between uhub at center, or u averaged
  !uturb(1:3) = xforc_material_nsi(  9:11) ! velocity at hub center.
  uturb(1:3) = xforc_material_nsi( 19:21)  ! averaged velocity in disk volume
  !  normal hub velocity u · n
  uhub1   = abs(dot_product(uturb(1:3), n_disk(1:3)))
  !
  !  axial induction factor from last Ct 
  !  using theoretical relation

  ! bisection method
     error = 10.0_rp
  !   U infty limits
     !lower limit
     uinfi = uhub1 ! a = 0.0 
     ! upper limit
     uinfj = 2.2_rp*uhub1 
  
     ! test before entering that the function changes sign
     veinf = uinfi
     ! manufacturer thrust curve Ct(veinf)
     call nsi_thrpow(veinf, thrco, powco, pmate)
     afact = 1.0_rp - uhub1/veinf 
     ct1   = 4.0_rp*afact*(1.0_rp-0.25_rp*afact*(5.0_rp - 3.0_rp*afact)) ! glauert
     zefun = thrco - ct1
     
     veinf = uinfj
     call nsi_thrpow(veinf, thrco, powco, pmate)
     
     afact = 1.0_rp - uhub1/veinf
     ct1   = 4.0_rp*afact*(1.0_rp-0.25_rp*afact*(5.0_rp - 3.0_rp*afact)) ! glauert
     if (zefun*( thrco - ct1)>1.0e-6_rp) then
        ! extend upper limit
        uinfj = uinfj*2.0_rp
        veinf = uinfj
        call nsi_thrpow(veinf, thrco, powco, pmate)
        
        afact = 1.0_rp - uhub1/veinf
        ct1   = 4.0_rp*afact*(1.0_rp-0.25_rp*afact*(5.0_rp - 3.0_rp*afact)) ! glauert
        if (zefun*( thrco - ct1)>1.0e-6_rp) then
           ! measurement of disalignment  (u.n)                                              
           if (uhub1 > 1.5_rp) then
              if (kfl_paral==1.or.kfl_paral==-1)  then
                 muhub = sqrt(dot_product( uturb(1:3),uturb(1:3)))
                 print *, 'hub velocity modulus:', muhub
                 print *, 'yawing angle (degrees) =', 180.0_rp/pi*acos(uhub1 /muhub)
                 write(messa, '(f10.3, i3)') uhub1, pmate-1
                 wmess = 'nsi_ctcalc: could not find Uinf,bad guess. Uhub, turbine numb= '//trim(messa)
                 call runend(wmess)
              end if
           else 
              write(messa, '(f10.3, i3)') uhub1, pmate-1
              wmess = 'nsi_ctcalc: could not find Uinf,bad guess. Uhub, turbine numb= '//trim(messa)
              ! assumes veinf = uhub ! a=0
              veinf =uhub1
              call nsi_thrpow(veinf, thrco, powco, pmate)            
              
              error=toler*0.1_rp 
           end if
        end if
     end if
     
     iters = 0
     do while (abs(error)>toler.and.iters<30)      
        veinf= 0.5_rp*(uinfi+uinfj)
        afact = 1.0_rp - uhub1/veinf
        ct1 = 4.0_rp*afact*(1.0_rp - 0.25_rp*afact*(5.0_rp - 3.0_rp*afact)) ! glauert
        !manufacturers
        call nsi_thrpow(veinf, thrco, powco, pmate)
        zefun = thrco - ct1
        if (zefun < 0.0_rp) then
           uinfj = veinf                     
        else
           uinfi = veinf            
        end if
        error = (ct1 -thrco) / thrco
        iters = iters +1 
     end do
     if (iters.gt.30) print *, 'warning nsi_ctcal:thrust coeff did not converge'
   
!!$     To load veave table
!!$     index = nint(veinf) -2
!!$     if (index.lt.1.or.index.gt.25) call runend('nsi_ctcalc:bad coefficient')
!!$     veinf= velta_nsi(index,2)
!!$     thrco= thrta_nsi(index,2) 
!!$     powco= powta_nsi(index,2)
 

  
!!weir east
   !!thrco = 0.63_rp
   !!powco = 0.44_rp
   !!veinf = 10.9_rp
   !!weir west
   !!veinf = 10.7
   
   !!low TI
!   thrco=0.79_rp
!   powco=0.47_rp
!   veinf=8.0_rp
 
  !!1T Sexbierum
!   thrco=0.75
!   powco=0.39
!   veinf=8.0_rp

!  call nsi_thrust(veinf, thrco, powco)                                                                                                                                                                                                   
  xforc_material_nsi(16) = thrco
  xforc_material_nsi(17) = powco
  xforc_material_nsi(18) = veinf



end subroutine nsi_ctcalc


  !-------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_ctcalc.f90
  !> @author  Matias
  !> @brief   Interpolates thrust and power coeffs for actuator disk in terms of veinf.
  !> @details 
  !> @} 
  !-------------------------------------
subroutine nsi_thrust(veinf, thrco, powco)  

  use def_kintyp, only     :  ip,rp

  implicit none
  real(rp), intent(in)   ::  veinf
  real(rp), intent(out)  ::  powco, thrco
  integer(ip)            ::  ntabl, iz, jz, kz
  real(rp)               ::  velin(17), thrta(17), powta(17), facto

  !
  !    loads table 
  !
  ntabl = 17 ! number of values
  
  !
  !    velocity values are supposed to be given in increasing order
  !
  
  velin(1)=5.0_rp
  velin(2)=6.0_rp
  velin(3)=7.0_rp
  velin(4)=8.0_rp
  velin(5)=9.0_rp
  velin(6)=10.0_rp
  velin(7)=11.0_rp
  velin(8)=12.0_rp
  velin(9)=12.4_rp
  velin(10)=13.0_rp
  velin(11)=14.0_rp
  velin(12)=15.0_rp
  velin(13)=16.0_rp
  velin(14)=17.0_rp
  velin(15)=18.0_rp
  velin(16)=19.0_rp
  velin(17)=20.0_rp

  thrta(1) = 0.876_rp
  thrta(2) = 0.810_rp
  thrta(3) = 0.750_rp
  thrta(4) = 0.750_rp
  thrta(5) = 0.750_rp
  thrta(6) = 0.750_rp
  thrta(7) = 0.700_rp
  thrta(8) = 0.660_rp
  thrta(9) = 0.644_rp
  thrta(10) = 0.500_rp
  thrta(11) = 0.370_rp
  thrta(12) = 0.300_rp
  thrta(13) = 0.240_rp
  thrta(14) = 0.200_rp
  thrta(15) = 0.170_rp  
  thrta(16) = 0.140_rp
  thrta(17) = 0.120_rp

  powta(1)  =0.360_rp
  powta(2)  =0.390_rp
  powta(3)  =0.390_rp
  powta(4)  =0.390_rp
  powta(5)  =0.390_rp
  powta(6)  =0.390_rp
  powta(7)  =0.380_rp
  powta(8)  =0.370_rp
  powta(9)  =0.366_rp
  powta(10) =0.323_rp
  powta(11) =0.257_rp
  powta(12) =0.210_rp
  powta(13) =0.172_rp
  powta(14) =0.144_rp
  powta(15) =0.121_rp
  powta(16) =0.102_rp
  powta(17) =0.087_rp

  if (veinf.lt.velin(1)) then  
     thrco = thrta(1)
     powco = powta(1)
  else if (veinf.gt.velin(ntabl)) then
     thrco =thrta(ntabl)
     powco =powta(ntabl)
  else !linear interpolation
     iz = 1
     jz = ntabl
     kz = ntabl/2            
     do while ((jz-iz).gt.1)                 
        if (veinf.lt.velin(kz)) then
           jz = kz                  
        else
           iz = kz
        end if
        kz = (iz+jz)/2
     end do

     facto = (veinf-velin(iz))/(velin(jz)-velin(iz))     
     thrco = thrta(iz) + facto*(thrta(jz)-thrta(iz))
     powco  =  powta(iz) + facto*(powta(jz)-powta(iz))
  end if


end subroutine nsi_thrust
!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_ctcalc.f90
!> @author  Matias
!> @brief   Calculates thrust coeff for actuator disk in terms of hub velocity
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_ctcalc_vmed(                     &
    xforc_material_nsi, pmate)

  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi
  use def_master, only     :  kfl_paral
  use def_nastin, only     :  velta_nsi, & ! velocity table
       thrta_nsi, &    ! thrust table
       powta_nsi   ! power table
  implicit none
  
  real(rp), intent(inout)  :: xforc_material_nsi(*)
  integer(ip), intent(in)  :: pmate
  integer(ip)              :: iters, index
  real(rp)                 :: thrco, toler, oldve, powco, uturb(3)
  real(rp)                 :: veinf, veave, u_ref(3), n_disk(3), afact, muhub
  real(rp)                 :: uinfi, uinfj, ct1, error, zefun, uhub(3)
  integer(ip)              :: method=2
  character(20)            :: messa
  character(100)           :: wmess
 
  iters = 0


  ! convergence tolerance between infinite velocity and Ct
  toler = 1.0e-7_rp
  !  disk normal 
  n_disk(1:3) = xforc_material_nsi( 3:5)
  ! choose between uhub at center, or u averaged
  !uturb(1:3) = xforc_material_nsi(  9:11) ! velocity at hub center.
  uturb(1:3) = xforc_material_nsi( 19:21)  ! averaged velocity in disk volume
  !  normal hub velocity u · n
  veave   = abs(dot_product(uturb(1:3), n_disk(1:3)))
  !
  !  axial induction factor from last Ct 
  !  using theoretical relation

  ! thrust curve Ct(veave) returns thrco powco and veinf in terms of veave
  call nsi_thrpow_veave(veave, veinf,thrco, powco, pmate)


!  call nsi_thrust(veinf, thrco, powco)                                                                                                                                                                                                   
  xforc_material_nsi(16) = thrco
  xforc_material_nsi(17) = powco
  xforc_material_nsi(18) = veinf



end subroutine nsi_ctcalc_vmed
