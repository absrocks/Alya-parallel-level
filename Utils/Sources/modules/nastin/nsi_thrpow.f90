!-------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_ctcalc.f90
!> @author  Matias
!> @brief   Interpolates thrust and power coeffs for actuator disk in terms of veinf.
!> @details 
!> @} 
!-------------------------------------
subroutine nsi_thrpow(veinf, thrco, powco, pmate)


  use def_kintyp, only     : ip,rp
  use def_nastin, only     : velta_nsi, & ! velocity table
       thrta_nsi, &    ! thrust table
       powta_nsi, &    ! power table
       ntabl_nsi       ! # data table

  implicit none
  real(rp), intent(in)   ::  veinf
  integer(ip), intent(in)::  pmate
  real(rp), intent(out)  ::  powco, thrco
  integer(ip)            ::  iz, jz, kz
  real(rp)               ::  facto

  !
  !  for interpolation  velocity values are supposed to be given in increasing order
  !
  if (veinf.lt.velta_nsi(1, pmate)) then  ! cuts production     
     thrco = thrta_nsi(1, pmate)      !0.0_rp  !
     powco = powta_nsi(1, pmate)      !0.0_rp  !   
  else if (veinf.gt.velta_nsi(ntabl_nsi(pmate), pmate)) then !cuts production
     thrco = thrta_nsi(ntabl_nsi(pmate), pmate)  !0.0_rp  !
     powco = 0.0_rp !powta_nsi(ntabl_nsi(pmate), pmate)  !0.0_rp  !
  else !linear interpolation
     iz = 1
     jz = ntabl_nsi(pmate)
     kz = ntabl_nsi(pmate)/2            
     do while ((jz-iz).gt.1)                 
        if (veinf.lt.velta_nsi(kz, pmate)) then
           jz = kz                  
        else
           iz = kz
        end if
        kz = (iz+jz)/2

     end do

     facto = (veinf-velta_nsi(iz, pmate))/(velta_nsi(jz, pmate)-velta_nsi(iz, pmate))          
     thrco = thrta_nsi(iz, pmate) + facto*(thrta_nsi(jz, pmate)-thrta_nsi(iz, pmate))
     powco = powta_nsi(iz, pmate) + facto*(powta_nsi(jz, pmate)-powta_nsi(iz, pmate))
  end if


end subroutine nsi_thrpow ! normal and tangential forces
!-------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_ctcalc.f90
!> @author  Matias
!> @brief   Interpolates thrust and power coeffs for actuator disk in terms of veinf.
!> @details 
!> @} 
!-------------------------------------
subroutine nsi_forcent(radiu, norma, tange, pmate)


  use def_kintyp, only     : ip,rp
  use def_nastin, only     : radiu_nsi, & ! velocity table
       forcn_nsi, &    ! normal force table
       forct_nsi, &    ! rotational force table
       ntabr_nsi       ! # data table

  implicit none
  real(rp), intent(in)   ::  radiu
  integer(ip), intent(in)::  pmate
  real(rp), intent(out)  ::  norma, tange
  integer(ip)            ::  iz, jz, kz
  real(rp)               ::  facto

  !
  !  for interpolation  velocity values are supposed to be given in increasing order
  !
  if (radiu.lt.radiu_nsi(1, pmate)) then !linear interp
     norma = forcn_nsi(1, pmate)*radiu/radiu_nsi(1, pmate)      !0.0_rp  !
     tange = forct_nsi(1, pmate)*radiu/radiu_nsi(1, pmate)    !0.0_rp  !
  else if (radiu.gt.radiu_nsi(ntabr_nsi(pmate), pmate)) then
     norma = forcn_nsi(ntabr_nsi(pmate), pmate)  !0.0_rp  !
     tange = forct_nsi(ntabr_nsi(pmate), pmate)  !0.0_rp  !
  else !linear interpolation
     iz = 1
     jz = ntabr_nsi(pmate)
     kz = ntabr_nsi(pmate)/2            
     do while ((jz-iz).gt.1)                 
        if (radiu.lt.radiu_nsi(kz, pmate)) then
           jz = kz                  
        else
           iz = kz
        end if
        kz = (iz+jz)/2

     end do

     facto = (radiu-radiu_nsi(iz, pmate))/(radiu_nsi(jz, pmate)-radiu_nsi(iz, pmate))          
     norma = forcn_nsi(iz, pmate) + facto*(forcn_nsi(jz, pmate)-forcn_nsi(iz, pmate))
     tange = forct_nsi(iz, pmate) + facto*(forct_nsi(jz, pmate)-forct_nsi(iz, pmate))
  end if


end subroutine nsi_forcent

subroutine nsi_thrpow_veave(veave, veinf, thrco, powco, pmate)


  use def_kintyp, only     : ip,rp
  use def_nastin, only     : velta_nsi, & ! velocity table
       thrta_nsi, &    ! thrust table
       powta_nsi, &    ! power table
       veave_nsi, &    ! averaged velocity
       ntabl_nsi       ! # data table

  implicit none
  real(rp), intent(in)   ::  veave
  integer(ip), intent(in)::  pmate
  real(rp), intent(out)  ::  powco, thrco, veinf
  integer(ip)            ::  iz, jz, kz
  real(rp)               ::  facto

  !
  !  for interpolation  velocity values are supposed to be given in increasing order
  !
  if (veave.lt.veave_nsi(1, pmate)) then  ! cuts production     
     thrco = thrta_nsi(1, pmate)      !0.0_rp  !
     powco = powta_nsi(1, pmate)      !0.0_rp  !
     veinf = veave * velta_nsi(1,pmate)/ veave_nsi(1, pmate) 
  else if (veave.gt.veave_nsi(ntabl_nsi(pmate), pmate)) then !cuts production
     thrco = thrta_nsi(ntabl_nsi(pmate), pmate)  !0.0_rp  !
     powco = 0.0_rp  !
     veinf = velta_nsi(ntabl_nsi(pmate), pmate)
  else !linear interpolation
     iz = 1
     jz = ntabl_nsi(pmate)
     kz = ntabl_nsi(pmate)/2            
     do while ((jz-iz).gt.1)                 
        if (veave.lt.veave_nsi(kz, pmate)) then
           jz = kz                  
        else
           iz = kz
        end if
        kz = (iz+jz)/2
     end do

     facto = (veave-veave_nsi(iz, pmate))/(veave_nsi(jz, pmate)-veave_nsi(iz, pmate))          
     thrco = thrta_nsi(iz, pmate) + facto*(thrta_nsi(jz, pmate)-thrta_nsi(iz, pmate))
     powco = powta_nsi(iz, pmate) + facto*(powta_nsi(jz, pmate)-powta_nsi(iz, pmate))
     veinf = velta_nsi(iz, pmate) + facto*(velta_nsi(jz, pmate)-velta_nsi(iz, pmate))
  end if


end subroutine nsi_thrpow_veave ! normal and tangential forces
