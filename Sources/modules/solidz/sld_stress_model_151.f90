!----------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_stress_model_151.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    November, 2015
!>          - Subroutine written
!> @author  Gerard Guillamet
!> @date    December, 2017
!>          - Adds strain measures corrections
!> @brief   Sant Venant - Kirchoff Orthotropic material model
!>
!> @details
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures\n
!>          E. J. Barbero. Introduction to Composite Materials Design\n
!>
!> @}
!----------------------------------------------------------------------------

subroutine sld_stress_model_151(pgaus,pmate,gpgdi,gpidg,gpdet,gpstr,ielem, &
     flagt,gpdds)

  use def_kintyp,                  only : ip,rp
  use def_domain,                  only : ndime
  use def_elmtyp
  use def_solidz,                  only : stiff0_sld
  use def_solidz,                  only : kfl_strai_sld, SLD_INFINITESIMAL, SLD_GREEN
  use mod_sld_stress_model_comput, only : sm_rotate_basis_creation
  use mod_sld_stress_model_comput, only : sm_rotate_voigt_second
  use mod_sld_stress_model_comput, only : sm_rotate_voigt_fourth
  use mod_sld_stress_model_comput, only : sm_strain_tensor
  use mod_sld_stress_model_comput, only : sm_stress_tensor
  use mod_sld_stress_model_comput, only : sm_tensor_to_voigt_fourth
  use mod_sld_stress_model_comput, only : sm_tensor_to_voigt_second
  use mod_sld_stress_model_comput, only : SM_stress_transport
  use mod_sld_stress_model_comput, only : SM_stiffness_transport
  use mod_sld_stress_model_comput, only : SM_PULLBACK

  implicit none

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  integer(ip), intent(in)  :: pgaus                                !< No. gauss points
  integer(ip), intent(in)  :: pmate                                !< Material code
  integer(ip), intent(in)  :: ielem                                !< Current element number
  integer(ip), intent(in)  :: flagt                                !< Integration scheme flag: 0: explicit / 1: implicit
  real(rp),    intent(in)  :: gpgdi(ndime,ndime,pgaus)             !< Displacement Gradient
  real(rp),    intent(in)  :: gpidg(ndime,ndime,pgaus)             !< Inverse of updated deformation gradient tensor
  real(rp),    intent(in)  :: gpdet(pgaus)                         !< Updated Jacobian
  real(rp),    intent(out) :: gpstr(ndime,ndime,pgaus)             !< 2nd Piola-Kirchoff stresses in tensor form
  real(rp),    intent(out) :: gpdds(ndime,ndime,ndime,ndime,pgaus) !< 2nd elasticity tensor in tensor form

  ! --------------------------------------------------------------------------------
  integer(ip)                             :: &
       igaus
  real(rp)                                :: &
       auxMA33(ndime,ndime),                 & ! Auxiliar 2nd Piola-Kirchoff stresses in tensor form
       gpgre(ndime,ndime,pgaus),             & ! Green-Lagrange strains in tensor form
       vogre(6),                             & ! Strains  in Voigt notation
       vostr(6),                             & ! Stresses in Voigt notation
       stiff(6,6),                           & ! Stiffness in Voigt form
       eltan(6,6),                           & ! Second elasticity tensor in Voig form
       traMa(3,3),                           & ! Transformation matrix
       auxV1(6),                             & ! Auxiliary stuff
       auxMA3333(3,3,3,3)

  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! INITIALIZE VARIABLES
  !
  ! Undamaged stiffness tensor
  !
  stiff(:,:) = stiff0_sld(:,:,pmate)
  !
  ! Initialise
  !
  gpstr = 0.0_rp
  gpdds = 0.0_rp
  ! --------------------------------------------------------------------------------
  ! LOOP OVER GAUSS POINTS
  !
  !...| Do gauss points |...........................................................
  do igaus = 1, pgaus
     !
     ! -----------------------------------------------------------------------------
     ! STRAIN TENSOR
     !
     ! Strain tensor
     !   0 - Infinitesimal tensor
     !   1 - Green-Lagrange tensor
     !
     if (kfl_strai_sld == SLD_INFINITESIMAL) then
        call SM_strain_tensor(0_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
     else if (kfl_strai_sld == SLD_GREEN) then
        call SM_strain_tensor(1_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
     end if
     !
     ! From tensor to Voigt notation
     !
     call SM_tensor_to_voigt_second(ndime, 0_ip, gpgre(:,:,igaus), vogre(:))
     !
     ! Rotate from the global to the material coordinate system
     !
     auxV1 = vogre
     call SM_rotate_basis_creation(ielem, gpgdi(:,:,igaus), traMa(:,:))
     call SM_rotate_voigt_second(1_ip, traMa(:,:), auxV1(:), vogre(:))
     !
     ! -----------------------------------------------------------------------------
     ! EFFECTIVE STRESS TENSOR
     !
     call SM_stress_tensor(0_ip, vogre(:), stiff(:,:), vostr(:))
     !
     ! Rotate from the local to the global coordinate system
     !
     auxV1 = vostr
     call SM_rotate_voigt_second(2_ip, traMa(:,:), vostr(:), auxV1(:))
     !
     ! From Voigt notation to tensorial notation
     !
     call SM_tensor_to_voigt_second(ndime, 1_ip, auxMA33(:,:), vostr(:))
     !
     ! Stress tensor due to strain measure
     !
     if (kfl_strai_sld == SLD_INFINITESIMAL) then
        !
        ! Assuming only Infinitesimal
        !
        call SM_stress_transport(SM_PULLBACK, ndime, gpdet(igaus), gpgdi(:,:,igaus), gpidg(:,:,igaus), &
             auxMA33(:,:), gpstr(:,:,igaus))

     else if (kfl_strai_sld == SLD_GREEN) then
        !
        ! Assuming full Green Lagrange
        !
        gpstr(:,:,igaus) = auxMA33(:,:)

     end if

     !
     ! ------------------------------------------------------------------------------
     ! SECOND ELASTICITY TENSOR (dS/dE = C^T)
     !
     ! [dS/dE] = [H]^-1 - [M]
     !
     ! ...| if implicit scheme |......................................................
     if (flagt .eq. 1_ip) then
        !
        ! Rotate from material to global coordinate system
        !
        call SM_rotate_voigt_fourth(2_ip, traMa(:,:), eltan(:,:), stiff(:,:))
        !
        ! From Voigt to tensor
        !
        call SM_tensor_to_voigt_fourth(ndime, 1_ip, auxMA3333(:,:,:,:), eltan(:,:))
        !
        ! Second elasticity tensor due to strain measure
        !
        if (kfl_strai_sld == SLD_INFINITESIMAL) then
           !
           ! Assuming only Infinitesimal
           !
           call SM_stiffness_transport(SM_PULLBACK, ndime, gpgdi(:,:,igaus),gpidg(:,:,igaus), &
                auxMA3333(:,:,:,:), gpdds(:,:,:,:,igaus))

        else if (kfl_strai_sld == SLD_GREEN) then
           !
           ! Assuming full Green Lagrange
           !
           gpdds(:,:,:,:,igaus) = auxMA3333(:,:,:,:)
        end if

     end if
     ! .....................................................| if implicit scheme |...
     !
  end do

end subroutine sld_stress_model_151

!----------------------------------------------------------------------------
!> @brief   Pre-calculus
!>
!> @details Get material properties and built stiffness matrix for material
!>          model 151
!----------------------------------------------------------------------------
subroutine sm151_precalculus(imate)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_parame, only                    :  &
       ip, rp
  use def_solidz, only                    :  &
       parco_sld, stiff0_sld
  ! --------------------------------------------------------------------------------
  implicit none
  ! --------------------------------------------------------------------------------
  integer(ip), intent(in) :: imate             !< Material code
  ! --------------------------------------------------------------------------------
  real(rp)                                :: &
       E11, E22, E33,                        & ! Material properties
       v12, v13, v23,                        &
       G12, G13, G23,                        &
       v21, v31, v32,                        & ! Not read, calculated
       delta, auxS1                            ! Auxiliary variables
  !
  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! MATERIAL PROPERTIES
  E11 = parco_sld(1,imate) ! Young moduli
  E22 = parco_sld(2,imate)
  E33 = parco_sld(3,imate)
  v12 = parco_sld(4,imate) ! Poison moduli
  v13 = parco_sld(5,imate)
  v23 = parco_sld(6,imate)
  G12 = parco_sld(7,imate) ! Shear moduli
  G13 = parco_sld(8,imate)
  G23 = parco_sld(9,imate)
  !
  v21 = (v12*E22)/E11
  v31 = (v13*E33)/E11
  v32 = (v23*E33)/E22
  !
  ! --------------------------------------------------------------------------------
  ! STIFFNESS TENSOR
  stiff0_sld(:,:,imate) = 0.0_rp
  delta = (1.0_rp - v12*v21 - v23*v32 - v31*v13 - 2.0_rp*v12*v23*v31)/(E11*E22*E33)
  !
  ! Stiff(1,:)
  !
  auxS1 = E22*E33*delta
  stiff0_sld(1,1,imate) = (1.0_rp - v23*v32)/auxS1
  stiff0_sld(1,2,imate) = (v21 + v31*v23)/auxS1
  stiff0_sld(1,3,imate) = (v31 + v21*v32)/auxS1
  !
  ! Stiff(2,:)
  !
  auxS1 = E33*E11*delta
  stiff0_sld(2,1,imate) = (v12 + v13*v32)/auxS1    ! It must be equal to stiff(1,2)
  stiff0_sld(2,2,imate) = (1.0_rp - v31*v13)/auxS1
  stiff0_sld(2,3,imate) = (v32 + v31*v12)/auxS1
  !
  ! Stiff(3,:)
  !
  auxS1 = E11*E22*delta
  stiff0_sld(3,1,imate) = (v13 + v12*v23)/auxS1    ! It must be equal to stiff(1,3)
  stiff0_sld(3,2,imate) = (v23 + v13*v21)/auxS1    ! It must be equal to stiff(2,3)
  stiff0_sld(3,3,imate) = (1.0_rp - v12*v21)/auxS1
  !
  ! Stiff(4,:), stiff(5,:) and stiff(6,:)
  !
  stiff0_sld(4,4,imate) = G23
  stiff0_sld(5,5,imate) = G13
  stiff0_sld(6,6,imate) = G12
  !
  ! ============================================================|    MAIN     |=====

end subroutine sm151_precalculus

