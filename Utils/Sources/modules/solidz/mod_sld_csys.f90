!------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_csys.f90
!> @author  Gerard Guillamet
!> @date    June, 2018
!>          - Subroutine written
!> @brief   ToolBox for coordinate system transformations and material axes
!>
!> @details ToolBox for coordinate system transformations and material axes
!>
!>          \verbatim
!>          This toolBox includes useful transformations using coordinate
!>          systems and material axes or fibers:
!>           - Basic rotations
!>           - Rotation matrix from axis and angle
!>           - Vector rotations
!>           - Vector transformation between Cylindrical and Cartesian basis
!>           - Tensor transformation between Cylindrical and Cartesian basis
!>           - Assign material system or axis (fibers) using different methods
!>           - Element normal based on midsurface method (HEX08 and PEN06)
!>           - Rotation vector of nodal variables due to local axes
!>           - Rotation matrix of nodal variables due to local axes
!>
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------------

module mod_sld_csys

  use def_kintyp, only : ip, rp, lg
  use mod_maths,  only : maths_vectorial_product
  use mod_maths,  only : maths_normalize_vector

  implicit none

  public  :: sld_csys_rotation_matrix
  public  :: sld_csys_rotation_matrix_from_axis_and_angle
  public  :: sld_csys_vector_rotation
  public  :: sld_csys_vector_cylindrical_to_cartesian
  public  :: sld_csys_tensor_cylindrical_to_cartesian
  public  :: sld_csys_assign_material_axes
  public  :: sld_csys_build_jacrot
  public  :: sld_csys_rotuni
  public  :: sld_csys_rotmat
  public  :: sld_csys_midsurface_element_normal

  private :: sld_midsurface_basis_rotation_element_normal
  private :: sld_midsurface_coordinates
  private :: sld_midsurface_shape_functions
  private :: sld_midsurface_gp_vectors
  private :: sld_midsurface_primary_axis

contains

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Basic rotations
  !> @details Basic rotations
  !----------------------------------------------------------------------------

  subroutine sld_csys_rotation_matrix(itask, pdime, theta, rotmat)

    use def_parame, only       : pi

    implicit none

    integer(ip), intent(in)   :: itask               !< Rotation axes
    integer(ip), intent(in)   :: pdime               !< Number of dimension
    real(rp),    intent(in)   :: theta               !< Rotation angle
    real(rp),    intent(out)  :: rotmat(pdime,pdime) !< Rotation matrix

    select case (itask)

    case(1_ip)
       !
       ! Rx(theta)
       !
       rotmat(:,:) =  0.0_rp
       if (pdime == 2_ip) then
          call runend('MOD_SLD_CSYS: ROTATION MATRIX NOT POSSIBLE')
       else if (pdime == 3_ip) then
          rotmat(1,1) =  1.0_rp
          rotmat(2,2) =  cos(theta*pi/180.0_rp)
          rotmat(2,3) = -sin(theta*pi/180.0_rp)
          rotmat(3,2) =  sin(theta*pi/180.0_rp)
          rotmat(3,3) =  cos(theta*pi/180.0_rp)
       end if

    case(2_ip)
       !
       ! Ry(theta)
       !
       rotmat(:,:) =  0.0_rp
       if (pdime == 2_ip) then
          call runend('MOD_SLD_CSYS: ROTATION MATRIX NOT POSSIBLE')
       else if (pdime == 3_ip) then
          rotmat(1,1) =  cos(theta*pi/180.0_rp)
          rotmat(1,3) =  sin(theta*pi/180.0_rp)
          rotmat(2,2) =  1.0_rp
          rotmat(3,1) = -sin(theta*pi/180.0_rp)
          rotmat(3,3) =  cos(theta*pi/180.0_rp)
       end if

    case(3_ip)
       !
       ! Rz(theta)
       !
       rotmat(:,:) =  0.0_rp
       rotmat(1,1) =  cos(theta*pi/180.0_rp)
       rotmat(1,2) = -sin(theta*pi/180.0_rp)
       rotmat(2,1) =  sin(theta*pi/180.0_rp)
       rotmat(2,2) =  cos(theta*pi/180.0_rp)
       if (pdime == 3_ip) then
          rotmat(3,3) =  1.0_rp
       end if

    case default

       call runend('MOD_SLD_CSYS: ROTATION MATRIX NOT DEFINED')

    end select

    return

  end subroutine sld_csys_rotation_matrix

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Rotation matrix from axis and angle
  !> @details Rotation matrix from axis and angle
  !----------------------------------------------------------------------------

  subroutine sld_csys_rotation_matrix_from_axis_and_angle( u, theta, rotmat)

    use def_parame, only       : pi

    implicit none

    real(rp),    intent(in)   :: theta       !< Rotation angle
    real(rp),    intent(in)   :: u(3)        !< Unit vector or rotation axis
    real(rp),    intent(out)  :: rotmat(3,3) !< Rotation matrix

    real(rp)                  :: c, s

    c = cos(theta*pi/180.0_rp)
    s = sin(theta*pi/180.0_rp)

    rotmat(:,:) = 0.0_rp
    rotmat(1,1) = c + (1.0_rp - c)*u(1)**2
    rotmat(1,2) = u(1)*u(2)*(1.0_rp - c) - u(3)*s
    rotmat(1,3) = u(1)*u(3)*(1.0_rp - c) + u(2)*s
    rotmat(2,1) = u(2)*u(1)*(1.0_rp - c) + u(3)*s
    rotmat(2,2) = c + (1.0_rp - c)*u(2)**2
    rotmat(2,3) = u(2)*u(3)*(1.0_rp - c) - u(1)*s
    rotmat(3,1) = u(3)*u(1)*(1.0_rp - c) - u(2)*s
    rotmat(3,2) = u(3)*u(2)*(1.0_rp - c) + u(1)*s
    rotmat(3,3) = c + (1.0_rp - c)*u(3)**2

  end subroutine sld_csys_rotation_matrix_from_axis_and_angle

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Rotate a vector
  !> @details Rotate a vector from xy-cartesian to an x'y'-cartesian CSYS and
  !>          viceversa
  !----------------------------------------------------------------------------

  subroutine sld_csys_vector_rotation(itask, pdime, rotmat, vecold, vecnew)

    implicit none

    integer(ip), intent(in)     :: itask               !< Mapping direction
    integer(ip), intent(in)     :: pdime               !< Number of dimensions
    real(rp),    intent(in)     :: rotmat(pdime,pdime) !< Rotation matrix
    real(rp),    intent(inout)  :: vecold(pdime)       !< Vector x - y - z
    real(rp),    intent(inout)  :: vecnew(pdime)       !< Vector x'- y'- z'

    integer(ip)                 :: idime, jdime

    select case (itask)

    case(1_ip)
       !
       ! Vector rotation from x-y-z to x'-y'-z'
       !
       vecnew = 0.0_rp
       do idime = 1, pdime
          do jdime = 1, pdime
             vecnew(idime) = vecnew(idime) + rotmat(idime,jdime)*vecold(jdime)
          end do
       end do

    case(2_ip)
       !
       ! Vector rotation x'-y'-z' to x-y-z
       !
       vecold = 0.0_rp
       do idime = 1, pdime
          do jdime = 1, pdime
             vecold(idime) = vecold(idime) + rotmat(jdime,idime)*vecnew(jdime)
          end do
       end do

    case default
       return

    end select

    return

  end subroutine sld_csys_vector_rotation

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Vector cylindrical to cartesian an viceversa
  !> @details
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_vector_cylindrical_to_cartesian(itask, pdime, veccar, veccyl)

    use def_parame, only         : pi

    implicit none

    integer(ip), intent(in)     :: itask         !< Transformation
    integer(ip), intent(in)     :: pdime         !< Number of dimensions
    real(rp),    intent(inout)  :: veccar(pdime) !< Vector in Cartesian coordinates
    real(rp),    intent(inout)  :: veccyl(pdime) !< Vector in Cylindrical coordinates

    real(rp)                    :: r, theta
    real(rp)                    :: x, y, z

    select case (itask)

    case(1_ip)
       !
       ! Cylindrical to Cartesian relationship
       !
       r     = veccyl(1)
       theta = veccyl(2)
       if (pdime == 3_ip) z = veccyl(3)
       !
       veccar(1) = r*cos(theta*pi/180_rp)
       veccar(2) = r*sin(theta*pi/180_rp)
       if (pdime == 3_ip)  veccyl(3) = veccar(3)

    case(2_ip)
       !
       ! Cartesian to Cylindrical relationship
       !
       x = veccar(1)
       y = veccar(2)
       if (pdime == 3_ip) z = veccar(3)
       veccyl(1) = sqrt(x**2 + y**2)
       veccyl(2) = atan2(y,x)
       if (pdime == 3_ip)  veccyl(3) = z

    case default
       return

    end select

    return

  end subroutine sld_csys_vector_cylindrical_to_cartesian

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Tensor transformation cylindrical to cartesian an viceversa
  !> @details
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_tensor_cylindrical_to_cartesian(itask, pdime, pocoo, csysc, tenca, tency)

    use def_parame, only : pi
    use def_solidz, only : csysl_sld

    implicit none

    integer(ip), intent(in)    :: itask              !< Transformation
    integer(ip), intent(in)    :: pdime              !< Number of dimensions
    real(rp),    intent(in)    :: pocoo(pdime)       !< Point coordinates
    real(rp),    intent(in)    :: csysc(pdime)       !< Center of coordinate system
    real(rp),    intent(inout) :: tenca(pdime,pdime) !< Tensor in Cartesian coordinates
    real(rp),    intent(inout) :: tency(pdime,pdime) !< Tensor in Cylindrical coordinates

    real(rp)                   :: vectr(pdime)       !< Vector r
    real(rp)                   :: tenau(pdime,pdime) !< Tensor auxiliar
    real(rp)                   :: rotma(pdime,pdime) !< Rotation matrix
    real(rp)                   :: rotmt(pdime,pdime) !< Transpose of rotation matrix
    real(rp)                   :: theta              !< Theta

    select case (itask)

    case(1_ip)
       !
       ! Cylindrical to Cartesian relationship
       !
       call runend('MOD_SLD_CSYS: CYLINDRICAL TO CARTESIAN TENSOR NOT IMPLMENTED')

    case(2_ip)
       !
       ! Cartesian to Cylindrical relationship
       !
       vectr(:) = pocoo(:) - csysc(:)
       if ( pdime == 2_ip ) then
             theta = atan(vectr(2)/vectr(1))*180.0_rp/pi
             call sld_csys_rotation_matrix(3_ip, pdime, theta, rotma)
       else if ( pdime == 3_ip )then
          if ( csysl_sld(3*pdime) /= 0_ip ) then
             !
             ! X-Y Plane (Axial=Z)
             theta = atan(vectr(2)/vectr(1))*180.0_rp/pi
             call sld_csys_rotation_matrix(3_ip, pdime, theta, rotma)
          else if ( csysl_sld(3*pdime-1) /= 0_ip  ) then
             !
             ! X-Z Plane (Axial=Y)
             theta = atan(vectr(3)/vectr(1))*180.0_rp/pi
             call sld_csys_rotation_matrix(2_ip, pdime, theta, rotma)
          else if ( csysl_sld(3*pdime-2) /= 0_ip ) then
             !
             ! Y-Z Plane (Axial=X)
             theta = atan(vectr(3)/vectr(2))*180.0_rp/pi
             call sld_csys_rotation_matrix(1_ip, pdime, theta, rotma)
          end if
       end if
       rotmt(:,:) = transpose(rotma(:,:))
       tenau(:,:) = matmul(rotma(:,:),tenca(:,:))
       tency(:,:) = matmul(tenau(:,:),rotmt(:,:))

    case default
       return

    end select

    return

  end subroutine sld_csys_tensor_cylindrical_to_cartesian

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Assign material coordinate system using different methods
  !> @details
  !>
  !>          The material orientation assignment methods are selected
  !>          using kfl_fiber_sld:\n
  !>          =  1 ... Fiber model (Nodal)\n
  !>          =  2 ... Orthotropic Fiber model (Nodal)\n
  !>          =  3 ... Isotropic Fiber model (Nodal)\n
  !>          =  4 ... Material CSYS using vector fields (Element)\n
  !>          =  5 ... Material CSYS using Global system (Element)\n
  !>          =  6 ... Material CSYS using Local coordinate system (Element)\n
  !>          =  7 ... Material CSYS using element normal method (Element)\n
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_assign_material_axes()

    use def_domain,    only : ndime, nelem
    use def_domain,    only : ltype, lmate, nnode
    use def_domain,    only : xfiel
    use def_solidz,    only : kfl_fiber_sld
    use def_solidz,    only : lawst_sld, lawco_sld
    use def_solidz,    only : axis1_sld, axis2_sld, axis3_sld
    use def_solidz,    only : csysm_sld, orien_sld, oripa_sld

    implicit none

    integer(ip)            :: idime, ielem, imate, pnode, pelty
    integer(ip)            :: nposi, refax_id, rotax_id
    real(rp)               :: refaxis(ndime), rotaxis(ndime), axis(ndime), useraxis(ndime)
    real(rp)               :: rotmat(ndime,ndime), theta
    real(rp)               :: P(ndime), S(ndime), N(ndime)
    real(rp)               :: globalbasis(ndime,ndime), localbasis(ndime,ndime), userbasis(ndime,ndime)

    if ( kfl_fiber_sld == 4_ip ) then

       !------------------------------------------------------------------
       !
       ! Using Vector Fields
       !
       !------------------------------------------------------------------

       !
       ! Assign material CSYS (Vector fields)
       !
       axis1_sld => xfiel(1_ip) % a(:,:,1)
       axis2_sld => xfiel(2_ip) % a(:,:,1)
       if ( ndime == 3_ip ) axis3_sld => xfiel(3_ip) % a(:,:,1)

    else if ( kfl_fiber_sld == 5_ip ) then

       !------------------------------------------------------------------
       !
       !  Using global coordinate system (Orientation fields)
       !
       !------------------------------------------------------------------

       call sld_memphy(2_ip)

       ! Get orientations from fields
       orien_sld => xfiel(oripa_sld(1)) % a(1,:,1)

       ! Get rotation axis id
       rotax_id = oripa_sld(2)

       ! Create global CSYS
       globalbasis(:,:) = 0.0_rp
       do idime = 1, ndime
          globalbasis(idime,idime) = 1.0_rp
       end do

       !
       ! Build material coordinate system
       !
       do ielem = 1, nelem
          imate = lmate(ielem)
          theta = orien_sld(ielem)

          if ( lawst_sld(imate) == 151_ip .or. &
               lawst_sld(imate) == 152_ip .or. &
               lawst_sld(imate) == 154_ip .or. &
               (lawst_sld(imate) == 200_ip .and. lawco_sld(imate) == 2_ip) ) then
             !
             ! Orthotropic material models
             !
             ! Rotation matrix (according to the rotation axis)
             call sld_csys_rotation_matrix(rotax_id, ndime, theta, rotmat(:,:))
             !
             ! Create material CSYS
             do idime = 1, ndime
                call sld_csys_vector_rotation(1_ip, ndime, rotmat(:,:), globalbasis(:,idime), rotaxis(:))
                localbasis(:,idime) = rotaxis(:)
             end do
          else
             !
             ! Other
             !
             ! Local basis as global coordinate system
             localbasis(:,:) = 0.0_rp
             do idime = 1, ndime
                localbasis(idime,idime) = 1.0_rp
             end do

          end if
          !
          ! Assign material CSYS vectors
          !
          axis1_sld(1:ndime,ielem) = localbasis(1:ndime,1)
          axis2_sld(1:ndime,ielem) = localbasis(1:ndime,2)
          if ( ndime == 3_ip ) axis3_sld(1:3,ielem) = localbasis(1:ndime,ndime)

       end do

    else if ( kfl_fiber_sld == 6_ip ) then

       !------------------------------------------------------------------
       !
       !  Using local coordinate system (Orientation fields)
       !
       !------------------------------------------------------------------

       call sld_memphy(2_ip)

       ! Get orientations from fields
       orien_sld => xfiel(oripa_sld(1)) % a(1,:,1)

       ! Get rotation axis id
       rotax_id = oripa_sld(2)

       ! Get and build user coordinate system
       nposi = 1_ip
       userbasis(:,:) = 0.0_rp
       do idime = 1, 2
          nposi = nposi + ndime
          useraxis = csysm_sld(nposi:nposi+2)
          call maths_normalize_vector(ndime,useraxis)
          userbasis(:,idime) = useraxis
       end do
       if ( ndime == 3_ip ) then
          P = userbasis(:,1)
          S = userbasis(:,2)
          call maths_vectorial_product(P,S,N,ndime)
          userbasis(:,3) = N
       end if

       !
       ! Build material coordinate system
       !
       do ielem = 1, nelem
          imate = lmate(ielem)
          theta = orien_sld(ielem)

          if ( lawst_sld(imate) == 151_ip .or. &
               lawst_sld(imate) == 152_ip .or. &
               lawst_sld(imate) == 154_ip .or. &
              (lawst_sld(imate) == 200_ip .and. lawco_sld(imate) == 2_ip)) then
             !
             ! Orthotropic material models
             !
             ! Rotation matrix (according to the rotation axis)
             call sld_csys_rotation_matrix(rotax_id, ndime, theta, rotmat(:,:))
             !
             ! Rotation material CSYS
             do idime = 1, ndime
                axis = userbasis(:,idime)
                call sld_csys_vector_rotation(1_ip, ndime, rotmat(:,:), axis(:), rotaxis(:))
                localbasis(:,idime) = rotaxis(:)
             end do
          else
             !
             ! Other
             !
             ! Local basis as global coordinate system
             localbasis(:,:) = 0.0_rp
             do idime = 1, ndime
                localbasis(idime,idime) = 1.0_rp
             end do

          end if
          !
          ! Assign material CSYS vectors
          !
          axis1_sld(1:ndime,ielem) = localbasis(1:ndime,1)
          axis2_sld(1:ndime,ielem) = localbasis(1:ndime,2)
          if ( ndime == 3_ip ) axis3_sld(1:3,ielem) = localbasis(1:ndime,ndime)

       end do

    else if ( kfl_fiber_sld == 7_ip ) then

       !------------------------------------------------------------------
       !
       !  Using element normal respect to mid-surface (Orientation fields)
       !
       !------------------------------------------------------------------

       call sld_memphy(2_ip)

       ! Get orientations from fields
       orien_sld => xfiel(oripa_sld(1)) % a(1,:,1)

       ! Get rotation axis id
       rotax_id = oripa_sld(2)

       ! Get reference axis id and build vector
       refax_id = oripa_sld(3)
       if ( refax_id == 1_ip ) then
          refaxis(1) = 1.0_rp
          refaxis(2) = 0.0_rp
          refaxis(3) = 0.0_rp
       else if ( refax_id == 2_ip ) then
          refaxis(1) = 0.0_rp
          refaxis(2) = 1.0_rp
          refaxis(3) = 0.0_rp
       else if ( refax_id == 3_ip ) then
          refaxis(1) = 0.0_rp
          refaxis(2) = 0.0_rp
          refaxis(3) = 1.0_rp
       end if

       !
       ! Build material coordinate system
       !
       do ielem = 1, nelem
          pelty = ltype(ielem)
          imate = lmate(ielem)
          pnode = nnode(abs(pelty))
          theta = orien_sld(ielem)

          if ( lawst_sld(imate) == 151_ip .or. &
               lawst_sld(imate) == 152_ip .or. &
               lawst_sld(imate) == 154_ip .or. &
              (lawst_sld(imate) == 200_ip .and. lawco_sld(imate) == 2_ip)) then
             !
             ! Orthotropic material models
             !
             ! Create material CSYS
             call sld_midsurface_basis_rotation_element_normal(ielem, &
                  pnode, theta, refaxis(:), localbasis(:,:))

          else
             !
             ! Other
             !
             ! Local basis as global coordinate system
             localbasis(:,:) = 0.0_rp
             do idime = 1, 3
                localbasis(idime,idime) = 1.0_rp
             end do

          end if
          !
          ! Assign material CSYS vectors
          !
          axis1_sld(1:3,ielem) = localbasis(1:3,1)
          axis2_sld(1:3,ielem) = localbasis(1:3,2)
          axis3_sld(1:3,ielem) = localbasis(1:3,3)

       end do

    end if

  end subroutine sld_csys_assign_material_axes

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Build rotation matrix by using different system of coordinates
  !> @details
  !>
  !>          \verbatim
  !>
  !>           Nomenclature:
  !>           X',Y',Z' : Local cartesian axes
  !>           R        : Radial axis
  !>           T        : Tangencial axis
  !>           A        : Axial axis
  !>           N        : Exterior normal
  !>           T1,T2    : First and sencond tangencial axes
  !>
  !>                                               (Rotation matrix
  !>                                                  colum order)
  !>           - Cartesian system   (User defined)    X'- Y'- Z'
  !>
  !>           - Cylindrical system (User defined)    R - T - A
  !>
  !>           - Spherical system   (User defined)
  !>
  !>           - Exterior normal    (Alya)            N - T1 - T2
  !>
  !>          \endverbatim
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_build_jacrot()

    use def_domain, only : ndime, npoin, lpoty, exnor, coord
    use mod_maths,  only : maths_vectorial_product
    use mod_maths,  only : maths_local_orthonormal_basis
    use mod_maths,  only : maths_normalize_vector
    use def_solidz, only : jacrot_dq_du_sld, jacrot_du_dq_sld
    use def_solidz, only : kfl_fixrs_sld, kfl_csysl_sld, csysl_sld
    use def_solidz, only : SLD_CSYS_CARTESIAN, SLD_CSYS_CYLINDRICAL
    use def_solidz, only : SLD_CSYS_SPHERICAL, SLD_CSYS_EXNOR

    implicit none

    integer(ip)         :: ipoin, ibopo
    real(rp)            :: rho
    real(rp)            :: N(ndime)
    real(rp)            :: R(ndime), A(ndime), T(ndime)
    real(rp)            :: localbasis(ndime,ndime)

    do ipoin = 1, npoin
       ibopo = lpoty(ipoin)
       localbasis = 0.0_rp

       if ( ibopo > 0 ) then

          if ( kfl_fixrs_sld(ibopo) /= 0_ip ) then
             !
             ! User-defined CSYS
             !
             if ( kfl_csysl_sld == SLD_CSYS_CARTESIAN ) then
                !
                ! Cartesian CSYS
                !
                localbasis(:,1) = csysl_sld(  ndime+1:ndime*2) - csysl_sld(1)
                localbasis(:,2) = csysl_sld(ndime*2+1:ndime*3) - csysl_sld(2)
                call maths_normalize_vector(ndime,localbasis(:,1),rho)
                call maths_normalize_vector(ndime,localbasis(:,2),rho)
                !
                ! N = P x S
                if ( ndime == 3_ip ) then
                   call maths_vectorial_product(localbasis(:,1), localbasis(:,2), N(:), 3_ip)
                   localbasis(:,3) = N(:)
                end if

             else if ( kfl_csysl_sld == SLD_CSYS_CYLINDRICAL ) then
                !
                ! Cylindrical CSYS (A = Local Axial axis)
                !
                if ( ndime == 2_ip ) then
                   !
                   ! Radial axis
                   R(1) = coord(1,ipoin) - csysl_sld(1)
                   R(2) = coord(2,ipoin) - csysl_sld(2)
                   call maths_normalize_vector(ndime,R(:),rho)
                   ! Tangencial axis
                   T(1) = -R(2)
                   T(2) =  R(1)

                else if ( ndime == 3_ip ) then
                   !
                   ! Radial axis
                   if ( csysl_sld(3*ndime) /= 0_ip ) then
                      !
                      ! X-Y Plane (Axial=Z)
                      !
                      R(1) = coord(1,ipoin) - csysl_sld(1)
                      R(2) = coord(2,ipoin) - csysl_sld(2)
                      R(3) = 0.0_rp
                      call maths_normalize_vector(ndime,R(:),rho)

                   else if ( csysl_sld(3*ndime-1) /= 0_ip  ) then
                      !
                      ! X-Z Plane (Axial=Y)
                      !
                      R(1) = coord(1,ipoin) - csysl_sld(1)
                      R(2) = 0.0_rp
                      R(3) = coord(3,ipoin) - csysl_sld(3)
                      call maths_normalize_vector(ndime,R(:),rho)

                   else if ( csysl_sld(3*ndime-2) /= 0_ip ) then
                      !
                      ! Y-Z Plane (Axial=X)
                      !
                      R(1) = 0.0_rp
                      R(2) = coord(2,ipoin) - csysl_sld(2)
                      R(3) = coord(3,ipoin) - csysl_sld(3)
                      call maths_normalize_vector(ndime,R(:),rho)

                   end if
                   !
                   ! Axial axis
                   A(1) = csysl_sld(3*ndime-2)
                   A(2) = csysl_sld(3*ndime-1)
                   A(3) = csysl_sld(3*ndime)

                end if
                !
                localbasis(:,1) = R(:)
                if ( ndime == 2_ip) then
                   localbasis(:,2) = T(:)
                else
                   !
                   ! T = R x A
                   call maths_vectorial_product(R(:), A(:), T(:), 3_ip)
                   !localbasis(:,2) = A(:)
                   !localbasis(:,3) = T(:)
                   !
                   ! new
                   localbasis(:,2) = T(:)
                   localbasis(:,3) = A(:)
                end if

             else if ( kfl_csysl_sld == SLD_CSYS_SPHERICAL ) then
                !
                ! Spherical CSYS
                !
                call runend('MOD_SLD_CSYS: NOT IMPLEMENTED')

             else if ( kfl_csysl_sld == SLD_CSYS_EXNOR ) then
                !
                ! Exnor based CSYS
                !
                localbasis(:,:) = exnor(:,:,ibopo)

             end if

          else
             !
             ! Alya Exnor based CSYS (Default)
             !
             localbasis(:,:) = exnor(:,:,ibopo)

          end if
          !
          ! Built rotation matrix
          !
          jacrot_du_dq_sld(:,:,ibopo) = localbasis(:,:)            ! From U (Global) to Q (Local)
          jacrot_dq_du_sld(:,:,ibopo) = transpose(localbasis(:,:)) ! From Q (Local)  to U (Global)

       end if
    end do

  end subroutine sld_csys_build_jacrot

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January, 2019
  !> @brief   This routine rotates nodal variables from a vector
  !> @details Rotation of nodal variables using the appropiate rotation
  !>          vector at point level.
  !>          Rotations are performed only for local axes.
  !>
  !>          \verbatim
  !>          ITASK = 1 ... From global to local
  !>          ITASK = 2 ... From local to global
  !>          \endverbatim
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_rotuni(itask,ndofn,ipoin,unrot)

    use def_master, only : INOTMASTER
    use def_domain, only : lpoty
    use def_solidz, only : jacrot_du_dq_sld

    implicit none

    integer(ip), intent(in)    :: itask          !< Transformation
    integer(ip), intent(in)    :: ndofn          !< Dimensions
    integer(ip), intent(in)    :: ipoin          !< Point id
    real(rp),    intent(inout) :: unrot(ndofn)   !< Vector transformed

    integer(ip)                :: ibopo
    real(rp)                   :: worma(ndofn), worve(ndofn),romatt(ndofn,ndofn)

    if ( INOTMASTER ) then

       select case (itask)

       case( 1_ip )
          !
          ! Global to local A = B^T * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the user local system
             !
             worve(1:ndofn) = unrot
             romatt = transpose(jacrot_du_dq_sld(:,:,ibopo))
             worma = matmul(romatt,worve)
             unrot(1:ndofn) = worma(1:ndofn)
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       case( 2_ip )
          !
          ! Local to global: A = B * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the local system
             !
             worve(1:ndofn) = unrot(1:ndofn)
             worma = matmul(jacrot_du_dq_sld(:,:,ibopo),worve)
             unrot(1:ndofn) = worma(1:ndofn)
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       end select

    end if

  end subroutine sld_csys_rotuni

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January, 2019
  !> @brief   This routine rotates nodal variables from a matrix
  !> @details Rotation of nodal variables using the appropiate rotation
  !>          matrix at point level.
  !>          Rotations are performed only for local axes.
  !>
  !>          \verbatim
  !>          ITASK = 1 ... From global to local
  !>          ITASK = 2 ... From local to global
  !>          \endverbatim
  !>
  !----------------------------------------------------------------------------

  subroutine sld_csys_rotmat(itask,ndofn,ipoin,unrot)

    use def_master, only : INOTMASTER
    use def_domain, only : lpoty
    use def_solidz, only : jacrot_du_dq_sld

    implicit none

    integer(ip), intent(in)    :: itask              !< Transformation
    integer(ip), intent(in)    :: ndofn              !< Dimensions
    integer(ip), intent(in)    :: ipoin              !< Point id
    real(rp),    intent(inout) :: unrot(ndofn,ndofn) !< Matrix transformed

    integer(ip)                :: ibopo
    real(rp)                   :: worma(ndofn,ndofn), romatt(ndofn,ndofn)

    if ( INOTMASTER ) then

       select case (itask)

       case( 1_ip )
          !
          ! Global to local A = B^T * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the user local system
             !
             worma(1:ndofn,1:ndofn) = unrot(1:ndofn,1:ndofn)
             romatt = transpose(jacrot_du_dq_sld(:,:,ibopo))
             worma  = matmul(worma,romatt)
             unrot  = matmul(worma,jacrot_du_dq_sld(:,:,ibopo))
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       case( 2_ip )
          !
          ! Local to global: A = B * C
          !
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             !
             ! Boundary conditions in the local system
             !
             call runend('MOD_SLD_CSYS: ROTATION MATRIX FROM LOCAL TO GLOBAL NOT IMPLEMENTED ')
          else
             call runend('MOD_SLD_CSYS: THIS POINT DOES NOT HAVE A ROTATION MATRIX FOR LOCAL AXES')
          end if

       end select

    end if

  end subroutine sld_csys_rotmat

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Basis creation using element normal and a reference axis
  !> @details
  !>
  !>          /verbatim
  !>          The basis creation is orthonormal. The steps are the
  !>          following:
  !>          Step-1: Calculation element normal N and primary axis P.
  !>                  P is orthonormal to N and belongs to the midsurface
  !>                  of the element.
  !>          Step-2: Rotation matrix calculated from N and user rotation angle
  !>          Step-3: Rotation of primary axis P
  !>          Step-4: Calculation secondary axis S as the vectorial product of
  !>                  S = N x P
  !>          /endverbatim
  !>
  !> @note    This method can only be used for HEX08 and PEN06 elements with
  !>          certain limitions:
  !>           - Element connectivity (bottom to top)
  !>           - A reference axis is required
  !>
  !----------------------------------------------------------------------------

  subroutine sld_midsurface_basis_rotation_element_normal(ielem, pnode, &
       theta, refvec, basis)

    implicit none

    integer(ip), intent(in)    :: ielem      !< Element number
    integer(ip), intent(in)    :: pnode      !< Number of nodes
    real(rp),    intent(in)    :: theta      !< Rotation angle by the user
    real(rp),    intent(in)    :: refvec(3)  !< Reference vector by the user
    real(rp),    intent(out)   :: basis(3,3) !< Orthonormal basis

    real(rp)                   :: rotmat(3,3), crefaxis(3)
    real(rp)                   :: P(3), S(3), N(3)

    !
    ! Normal element vector (n) and primary axis or reference (p)
    !
    call sld_csys_midsurface_element_normal(ielem, 3_ip, pnode, N(:), refvec(:), crefaxis(:))
    !
    ! Rotation matrix calculated from normal axis and angle
    !
    call sld_csys_rotation_matrix_from_axis_and_angle(N(:), theta, rotmat(:,:))
    !
    ! Rotate reference vector (Global to Local)
    !
    call sld_csys_vector_rotation(1_ip, 3_ip, rotmat(:,:), crefaxis(:), P(:))
    !
    ! S = N x P
    !
    call maths_vectorial_product(N(:), P(:), S(:), 3_ip)
    !
    ! Orthonormal CSYS
    !
    basis(:,:) = 0.0_rp
    basis(:,1) = P(:)
    basis(:,2) = S(:)
    basis(:,3) = N(:)

  end subroutine sld_midsurface_basis_rotation_element_normal

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Element normal vector or thickness direction of the element
  !> @details
  !>          /verbatim
  !>          To compute the thickness direction or normal direction for
  !>          2-d and 3-d elements, Alya forms a midsurface by averaging
  !>          the coordinates of the node pairs forming the bottom and top
  !>          surfaces of the element. This midsurface passes through the
  !>          integration points of the element using a Newton Cotes scheme.
  !>          For each integration point Alya computes a tangent whose direction
  !>          is defined by the sequence of nodes given by the midsurface
  !>
  !>          2-D
  !>          For 2-d elements the normal direction is obtained as the cross
  !>          product between the out-of-plane and tangent directions.
  !>
  !>             (4) _ n1
  !>              o  /|     (3) _ n2
  !>              | /        o  /|
  !>              |/         | /
  !>              x--> t1    |/
  !>              |          x --> t2
  !>              |          |
  !>              o----------o
  !>             (1)        (2)
  !>
  !>          3-D
  !>          For 3-d elements the normal direction to the midsurface
  !>          corresponds to the positive direction which is obtained with the
  !>          right-hand rule going around the nodes of the element on the
  !>          bottom or top surface.
  !>
  !>             (5) _ n1               (8) _ n4   (7)
  !>              o  /|     (6) _ n2     o  /|      o  _ n3
  !>              | /        o  /|       | /        |  /|
  !>              |/         | /         |/         | /
  !>              x--> t1    |/          x--> t4    x --> t3
  !>              |          x --> t2    |          |
  !>              |          |           |          |
  !>              o----------o           o----------o
  !>             (1)        (2)         (4)        (3)
  !>
  !>           /endverbatim
  !>
  !------------------------------------------------------------------------------

  subroutine sld_csys_midsurface_element_normal(ielem, pdime, pnode, &
       N, refvec, P)

    use def_elmtyp,       only : HEX08, QUA04
    use def_domain,       only : lnods, ltype
    use def_domain,       only : coord

    implicit none

    integer(ip), intent(in)            :: ielem         !< Element number
    integer(ip), intent(in)            :: pdime         !< Number of dimension
    integer(ip), intent(in)            :: pnode         !< Number of nodes
    real(rp),    intent(out)           :: N(pdime)      !< Normal unit vector
    real(rp),    intent(in),  optional :: refvec(pdime) !< Reference unit vector
    real(rp),    intent(out), optional :: P(pdime)      !< Primary unit  vector

    integer(ip)               :: pelty
    integer(ip)               :: inode, igaus, ipoin
    integer(ip)               :: sdime, snode, sgaus
    logical(lg)               :: kflag_axisalig
    real(rp)                  :: elcoo(pdime,pnode)     ! Nodal element coordinates
    real(rp)                  :: sfcoo(pdime,pnode/2)   ! Midsurface coordinates
    real(rp)                  :: gpnorm(pdime)          ! GP Unit midsurface normal
    real(rp)                  :: gptan1(pdime)          ! GP Unit midsurface first tangent
    real(rp)                  :: gptan2(pdime)          ! GP Unit midsurface second tangent
    real(rp)                  :: gpder(pdime-1,pnode/2) ! GP midsurface derivative of the shape function dN/dn
    real(rp)                  :: gprefv(pdime)          ! GP midsurface tangent correction

    kflag_axisalig = .false.  !<GGU> Testing
    !
    ! Compatibility checks
    !
    pelty = ltype(ielem)
    if ( pelty /= HEX08 .and. pelty /= QUA04 ) then
       call runend('SLD_MID_SURFACE_ELEMENT_NORMAL: WRONG ELEMENT TYPE')
    end if
    !
    ! Auxiliar dimensions for mid-surface
    !
    sdime = int(pdime-1) ! Mid-surface number of dimensions
    snode = int(pnode/2) ! Mid-surface number of nodes
    sgaus = int(pnode/2) ! Mid-surface gauss points (It has to be equal to half o t)
    !
    ! Element coordinates
    !
    do inode = 1, pnode
       ipoin = lnods(inode,ielem)
       elcoo(1:pdime,inode) = coord(1:pdime,ipoin)
    end do
    !
    ! Mid-surface coordinates
    !
    call sld_midsurface_coordinates(pdime, pnode, snode, elcoo, sfcoo)
    !
    ! Compute averaged normal vector and primary axis vectors based on the midsurface
    !
    N(:) = 0.0_rp
    P(:) = 0.0_rp
    gausspoints: do igaus = 1, sgaus
       !
       ! Interpolating function and its derivative (N and dN/dn)
       !
       call sld_midsurface_shape_functions(igaus, pdime, sdime, snode, gpder(:,:))
       !
       ! Vectors from midsurface ( Normal and tangents)
       !
       call sld_midsurface_gp_vectors(pdime, sdime, snode, sfcoo(:,:), gpder(:,:), &
            gpnorm(:), gptan1(:), gptan2(:))
       !
       ! Build primary axis based on refaxis and tangents
       !
       if ( present(refvec) ) then
          call sld_midsurface_primary_axis(pdime, kflag_axisalig, refvec(:),       &
               gpnorm(:), gptan1(:), gptan2(:), gprefv(:))
       end if
       !
       ! Averaged normal and primary/reference vectors
       !
       N(:) = N(:) + gpnorm(:)
       if ( present(P) ) P(:) = P(:) + gprefv(:)

    end do gausspoints
    !
    ! Normalize vectors
    !
    N(:) = N(:)/real(sgaus,rp)
    call maths_normalize_vector(pdime,N(:))
    if ( present(P) ) call maths_normalize_vector(pdime,P(:))

  end subroutine sld_csys_midsurface_element_normal

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Midsurface coordinates for hexahedron and pentahedron element
  !>          types using the reference configuration (undeformed mesh)
  !> @details
  !>          2-D elements:
  !>            - QUA04
  !>          3-D elements:
  !>            - HEX08 and PEN06
  !------------------------------------------------------------------------------

  subroutine sld_midsurface_coordinates(pdime, pnode, snode, elcoo, sfcoo)

    implicit none

    integer(ip), intent(in)   :: pdime              !< Dimension
    integer(ip), intent(in)   :: pnode              !< Number of nodes
    integer(ip), intent(in)   :: snode              !< Mid surface number of nodes
    real(rp),    intent(in)   :: elcoo(pdime,pnode) !< Nodal element coordinates
    real(rp),    intent(out)  :: sfcoo(pdime,snode) !< Midsurface coordinates

    integer(ip)               :: inode

    !
    ! Mid-surface coordinates
    !
    sfcoo(:,:) = 0.0_rp
    if ( pdime == 2_ip ) then       ! 2-d
       sfcoo(1:pdime,1) = 0.5_rp*(elcoo(1:pdime,1) + elcoo(1:pdime,4))
       sfcoo(1:pdime,2) = 0.5_rp*(elcoo(1:pdime,2) + elcoo(1:pdime,3))
    else                            ! 3-d
       do inode = 1, snode
          sfcoo(1:pdime,inode) = 0.5_rp*(elcoo(1:pdime,inode) + elcoo(1:pdime,inode+snode))
       end do
    end if

  end subroutine sld_midsurface_coordinates

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Compute shape function based on Newton Cotes
  !> @details
  !>          Natural coordiantes of the integration points position
  !>
  !>            s                         s
  !>  2-d  o----|----o          3-d       /
  !>       |    |    |            4 X ---/--- X 3
  !>       |    |    |             /    /    /
  !>     1 X---------X 2 --t      /    /_ _ /_ _ t
  !>       |         |           /         /
  !>       |         |        1 X ------- X 2
  !>       o---------o
  !------------------------------------------------------------------------------

  subroutine sld_midsurface_shape_functions( igaus, pdime, sdime, snode, gpderi )

    implicit none

    integer(ip), intent(in)  :: igaus               !< Integration point identification
    integer(ip), intent(in)  :: pdime               !< Dimension
    integer(ip), intent(in)  :: sdime               !< Dimension number of the mid-surface
    integer(ip), intent(in)  :: snode               !< Number of nodes of the mid-surface
    real(rp),    intent(out) :: gpderi(sdime,snode) !< GP Derivative shape function

    real(rp)                 :: t, s                ! Natural coordinates (s, t)
    real(rp)                 :: gpposi(sdime,snode) ! Natural coordintaes of the integration points

    ! Init
    gpposi(:,:) = 0.0_rp
    gpderi(:,:) = 0.0_rp

    if ( pdime == 2_ip ) then             ! 2-d
       gpposi(1,1) = -1.0_rp
       gpposi(1,2) =  1.0_rp
       t = gpposi(1,igaus)
       ! GP Derivative Shape functions w.r.t. natural coordinates (dN/dn)
       gpderi(1,1) = -0.5_rp
       gpderi(1,2) =  0.5_rp

    else                                  ! 3-d
       gpposi(1,1) = -1.0_rp
       gpposi(2,1) = -1.0_rp
       gpposi(1,2) =  1.0_rp
       gpposi(2,2) = -1.0_rp
       gpposi(1,3) =  1.0_rp
       gpposi(2,3) =  1.0_rp
       gpposi(1,4) = -1.0_rp
       gpposi(2,4) =  1.0_rp
       t = gpposi(1,igaus)
       s = gpposi(2,igaus)
       ! GP Derivative Shape function w.r.t. natural coordinates (dN/dn)
       gpderi(1,1) = 0.25_rp*(-1.0_rp + s)
       gpderi(1,2) = 0.25_rp*( 1.0_rp - s)
       gpderi(1,3) = 0.25_rp*( 1.0_rp + s)
       gpderi(1,4) = 0.25_rp*(-1.0_rp - s)
       !
       gpderi(2,1) = 0.25_rp*(-1.0_rp + t)
       gpderi(2,2) = 0.25_rp*(-1.0_rp - t)
       gpderi(2,3) = 0.25_rp*( 1.0_rp + t)
       gpderi(2,4) = 0.25_rp*( 1.0_rp - t)
    end if

  end subroutine sld_midsurface_shape_functions

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Midsurface vectors
  !> @details Maths:
  !>            x  = (X^{top} + X^{bot} + q^{top} + q^{bot})
  !>            t1 = 1/2*(dN/dn1)*x
  !>            t2 = 1/2*(dN/dn2)*x
  !>            n = t1 x t2
  !>          where:
  !>            X are the nodal coordinates at the reference configuration
  !>            q are the nodal displacements at the reference configuration
  !>            N are the shape functions of the interface element, i.e. for
  !>              - 3-D 8-nodes element -> 4-nodes element.
  !>                (*)^{bot} are the bottom nodes, i.e. 1 - 2 - 3 - 4 (3D 8-nodes)
  !>                (*)^{top} are the top nodes,    i.e. 5 - 6 - 7 - 8 (3D 8-nodes)
  !>              - 2-d 4-nodes element -> 2-nodes element.
  !>                (*)^{bot} are the bottom nodes, i.e. 1 - 2 (2D 4-nodes)
  !>                (*)^{top} are the top nodes,    i.e. 3 - 4 (2D 4-nodes)
  !------------------------------------------------------------------------------

  subroutine sld_midsurface_gp_vectors( &
       pdime, sdime, snode, sfcoo, gpder, gpnorm, gptan1, gptan2, gparea)

    implicit none

    integer(ip), intent(in)            :: pdime              !< Number of dimensions
    integer(ip), intent(in)            :: sdime              !< Number of dimension midsurface
    integer(ip), intent(in)            :: snode              !< Mid-surface number of nodes
    real(rp),    intent(in)            :: sfcoo(pdime,snode) !< Mid-surface coordinates
    real(rp),    intent(in)            :: gpder(sdime,snode) !< GP derivative shape function
    real(rp),    intent(out)           :: gpnorm(pdime)      !< GP Unit normal vector
    real(rp),    intent(out)           :: gptan1(pdime)      !< GP Unit first tangent
    real(rp),    intent(out)           :: gptan2(pdime)      !< GP Unit second tangent
    real(rp),    intent(out), optional :: gparea             !< GP area

    integer(ip)                        :: idime, inode
    real(rp)                           :: mod

    !
    ! Initialize variables
    !
    if ( present(gparea) ) gparea = 0.0_rp
    gpnorm(:) = 0.0_rp
    gptan1(:) = 0.0_rp
    gptan2(:) = 0.0_rp
    !
    if (pdime == 2_ip ) then      ! 2-D

       ! Tangent vector midsurface
       do inode = 1, snode
          do idime = 1, pdime
             gptan1(idime) = gptan1(idime) + gpder(1,inode)*sfcoo(idime,inode)
          end do
       end do
       call maths_normalize_vector(pdime,gptan1,mod)
       ! Normal unit vector
       gpnorm(1) = -gptan1(2)
       gpnorm(2) =  gptan1(1)

       if ( present(gparea) ) gparea = mod

    else                          ! 3-D

       ! Tangent vectors midsurface
       do inode = 1, snode
          do idime = 1, pdime
             gptan1(idime) = gptan1(idime) + gpder(1,inode)*sfcoo(idime,inode)
             gptan2(idime) = gptan2(idime) + gpder(2,inode)*sfcoo(idime,inode)
          end do
       end do
       ! Compute normal director vector n = t1 X t2
       call maths_vectorial_product(gptan1,gptan2,gpnorm,pdime)

       ! Normalize vectors
       call maths_normalize_vector(pdime,gptan1)
       call maths_normalize_vector(pdime,gptan2)
       call maths_normalize_vector(pdime,gpnorm,mod)

       if ( present(gparea) ) gparea = mod

    end if

  end subroutine sld_midsurface_gp_vectors

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Reference axis correction
  !> @details Reference axis correction
  !>
  !>          /verbatim
  !>          The primary axis is created based on the reference axis and the
  !>          midsurface plane of the element. The steps are the following:
  !>          Step-1:            Reference axis is projected to the midsurface
  !>                             plane in order to be an orthonomal vector.
  !>          Step-2 (Optional): Once the reference axis is projected to the
  !>                             the midsurface, a correction is applied which
  !>                             forces the primary axis to midsurface tangent
  !>                             vector director.
  !>          /endverbatim
  !>
  !------------------------------------------------------------------------------

  subroutine sld_midsurface_primary_axis(pdime, kflag_axisalig, refvec, gpnorm, &
       gptan1, gptan2, gprefc)

    implicit none

    integer(ip), intent(in)  :: pdime           !< Number of dimensions
    logical(lg), intent(in)  :: kflag_axisalig  !< Flag axis element edge alignement
    real(rp),    intent(in)  :: refvec(pdime)   !< Unit reference vector
    real(rp),    intent(in)  :: gpnorm(pdime)   !< GP Unit normal vector
    real(rp),    intent(in)  :: gptan1(pdime)   !< GP Unit first tangent
    real(rp),    intent(in)  :: gptan2(pdime)   !< GP Unit second tangent
    real(rp),    intent(out) :: gprefc(pdime)   !< GP Unit reference vector corrected

    real(rp)                 :: aux, aux1, aux2

    !
    ! Vector projection to the midsurface
    !
    aux = dot_product(refvec, gpnorm)
    gprefc(:) = refvec(:) -  aux*gpnorm(:)

    !
    ! Force alignment with element midsurface edge (most parallel with refaxis)
    !
    if ( kflag_axisalig ) then
       !
       ! Reference vector correction to the element midsurface vectors
       !
       aux1 = dot_product(gprefc,gptan1)
       aux2 = dot_product(gprefc,gptan2)

       if (abs(aux1) /= 1.0_rp) then      ! t1 not parallel to refaxis

          if (abs(aux1) > abs(aux2)) then ! Choose most parallel vector
             gprefc(:) = gptan1(:)
          else
             gprefc(:) = gptan2(:)
          end if

       end if

       if (abs(aux2) /= 1.0_rp) then      ! t2 not parallel to refaxis

          if (abs(aux1) < abs(aux2)) then
             gprefc(:) = gptan2(:)
          else
             gprefc(:) = gptan1(:)
          end if

       end if

    end if

  end subroutine sld_midsurface_primary_axis

end module mod_sld_csys
