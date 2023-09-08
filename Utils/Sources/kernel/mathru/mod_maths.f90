!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical functions and subroutines
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for maths
!> @details ToolBox for maths
!
!-----------------------------------------------------------------------

module mod_maths

  use def_kintyp, only : ip,rp,lg
  use def_parame, only : pi
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_memory, only : memory_resize
  use mod_memory, only : memory_size
  use mod_std                           ! defintion of qp and count
  implicit none
  integer(8)                            :: memor(2)
#ifndef __PGI 
  integer,     parameter                :: qp = 16  
#endif
  real(rp),    parameter                :: epsil = epsilon(1.0_rp)
  integer(ip), parameter,dimension(192) :: postb3 = (/ &
       0,2,3,1,5,7,6,4,0,4,5,1,3,7,6,2,0,4,6,2,3,7,5,1,5,4,0,1,3,2,6,7,&
       3,7,5,1,0,4,6,2,6,4,5,7,3,1,0,2,0,2,6,4,5,7,3,1,0,1,5,4,6,7,3,2,&
       6,2,0,4,5,1,3,7,5,7,6,4,0,2,3,1,3,2,6,7,5,4,0,1,3,2,0,1,5,4,6,7,&
       5,7,3,1,0,2,6,4,6,7,3,2,0,1,5,4,6,2,3,7,5,1,0,4,5,1,3,7,6,2,0,4,&
       0,1,3,2,6,7,5,4,3,1,0,2,6,4,5,7,6,4,0,2,3,1,5,7,3,7,6,2,0,4,5,1,&
       5,4,6,7,3,2,0,1,6,7,5,4,0,1,3,2,3,1,5,7,6,4,0,2,5,1,0,4,6,2,3,7/)
  integer(2), parameter,dimension(192) :: typtb3 = (/ &
       1_2, 6_2, 6_2,11_2,11_2,12_2,12_2,14_2, 0_2, 2_2, 2_2, 3_2, 3_2, 4_2, 4_2, 5_2,&
       16_2, 1_2, 1_2,18_2,18_2,19_2,19_2,20_2,12_2,20_2,20_2, 1_2, 1_2,11_2,11_2,18_2,&
       11_2,19_2,19_2,12_2,12_2, 1_2, 1_2,21_2,14_2,18_2,18_2,20_2,20_2,22_2,22_2, 1_2,&
       7_2, 0_2, 0_2, 8_2, 8_2, 9_2, 9_2,10_2, 6_2,16_2,16_2,23_2,23_2,21_2,21_2,22_2,&
       21_2,14_2,14_2, 6_2, 6_2,23_2,23_2,11_2,23_2,12_2,12_2,21_2,21_2, 6_2, 6_2,19_2,&
       22_2,11_2,11_2,14_2,14_2,20_2,20_2, 6_2, 4_2,10_2,10_2, 0_2, 0_2, 3_2, 3_2, 8_2,&
       3_2, 9_2, 9_2, 4_2, 4_2, 0_2, 0_2,13_2,18_2,21_2,21_2,19_2,19_2,16_2,16_2,12_2,&
       5_2, 8_2, 8_2,10_2,10_2,15_2,15_2, 0_2,20_2,23_2,23_2,22_2,22_2,14_2,14_2,16_2,&
       2_2, 7_2, 7_2,17_2,17_2,13_2,13_2,15_2,19_2,22_2,22_2,16_2,16_2,18_2,18_2,23_2,&
       13_2, 5_2, 5_2, 2_2, 2_2,17_2,17_2, 3_2,17_2, 4_2, 4_2,13_2,13_2, 2_2, 2_2, 9_2,&
       15_2, 3_2, 3_2, 5_2, 5_2,10_2,10_2, 2_2, 8_2,13_2,13_2, 9_2, 9_2, 7_2, 7_2, 4_2,&
       10_2,17_2,17_2,15_2,15_2, 5_2, 5_2, 7_2, 9_2,15_2,15_2, 7_2, 7_2, 8_2, 8_2,17_2/)


  private

  interface maths_equalize_arrays
     module procedure maths_equalize_arrays_RP_11,&
          &           maths_equalize_arrays_RP_12,&
          &           maths_equalize_arrays_RP_22,&
          &           maths_equalize_arrays_RP_ndim1,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_22,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_32,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_23,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_33
  end interface maths_equalize_arrays

  interface maths_array_permutation
     module procedure maths_array_permutation_IP_1,&
          &           maths_array_permutation_IP_2
  end interface maths_array_permutation

  interface maths_sfc_1d_to2d3d_tab
     module procedure maths_sfc_d2xy_tab,&
          &           maths_sfc_d2xyz_tab
  end interface maths_sfc_1d_to2d3d_tab

  interface maths_heap_sort
     module procedure maths_heap_sort_I1,&
          &           maths_heap_sort_I2,&
          &           maths_heap_sort_RP1
  end interface maths_heap_sort
  interface maths_quick_sort
     module procedure maths_quick_sort_rp,&
          &           maths_quick_sort_rp_1,&
          &           maths_quick_sort_ip,&
          &           maths_quick_sort_ip_1
  end interface maths_quick_sort

  public :: maths_sfc_d2xy_tab
  public :: maths_mapping_3d_to_1d
  public :: maths_mapping_1d_to_3d
  public :: maths_mapping_1d_to_3d_x
  public :: maths_mapping_1d_to_3d_y
  public :: maths_mapping_1d_to_3d_z
  public :: maths_mapping_coord_to_3d
  public :: maths_in_box
  public :: maths_invert_matrix
  public :: maths_equalize_arrays
  public :: maths_heap_sort
  public :: maths_quick_sort
  public :: maths_geometrical_sort_using_sfc
  public :: maths_geometrical_sort_using_coordinates
  public :: maths_rank_table
  public :: maths_local_orthonormal_basis
  public :: maths_vector_to_new_basis
  public :: maths_vector_from_new_basis
  public :: maths_array_permutation
  public :: maths_time_unit
  public :: maths_schur_complement
  public :: maths_angle_of_a_vector
  public :: maths_outer_product
  public :: maths_norm2                         ! L2 norm of an array
  public :: maths_backsu
  public :: maths_matrix_multiplication
  public :: maths_matrix_vector_multiplication
  public :: maths_economy_qrhousehold
  public :: maths_sfc_part
  public :: maths_sfc_par_part
  public :: maths_sfc_par_part_2
  public :: maths_sfc_1d_to2d3d_tab
  public :: maths_day_of_week                   ! Day of the week
  public :: maths_list_different_elements       ! List the different elements of an array
  public :: maths_unit                          ! Scaling and basic unit of a value
  public :: maths_vectorial_product             ! Vectorial product
  public :: maths_normalize_vector              ! Normalize a vector
  public :: maths_merge_ordered_lists           ! Merge ordered list
  public :: maths_solve_overdetermined_system
  public :: maths_quadratic_equation            ! Solve a quadratic equation
  public :: maths_Szudzik_pairing_function
  public :: maths_maxloc_nonzero                ! Last non-zero position
  
contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  subroutine maths_equalize_arrays_RP_ndim1(ndim1,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    real(rp),          intent(in)  :: x_in(*)
    real(rp),          intent(out) :: x_out(*)
    integer(ip)                    :: idim1
    do idim1 = 1,ndim1
       x_out(idim1) = x_in(idim1)
    end do
  end subroutine maths_equalize_arrays_RP_ndim1
  subroutine maths_equalize_arrays_RP_ndim1_ndim2_22(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2)
    real(rp),          intent(out) :: x_out(ndim1,ndim2)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2) = x_in(idim1,idim2)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_22
  subroutine maths_equalize_arrays_RP_ndim1_ndim2_32(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2,*)
    real(rp),          intent(out) :: x_out(ndim1,ndim2)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2) = x_in(idim1,idim2,1)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_32
  subroutine maths_equalize_arrays_RP_ndim1_ndim2_23(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2)
    real(rp),          intent(out) :: x_out(ndim1,ndim2,*)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2,1) = x_in(idim1,idim2)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_23
  subroutine maths_equalize_arrays_RP_ndim1_ndim2_33(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2,*)
    real(rp),          intent(out) :: x_out(ndim1,ndim2,*)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2,1) = x_in(idim1,idim2,1)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_33

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  subroutine maths_equalize_arrays_RP_11(x_in,x_out)
    real(rp), pointer, intent(in)    :: x_in(:)
    real(rp), pointer, intent(inout) :: x_out(:)
    integer(ip)                      :: idim1
    integer(ip)                      :: ndim1

    ndim1 = min(size(x_in,KIND=ip),size(x_out,KIND=ip))

    do idim1 = 1,ndim1
       x_out(idim1) = x_in(idim1)
    end do

  end subroutine maths_equalize_arrays_RP_11

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  subroutine maths_equalize_arrays_RP_12(x_in,x_out)
    real(rp), pointer, intent(in)    :: x_in(:)
    real(rp), pointer, intent(inout) :: x_out(:,:)
    integer(ip)                      :: idim1,idim2,idime
    integer(ip)                      :: ndim1,ndim2

    ndim1 = size(x_out,1,KIND=ip)
    ndim2 = size(x_out,2,KIND=ip)

    if( size(x_in,KIND=ip) /= ndim1*ndim2 ) then
       call runend('maths_equalize_arrays: ARRAYS ARE NOT OF SAME DIMENSIONS')
    else
       idime = 0
       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             idime = idime + 1
             x_out(idim1,idim2) = x_in(idime)
          end do
       end do
    end if

  end subroutine maths_equalize_arrays_RP_12

  subroutine maths_equalize_arrays_RP_22(x_in,x_out)
    real(rp), pointer, intent(in)    :: x_in(:,:)
    real(rp), pointer, intent(inout) :: x_out(:,:)
    integer(ip)                      :: idim1,idim2
    integer(ip)                      :: ndim1,ndim2

    ndim1 = size(x_out,1,KIND=ip)
    ndim2 = size(x_out,2,KIND=ip)

    if( size(x_in,1,KIND=ip) /= ndim1 .and. size(x_in,2,KIND=ip) /= ndim2 ) then
       call runend('maths_equalize_arrays: ARRAYS ARE NOT OF SAME DIMENSIONS')
    else
       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             x_out(idim1,idim2) = x_in(idim1,idim2)
          end do
       end do
    end if

  end subroutine maths_equalize_arrays_RP_22

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping (x,y,z) => s
  !> @details Mapping 3D array <=> 1D array
  !>          x: 1 to mx
  !>          y: 1 to my
  !>          z: 1 to mz
  !>          s: 1 to mx*my*mx
  !>
  !>          3D to 1D: s = x + mx*[ (y-1) + (z-1)*my ]
  !>          1D to 3D: x = mod(s-1,mx) + 1
  !>                    y = mod((s-1)/mx,my) + 1
  !>                    z = (s-1)/(mx*my) + 1
  !
  !----------------------------------------------------------------------

  function maths_mapping_3d_to_1d(mx,my,mz,xx,yy,zz)
    integer(ip), intent(in) :: mz,my,mx
    integer(ip), intent(in) :: zz,yy,xx
    integer(ip)             :: maths_mapping_3d_to_1d
    maths_mapping_3d_to_1d = xx+mx*( yy-1 + (zz-1)*my)
  end function maths_mapping_3d_to_1d

  function maths_mapping_1d_to_3d_x(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_x
    maths_mapping_1d_to_3d_x = modulo(ss-1,mx)+1
  end function maths_mapping_1d_to_3d_x

  function maths_mapping_1d_to_3d_y(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_y
    maths_mapping_1d_to_3d_y = modulo((ss-1)/mx,my)+1
  end function maths_mapping_1d_to_3d_y

  function maths_mapping_1d_to_3d_z(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_z
    maths_mapping_1d_to_3d_z = (ss-1) / (mx*my)+1
  end function maths_mapping_1d_to_3d_z

  subroutine maths_mapping_1d_to_3d(ss,mx,my,mz,xx,yy,zz)
    integer(ip), intent(in)  :: ss,mz,my,mx
    integer(ip), intent(out) :: zz,yy,xx

    xx = modulo(ss-1,mx)+1
    yy = modulo((ss-1)/mx,my)+1
    zz = (ss-1) / (mx*my)+1

  end subroutine maths_mapping_1d_to_3d

  !-----------------------------------------------------------------------
  !
  !> @date    27/02/2014
  !> @author  Guillaume Houzeaux
  !> @brief   Find the box in a bin structure
  !> @details Find the box (II,JJ,KK) the point COORD is located in.
  !>          To detect if the point is outside, put the message
  !>          'DETECT OUTSIDE' as an argument. Then,
  !>          II = 0 if points is out of the bin
  !>
  !
  !-----------------------------------------------------------------------

  subroutine maths_mapping_coord_to_3d(ndime,boxes,comin,comax,coord,ii,jj,kk,message)
    integer(ip), intent(in)           :: ndime          !< Dimension of problem
    integer(ip), intent(in)           :: boxes(ndime)   !< # boxes in each dimension
    real(rp),    intent(in)           :: comin(ndime)   !< Minimum bin coordinates
    real(rp),    intent(in)           :: comax(ndime)   !< maximum bin coordinates
    real(rp),    intent(in)           :: coord(ndime)   !< Coordinate of the test point
    integer(ip), intent(out)          :: ii             !< Box in x direction
    integer(ip), intent(out)          :: jj             !< Box in y direction
    integer(ip), intent(out)          :: kk             !< Box in z direction
    character(*),intent(in), optional :: message        !< Box in z direction
    logical(lg)                       :: detect_outside

    ii = int( ( (coord(1)-comin(1)-epsil) / (comax(1)-comin(1)) )*real(boxes(1),rp), ip ) + 1
    if( ndime >= 2 ) then
       jj = int( ( (coord(2)-comin(2)-epsil) / (comax(2)-comin(2)) )*real(boxes(2),rp), ip ) + 1
    else
       jj = 1
    end if
    if( ndime >= 3 ) then
       kk = int( ( (coord(3)-comin(3)-epsil) / (comax(3)-comin(3)) )*real(boxes(3),rp), ip ) + 1
    else
       kk = 1
    end if
    !
    ! Out of bin?
    !
    detect_outside = .false.
    if( present(message) ) then
       if( trim(message) == 'DETECT OUTSIDE') detect_outside = .true.
    end if

    if( detect_outside ) then
       if( ii < 1 .or. ii > boxes(1) ) ii = 0
       if( ndime >= 2 .and. ( jj < 1 .or. jj > boxes(2) ) ) ii = 0
       if( ndime == 3 .and. ( kk < 1 .or. kk > boxes(3) ) ) ii = 0
    else
       ii = min(max(1_ip,ii),boxes(1))
       if( ndime >= 2 ) jj = min(max(1_ip,jj),boxes(2))
       if( ndime == 3 ) kk = min(max(1_ip,kk),boxes(3))
    end if

  end subroutine maths_mapping_coord_to_3d

  !-----------------------------------------------------------------------
  !
  !> @date    23/05/2014
  !> @author  Guillaume Houzeaux
  !> @brief   If in a box
  !> @details If in a box
  !
  !-----------------------------------------------------------------------

  function maths_in_box(ndime,xx,box_comin,box_comax)
    integer(ip), intent(in) :: ndime
    real(rp),    intent(in) :: xx(ndime)
    real(rp),    intent(in) :: box_comin(ndime)
    real(rp),    intent(in) :: box_comax(ndime)
    logical(lg)             :: maths_in_box
    integer(ip)             :: idime

    maths_in_box = .true.
    do idime = 1,ndime
       if( xx(idime) < box_comin(idime) .or. xx(idime) > box_comax(idime) ) then
          maths_in_box = .false.
          return
       end if
    end do

  end function maths_in_box

  !-----------------------------------------------------------------------
  !
  !> @date    27/06/2014
  !> @author  Guillaume Houzeaux
  !> @brief   Inverse a matrix
  !> @details Invert
  !
  !-----------------------------------------------------------------------

  subroutine maths_invert_matrix(nsize,a,deter,inva,ACCURACY)
    implicit none
    integer(ip), intent(in)             :: nsize
    real(rp),    intent(inout)          :: a(nsize,nsize)
    real(rp),    intent(out),  optional :: inva(nsize,nsize)
    real(rp),    intent(out),  optional :: deter
    logical(lg), intent(in),   optional :: accuracy
    logical(lg)                         :: accurate_method
    real(rp)                            :: denom,t1,t2,t3,t4,det
    real(rp)                            :: tt(6),tp,tn
    real(rp),    allocatable            :: b(:,:)
    integer(ip)                         :: ii,jj,n

    accurate_method = .false.
    if( present(ACCURACY) ) then
       accurate_method = ACCURACY
    end if

    if( present(inva) ) then

       select case ( nsize )

       case ( 1_ip )

          det = a(1,1)
          if( abs(det) == 0.0_rp ) goto 10
          inva(1,1) = 1.0_rp/a(1,1)

       case( 2_ip )

          det = a(1,1)*a(2,2)-a(2,1)*a(1,2)
          if( abs(det) == 0.0_rp ) goto 10
          denom  = 1.0_rp/det
          inva(1,1) = a(2,2)*denom
          inva(2,2) = a(1,1)*denom
          inva(2,1) =-a(2,1)*denom
          inva(1,2) =-a(1,2)*denom

       case ( 3_ip )

          if( accurate_method ) then
             tt(1) =  a(1,1)*a(2,2)*a(3,3)
             tt(2) = -a(1,1)*a(3,2)*a(2,3)
             tt(3) = -a(1,2)*a(2,1)*a(3,3)
             tt(4) =  a(1,2)*a(3,1)*a(2,3)
             tt(5) =  a(1,3)*a(2,1)*a(3,2)
             tt(6) = -a(1,3)*a(3,1)*a(2,2)
             tp    =  sum(tt,mask=tt >= 0.0_rp)
             tn    =  sum(tt,mask=tt <  0.0_rp)
             det   =  tp-tn
          else
             t1    = a(2,2)*a(3,3) - a(3,2)*a(2,3)
             t2    =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
             t3    = a(2,1)*a(3,2) - a(3,1)*a(2,2)
             det   = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
          end if

          if( abs(det) == 0.0_rp ) goto 10
          denom     = 1.0_rp/det
          inva(1,1) = t1*denom
          inva(2,1) = t2*denom
          inva(3,1) = t3*denom
          inva(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
          inva(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
          inva(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
          inva(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
          inva(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
          inva(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

       case ( 4_ip )

          t1=    a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
               + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
               - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
          t2=  - a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
               - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
               + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
          t3=    a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
               + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
               - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
          t4=  - a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
               - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
               + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
          det= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
          if( abs(det) == 0.0_rp ) goto 10
          denom     = 1.0_rp/det
          inva(1,1) = t1*denom
          inva(2,1) = t2*denom
          inva(3,1) = t3*denom
          inva(4,1) = t4*denom
          inva(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
               &      - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
               &      + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
          inva(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
               &      + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
               &      - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
          inva(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
               &      - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
               &      + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
          inva(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
               &      + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
               &      - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
          inva(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
               &      + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
               &      - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
          inva(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
               &      - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
               &      + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
          inva(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
               &      + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
               &      - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
          inva(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
               &      - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
               &      + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
          inva(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
               &      - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
               &      + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
          inva(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
               &      + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
               &      - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
          inva(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
               &      - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
               &      + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
          inva(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
               &      + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
               &      - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

       case default

          inva(1:nsize,1:nsize) = a(1:nsize,1:nsize)
          det = 1.0_rp

          do n = 1,nsize
             t1 = inva(n,n)
             if( t1 == 0.0_rp ) then
                det = 0.0_rp
                goto 10
             end if
             do jj = 1,nsize
                inva(n,jj) = -inva(n,jj)/t1
             end do
             do ii = 1,nsize
                if ( n /= ii ) then
                   do jj = 1,nsize
                      if( n /= jj ) inva(ii,jj) = inva(ii,jj) + inva(ii,n) * inva(n,jj)
                   end do
                end if
                inva(ii,n) = inva(ii,n)/t1
             end do
             inva(n,n) = 1.0_rp/t1
          end do

       end select

    else

       if( nsize <= 4 ) then
          allocate( b(nsize,nsize) )
          b(1:nsize,1:nsize) = a(1:nsize,1:nsize)
       end if

       select case ( nsize )

       case ( 1_ip )

          det = b(1,1)
          if( abs(det) == 0.0_rp ) goto 20
          a(1,1) = 1.0_rp/b(1,1)

       case( 2_ip )

          det = b(1,1)*b(2,2)-b(2,1)*b(1,2)
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = b(2,2)*denom
          a(2,2) = b(1,1)*denom
          a(2,1) =-b(2,1)*denom
          a(1,2) =-b(1,2)*denom

       case ( 3_ip )

          t1  = b(2,2)*b(3,3) - b(3,2)*b(2,3)
          t2  =-b(2,1)*b(3,3) + b(3,1)*b(2,3)
          t3  = b(2,1)*b(3,2) - b(3,1)*b(2,2)
          det = b(1,1)*t1 + b(1,2)*t2 + b(1,3)*t3
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = t1*denom
          a(2,1) = t2*denom
          a(3,1) = t3*denom
          a(2,2) = ( b(1,1)*b(3,3) - b(3,1)*b(1,3))*denom
          a(3,2) = (-b(1,1)*b(3,2) + b(1,2)*b(3,1))*denom
          a(3,3) = ( b(1,1)*b(2,2) - b(2,1)*b(1,2))*denom
          a(1,2) = (-b(1,2)*b(3,3) + b(3,2)*b(1,3))*denom
          a(1,3) = ( b(1,2)*b(2,3) - b(2,2)*b(1,3))*denom
          a(2,3) = (-b(1,1)*b(2,3) + b(2,1)*b(1,3))*denom

       case ( 4_ip )

          t1=    b(2,2)*b(3,3)*b(4,4) + b(2,3)*b(3,4)*b(4,2)&
               + b(2,4)*b(3,2)*b(4,3) - b(2,3)*b(3,2)*b(4,4)&
               - b(2,2)*b(3,4)*b(4,3) - b(2,4)*b(3,3)*b(4,2)
          t2=  - b(2,1)*b(3,3)*b(4,4) - b(2,3)*b(3,4)*b(4,1)&
               - b(2,4)*b(3,1)*b(4,3) + b(2,4)*b(3,3)*b(4,1)&
               + b(2,3)*b(3,1)*b(4,4) + b(2,1)*b(3,4)*b(4,3)
          t3=    b(2,1)*b(3,2)*b(4,4) + b(2,2)*b(3,4)*b(4,1)&
               + b(2,4)*b(3,1)*b(4,2) - b(2,4)*b(3,2)*b(4,1)&
               - b(2,2)*b(3,1)*b(4,4) - b(2,1)*b(3,4)*b(4,2)
          t4=  - b(2,1)*b(3,2)*b(4,3) - b(2,2)*b(3,3)*b(4,1)&
               - b(2,3)*b(3,1)*b(4,2) + b(2,3)*b(3,2)*b(4,1)&
               + b(2,2)*b(3,1)*b(4,3) + b(2,1)*b(3,3)*b(4,2)
          det= b(1,1)*t1 + b(1,2)*t2 + b(1,3)*t3 + b(1,4)*t4
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = t1*denom
          a(2,1) = t2*denom
          a(3,1) = t3*denom
          a(4,1) = t4*denom
          a(1,2) =(- b(1,2)*b(3,3)*b(4,4) - b(1,3)*b(3,4)*b(4,2)&
               &   - b(1,4)*b(3,2)*b(4,3) + b(1,3)*b(3,2)*b(4,4)&
               &   + b(1,2)*b(3,4)*b(4,3) + b(1,4)*b(3,3)*b(4,2))*denom
          a(2,2) =(  b(1,1)*b(3,3)*b(4,4) + b(1,3)*b(3,4)*b(4,1)&
               &   + b(1,4)*b(3,1)*b(4,3) - b(1,4)*b(3,3)*b(4,1)&
               &   - b(1,3)*b(3,1)*b(4,4) - b(1,1)*b(3,4)*b(4,3))*denom
          a(3,2) =(- b(1,1)*b(3,2)*b(4,4) - b(1,2)*b(3,4)*b(4,1)&
               &   - b(1,4)*b(3,1)*b(4,2) + b(1,4)*b(3,2)*b(4,1)&
               &   + b(1,2)*b(3,1)*b(4,4) + b(1,1)*b(3,4)*b(4,2))*denom
          a(4,2) =(  b(1,1)*b(3,2)*b(4,3) + b(1,2)*b(3,3)*b(4,1)&
               &   + b(1,3)*b(3,1)*b(4,2) - b(1,3)*b(3,2)*b(4,1)&
               &   - b(1,2)*b(3,1)*b(4,3) - b(1,1)*b(3,3)*b(4,2))*denom
          a(1,3) =(  b(1,2)*b(2,3)*b(4,4) + b(1,3)*b(2,4)*b(4,2)&
               &   + b(1,4)*b(2,2)*b(4,3) - b(1,3)*b(2,2)*b(4,4)&
               &   - b(1,2)*b(2,4)*b(4,3) - b(1,4)*b(2,3)*b(4,2))*denom
          a(2,3) =(- b(1,1)*b(2,3)*b(4,4) - b(1,3)*b(2,4)*b(4,1)&
               &   - b(1,4)*b(2,1)*b(4,3) + b(1,4)*b(2,3)*b(4,1)&
               &   + b(1,3)*b(2,1)*b(4,4) + b(1,1)*b(2,4)*b(4,3))*denom
          a(3,3) =(  b(1,1)*b(2,2)*b(4,4) + b(1,2)*b(2,4)*b(4,1)&
               &   + b(1,4)*b(2,1)*b(4,2) - b(1,4)*b(2,2)*b(4,1)&
               &   - b(1,2)*b(2,1)*b(4,4) - b(1,1)*b(2,4)*b(4,2))*denom
          a(4,3) =(- b(1,1)*b(2,2)*b(4,3) - b(1,2)*b(2,3)*b(4,1)&
               &   - b(1,3)*b(2,1)*b(4,2) + b(1,3)*b(2,2)*b(4,1)&
               &   + b(1,2)*b(2,1)*b(4,3) + b(1,1)*b(2,3)*b(4,2))*denom
          a(1,4) =(- b(1,2)*b(2,3)*b(3,4) - b(1,3)*b(2,4)*b(3,2)&
               &   - b(1,4)*b(2,2)*b(3,3) + b(1,4)*b(2,3)*b(3,2)&
               &   + b(1,3)*b(2,2)*b(3,4) + b(1,2)*b(2,4)*b(3,3))*denom
          a(2,4) =(  b(1,1)*b(2,3)*b(3,4) + b(1,3)*b(2,4)*b(3,1)&
               &   + b(1,4)*b(2,1)*b(3,3) - b(1,4)*b(2,3)*b(3,1)&
               &   - b(1,3)*b(2,1)*b(3,4) - b(1,1)*b(2,4)*b(3,3))*denom
          a(3,4) =(- b(1,1)*b(2,2)*b(3,4) - b(1,2)*b(2,4)*b(3,1)&
               &   - b(1,4)*b(2,1)*b(3,2) + b(1,4)*b(2,2)*b(3,1)&
               &   + b(1,2)*b(2,1)*b(3,4) + b(1,1)*b(2,4)*b(3,2))*denom
          a(4,4) =(  b(1,1)*b(2,2)*b(3,3) + b(1,2)*b(2,3)*b(3,1)&
               &   + b(1,3)*b(2,1)*b(3,2) - b(1,3)*b(2,2)*b(3,1)&
               &   - b(1,2)*b(2,1)*b(3,3) - b(1,1)*b(2,3)*b(3,2))*denom

       case default

          det = 1.0_rp
          do n = 1,nsize
             t1 = a(n,n)
             if( t1 == 0.0_rp ) then
                det = 0.0_rp
                goto 20
             end if
             do jj = 1,nsize
                a(n,jj) = -a(n,jj)/t1
             end do
             do ii = 1,nsize
                if ( n /= ii ) then
                   do jj = 1,nsize
                      if( n /= jj ) a(ii,jj) = a(ii,jj) + a(ii,n) * a(n,jj)
                   end do
                end if
                a(ii,n) = a(ii,n)/t1
             end do
             a(n,n) = 1.0_rp/t1
          end do

       end select

20     if( nsize <= 4 ) deallocate(b)

    end if

10  if( present(deter) ) deter = det

  end subroutine maths_invert_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Geometrical ordering
  !> @details Uses a geometrical ordering based on coordinates
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_geometrical_sort_using_coordinates(itask,pdime,nrows,xcoor,ivin)

    integer(ip), intent(in)           :: itask
    integer(ip), intent(in)           :: pdime
    integer(ip), intent(inout)        :: nrows
    real(rp),    intent(inout)        :: xcoor(pdime,nrows)
    integer(ip), intent(inout)        :: ivin(*)
    integer(ip)                       :: boxes(3),ii,jj,kk,ll,mm
    integer(ip)                       :: nsfc,nd,d,isame,irows,nsize,jrows,jdime
    real(rp)                          :: xx
    integer(ip)                       :: idime
    integer(ip), allocatable          :: permu(:)
    integer(ip), allocatable          :: ivin_cpy(:)
    real(rp),    allocatable          :: xcoor_cpy(:,:)
    real(rp),    allocatable          :: xtmp(:)

    if( nrows < 2 ) return

    allocate(permu(nrows))
    allocate(xtmp(nrows))    
    do ii = 1,nrows
       permu(ii) = ii
    end do

    xtmp(1:nrows) = xcoor(1,1:nrows)
    call maths_heap_sort_RP1(itask,nrows,xtmp,permu)

    do idime = 1,pdime-1
       do ii = 1,nrows-1
          kk = permu(ii)
          xx = xcoor(idime,kk)
          jj = ii + 1 
          loop_jj: do while( jj <= nrows )
             if( abs(xx-xcoor(idime,permu(jj))) > epsil ) exit loop_jj
             jj = jj + 1
          end do loop_jj
          nsize = jj-ii
          if( nsize > 1 ) then
             mm = 0
             do kk = ii,jj-1
                mm = mm + 1
                ll = permu(kk)
                xtmp(mm) = xcoor(idime+1,ll)
             end do
             call maths_heap_sort_RP1(itask,nsize,xtmp,permu(ii))
          end if
       end do
    end do

    allocate(xcoor_cpy(pdime,nrows))
    allocate(ivin_cpy(nrows))
    xcoor_cpy(1:pdime,1:nrows) = xcoor(1:pdime,1:nrows)
    ivin_cpy(1:nrows)          = ivin(1:nrows)
    do ii = 1,nrows
       jj          = permu(ii)
       ivin(ii)    = ivin_cpy(jj)
       xcoor(:,ii) = xcoor_cpy(:,jj)
    end do
    deallocate(xcoor_cpy)
    deallocate(ivin_cpy)
    deallocate(permu)
    deallocate(xtmp)

  end subroutine maths_geometrical_sort_using_coordinates

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Geometrical ordering
  !> @details Uses a geometrical ordering based on SFC
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_geometrical_sort_using_sfc(itask,pdime,nrows,xcoor,ivin,NUMBER_BOXES)

    integer(ip), intent(in)           :: itask
    integer(ip), intent(in)           :: pdime
    integer(ip), intent(inout)        :: nrows
    real(rp),    intent(inout)        :: xcoor(pdime,nrows)
    integer(ip), intent(inout)        :: ivin(*)
    integer(ip), intent(in), optional :: NUMBER_BOXES
    integer(ip)                       :: boxes(3),ii,jj,kk,idime
    integer(ip)                       :: nsfc,nd,d,isame,irows
    real(rp)                          :: xleng(3)
    real(rp)                          :: xnorm,rsfc
    real(rp)                          :: xmini(3)
    real(rp)                          :: xmaxi(3)
    real(rp)                          :: toler
    integer(ip), allocatable          :: i2dto1d(:,:,:)
    integer(ip), allocatable          :: coord1d(:)

    if( nrows <= 1 ) then

       return

    else

       if( present(NUMBER_BOXES) ) then
          nsfc = NUMBER_BOXES
       else
          rsfc = real(nrows,rp)**(1.0_rp/real(pdime,rp))  ! Number of nodes in each direction
          nsfc = min(64_ip,2_ip**int(rsfc+1.0_rp,ip))
       end if
       !
       ! Sizes
       !
       toler = 1.0e-6_rp
       xleng = 0.0_rp
       xmini = 0.0_rp
       xmaxi = 0.0_rp    
       do idime = 1,pdime       
          xmini(idime) = minval(xcoor(idime,:))
          xmaxi(idime) = maxval(xcoor(idime,:))
       end do
       xleng(1:pdime) = xmaxi(1:pdime)-xmini(1:pdime)
       xnorm          = maxval(xleng)
       !
       ! Assign minimum size in each direction
       !
       xmini = xmini - 0.01_rp*xnorm
       xmaxi = xmaxi + 0.01_rp*xnorm
       !
       ! Look for 1D coordinates
       !
       allocate(coord1d(nrows))
       isame = 1
       do while( isame == 1 )

          boxes = nsfc
          if( pdime == 2 ) boxes(3) = 1
          nd = nsfc**pdime
          allocate(i2dto1d(boxes(1),boxes(2),boxes(3)))

          if( pdime == 2 ) then
             do d = 1,nd
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj)
                i2dto1d(ii,jj,1) = d
             end do
             do irows = 1,nrows
                call maths_mapping_coord_to_3d(pdime,boxes,xmaxi,xmini,xcoor(:,irows),ii,jj,kk)
                coord1d(irows) =  i2dto1d(ii,jj,1)
             end do
          else
             do d = 1,nd
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj,kk)
                i2dto1d(ii,jj,kk) = d
             end do
             do irows = 1,nrows
                call maths_mapping_coord_to_3d(pdime,boxes,xmaxi,xmini,xcoor(:,irows),ii,jj,kk)
                coord1d(irows) =  i2dto1d(ii,jj,kk)
             end do
          end if
          !
          ! Check if some nodes have the same 1D coordinate
          !
          isame = 0
          loop_irows: do irows = 1,nrows-1
             d     = coord1d(irows)
             isame = count(coord1d(irows+1:nrows)==d,KIND=ip)
             if( isame > 0 ) exit loop_irows
          end do loop_irows
          deallocate(i2dto1d)
          if( isame /= 0 ) nsfc = 2_ip * nsfc

       end do
       !
       ! Order according to 1D coordinate
       !
       call maths_heap_sort(itask,nrows,coord1d,' ',ivin) 
       deallocate(coord1d)

    end if

  end subroutine maths_geometrical_sort_using_sfc

  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    14/06/2016
  !> @brief   obtain rank table from index table
  !> @details See Figure 8.4.1 from numerical recipes (sorting)
  !>          adapted from NR - Given indx(1:n) as output from the routine indexx, this routine returns an array irank(1:n),
  !>          the corresponding table of ranks.
  !>          In order to obtain indx you have to give maths_heap_sort   on input the vector 1,2,3,....n to ivo1
  !>
  !----------------------------------------------------------------------

  subroutine maths_rank_table(n,indx,irank)
    integer(ip), intent(in)            :: n               !< Dimension
    integer(ip), intent(in)            :: indx(n)         !< index table
    integer(ip), intent(out)           :: irank(n)        !< rank table

    integer(ip)        :: j

    do  j=1,n
       irank(indx(j))=j
    end do
  end subroutine maths_rank_table


  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Computes tangent vectors
  !> @details Computes tangent vectors given a normal vector which form
  !>          an orthogonal basis
  !>
  !----------------------------------------------------------------------

  subroutine maths_local_orthonormal_basis(ndime,basis,ierro)

    integer(ip), intent(in)            :: ndime               !< Dimension
    real(rp),    intent(inout)         :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    integer(ip), intent(out), optional :: ierro               !< Error counter
    integer(ip)                        :: idime,kdime
    real(rp)                           :: xnor1,xnor2
    real(rp)                           :: xnor3,xnoma
    real(rp)                           :: exwor(3)

    if( ndime == 2 ) then

       basis(1,2) = -basis(2,1)
       basis(2,2) =  basis(1,1)

    else if( ndime == 3 ) then
       !
       ! Look for e_k such that n x e_k is maximum
       !
       xnor1    = basis(1,1) * basis(1,1)
       xnor2    = basis(2,1) * basis(2,1)
       xnor3    = basis(3,1) * basis(3,1)
       exwor(1) = xnor2 + xnor3
       exwor(2) = xnor1 + xnor3
       exwor(3) = xnor1 + xnor2
       xnoma    = 0.0_rp
       kdime    = 0
       do idime = 1,3
          if( exwor(idime) > xnoma ) then
             xnoma = exwor(idime)
             kdime = idime
          end if
       end do
       xnoma = 1.0_rp / sqrt(xnoma)
       !
       ! Set t_1 = e_k x n, first tangent vector
       !
       if( kdime == 1 ) then
          basis(1,2) =  0.0_rp
          basis(2,2) = -basis(3,1) * xnoma
          basis(3,2) =  basis(2,1) * xnoma
       else if( kdime == 2 ) then
          basis(1,2) =  basis(3,1) * xnoma
          basis(2,2) =  0.0_rp
          basis(3,2) = -basis(1,1) * xnoma
       else if( kdime == 3 ) then
          basis(1,2) = -basis(2,1) * xnoma
          basis(2,2) =  basis(1,1) * xnoma
          basis(3,2) =  0.0_rp
       else
          if( present(ierro) ) then
             ierro = ierro + 1
          else
             call runend('MATHS_LOCAL_BASIS: WRONG NORMAL')
          end if
       end if
       !
       ! Set t_2 = n x t_1, second tangent vector
       !
       basis(1,3) = basis(2,1) * basis(3,2) - basis(3,1) * basis(2,2)
       basis(2,3) = basis(3,1) * basis(1,2) - basis(1,1) * basis(3,2)
       basis(3,3) = basis(1,1) * basis(2,2) - basis(2,1) * basis(1,2)

    end if

  end subroutine maths_local_orthonormal_basis

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Rotate a vector
  !> @details Rotate a vector to a new basis, global -> local
  !>
  !>          Let ei' be the cartesian basis
  !>          e1' = [ 1 0 0 ]
  !>          e2' = [ 0 1 0 ]
  !>          e3' = [ 0 0 1 ]
  !>          and ei the new basis where
  !>          e1  = [ basis(1,1) basis(2,1) basis(3,1) ]
  !>          e2  = [ basis(1,2) basis(2,2) basis(3,2) ]
  !>          e3  = [ basis(1,3) basis(2,3) basis(3,3) ]
  !>
  !>          The base change is the following:
  !>
  !>          +-  -+   +-                                -+ +-   -+
  !>          | x1 |   | basis(1,1) basis(2,1) basis(3,1) | | x1' |
  !>          | x2 | = | basis(1,2) basis(2,2) basis(3,2) | | x2' |
  !>          | x3 |   | basis(1,3) basis(2,3) basis(3,3) | | x3' |
  !>          +-  -+   +-                                -+ +-   -+
  !>
  !----------------------------------------------------------------------

  subroutine maths_vector_to_new_basis(ndime,basis,vv,ww)

    integer(ip), intent(in)              :: ndime               !< Dimension
    real(rp),    intent(in)              :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    real(rp),    intent(inout)           :: vv(ndime)           !< Vector 1
    real(rp),    intent(inout), optional :: ww(ndime)           !< Vector 2
    real(rp)                             :: vv_tmp(ndime)
    real(rp)                             :: ww_tmp(ndime)
    integer(ip)                          :: ii,jj

    if( present(ww) ) then
       vv_tmp(1:ndime) = vv(1:ndime)
       ww_tmp(1:ndime) = ww(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          ww(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(jj,ii) * vv_tmp(jj)
             ww(ii) = ww(ii) + basis(jj,ii) * ww_tmp(jj)
          end do
       end do
    else
       vv_tmp(1:ndime) = vv(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(jj,ii) * vv_tmp(jj)
          end do
       end do
    end if

  end subroutine maths_vector_to_new_basis

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Rotate a vector
  !> @details Rotate a vector from a new basis: local -> global
  !>
  !>                         n1         t1         t2
  !>          +-   -+   +-                                -+ +-  -+
  !>          | x1' |   | basis(1,1) basis(1,2) basis(1,3) | | x1 |
  !>          | x2' | = | basis(2,1) basis(2,2) basis(2,3) | | x2 |
  !>          | x3' |   | basis(3,1) basis(3,2) basis(3,3) | | x3 |
  !>          +-   -+   +-                                -+ +-  -+
  !>
  !----------------------------------------------------------------------

  subroutine maths_vector_from_new_basis(ndime,basis,vv,ww)

    integer(ip), intent(in)              :: ndime               !< Dimension
    real(rp),    intent(in)              :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    real(rp),    intent(inout)           :: vv(ndime)           !< Vector 1
    real(rp),    intent(inout), optional :: ww(ndime)           !< Vector 2
    real(rp)                             :: vv_tmp(ndime)
    real(rp)                             :: ww_tmp(ndime)
    integer(ip)                          :: ii,jj

    if( present(ww) ) then
       vv_tmp(1:ndime) = vv(1:ndime)
       ww_tmp(1:ndime) = ww(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          ww(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(ii,jj) * vv_tmp(jj)
             ww(ii) = ww(ii) + basis(ii,jj) * ww_tmp(jj)
          end do
       end do
    else
       vv_tmp(1:ndime) = vv(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(ii,jj) * vv_tmp(jj)
          end do
       end do
    end if

  end subroutine maths_vector_from_new_basis


  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    15/02/2016
  !> @brief   Permute arrays
  !> @details Permute arrays
  !-----------------------------------------------------------------------

  subroutine maths_array_permutation_IP_1(invpr,array)
    integer(ip), intent(in) ,   pointer :: invpr(:)   !< NEW  = INVPR(OLD)
    integer(ip), intent(inout), pointer :: array(:)   !< Array
    integer(ip)                         :: ndofn
    integer(ip)                         :: ii,jj,nn
    integer(ip), allocatable            :: array_tmp(:)

    if( .not. associated(invpr) ) return
    nn = size(invpr,KIND=ip)
    if( nn <= 0 ) return
    if( nn /= size(array,KIND=ip) ) call runend('MATHS_ARRAY_PERMUTATION: WROND DIMENSIONS')

    allocate(array_tmp(nn))
    array_tmp = array
    do ii = 1,nn
       jj        = invpr(ii)
       array(jj) = array_tmp(ii)
    end do
    deallocate(array_tmp)

  end subroutine maths_array_permutation_IP_1

  subroutine maths_array_permutation_IP_2(invpr,array)
    integer(ip), intent(in) ,   pointer :: invpr(:)   !< NEW  = INVPR(OLD)
    integer(ip), intent(inout), pointer :: array(:,:) !< Array
    integer(ip)                         :: ndofn
    integer(ip)                         :: ii,jj,nn
    integer(ip), allocatable            :: array_tmp(:,:)

    if( .not. associated(invpr) ) return
    nn    = size(invpr,KIND=ip)
    ndofn = size(array,1,KIND=ip)
    if( nn /= size(array,2,KIND=ip) ) call runend('MATHS_ARRAY_PERMUTATION: WROND DIMENSIONS')

    allocate(array_tmp(ndofn,nn))
    array_tmp = array
    do ii = 1,nn
       jj          = invpr(ii)
       array(:,jj) = array_tmp(:,ii)
    end do
    deallocate(array_tmp)

  end subroutine maths_array_permutation_IP_2

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Time units
  !> @details Returns time scaling info according to a given value
  !>
  !-----------------------------------------------------------------------

  subroutine maths_time_unit(time_value,time_char,time_factor)
    real(rp),      intent(in)  :: time_value  !< Time in seconds
    character(*),  intent(out) :: time_char   !< Time unit character
    real(rp),      intent(out) :: time_factor !< Time scaling

    time_char = ' '

    if(      time_value < 1.0e-12_rp ) then  ! femto seconds
       time_factor = 1.0e15_rp
       time_char   = 'fs'
    else if( time_value < 1.0e-9_rp  ) then  ! pico seconds
       time_factor = 1.0e12_rp
       time_char   = 'ps'
    else if( time_value < 1.0e-6_rp  ) then  ! nano seconds
       time_factor = 1.0e9_rp
       time_char   = 'ns'
    else if( time_value < 1.0e-3_rp  ) then  ! micro seconds
       time_factor = 1.0e6_rp
       time_char   = 'mus'
    else if( time_value < 1.0_rp     ) then  ! milli seconds
       time_factor = 1.0e3_rp
       time_char   = 'ms '
    else if( time_value < 3600.0_rp  ) then  ! seconds
       time_factor = 1.0_rp
       time_char   = 's  '
    else                                     ! hours
       time_factor = 1.0_rp/3600.0_rp
       time_char   = 'h  '
    end if

  end subroutine maths_time_unit

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Schur complement
  !> @details Compute the Schur complement of dense matrices
  !>          A11 = A11 - A12.A22^-1.A21
  !>          A11(n1,n2), A12(n1,n3), A22(n3,n3), A21(n3,n2)
  !>
  !-----------------------------------------------------------------------

  subroutine maths_schur_complement(n1,n2,n3,A11,A12,A22,A21)

    integer(ip), intent(in)    :: n1
    integer(ip), intent(in)    :: n2
    integer(ip), intent(in)    :: n3
    real(rp),    intent(inout) :: A11(n1,n2)
    real(rp),    intent(in)    :: A12(n1,n3)
    real(rp),    intent(in)    :: A22(n3,n3)
    real(rp),    intent(in)    :: A21(n3,n2)
    real(rp)                   :: B21(n3,n2)
    integer(ip)                :: i,j,k
    real(rp)                   :: deter
    real(rp)                   :: A22inv(n3,n3)
    !
    ! Using the fact that: A_{n x m} X B_{m x p} => (AB)_ij = sum_k=1^m Aik*Bkj
    !
    if( n3 == 1 ) then
       !
       ! B21_{n3 x n2} = A22^-1_{n3 x n3} X A21_{n3 X n2}
       !
       B21 = 0.0_rp
       do j = 1,n2
          B21(1,j) = B21(1,j) + A21(1,j) / A22(1,1)
       end do
       !
       ! A11_{n1 x n2} = A11_{n1 x n2} - A12_{n1 x n3} X B21_{n3 x n2}
       !
       do i = 1,n1
          do j = 1,n2
             A11(i,j) = A11(i,j) - A12(i,1) * B21(1,j)
          end do
       end do

    else
       !
       ! Invert A22
       !
       A22inv = A22
       call maths_invert_matrix(n3,A22inv,deter)
       !
       ! B21_{n3 x n2} = A22^-1_{n3 x n3} X A21_{n3 X n2}
       !
       B21 = 0.0_rp
       do i = 1,n3
          do j = 1,n2
             do k = 1,n3
                B21(i,j) = B21(i,j) + A21(k,j) * A22inv(i,k)
             end do
          end do
       end do
       !
       ! A11_{n1 x n2} = A11_{n1 x n2} - A12_{n1 x n3} X B21_{n3 x n2}
       !
       do i = 1,n1
          do j = 1,n2
             do k = 1,n3
                A11(i,j) = A11(i,j) - A12(i,k) * B21(k,j)
             end do
          end do
       end do
    end if

  end subroutine maths_schur_complement

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Angle
  !> @details Compute the angle formed by a vecfir
  !>
  !-----------------------------------------------------------------------

  subroutine maths_angle_of_a_vector(ndime,vecto,angl1,angl2)

    integer(ip), intent(in)  :: ndime        !< Dimension
    real(rp),    intent(in)  :: vecto(ndime) !< Vector
    real(rp),    intent(out) :: angl1        !< Angle 1
    real(rp),    intent(out) :: angl2        !< Angle 2
    real(rp)                 :: vect2(3)
    real(rp)                 :: cosa,sina

    if( ndime == 2 ) then
       vect2(1:ndime) = vecto(1:ndime) / sqrt(dot_product(vecto,vecto)+epsil)
       cosa  = vect2(1)
       sina  = vect2(2)
       angl1 = acos(cosa)*180.0_rp/pi
       if( sina < 0.0_rp ) angl1 = -angl1
    else
       angl1 = 0.0_rp
       angl2 = 0.0_rp
    end if

  end subroutine maths_angle_of_a_vector

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    15/06/2016
  !> @brief   economy size qrdecomposition with household reflexions for rectangular matrices
  !> @details compute a QR=A decomposition were Q is orthogonal and R upper triangular
  !>
  ! This subroutine decomposes the matrix A(m,n)
  ! in two matrices Q and R, where Q is orthogonal
  ! and R is upper triangular.
  !
  ! as Q is orthogonal Q^T*Q=I so  the inversion of
  ! this matrix is trivial.
  !
  ! Also:
  !       Q*R=A
  !
  ! In the original version of the algorithm
  ! Q(m,m) as it is orthogonal and R(m,n) for
  ! consistency in the size of the matrix.
  !
  ! But in the economy-size version, as R is upper
  ! triangular, you can build an R(n,n) matrix
  ! (because you know that the rest of the elements
  ! are zero and an Q(m,n), because the columns
  ! n->m are going to be multiplied by zero. In this
  ! algorithm A=Q*R holds.
  !                 n
  !            +--+--+--+
  !            |  |  |  |
  !            +--+--+--+
  !            |  |  |  |
  !  A(m,n)=   +--+--+--+ m
  !            |  |  |  |
  !            +--+--+--+
  !            |  |  |  |
  !            +--+--+--+
  !
  !
  !                n
  !            +--+--+--+
  !            |  |  |  |                   n
  !            +--+--+--+               +--+--+--+
  !            |  |  |  |               |  |  |  |
  !  Q(m,n)=   +--+--+--+ m             +--+--+--+
  !            |  |  |  |      R(m,n)=  |  |  |  | n
  !            +--+--+--+               +--+--+--+
  !            |  |  |  |               |  |  |  |
  !            +--+--+--+               +--+--+--+
  !
  !                       Q(m,n)*R(n,n) = A(m,n)
  !
  !
  !
  !
  !----------------------------------------------------
  subroutine maths_economy_qrhousehold(A,Q,R,mmax,nmax)

    implicit none

    real(rp), intent (in)               :: A(:,:)
    real(rp), intent (out)              :: Q(:,:)
    real(rp), intent (out)              :: R(:,:)
    integer(ip), intent(in)             :: mmax, nmax !! maximum values of the matrix to decompose
    !  integer(ip), intent(in)             :: m, n       !! actuall full values of the matrix
    integer(ip)                         :: i_col, i, j, k
    real(rp)                            :: norm_column, re_aux
    real(rp)                            :: v(mmax,nmax)
    real(rp)                            :: A_aux(mmax,nmax)

    !  m=size(A,1,KIND=ip)
    !  n=size(A,2,KIND=ip)
    !  if( ( m .lt. 0)             .or. &       ! A is pointing to memory
    !      ( n .lt. 0)             .or. &       !
    !      ( size(Q,1,KIND=ip) .lt. mmax)  .or. &       ! Q should be as big as mmax*mmax
    !      ( size(Q,2,KIND=ip) .lt. mmax)  .or. &       !
    !      ( size(R,1,KIND=ip) .lt. mmax)  .or. &       ! R shold be as big as mmax*nmax
    !      ( size(R,2,KIND=ip) .lt. nmax)  .or. &       !
    !      ( mmax .gt. size(A,1,KIND=ip) ) .or. &       ! mmax is always smaller than m
    !      ( nmax .gt. size(A,2,KIND=ip) ) .or. ) then  ! nmax is always smaller than n
    !      call runend('MATRIX DECOMPOSITION WRONG DIMENSIONS')
    !  endif

    Q=0.0_rp
    A_aux=A


    columns: do i_col=1, nmax
       ! a_i=A(i_col,i_col:mmax)
       v(:,i_col)=0.0_rp
       v(i_col:mmax,i_col) = A_aux(i_col:mmax,i_col)


       ! ||alpha|| = sqrt(sum(A(i)^2))
       !
       norm_column = 0.0_rp
       do i=i_col,mmax
          norm_column =norm_column + v(i,i_col) * v(i,i_col)
       enddo
       norm_column = sqrt(norm_column)

       ! u=a_i-alpha*e_1
       !
       !! IMPORTANT THING HERE. Some algoritmhs (i.e. matlab)
       !! uses a_i + alpha* e_i. This algorithm also converges
       !! but slightly different. The minus sign has been
       !! left here trivially.
       v(i_col,i_col) = v(i_col,i_col) - norm_column

       ! ||u|| = sqrt(sum(u(i)^2))
       !
       norm_column = 0.0_rp
       do i=i_col,mmax
          norm_column =norm_column + v(i,i_col) * v(i,i_col)
       enddo
       norm_column = sqrt(norm_column)

       ! v=u/||u||
       !
       do i=i_col,mmax
          v(i,i_col)=v(i,i_col)/norm_column
       enddo

       ! Q_i=I-2*v*v^T
       !
       ! obtained in function Q_mat(v,i_ini,i,j)

       ! A_i+1=Q_i * A
       !
       if(.not.(i_col.eq.nmax))then

          ! Tener la matrix A_i es indispensable no se puede tener solo
          ! el vector a_i ya que para el paso siguiente necesitaremos
          ! el a_ip1 y para el siguiente el a_ip2 y asi sucesivamente
          ! Por lo que se ha de tener toda la matrix A

          call vQ_times_Qaux(v(:,i_col),i_col,A_aux)

       endif

    enddo columns

    ! Q=Q_1*Q_2*Q_3*....Q_n
    !

    Q=0.0_rp
    Q=obtain_Q_mat(v(:,nmax),nmax,mmax)

    do i_col=nmax-1,1,-1
       call vQ_times_Qaux(v(:,i_col),i_col,Q)
    enddo

    ! A=Q*R  -> R=Q^T*A
    !
    !

    do i=1,nmax !swipe in rows
       do j=1,nmax
          re_aux=0.0_rp
          do k=1,mmax !swipe in columns
             re_aux = re_aux + Q(k,i) * A(k,j) !contract columns. Q has inversed indexes because it is the transpost
          enddo
          R(i,j)=re_aux
       enddo
    enddo


  contains
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    function Q_mat(v,i_ini,i,j) result(q_ij)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini,i,j
      real(rp)                    :: q_ij

      if (i .lt. i_ini) then
         if(i .eq. j) then
            q_ij = 1.0_rp
         else
            q_ij = 0.0_rp
         endif
      else
         if(i .eq. j) then
            q_ij = 1.0_rp - 2.0_rp*v(i)*v(j)
         else
            q_ij = -2.0_rp*v(i)*v(j)
         endif
      endif

    end function Q_mat
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    function obtain_Q_mat(v,i_ini,max_col) result(Q)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini, max_col
      integer(ip)                 :: i, j, m
      real(rp)                    :: Q(size(v,1,KIND=ip),size(v,1,KIND=ip))

      m=size(v,KIND=ip)

      if ( (max_col .gt. m)     .or. &
           (max_col .lt. i_ini) ) then
         ! call runend('mod_maths: problem with obtain_Q_mat. Wrong max column number')
         write(6,*) 'OBTAIN_Q_MAT ERROR: MAXIMUM COLUMN NUMBER PROBLEM '
      endif
      do i=1,m
         do j=1,max_col
            Q(i,j) = Q_mat(v,i_ini,i,j)
         enddo
      enddo

    end function obtain_Q_mat

    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    subroutine vQ_times_Qaux(v,i_ini,Qaux)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini
      real(rp), intent(inout)     :: Qaux(:,:)
      real(rp), allocatable       :: mat_aux(:,:)
      real(rp)                    :: re_aux
      integer(ip)                 :: i, j, mv,mQ,nQ

      mv=size(v,KIND=ip)
      mQ=size(Qaux,1,KIND=ip)
      nQ=size(Qaux,2,KIND=ip)

      if ( (mv .ne. mQ ) ) write(6,*) 'vQ_times_Qaux: wrong dimensions' !call runend(mod_maths: vQ_times_Qaux: wrong dimensions')

      allocate(mat_aux(mQ,nQ))


      do i=1,mv ! NWET swipe in rows
         do j=1,nQ
            re_aux=0.0_rp
            do k=1,mv !swipe in columns
               re_aux = re_aux + Q_mat(v, i_ini, i, k) * Qaux(k,j) !contract columns
            enddo
            mat_aux(i, j) = re_aux
         enddo
      enddo


      Qaux=mat_aux

      deallocate(mat_aux)

    end subroutine vQ_times_Qaux

    !--------------------------------------
    ! Function that multiplies two Q
    ! special matrices
    !--------------------------------------
    function vQ_times_vQ(v1, i_ini_1, v2, i_ini_2) result(Q)
      implicit none

      real(rp), intent(in)        :: v1(:), v2(:)
      integer(ip), intent(in)     :: i_ini_1, i_ini_2
      real(rp)                    :: Q(size(v1,1,KIND=ip),size(v1,1,KIND=ip))
      real(rp)                    :: rea_aux,part_V1, part_v2
      integer(ip)                 :: m,i,j,k

      m=size(v1,1,KIND=ip)
      Q=0.0_rp



      do i=1,m !swipe in rows
         do j=1,m
            rea_aux=0.0_rp
            do k=1,m !swipe in columns
               rea_aux = rea_aux + Q_mat(v1, i_ini_1, i, k) * Q_mat(v2, i_ini_2, k, j)!contract columns
            enddo
            Q(i,j)=rea_aux
         enddo
      enddo


    end function vQ_times_vQ
  end subroutine maths_economy_qrhousehold

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   matrix multiplication
  !> @details Multiplies two rectangular matrices
  !>
  !
  !
  !----------------------------------------------------
  subroutine maths_matrix_multiplication(A,B,C)
    !
    ! A(m,n)*B(n,p)=C(m,p)
    !
    real(rp), intent(in)            :: A(:,:), B(:,:)
    real(rp), intent(out)           :: C(:,:)
    integer(ip)                     :: m,n,r
    integer(ip)                     :: i,j,k
    real(rp)                        :: acum

    m=size(A,1,KIND=ip) ! Rows
    n=size(B,2,KIND=ip) ! Columns
    r=size(A,2,KIND=ip) ! reduce

    if( r .ne. size(B,1,KIND=ip) ) then
       CALL RUNEND('MATRIX_MULTIPLIUCATION: WRONG DIMENSIONS IN MATRIX')
    endif

    C=0.0_rp

    do i=1,m !swipe in rows
       do j=1,n
          acum=0.0_rp
          do k=1,r !swipe in columns
             acum=acum+A(i,k)*B(k,j) !contract columns
          enddo
          C(i,j)=acum
       enddo
    enddo

  end subroutine maths_matrix_multiplication

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   matrix multiplication
  !> @details Multiplies two rectangular matrices
  !>
  !
  !
  !----------------------------------------------------

  subroutine maths_matrix_vector_multiplication(A,B,C,maxn)
    !
    ! A(m,n)*B(n,1)=C(m,1)
    !
    real(rp),    intent(in)  :: A(:,:)
    real(rp),    intent(in)  :: B(:)
    real(rp),    intent(out) :: C(:)
    integer(ip), intent(in)  :: maxn
    integer(ip)              :: m,n
    integer(ip)              :: i,j
    real(rp)                 :: acum

    m=size(A,1,KIND=ip) ! Rows
    n=size(A,2,KIND=ip) ! Columns

    if(size(C,1,KIND=ip) == 0_ip) return

    if( (n /= size(B,1,KIND=ip)) .or. (maxn > n) ) then
       CALL RUNEND('MATRIX_VECTOR_MULTIPLIUCATION: WRONG DIMENSIONS IN MATRIX')
    endif

    C = 0.0_rp

    do i = 1,m !swipe in rows
       acum = 0.0_rp
       do j = 1,maxn
          acum = acum+A(i,j)*B(j) !contract columns
       enddo
       C(i) = acum
    enddo

  end subroutine maths_matrix_vector_multiplication

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   Back substitution
  !> @details backsubtitutes the sistem Ax=B were A is upper triangular
  !>
  !
  !
  !----------------------------------------------------
  subroutine maths_backsu(A,x,b, nmax)

    real(rp), intent(in)        :: A(:,:)
    real(rp), intent(in)        :: b(:)
    integer(ip), intent(in)     :: nmax !bounds of the operation
    real(rp), intent(out)       :: x(:)
    integer(ip)                 :: i,j
    real(rp)                    :: aux


    x=0.0_rp


    if (size(b,1,KIND=ip).eq.0_ip) return

    do i=nmax,1,-1

       aux=0.0_rp

       do j=nmax,i+1,-1
          aux=aux+A(i,j)*x(j)
       enddo

       if(A(i,i).eq. 0.0_rp) call runend('mod_maths backsubstitution: diagonal equal to zero in backsubstitution')

       x(i)=(b(i)-aux)/A(i,i)

    enddo

  end subroutine maths_backsu

  !-----------------------------------------------------------
  !>
  !> @author  J.C. Cajas
  !> @date    22/06/2016
  !> @brief   Outer product
  !> @details Outer product between two vectors A(m) and B(n)
  !> @        result is C(mxn)
  !>
  !>       _    _
  !>      |  a1  |
  !> A =  |  a2  |
  !>      |   .  |
  !>      |   .  |
  !>      |   .  |
  !>      |_ am _|
  !>
  !> B = [b1,b2,...,bn]
  !>
  !>           --                      --
  !>          | a1 b1   a1 b2  ... a1 bn |
  !>          | a2 b1   a2 b2      a2 bn |
  !>          |  .                       |
  !> A xo B = |  .                       |
  !>          |  .                       |
  !>          | am b1     ...      am bn |
  !>           --                      --
  !------------------------------------------------------------
  subroutine maths_outer_product(A,B,C)
    real(rp), intent(in) ,   pointer :: A(:),B(:)   !< Input vectors
    real(rp), intent(inout), pointer :: C(:,:)      !< Output matrix
    integer(ip)                      :: mm,nn       !< Dimensions of the vectors
    integer(ip)                      :: ii,jj

    mm = size(A,KIND=ip)
    nn = size(B,KIND=ip)

    if( mm /= size(C,1_ip,KIND=ip) .or. nn /= size(C,2_ip,KIND=ip) )call runend('mod_maths, outer product: Wrong dimensions')

    do jj = 1_ip, nn
       do ii = 1_ip, mm
          C(ii,jj) = A(ii) * B(jj)
       end do
    end do


  end subroutine maths_outer_product

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   L2 norm of an array
  !> @details L2 norm of an array
  !
  !------------------------------------------------------------------------

  function maths_norm2(ndime,array)
    integer(ip), intent(in) :: ndime
    real(rp),    intent(in) :: array(ndime)
    real(rp)                :: maths_norm2

    maths_norm2 = sqrt(dot_product(array,array))

  end function maths_norm2

  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Partition of bin boxes by means of space filling curves
  !> @details npart can be a real number
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_part(ndime,boxes,weigh,npart,parts,fronw,irank,nrank,ityp)

    integer(ip),          intent(in)     :: ndime         !< Dimension of problem
    integer(ip),          intent(in)     :: boxes(ndime)  !< # boxes in each dimension
    real(rp),    pointer, intent(in)     :: weigh(:)      !< boxes weights
    real(rp),             intent(in)     :: npart         !< # partitions
    integer(ip), pointer, intent(inout)  :: parts(:)      !< Partition to which each box is assigned
    real(rp),    optional,intent(in)     :: fronw         !< weigh of first partition
    integer(ip), optional,intent(in)     :: irank         !> first output partitions rank
    integer(ip), optional,intent(in)     :: nrank         !< # ranks used in the partition
    integer(2), optional,intent(inout)   :: ityp          ! orientation type

    integer(ip)                       :: irank_,nrank_
    integer(ip)                       :: aux,icont,ipart
    integer(ip)                       :: nbox,ibox,frank
    integer(ip)                       :: hcoor(ndime),coo1d
    integer(2)                        :: auxtyp,auxtyp2
    real(rp)                          :: fronw_
    real(rp)                          :: sump,iweig,wsum
    real(rp)                          :: wprev,waux,wobj

    !
    ! Check arguments
    !

    if(ndime/= 2_ip .and. ndime/= 3_ip)&
         call runend("maths_sfc_part: wrong ndime value")
    if(boxes(1) <= 0_ip) call runend("maths_sfc_part: boxes(1)<=0")
    nbox = boxes(1)
    do icont =2_ip,ndime
       if(boxes(icont) /= boxes(icont-1_ip)) call runend("maths_sfc_part: different boxes sizes")
       nbox = nbox * boxes(1)
    enddo
    if(iand(nbox,nbox-1_ip)/=0)&
         call runend("math_sfc_part: nbox is not power of two")
    if( size(weigh,KIND=ip) /= nbox )&
         call runend("math_sfc_part: Incongruent number of boxes and weigh dim")

    !
    ! Check or optional inputs
    !

    if(present(fronw))then
       if(fronw <= 0_rp .or. fronw> 1_rp) then
          call runend("maths_sfc_part: fronw <= 0 or fronw > 1.0!")
       endif
       fronw_=fronw
    else
       fronw_ = 1_rp
    endif

    if(present(irank)) then
       if(irank<=0_ip) call runend("maths_sfc_part: irank<=1")
       irank_=irank
    else
       irank_ = 1_ip
    endif

    aux = npart-int(fronw_,KIND=ip)
    if( npart-int(fronw_,KIND=ip)-aux > 0_ip) aux = aux + 1_ip
    aux = aux + 1
    if(present(nrank)) then
       if(nrank > aux .or. nrank < aux-1) then
          call runend("maths_sfc_part: nrank out of rank")
       endif
       nrank_=nrank
    else
       nrank_ = aux
    endif

    if(present(ityp)) then
       auxtyp = ityp
    else
       auxtyp = 0_2
    endif

    !
    ! Allocate parts
    !

    nullify(parts)
    call memory_alloca(memor,'parts',"maths_sfc_part",parts,nbox)

    !
    ! calculate sum of weights and initial obj. weight
    !
    frank = nrank_ + irank_ - 1_ip
    wsum = 0_rp
    do ibox=1,nbox
       wsum=wsum+weigh(ibox)
    enddo
    wobj = (wsum/npart)*fronw_
    sump=fronw_
    !
    ! assign a partition to each box
    !
    iweig = 0_rp
    wprev = 0_rp
    ipart = irank_

    do ibox=1_ip,nbox

       if(ndime==2_ip) then
          auxtyp2=auxtyp
          call maths_sfc_d2xy_tab(boxes(1),ibox,hcoor(1),hcoor(2),auxtyp2)
          coo1d = (hcoor(2)-1)*boxes(1)+hcoor(1)
       else
          auxtyp2=auxtyp
          call maths_sfc_d2xyz_tab(boxes(1),ibox,hcoor(1),hcoor(2),hcoor(3),auxtyp2)
          coo1d = (hcoor(3)-1)*boxes(2)*boxes(1)+(hcoor(2)-1)*boxes(1)+hcoor(1)
       endif
       waux  = weigh(coo1d)

       if( abs(iweig+waux-wobj) <= abs(iweig-wobj) .or. frank == ipart ) then
          parts(coo1d)=ipart
          iweig = iweig+waux
          wprev = wprev+waux
       else
          wobj  = (wsum-wprev)/real(npart-sump,rp)
          ipart = ipart + 1_ip
          parts(coo1d) = ipart
          iweig = waux
          wprev = wprev + waux
          sump = sump + 1_ip
       endif

    enddo

  end subroutine maths_sfc_part
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Partition of bin boxes by means of space filling curves
  !> @details npart can be a real number
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_part_2(ndime,boxes,weigh,parts,ityp,irank,nrank,wrank)

    integer(ip),          intent(in)     :: ndime         !< Dimension of problem
    integer(ip),          intent(in)     :: boxes(ndime)  !< # boxes in each dimension
    real(rp),    pointer, intent(in)     :: weigh(:)      !< boxes weights
    integer(ip), pointer, intent(inout)  :: parts(:)      !< Partition to which each box is assigned

    integer(2),  optional,intent(inout)  :: ityp          ! orientation type
    integer(ip), optional,intent(in)     :: irank         !> first output partitions rank
    integer(ip), optional,intent(in)     :: nrank         !< # ranks used in the partition
    real(rp),    pointer,optional,intent(in)  :: wrank(:) !< relative weight per rank  

    integer(ip)                       :: irank_,nrank_
    integer(ip)                       :: aux,icont,ipart,ilpart
    integer(ip)                       :: nbox,ibox,frank
    integer(ip)                       :: hcoor(ndime),coo1d
    integer(2)                        :: auxtyp,auxtyp2
    real(rp)                          :: sump,iweig,wsum
    real(rp)                          :: wprev,waux,wobj
    real(rp), pointer                 :: wrank_(:)
    real(rp)                          :: npart         !< # partitions

    character(100), PARAMETER :: vacal = "maths_sfc_part_2"
    !
    ! TODO: I supose by now that all optional arguments are present...
    !

    !
    ! Check arguments
    !

    if(ndime/= 2_ip .and. ndime/= 3_ip)&
         call runend("maths_sfc_part: wrong ndime value")
    if(boxes(1) <= 0_ip) call runend("maths_sfc_part: boxes(1)<=0")
    nbox = boxes(1)
    do icont =2_ip,ndime
       if(boxes(icont) /= boxes(icont-1_ip)) call runend("maths_sfc_part: different boxes sizes")
       nbox = nbox * boxes(1)
    enddo
    if(iand(nbox,nbox-1_ip)/=0)&
         call runend("math_sfc_part: nbox is not power of two")
    if( size(weigh,KIND=ip) /= nbox )&
         call runend("math_sfc_part: Incongruent number of boxes and weigh dim")

    ! rick: check combinations of elements wich have to be present...

    !
    ! Check optional inputs
    !

    if(present(ityp)) then
       auxtyp = ityp
    else
       auxtyp = 0_2
    endif

    if(present(irank)) then
       if(irank<=0_ip) call runend("maths_sfc_part: irank<=1")
       irank_=irank
    else
       irank_ = 1_ip
    endif

    if(present(nrank)) then
       nrank_=nrank
       nullify(wrank_) 
       call memory_alloca(memor,'wrank_',vacal,wrank_,nrank_)
       wrank_(1:nrank)=wrank(1:nrank)
    endif

    !
    ! Allocate parts
    !

    nullify(parts)
    call memory_alloca(memor,'parts',"maths_sfc_part",parts,nbox)

    !
    ! calculate sum of weights and initial obj. weight
    !
    npart = sum(wrank_)
    frank = nrank_ + irank_ - 1_ip  !! final rank?
    wsum = 0_rp
    do ibox=1,nbox
       wsum=wsum+weigh(ibox)
    enddo

    wobj = (wsum/npart)*wrank_(1)
    sump=wrank_(1)

    ! assign a partition to each box
    !
    iweig = 0_rp
    wprev = 0_rp
    ipart = irank_
    ilpart = 1

    do ibox=1_ip,nbox

       if(ndime==2_ip) then
          auxtyp2=auxtyp
          call maths_sfc_d2xy_tab(boxes(1),ibox,hcoor(1),hcoor(2),auxtyp2)
          coo1d = (hcoor(2)-1)*boxes(1)+hcoor(1)
       else
          auxtyp2=auxtyp
          call maths_sfc_d2xyz_tab(boxes(1),ibox,hcoor(1),hcoor(2),hcoor(3),auxtyp2)
          coo1d = (hcoor(3)-1)*boxes(2)*boxes(1)+(hcoor(2)-1)*boxes(1)+hcoor(1)
       endif
       waux  = weigh(coo1d)

       if( abs(iweig+waux-wobj) <= abs(iweig-wobj) .or. frank == ipart ) then
          parts(coo1d)=ipart
          iweig = iweig+waux
          wprev = wprev+waux
       else
          ilpart = ilpart + 1_ip
          wobj  = ((wsum-wprev)/real(npart-sump,rp))*wrank_(ilpart)
          ipart = ipart + 1_ip
          parts(coo1d) = ipart
          iweig = waux
          wprev = wprev + waux
          sump = sump + wrank_(ilpart)
       endif

    enddo

    !rick: see what happends when no present
    call memory_deallo(memor,'wrank_',vacal,wrank_)


  end subroutine maths_sfc_part_2
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Parallel Partition of bin boxes by means of space filling curves
  !> @details ...
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_par_part(ndime,loc_boxes,loc_weigh,npar,proc_boxes,mypro,proc_wdist,loc_part)

    integer(ip),          intent(in)     :: ndime             !< Dimension of problem
    integer(ip),          intent(in)     :: loc_boxes(ndime)  !< # boxes in each dimension (local)
    real(rp),    pointer, intent(in)     :: loc_weigh(:)      !< boxes weights (local)
    integer(ip),          intent(in)     :: npar              !< # partitions (global)
    integer(ip),          intent(in)     :: proc_boxes(ndime) !< # processes performing partition
    integer(ip),          intent(in)     :: mypro             !< # my process rank
    real(rp), pointer,    intent(in)     :: proc_wdist(:)     !< weig accumulated by each process
    integer(ip), pointer, intent(inout)  :: loc_part(:)       !< Partition to which each box is assigned (local)

    real(rp),    pointer                 :: loc_npar(:)
    real(rp),    pointer                 :: fronw(:)
    integer(ip), pointer                 :: irank(:)
    integer(ip), pointer                 :: reord(:)
    integer(ip), pointer                 :: invor(:)
    integer(2),  pointer                 :: types(:)
    integer(ip)                          :: coo1d
    integer(2)                           :: auxty
    real(rp)                             :: wsum,npsum
    integer(ip)                          :: iproc,x,y,z
    integer(ip)                          :: nproc,idime
    integer(ip)                          :: invmy

    character(100), PARAMETER :: vacal = "maths_sfc_par_part"

    !
    ! Check proc_boxes input
    !
    do idime=2_ip,ndime
       if(proc_boxes(idime)/=proc_boxes(idime-1_ip))&
            call runend("maths_sfc_par_part: wrong proc_boxes")
    enddo
    nproc = proc_boxes(1_ip)**(ndime)
    if(iand(nproc,nproc-1)/=0)&
         call runend("math_sfc_part: boxes(1) is not a power of 2")
    !
    ! Eval the new order of proc. boxes and its orientation type
    !
    nullify(loc_npar,fronw,irank,reord,invor,types)
    call memory_alloca(memor,'loc_npar',vacal,loc_npar,nproc)
    call memory_alloca(memor,'fronw',vacal,fronw,nproc)
    call memory_alloca(memor,'irank',vacal,irank,nproc)
    call memory_alloca(memor,'reord',vacal,reord,nproc)
    call memory_alloca(memor,'invor',vacal,invor,nproc)
    allocate(types(nproc))
    wsum  = 0_rp
    do iproc=1,nproc
       wsum = wsum + proc_wdist(iproc)
       auxty = 0_ip
       if(ndime==2_ip) then
          call maths_sfc_d2xy_tab(proc_boxes(1),iproc,x,y,auxty)
          coo1d=(y-1)*proc_boxes(1)+x
       else
          call maths_sfc_d2xyz_tab(proc_boxes(1),iproc,x,y,z,auxty)
          coo1d=(z-1)*proc_boxes(1)*proc_boxes(1)+(y-1)*proc_boxes(1)+x
       endif
       reord(iproc) = coo1d
       types(iproc) = auxty
       invor(coo1d) = iproc
    enddo

    if(wsum == 0_rp) return
    !
    ! Perform partition
    !
    npsum = 0_rp
    do iproc=1,nproc
       loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
       irank(iproc) = int(npsum) + 1_ip
       fronw(iproc) = real(irank(iproc),rp) - npsum
       npsum = npsum + loc_npar(iproc)
    enddo

    invmy = invor(mypro)
    if(loc_npar(invmy) /= 0_rp) then
       call maths_sfc_part(ndime,loc_boxes,loc_weigh,loc_npar(invmy),loc_part,&
            fronw=fronw(invmy),irank=irank(invmy),ityp=types(invmy))
    endif
    call memory_deallo(memor,'loc_npar',vacal,loc_npar)
    call memory_deallo(memor,'fronw',vacal,fronw)
    call memory_deallo(memor,'irank',vacal,irank)
    call memory_deallo(memor,'reord',vacal,reord)
    call memory_deallo(memor,'invor',vacal,invor)
    deallocate(types)

  end subroutine maths_sfc_par_part
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Parallel Partition of bin boxes by means of space filling curves
  !> @details ...
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_par_part_2(ndime,loc_boxes,loc_weigh,npar,proc_boxes,&
       mypro,proc_wdist,loc_part,corr_dist,totalw)

    integer(ip),          intent(in)     :: ndime                 !< Dimension of problem
    integer(ip),          intent(in)     :: loc_boxes(ndime)      !< # boxes in each dimension (local)
    real(rp),    pointer, intent(in)     :: loc_weigh(:)          !< boxes weights (local)
    integer(ip),          intent(in)     :: npar                  !< # partitions (global)
    integer(ip),          intent(in)     :: proc_boxes(ndime)     !< # processes performing partition
    integer(ip),          intent(in)     :: mypro                 !< # my process rank
    real(rp), pointer,    intent(in)     :: proc_wdist(:)         !< weig accumulated by each process
    integer(ip), pointer, intent(inout)  :: loc_part(:)           !< Partition to which each box is assigned (local)
    real(rp), pointer,    intent(inout), optional :: corr_dist(:) !< correction factors for the distribution
    real(rp),             intent(in),    optional :: totalw

    real(rp),    pointer                 :: loc_npar(:)
    real(rp),    pointer                 :: dist_corr(:)
    integer(ip), pointer                 :: irank(:)
    integer(ip), pointer                 :: reord(:)
    integer(ip), pointer                 :: invor(:)
    integer(2),  pointer                 :: types(:)
    integer(ip)                          :: coo1d
    integer(2)                           :: auxty
    real(rp)                             :: wsum
    integer(ip)                          :: iproc,x,y,z
    integer(ip)                          :: nproc,idime
    integer(ip)                          :: nrank
    real(rp),    pointer                 :: wrank(:)
    logical                              :: auxlg    
    integer(ip)                          :: iaux

    character(100), PARAMETER :: vacal = "maths_sfc_par_part_2"

    !
    ! Check proc_boxes input
    !
    do idime=2_ip,ndime
       if(proc_boxes(idime)/=proc_boxes(idime-1_ip))&
            call runend("maths_sfc_par_part: wrong proc_boxes")
    enddo
    nproc = proc_boxes(1_ip)**(ndime)
    if(iand(nproc,nproc-1)/=0)&
         call runend("math_sfc_part: boxes(1) is not a power of 2")
    !
    ! Eval the new order of proc. boxes and its orientation type
    !
    nullify(loc_npar,irank,reord,invor,types,dist_corr)
    call memory_alloca(memor,'loc_npar' ,vacal,loc_npar,nproc)
    call memory_alloca(memor,'dist_corr',vacal,dist_corr,npar)
    call memory_alloca(memor,'irank'    ,vacal,irank,nproc)
    call memory_alloca(memor,'reord'    ,vacal,reord,nproc)
    call memory_alloca(memor,'invor'    ,vacal,invor,nproc)
    allocate(types(nproc))
    if(.NOT. present(corr_dist)) then
       dist_corr(:) = 1_rp
    else
       iaux = min(size(dist_corr,KIND=ip),size(corr_dist,KIND=ip))
       dist_corr(1:iaux) = corr_dist(1:iaux)
    endif
    !
    ! Ordering of partition processess acording to coarse SFC
    !
    wsum  = 0_rp
    do iproc=1,nproc
       wsum = wsum + proc_wdist(iproc)
       auxty = 0_ip
       if(ndime==2_ip) then
          call maths_sfc_d2xy_tab(proc_boxes(1),iproc,x,y,auxty)
          coo1d=(y-1)*proc_boxes(1)+x
       else
          call maths_sfc_d2xyz_tab(proc_boxes(1),iproc,x,y,z,auxty)
          coo1d=(z-1)*proc_boxes(1)*proc_boxes(1)+(y-1)*proc_boxes(1)+x
       endif
       reord(iproc) = coo1d
       types(iproc) = auxty
       invor(coo1d) = iproc
    enddo
    if(present(totalw)) wsum=totalw;
    if(wsum == 0_rp) return
    !
    ! Perform partition
    !
    iaux = 1_ip
    nullify(wrank)
    call memory_alloca(memor,'wrank',vacal,wrank,npar)

    do iproc=1,nproc

       nrank = 0
       !loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
       loc_npar(iproc) = (proc_wdist(iproc)/wsum) * real(npar,rp)
       if(loc_npar(iproc) > 0_rp) then

          irank(iproc) = iaux
          auxlg = .true.

          do while(auxlg .and. nrank < npar .and. iaux <= npar)
             loc_npar(iproc) = loc_npar(iproc) - dist_corr(iaux)            
             if(loc_npar(iproc) < 0) then 
                auxlg = .false.
                wrank(nrank+1)= loc_npar(iproc) + dist_corr(iaux)
                dist_corr(iaux) = -loc_npar(iproc)
             else
                wrank(nrank+1)=dist_corr(iaux)
                iaux = iaux + 1
             endif
             nrank = nrank + 1    ! si queda 0 el resido petara - repassar aquest cas
             if(loc_npar(iproc) == 0) auxlg = .false.
          end do
          if(iproc == nproc) then
             nrank = npar-irank(iproc)+1
          endif

          !loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
          loc_npar(iproc) = (proc_wdist(iproc)/wsum) * real(npar,rp)

          !if(invor(mypro) == iproc) then
          if(mypro == iproc) then
             call maths_sfc_part_2(ndime,loc_boxes,loc_weigh,loc_part,types(iproc),&
                  irank(iproc),nrank,wrank)
             exit
          endif
       endif
    enddo
    !
    ! deallocate pointers
    !
    call memory_deallo(memor,'wrank',vacal,wrank)
    call memory_deallo(memor,'loc_npar',vacal,loc_npar)
    call memory_deallo(memor,'dist_corr',vacal,dist_corr)
    call memory_deallo(memor,'irank',vacal,irank)
    call memory_deallo(memor,'reord',vacal,reord)
    call memory_deallo(memor,'invor',vacal,invor)
    deallocate(types)

  end subroutine maths_sfc_par_part_2
  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 2D cartesian coordinates
  !> @details Kernel based on bit operations
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n)
  !           x,y: 1<=x,y<=n
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_d2xy_bit(n,d,x,y)

    implicit none
    integer(ip), intent(in)  :: n,d
    integer(ip), intent(out) :: x,y
    integer(ip)              :: rx,ry,s,t

    x = 0_ip
    y = 0_ip
    t = d-1_ip
    s = 1_ip

    do while( s<n )
       rx = mod(t/2,2_ip)
       if( rx==0 ) then
          ry = mod(t,2_ip)
       else
          ry = mod(ieor(t,rx),2_ip)
       endif
       call maths_sfc_rot_bit(s,x,y,rx,ry)
       x = x+s*rx
       y = y+s*ry
       t = t/4
       s=s*2
    end do

    x=x+1_ip
    y=y+1_ip

  end subroutine maths_sfc_d2xy_bit
  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Rotates and flips a quadrant appropritely
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_rot_bit(n,x,y,rx,ry)

    implicit none
    integer(ip), intent(in)    :: n,rx,ry
    integer(ip), intent(inout) :: x,y
    integer(ip)                :: t

    if(ry==0) then

       !Reflect
       if( rx == 1_ip ) then
          x = n-1_ip-x
          y = n-1_ip-y
       endif

       !Flip
       t = x
       x = y
       y = t
    endif

  end subroutine maths_sfc_rot_bit
  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 2D cartesian coordinates
  !> @details Kernel based on the table of recursive divisions
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n)
  !           x,y: 1<=x,y<=n
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_d2xy_tab(n,d,x,y,typ)

    implicit none
    integer(ip), intent(in)             :: n,d
    integer(ip), intent(out)            :: x,y
    integer(2),  optional,intent(inout) :: typ
    integer(ip)                         :: ntot,naxi
    integer(ip)                         :: pos1,pos2,t
    integer(2)                          :: typ1

    integer(ip), PARAMETER,DIMENSION(16) :: postb2 = (/ 0,2,3,1,0,1,3,2,3,2,0,1,3,1,0,2/)
    integer(2),  PARAMETER,DIMENSION(16) :: typtb2 = (/ 1_2,0_2,0_2,2_2,0_2,1_2,1_2,3_2,3_2,2_2,2_2,0_2,2_2,3_2,3_2,1_2/)

    x    = 0_ip
    y    = 0_ip
    t    = d-1_ip
    typ1 = 0_2
    if(present(typ)) then
       typ1 = typ
    endif

    ntot = n*n
    naxi = n
    do while( ntot> 1_ip )
       naxi = naxi/2_ip
       pos1 = (t*4_ip)/ntot
       pos2  = int(postb2(int(typ1,ip)*4_ip+pos1+1_ip),KIND=ip)
       typ1  = typtb2(typ1*4+pos1+1)
       if(pos2==1_ip .or. pos2==3_ip) then
          x = x + naxi
       endif
       if(pos2>1_ip) then
          y = y + naxi
       endif
       t = t-(pos1*ntot)/4_ip
       ntot = ntot/4_ip
    end do

    x=x+1_ip
    y=y+1_ip
    if(present(typ)) typ = typ1

  end subroutine maths_sfc_d2xy_tab
  !-----------------------------------------------------------------------
  !
  !> @date    3/11/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 3D cartesian coordinates
  !> @details Kernel based on the table of recursive divisions
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n*n)
  !           x,y,z: 1<=x,y,z<=n
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_d2xyz_tab(n,d,x,y,z,typ)

    implicit none
    integer(ip), intent(in)             :: n,d
    integer(ip), intent(out)            :: x,y,z
    integer(2), optional,intent(inout)  :: typ
    integer(ip)                         :: ntot,naxi,t
    integer(ip)                         :: pos1,pos2,auxp
    integer(2)                          :: typ1

    x   = 0_ip
    y   = 0_ip
    z   = 0_ip
    t   = d-1_ip
    typ1 = 0_2
    if(present(typ)) typ1 = typ

    ntot = n*n*n
    naxi = n
    do while( ntot> 1_ip )
       naxi = naxi/2_ip
       pos1 = (t*8_ip)/ntot
       pos2  = postb3(int(typ1,ip)*8_ip+pos1+1_ip)
       typ1  = typtb3(int(typ1,ip)*8_ip+pos1+1_ip)
       auxp=mod(pos2,4_ip)
       if(auxp == 1_ip .or. auxp == 3_ip) x = x + naxi
       if(auxp > 1_ip) y = y + naxi
       if(pos2 > 3_ip) z = z + naxi
       t = t-(pos1*ntot)/8_ip
       ntot = ntot/8_ip
    end do

    x=x+1_ip
    y=y+1_ip
    z=z+1_ip
    if(present(typ)) typ = typ1

  end subroutine maths_sfc_d2xyz_tab
  !-----------------------------------------------------------------------
  !
  !> @date    3/11/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D lexicografical order to 3D cartesian coordinates
  !> @details n: number of boxes in each direction
  !>          d: position in lexicografical order (1<=d<=n*n*n)
  !>          x,y: 1<=x,y,z<=n
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_d2xyz_lex(n,d,x,y,z)
    implicit none
    integer(ip), intent(in)  :: n,d
    integer(ip), intent(out) :: x,y,z
    x=n
    y=n
    z=(d-1_ip)/(n*n)
    if(d-z*n*n >0_ip) then
       y=((d-z*n*n)-1_ip)/n
       x=d-z*n*n-y*n
    endif
    z=z+1_ip
    y=y+1_ip

  end subroutine maths_sfc_d2xyz_lex

  !-----------------------------------------------------------------------
  !
  !> @date    20/12/2017
  !> @author  Guillaume
  !> @brief   Computes the day of the week
  !> @details Computes the day of the week using Formula of Tomohiko
  !>          Sakamoto (1993)
  !
  !-----------------------------------------------------------------------

  function maths_day_of_week(d,m,y)

    integer(8), intent(in)    :: d
    integer(8), intent(in)    :: m
    integer(8), intent(in)    :: y
    character(9)              :: maths_day_of_week
    integer(8)                :: day_of_week_i
    integer(8), dimension(12) :: t = (/ 0_8, 3_8, 2_8, 5_8, 0_8, 3_8, 5_8, 1_8, 4_8, 6_8, 2_8, 4_8 /)
    integer(8)                :: y_new

    y_new = y
    if( m < 3 ) y_new = y_new - 1
    day_of_week_i = mod(y_new + y_new/4_8 - y_new*3_8/400_8  - y_new/4000_8 + t(m) + d,7_8)
    if( day_of_week_i == 0 ) day_of_week_i=7

    select case ( day_of_week_i )
    case ( 1_8 )    ; maths_day_of_week = 'Monday'
    case ( 2_8 )    ; maths_day_of_week = 'Tuesday'
    case ( 3_8 )    ; maths_day_of_week = 'Wednesday'
    case ( 4_8 )    ; maths_day_of_week = 'Thursday'
    case ( 5_8 )    ; maths_day_of_week = 'Friday'
    case ( 6_8 )    ; maths_day_of_week = 'Saturday'
    case ( 7_8 )    ; maths_day_of_week = 'Sunday'
    case default    ; maths_day_of_week = 'Unknown'
    end select

  end function maths_day_of_week

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2019
  !> @author  Guillaume Houzeaux
  !> @brief   Different elements in an array
  !> @details Extract the different elements in an array
  !>
  !-----------------------------------------------------------------------

  subroutine maths_list_different_elements(list_in,num_diff,list_diff,OPTION)

    integer(ip), intent(in),    pointer           :: list_in(:)
    integer(ip), intent(out)                      :: num_diff
    integer(ip), intent(inout), pointer, optional :: list_diff(:)
    character(*),                        optional :: OPTION
    logical(lg),                pointer           :: mask(:)
    integer(ip),                pointer           :: index_vector(:)
    integer(ip)                                   :: i,j,nn

    nn = size(list_in,KIND=ip)
    if( present(list_diff) ) nullify(list_diff)

    if( nn <= 0 .or. .not. associated(list_in) ) then

       num_diff = 0

    else

       allocate(mask(nn))
       mask = .true.

       do i = nn,2,-1
          mask(i) = .not.(any(list_in(:i-1)==list_in(i)))
       end do

       if( present(OPTION) ) then
          if( trim(OPTION) == 'STRICTLY POSITIVE' ) then
             where( list_in <= 0 ) mask = .false.
          end if
       end if

       ! Count

       num_diff = size(pack([(i,i=1,nn)],mask),KIND=ip)

       if( present(list_diff) .and. num_diff > 0 ) then
          !
          ! Make an index vector
          !
          allocate(index_vector(num_diff))
          index_vector=pack([(i,i=1,nn)],mask)
          !
          ! Now copy the unique elements of list_in to list_diff
          !
          allocate(list_diff(num_diff))

          list_diff=list_in(index_vector)
          deallocate(index_vector)

       end if

       deallocate(mask)

    end if

  end subroutine maths_list_different_elements

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2018
  !> @brief   Units
  !> @details Compute unit to scale a variable
  !>
  !-----------------------------------------------------------------------

  subroutine maths_unit(unit_value,unit_char,unit_factor)
    implicit none
    real(rp),      intent(in)  :: unit_value  !< Memory in bytes
    character(1),  intent(out) :: unit_char   !< Memory unit character
    real(rp),      intent(out) :: unit_factor !< Memory scaling

    if(      unit_value >= 1.0e12_rp ) then
       unit_factor = 1.0e-12_rp
       unit_char   = 'P'
    else if( unit_value >= 1.0e9_rp  ) then
       unit_factor = 1.0e-9_rp
       unit_char   = 'G'
    else if( unit_value >= 1.0e6_rp  ) then 
       unit_factor = 1.0e-6_rp
       unit_char   = 'M'     
    else if( unit_value >= 1.0e3_rp  ) then 
       unit_factor = 1.0e3_rp
       unit_char   = 'k'          
    else if( unit_value >= 1.0_rp ) then 
       unit_factor = 1.0_rp
       unit_char   = 'm'          
    else if( unit_value >= 1.0e-3_rp ) then 
       unit_factor = 1.0e3_rp
       unit_char   = 'm'          
    else if( unit_value >= 1.0e-6_rp ) then 
       unit_factor = 1.0e6_rp
       unit_char   = 'u'          
    else if( unit_value >= 1.0e-9_rp ) then 
       unit_factor = 1.0e9_rp
       unit_char   = 'n'          
    else if( unit_value >= 1.0e-12_rp ) then 
       unit_factor = 1.0e12_rp
       unit_char   = 'p'          
    else  
       unit_factor = 1.0_rp
       unit_char   = ' '     
    end if

  end subroutine maths_unit

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2018
  !> @author  Guillaume Houzeaux
  !> @brief   Different elements in an array
  !> @details Two and three-dimensional vectorial product of two vectors  v3 = v1 x v2.
  !>          The same pointer as for v1 or v2 may be used for v3. If N = 2, it is
  !>          assumed that v1 = (0,0,v1_3) and v2 = (v2_1,v2_2,0).      
  !>
  !-----------------------------------------------------------------------

  subroutine maths_vectorial_product(v1,v2,v3,n)

    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v2(n),v1(3)
    real(rp),    intent(out) :: v3(n)
    real(rp)                 :: c1,c2,c3

    if(n==2) then
       c1=-v1(3)*v2(2)
       c2= v1(3)*v2(1)
       v3(1)=c1
       v3(2)=c2
    else if(n==3) then
       c1=v1(2)*v2(3)-v1(3)*v2(2)
       c2=v1(3)*v2(1)-v1(1)*v2(3)
       c3=v1(1)*v2(2)-v1(2)*v2(1)
       v3(1)=c1
       v3(2)=c2
       v3(3)=c3
    end if

  end subroutine maths_vectorial_product

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   ???
  !> @details This routine computes the length of vector V and converts it to  
  !> a unit one 
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_normalize_vector(n,v,norm)

    integer(ip), intent(in)            :: n
    real(rp),    intent(inout)         :: v(n)
    real(rp),    intent(out), optional :: norm
    real(rp)                           :: rmod

    rmod = sqrt(dot_product(v(1:n),v(1:n)))
    if( rmod > epsilon(1.0_rp) ) v = v / rmod
    if( present(norm) ) norm = rmod

  end subroutine maths_normalize_vector

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-02
  !> @brief   Merge arrays
  !> @details Merge ordered arrays:
  !>          JA_IN  = [ 2 6 9 ]
  !>          JA_OUT = [ 3 6 8 10 ]
  !>          Result:
  !>          JA_OUT = [ 2 3 6 8 9 10 ]
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_merge_ordered_lists(nz_in,ja_in,nz_out,ja_out,MEMORY_COUNTER,RESIZE_GRAPH)

    integer(ip), intent(in)              :: nz_in
    integer(ip), intent(in)              :: ja_in(*)
    integer(ip), intent(out)             :: nz_out
    integer(ip), intent(inout), pointer  :: ja_out(:)
    integer(8),  intent(inout), optional :: MEMORY_COUNTER(2)
    logical(lg), intent(in),    optional :: RESIZE_GRAPH
    integer(ip)                          :: iz,ipoin,iz_out,last_position,iz_in
    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = 0_8
    end if
    !
    ! Look for last first available position 
    !    
    nz_out = size(ja_out)
    last_position = nz_out
    if( ja_out(last_position) /= 0 ) then
       nz_out = int(1.5_rp*real(size(ja_out,KIND=ip),rp))
       call memory_resize(memor,'JA_OUT','maths_merge_ordered_lists',ja_out,nz_out)
    end if
    do while( ja_out(last_position) == 0 )
       last_position = last_position - 1
    end do
    last_position = last_position + 1
    iz_out = 1
    !
    ! Scan entries of JA_IN and possibly add it to JA_OUT
    !
    do iz_in = 1,nz_in
       ipoin = ja_in(iz_in)
       if( ipoin > 0 ) then
          !
          ! Cjeck if entry already exists
          !
          loop_ipoin: do while( ipoin > ja_out(iz_out) .and. ja_out(iz_out) /= 0 )
             iz_out = iz_out + 1
             if( iz_out > nz_out ) exit loop_ipoin          
          end do loop_ipoin
          !
          ! Reallocate if necessay
          !
          if( iz_out > nz_out .or. last_position > nz_out ) then
             nz_out = int(1.5_rp*real(size(ja_out,KIND=ip),rp))
             call memory_resize(memor,'JA_OUT','maths_merge_ordered_lists',ja_out,nz_out)
          end if
          !
          ! Add new entry at position LAST_POSITION
          !
          if( ipoin /= ja_out(iz_out) ) then
             ja_out(last_position) = ipoin
             last_position = last_position + 1
          end if

       end if
    end do
    !
    ! Order JA_OUT
    !
    last_position = last_position - 1
    nz_out = last_position
    call maths_heap_sort(2_ip,last_position,ja_out)
    !
    ! Resize graph
    !
    if( present(RESIZE_GRAPH) ) then
       if( RESIZE_GRAPH ) & 
            call memory_resize(memor,'JA_OUT','maths_merge_ordered_lists',ja_out,nz_out)
    end if

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

  end subroutine maths_merge_ordered_lists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-02-25
  !> @brief   Solve an overdetermined system
  !> @details Solve Ax=b with A(n,m) n > m
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_solve_overdetermined_system(n,m,A,b,x,error)

    integer(ip), intent(in)  :: n
    integer(ip), intent(in)  :: m
    real(rp),    intent(in)  :: A(n,m)
    real(rp),    intent(in)  :: b(n,1)
    real(rp),    intent(out) :: x(m,1)
    real(rp),    intent(out) :: error
    integer(ip)              :: ii 
    real(rp),    allocatable :: At(:,:)
    real(rp),    allocatable :: At_A(:,:)
    real(rp),    allocatable :: At_A_inv(:,:)
    real(rp),    allocatable :: At_b(:,:)
    real(rp)                 :: deter,numer
    real(rp)                 :: dummr,denom

    allocate(At(m,n))
    allocate(At_A(m,m))
    allocate(At_A_inv(m,m))
    allocate(At_b(m,1))
    !
    ! x = (A^t A)^-1 (A^t b)
    !
    At   = transpose(A)
    At_A = matmul(At,A)
    call maths_invert_matrix(m,At_A,deter,At_A_inv)
    if( abs(deter) < epsil ) then
       error = -1.0_rp
    else
       At_b = matmul(At,b)
       x    = matmul(At_A_inv,At_b)
       !
       ! Error
       !
       error = 0.0_rp
       denom = 0.0_rp
       numer = 0.0_rp
       do ii = 1,n
          dummr = dot_product(x(:,1),A(ii,:))
          numer = numer + (dummr-b(ii,1))**2
          denom = denom + (b(ii,1))**2
       end do
       error = 100.0_rp*sqrt(numer/(denom+epsilon(1.0_rp)))
    end if

    deallocate(At,At_A,At_A_inv,At_b)

  end subroutine maths_solve_overdetermined_system

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-11
  !> @brief   Solve a quedratic equation
  !> @details Find roots of a*x^2+b*x+c=0
  !>          All the computations can be done in 16-byte precision:
  !>          Solutions are nevertheless converted into 8-byte precision
  !>          at the end.
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_quadratic_equation(aa,bb,cc,x1,x2,num_solutions,QUAD_REAL)

    real(rp),              intent(in)  :: aa             !< Polynomial coef
    real(rp),              intent(in)  :: bb             !< Polynomial coef
    real(rp),              intent(in)  :: cc             !< Polynomial coef
    real(rp),              intent(out) :: x1             !< First root
    real(rp),              intent(out) :: x2             !< Second root
    integer(ip),           intent(out) :: num_solutions  !< Number of solutions
    logical(lg), optional, intent(in)  :: QUAD_REAL      !< If 16-bytes precision is required
    real(rp)                           :: delta

    real(qp)                           :: aa_qp
    real(qp)                           :: bb_qp
    real(qp)                           :: cc_qp
    real(qp)                           :: x1_qp
    real(qp)                           :: delta_qp
    logical(lg)                        :: if_quad_real

    if( present(QUAD_REAL) ) then
       if_quad_real = QUAD_REAL
    else
       if_quad_real = .false.
    end if

    if( if_quad_real ) then

       !----------------------------------------------------------------
       !
       ! Quad-real precision
       !
       !----------------------------------------------------------------

       aa_qp = real(aa,qp)
       bb_qp = real(bb,qp)
       cc_qp = real(cc,qp)
       if( abs(aa_qp) > 0.0_qp ) then  
          delta_qp  = ((bb_qp/aa_qp)*bb_qp - 4.0_qp*cc_qp)/aa_qp
          if( delta_qp == 0.0_qp ) then
             !
             ! One real solution:
             !
             !      -b
             ! x1 = -- 
             !      2a
             !
             x1 = real(-bb_qp/(2.0_qp*aa_qp),rp)
             num_solutions = 1_ip
          else if( delta_qp < 0.0_qp ) then
             !
             ! No real solution
             !
             num_solutions = 0_ip        
          else
             !
             ! Two real solutions
             !
             delta_qp = sqrt(delta_qp)
             if( bb_qp >= 0.0_qp ) then
                x1_qp = 0.5_qp * (-bb_qp/aa_qp-delta_qp) 
                x1    = real(x1_qp,rp)
                if( x1_qp /= 0.0_qp ) then
                   x2 = real(cc_qp / (aa_qp*x1_qp),rp)
                else
                   x2 = real(0.5_qp * (-bb_qp/aa_qp+delta_qp),rp)
                end if
             else
                x1_qp = 0.5_qp * (-bb_qp/aa_qp+delta_qp)
                x1    = real(x1_qp,rp)
                if( x1_qp /= 0.0_qp ) then
                   x2 = real(cc_qp / (aa_qp*x1_qp),rp)
                else
                   x2 = real(0.5_qp * (-bb_qp/aa_qp-delta_qp),rp)
                end if
             end if
             num_solutions = 2_ip
          end if

       else
          !
          ! a=0
          !
          if( abs(bb_qp) > 0.0_qp ) then
             x1 = real(-bb_qp/cc_qp,rp)
             num_solutions = 1_ip
          else
             num_solutions = 0_ip
          end if

       end if

    else

       !----------------------------------------------------------------
       !
       ! RP real precision
       !
       !----------------------------------------------------------------

       if( abs(aa) > 0.0_rp ) then  
          delta  = ((bb/aa)*bb - 4.0_rp*cc)/aa
          if( delta == 0.0_rp ) then
             !
             ! One real solution:
             !
             !      -b
             ! x1 = -- 
             !      2a
             !
             x1 = -bb/(2.0_rp*aa)
             num_solutions = 1_ip
          else if( delta < 0.0_rp ) then
             !
             ! No real solution
             !
             num_solutions = 0_ip        
          else
             !
             ! Two real solutions
             !
             delta = sqrt(delta)
             if( bb >= 0.0_rp ) then
                x1 = 0.5_rp * (-bb/aa-delta) 
                if( x1 /= 0.0_rp ) then
                   x2 = cc / (aa*x1)
                else
                   x2 = 0.5_rp * (-bb/aa+delta)
                end if
             else
                x1 = 0.5_rp * (-bb/aa+delta) 
                if( x1 /= 0.0_rp ) then
                   x2 = cc / (aa*x1)
                else
                   x2 = 0.5_rp * (-bb/aa-delta) 
                end if
             end if
             num_solutions = 2_ip
          end if

       else
          !
          ! a=0
          !
          if( abs(bb) > 0.0_rp ) then
             x1 = -bb/cc
             num_solutions = 1_ip
          else
             num_solutions = 0_ip
          end if

       end if
    end if

  end subroutine maths_quadratic_equation

  subroutine maths_quadratic_equation_unity_test()

    real(rp)    :: aa,bb,cc,x1,x2,x3,x4,xtmp,delta
    real(qp)    :: xe1,xe2
    integer(ip) :: num1,num2,solution

    solution = 4

    select case ( solution )

    case ( 1 )

       aa  =  2.0_rp
       bb  = -3.0_rp
       cc  =  1.0_rp
       xe1 =  0.5_qp
       xe2 =  1.0_qp

    case ( 2 )

       aa  =  1.0_rp 
       bb  =  200.0_rp
       cc  = -0.000015_rp
       xe1 =-200.000000075_qp
       xe2 =   0.000000075_qp

    case ( 3 )

       xe1 = 1.000000000000000_qp
       xe2 = 1.000000028975958_qp
       aa  = 94906265.625_rp
       bb  = -189812534.0_rp
       cc  = 94906268.375_rp

    case ( 4 )

       aa  = 1.0_rp
       bb  = 1.0e8_rp
       cc  = 1.0_rp
       xe1 = -99999999.9999999899999999999999990024_qp
       xe2 = -1.0000000000000001000000000000000e-0008_qp

    end select
    
    x1 = 0.0_rp
    x2 = 0.0_rp
    x3 = 0.0_rp
    x4 = 0.0_rp
    !
    ! Our algorithm
    !
    call maths_quadratic_equation(aa,bb,cc,x1,x2,num1) !,QUAD_REAL=.true.) 
    !
    ! Original algorithm
    !
    delta = bb**2-4.0_rp*aa*cc
    if( delta > 0.0_rp ) then
       x3   = (-bb-sqrt(delta))/(2.0_rp*aa)
       x4   = (-bb+sqrt(delta))/(2.0_rp*aa)
       num2 = 2
    else if( delta == 0.0_rp ) then
       x3   = -bb/(2.0_rp*aa)
       x4   = 0.0_rp
       num2 = 1
    end if
    !
    ! Order solution
    !
    if( num2 == 2 ) then
       if( x4 < x3 ) then
          xtmp = x3
          x3   = x4
          x4   = xtmp
       end if
    end if
    if( num1 == 2 ) then
       if( x2 < x1 ) then
          xtmp = x1
          x1   = x2
          x2   = xtmp
       end if
    end if
    
    print*,' '
    print*,'Modified algorithm:'
    print*,'  - Number solution= ',num1
    print*,'  - Solution 1=      ',x1, abs(real(x1,qp)-xe1)/abs(xe1)*100.0_qp
    print*,'  - Solution 2=      ',x2, abs(real(x2,qp)-xe2)/abs(xe2)*100.0_qp
    print*,' '
    print*,'Original algorithm:'
    print*,'  - Number solution= ',num2
    print*,'  - Solution 1=      ',x3, abs(real(x3,qp)-xe1)/abs(xe1)*100.0_qp
    print*,'  - Solution 2=      ',x4, abs(real(x4,qp)-xe2)/abs(xe2)*100.0_qp

  end subroutine maths_quadratic_equation_unity_test

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   heap sort
  !> @details Equalize some arrays
  !>          Quick sorting. The element in ivin are sorting in:
  !>          ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !>          ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  !
  !----------------------------------------------------------------------

  subroutine maths_heap_sort_I1(itask,nrows,ivin,message,ivo1,ivo2,PERMUTATION,ELIMINATE_REPLICATES)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(inout)           :: nrows
    integer(ip),  intent(inout)           :: ivin(*)
    character(*), intent(in),    optional :: message
    integer(ip),  intent(inout), optional :: ivo1(*)
    integer(ip),  intent(inout), optional :: ivo2(*)
    integer(ip),  intent(inout), optional :: PERMUTATION(*)
    logical(lg),  intent(in),    optional :: ELIMINATE_REPLICATES
    integer(ip)                           :: leni,ir,ii,jj,iaux,krows
    integer(ip)                           :: iau1,iau2,kk
    logical(lg)                           :: if_eliminate_replicates

    if_eliminate_replicates = .false.
    if( present(ELIMINATE_REPLICATES) ) if_eliminate_replicates = ELIMINATE_REPLICATES
    if( present(message) ) then
       if( trim(message) == 'ELIMINATE DUPLICATES') if_eliminate_replicates = .true.
    end if

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux = ivin(ir)
          if( present(ivo1) ) iau1 = ivo1(ir)
          if( present(ivo2) ) iau2 = ivo2(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) ivo1(ir) = ivo1(1)
          if( present(ivo2) ) ivo2(ir) = ivo2(1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1     = ivo1(ir)
             ivo1(ir) = ivo1(1)
          end if
          if( present(ivo2) ) then
             iau2     = ivo2(ir)
             ivo2(ir) = ivo2(1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

       if( present(PERMUTATION) ) then
401       if( jj <= ir ) then
             if( jj < ir ) then
                if ( PERMUTATION(ivin(jj)) < PERMUTATION(ivin(jj+1)) ) then
                   jj = jj + 1
                end if
             end if

             if( PERMUTATION(iaux) < PERMUTATION(ivin(jj)) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 401
          end if
       else
400       if( jj <= ir ) then
             if( jj < ir ) then
                if ( ivin(jj) < ivin(jj+1) ) then
                   jj = jj + 1
                end if
             end if

             if( iaux < ivin(jj) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 400
          end if
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 300

    end select
    !
    ! Eliminate duplicates
    !
500 if( if_eliminate_replicates ) then
       if( itask == 1 ) stop
       jj = 0
       kk = 0
       do ii = 1,nrows
          if( ivin(ii) <= jj ) then
             ivin(ii) = 0
             if( present(ivo1) ) ivo1(ii) = 0
             if( present(ivo2) ) ivo2(ii) = 0
          else
             kk       = kk + 1
             ivin(kk) = ivin(ii)
             jj       = ivin(ii)
             if( present(ivo1) ) ivo1(kk) = ivo1(ii)
             if( present(ivo2) ) ivo2(kk) = ivo2(ii)
          end if
       end do
       if( kk < nrows ) then
          ivin(kk+1:nrows) = 0
          if( present(ivo1) ) ivo1(kk+1:nrows) = 0
          if( present(ivo2) ) ivo2(kk+1:nrows) = 0
       end if
       nrows = kk
    end if


  end subroutine maths_heap_sort_I1

  subroutine maths_heap_sort_I2(itask,nrows,ivin,ivo1,ivo2,message)

    integer(ip),  intent(in)                        :: itask
    integer(ip),  intent(inout)                     :: nrows
    integer(ip),  intent(inout)                     :: ivin(*)
    integer(ip),  intent(inout),  pointer           :: ivo1(:,:)
    integer(ip),  intent(inout),  pointer, optional :: ivo2(:,:)
    character(*), intent(in),              optional :: message
    integer(ip)                                     :: leni,ir,ii,jj,iaux,krows,ncol1,ncol2
    integer(ip),                  pointer           :: iau1(:)
    integer(ip),                  pointer           :: iau2(:)

    nullify(iau1)
    nullify(iau2)

    ncol1 = size(ivo1,2,KIND=ip)
    if( ncol1 < 1 ) return
    if( size(ivo1,1,KIND=ip) < nrows ) stop
    allocate(iau1(ncol1))

    if( present(ivo2) ) then
       ncol2 = size(ivo2,2,KIND=ip)
       if( ncol2 < 1 ) return
       if( size(ivo2,1,KIND=ip) < nrows ) stop
       allocate(iau2(ncol2))
    end if

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          iau1(1:ncol1) = ivo1(leni,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(leni,1:ncol2)
          end if
       else
          iaux = ivin(ir)
          iau1(1:ncol1) = ivo1(ir,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(ir,1:ncol2)
          end if
          ivin(ir) = ivin(1)
          ivo1(ir,1:ncol1) = ivo1(1,1:ncol1)
          if( present(ivo2) ) then
             ivo2(ir,1:ncol2) = ivo2(1,1:ncol2)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(1,1:ncol1) = iau1(1:ncol1)
             if( present(ivo2) ) then
                ivo2(1,1:ncol2) = iau2(1:ncol2)
             end if
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(ii,1:ncol1) = ivo1(jj,1:ncol1)
             if( present(ivo2) ) then
                ivo2(ii,1:ncol2) = ivo2(jj,1:ncol2)
             end if
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       ivo1(ii,1:ncol1) = iau1(1:ncol1)
       if( present(ivo2) ) then
          ivo2(ii,1:ncol2) = iau2(1:ncol2)
       end if

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          iau1(1:ncol1) = ivo1(leni,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2) = ivo2(leni,1:ncol2)
          end if
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          iau1(1:ncol1)    = ivo1(ir,1:ncol1)
          ivo1(ir,1:ncol1) = ivo1(1,1:ncol1)
          if( present(ivo2) ) then
             iau2(1:ncol2)    = ivo2(ir,1:ncol2)
             ivo2(ir,1:ncol2) = ivo2(1,1:ncol2)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(1,1:ncol1) = iau1(1:ncol1)
             if( present(ivo2) ) then
                ivo2(1,1:ncol2) = iau2(1:ncol2)
             end if
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(ii,1:ncol1) = ivo1(jj,1:ncol1)
             if( present(ivo2) ) then
                ivo2(ii,1:ncol2) = ivo2(jj,1:ncol2)
             end if
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ivin(ii) = iaux
       ivo1(ii,1:ncol1) = iau1(1:ncol1)
       if( present(ivo2) ) then
          ivo2(ii,1:ncol2) = iau2(1:ncol2)
       end if

       goto 300

    end select

    if( associated(iau1) ) deallocate(iau1)
    if( associated(iau2) ) deallocate(iau2)

    return
    !
    ! Eliminate duplicates
    !
500 continue
    if( present(message) ) then
       if( trim(message) == 'ELIMINATE DUPLICATES') then
          stop
       end if
    end if


  end subroutine maths_heap_sort_I2

  subroutine maths_heap_sort_RP1(itask,nrows,ivin,ivo1,ivo2)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(inout)           :: nrows
    real(rp),     intent(inout)           :: ivin(*)
    integer(ip),  intent(inout), optional :: ivo1(*)
    integer(ip),  intent(inout), optional :: ivo2(*)
    real(rp)                              :: iaux
    integer(ip)                           :: leni,ir,ii,jj
    integer(ip)                           :: iau1,iau2

    select case ( itask )
 
    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

100    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)
          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux = ivin(ir)
          if( present(ivo1) ) iau1 = ivo1(ir)
          if( present(ivo2) ) iau2 = ivo2(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) ivo1(ir) = ivo1(1)
          if( present(ivo2) ) ivo2(ir) = ivo2(1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 100      

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)

          if( present(ivo1) ) iau1 = ivo1(leni)
          if( present(ivo2) ) iau2 = ivo2(leni)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1     = ivo1(ir)
             ivo1(ir) = ivo1(1)
          end if
          if( present(ivo2) ) then
             iau2     = ivo2(ir)
             ivo2(ir) = ivo2(1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 300

    end select

500 continue

  end subroutine maths_heap_sort_RP1

  subroutine maths_heap_sort_RP2(itask,ndofn,nrows,ivin,ivo1,ivo2)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(inout)           :: nrows
    integer(ip),  intent(inout)           :: ivin(*)
    real(rp),     intent(inout), optional :: ivo1(ndofn,*)
    real(rp),     intent(inout), optional :: ivo2(ndofn,*)
    integer(ip)                           :: iaux
    integer(ip)                           :: leni,ir,ii,jj
    real(rp)                              :: iau1(ndofn),iau2(ndofn)

    select case ( itask )

    case ( 1_ip )
       !
       ! Oups
       !
       stop

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       leni = (nrows/2) + 1
       ir  = nrows

300    continue

       if( leni > 1 ) then
          leni = leni - 1
          iaux = ivin(leni)

          if( present(ivo1) ) iau1 = ivo1(:,leni)
          if( present(ivo2) ) iau2 = ivo2(:,leni)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1       = ivo1(:,ir)
             ivo1(:,ir) = ivo1(:,1)
          end if
          if( present(ivo2) ) then
             iau2       = ivo2(:,ir)
             ivo2(:,ir) = ivo2(:,1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(:,1) = iau1
             if( present(ivo2) ) ivo2(:,1) = iau2
             goto 500
          end if
       end if

       ii = leni
       jj = leni + leni

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(:,ii) = ivo1(:,jj)
             if( present(ivo2) ) ivo2(:,ii) = ivo2(:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(:,ii) = iau1
       if( present(ivo2) ) ivo2(:,ii) = iau2

       goto 300

    end select

500 continue

  end subroutine maths_heap_sort_RP2

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   heap sort
  !> @details Quick sort
  !>          Quick sorting. The element in ivin are sorting in:
  !>          ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !>          ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  !
  !----------------------------------------------------------------------
    
  recursive subroutine maths_quick_sort_rp_1(itask,ivin)

    integer(ip),          intent(in)    :: itask
    real(rp),    pointer, intent(inout) :: ivin(:)
    integer(ip)                         :: iq,nrows

    nrows = size(ivin)
    if( nrows > 1) then
       call maths_partition_rp(itask,nrows,ivin,iq)
       call maths_quick_sort_rp(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_rp(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_rp_1

 recursive subroutine maths_quick_sort_rp(itask,nrows,ivin)

    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    real(rp),    intent(inout) :: ivin(*)
    integer(ip)                :: iq

    if( nrows > 1) then
       call maths_partition_rp(itask,nrows,ivin,iq)
       call maths_quick_sort_rp(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_rp(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_rp

  recursive subroutine maths_quick_sort_ip_1(itask,ivin)

    integer(ip),          intent(in)    :: itask
    integer(ip), pointer, intent(inout) :: ivin(:)
    integer(ip)                         :: iq,nrows

    nrows = size(ivin)
    if( nrows > 1) then
       call maths_partition_ip(itask,nrows,ivin,iq)
       call maths_quick_sort_ip(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_ip(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_ip_1

 recursive subroutine maths_quick_sort_ip(itask,nrows,ivin)

    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    integer(ip), intent(inout) :: ivin(*)
    integer(ip)                :: iq

    if( nrows > 1) then
       call maths_partition_ip(itask,nrows,ivin,iq)
       call maths_quick_sort_ip(itask,iq-1_ip,ivin(1:iq-1))
       call maths_quick_sort_ip(itask,nrows-iq+1_ip,ivin(iq:nrows))
    endif

  end subroutine maths_quick_sort_ip

  subroutine maths_partition_rp(itask,nrows,ivin,marker)
    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    real(rp),    intent(inout) :: ivin(*)
    integer(ip), intent(out)   :: marker
    integer(ip)                :: i,j
    real(rp)                   :: temp
    real(rp)                   :: x      ! pivot point

    x = ivin(1)
    i= 0
    j= nrows + 1

    if( itask == 1 ) then
       do
          j = j-1
          do
             if( ivin(j) >= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) <= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    else
       do
          j = j-1
          do
             if( ivin(j) <= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) >= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    end if
    
  end subroutine maths_partition_rp

  subroutine maths_partition_ip(itask,nrows,ivin,marker)
    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: nrows
    integer(ip), intent(inout) :: ivin(*)
    integer(ip), intent(out)   :: marker
    integer(ip)                :: i,j
    integer(ip)                :: temp
    integer(ip)                :: x      ! pivot point

    x = ivin(1)
    i= 0
    j= nrows + 1

    if( itask == 1 ) then
       do
          j = j-1
          do
             if( ivin(j) >= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) <= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    else
       do
          j = j-1
          do
             if( ivin(j) <= x ) exit
             j = j-1
          end do
          i = i+1
          do
             if( ivin(i) >= x ) exit
             i = i+1
          end do
          if (i < j) then
             ! exchange ivin(i) and ivin(j)
             temp    = ivin(i)
             ivin(i) = ivin(j)
             ivin(j) = temp
          elseif (i == j) then
             marker = i+1
             return
          else
             marker = i
             return
          endif
       end do
    end if
    
  end subroutine maths_partition_ip

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-01
  !> @brief   Sorting unity test
  !> @details Unity tests of sorting algorithms for real and integer
  !> 
  !-----------------------------------------------------------------------

!!$  subroutine maths_sort_unity_test()
!!$
!!$    integer(ip)           :: nn,ii,itask,imeth,ierro,jj
!!$    integer(ip),parameter :: seed = 86456
!!$    real(rp),   pointer   :: aa(:)
!!$    integer(ip),pointer   :: kk(:)
!!$    real(rp)              :: time1,time2
!!$    character(10)         :: method(2)
!!$    character(10)         :: realint(2)
!!$    character(10)         :: order(2)  
!!$
!!$    order(1)   = 'DECREASING'
!!$    order(2)   = 'INCREASING'
!!$    method(1)  = 'QUICK SORT'
!!$    method(2)  = 'HEAP  SORT'
!!$    realint(1) = 'REAL'
!!$    realint(2) = 'INT.'
!!$
!!$    nn = 600000
!!$    allocate(aa(nn))
!!$    allocate(kk(nn))
!!$    aa = 0.0_rp
!!$    kk = 0_ip
!!$    call srand(seed)  
!!$
!!$    do itask = 1,2
!!$       !
!!$       ! ITASK = 1 ... Decreasing order
!!$       ! ITASK = 2 ... Increasing order
!!$       !     
!!$       ierro = 0
!!$       do imeth = 1,2
!!$          do jj = 1,2
!!$             call cputim(time1)
!!$             do ii = 1,nn
!!$                aa(ii) = rand()
!!$                kk(ii) = -nn/2+FLOOR(nn*aa(ii)) 
!!$             end do
!!$             if( imeth == 1 ) then
!!$                if( jj == 1 ) then
!!$                   call maths_quick_sort(itask,nn,aa)
!!$                else
!!$                   call maths_quick_sort(itask,nn,kk)
!!$                end if
!!$             else
!!$                if( jj == 1 ) then
!!$                   call maths_heap_sort(itask,nn,aa)
!!$                else
!!$                   call maths_heap_sort(itask,nn,kk)
!!$                end if
!!$             end if
!!$             call cputim(time2)
!!$             print*,trim(method(imeth))//', '//trim(order(ITASK))//', '//trim(realint(jj))//', ','TIME= ',time2-time1
!!$             if( jj == 1 ) then
!!$                if( itask == 2 ) then
!!$                   do ii = 1,nn-1
!!$                      if( aa(ii) > aa(ii+1) ) then
!!$                         ierro = 1
!!$                         goto 1
!!$                      end if
!!$                   end do
!!$                else
!!$                   do ii = 1,nn-1
!!$                      if( aa(ii) < aa(ii+1) ) then
!!$                         ierro = 1
!!$                         goto 1
!!$                      end if
!!$                   end do
!!$                end if
!!$             else
!!$                if( itask == 2 ) then
!!$                   do ii = 1,nn-1
!!$                      if( kk(ii) > kk(ii+1) ) then
!!$                         ierro = 1
!!$                         goto 1
!!$                      end if
!!$                   end do
!!$                else
!!$                   do ii = 1,nn-1
!!$                      if( kk(ii) < kk(ii+1) ) then
!!$                         ierro = 1
!!$                         goto 1
!!$                      end if
!!$                   end do
!!$                end if
!!$             end if
!!$          end do
!!$       end do
!!$
!!$    end do
!!$
!!$    deallocate(aa,kk)
!!$    stop
!!$1   print*,'error: ',imeth,itask
!!$    stop
!!$    
!!$  end subroutine maths_sort_unity_test

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-04
  !> @brief   Szudzik pairing function
  !> @details Provides a bijective 2 to 1 mapping
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function maths_Szudzik_pairing_function(ii,jj)

    integer(ip), intent(in) :: ii
    integer(ip), intent(in) :: jj
    
    if( ii /= max(ii,jj) ) then
       maths_Szudzik_pairing_function = jj*jj+ii
    else
       maths_Szudzik_pairing_function = ii*(ii+1)+jj
    end if

  end function maths_Szudzik_pairing_function

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-19
  !> @brief   Return position of last non-zero value
  !> @details Return position of last non-zero value
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function maths_maxloc_nonzero(ll)
    
    integer(ip), intent(in) :: ll(:)
    integer(ip)             :: ii(1)

    ii = maxloc(ll(:),ll(:)>0)
    maths_maxloc_nonzero = ii(1)
    
  end function maths_maxloc_nonzero

end module mod_maths
!> @}
