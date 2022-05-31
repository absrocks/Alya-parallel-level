!-----------------------------------------------------------------------
!> @addtogroup AMR
!> @{
!> @file    mod_error_estimate.f90
!> @author  houzeaux
!> @date    2020-05-05
!> @brief   Error estimate
!> @details Some error estimation techniques
!-----------------------------------------------------------------------

module mod_error_estimate

  use def_kintyp_basic,        only : ip,rp
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use def_master,              only : optional_argument
  use def_master,              only : zeror
  use def_master,              only : npart
  use def_elmtyp
  use def_domain
  use mod_elmgeo
  use mod_communications
  use mod_memory
  use mod_maths
  use mod_gradie
  use mod_postpr
  implicit none

  public :: error_estimate_mesh_size
  private
  
contains
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Mesh size
  !> @details Compute the error and then the optimum mesh size
  !> 
  !-----------------------------------------------------------------------

   subroutine error_estimate_mesh_size(solut,hh_opt,nelem_opt)

    real(rp),    pointer,  intent(in)    :: solut(:)
    real(rp),    pointer,  intent(inout) :: hh_opt(:)
    integer(ip), optional, intent(in)    :: NELEM_OPT
    real(rp),    pointer                 :: hh(:)
    real(rp),    pointer                 :: ee(:)

    nullify(hh,ee)

    call error_estime_circumradius(hh)
    call error_estime_hessian(solut,ee)
    call error_estime_e_to_h(ee,hh,hh_opt,FACTOR=1.0_rp,NELEM_OPT=nelem_opt)
    
  end subroutine error_estimate_mesh_size
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Compute the circumradius of the elements
  !>          Circumradius is only coded for TRI03, QUA04, TET04. For
  !>          other elements, the average edge length is taken
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_circumradius(hh)

    real(rp), pointer, intent(inout) :: hh(:)
    integer(ip)                      :: ielem
    integer(ip)                      :: pelty,pnode
    
    if( .not. associated(hh) ) & 
         call memory_alloca(memor_dom,'HH'    ,'mod_error_estimate',hh,nelem)
    
    do ielem = 1,nelem
       pelty     = ltype(ielem)
       pnode     = lnnod(ielem)
       hh(ielem) = elmgeo_circumradius(pelty,coord(1:ndime,lnods(1:pnode,ielem)))
   end do

  end subroutine error_estime_circumradius
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Convert error
  !> @details Convert the error into a mesh size
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_hessian(uu,ee)

    real(rp), pointer, intent(in)    :: uu(:)
    real(rp), pointer, intent(inout) :: ee(:)
    real(rp), pointer                :: gradu(:,:)
    integer(ip)                      :: ipoin,idime,inode,ielem
    integer(ip)                      :: pnode,pgaus,plapl,pelty
    real(rp)                         :: gpcar(ndime,mnode,1)              ! dN/dxi
    real(rp)                         :: gpvol(1)                          ! w*|J|, |J|
    real(rp)                         :: elgra(ndime,mnode)
    real(rp)                         :: elcod(ndime,mnode)
    real(rp)                         :: gplap,x

    nullify(gradu)
    call memory_alloca(memor_dom,'GRADU','mod_error_estimate',gradu,ndime,npoin)
    if( npoin > 0 ) call gradie(uu,gradu)

    if( .not. associated(ee) ) & 
         call memory_alloca(memor_dom,'EE','mod_error_estimate',ee,nelem)

    do ielem = 1,nelem
       pelty = ltype(ielem)
       pnode = lnnod(ielem)
       do inode = 1,pnode
          ipoin                = lnods(inode,ielem)
          elgra(1:ndime,inode) = gradu(1:ndime,ipoin)
          elcod(1:ndime,inode) = coord(1:ndime,ipoin)
       end do

       call elmgeo_cartesian_derivatives(ndime,pnode,elcod,elmar(pelty) % dercg,gpcar)
      
       gplap = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             gplap = gplap + gpcar(idime,inode,1) * elgra(idime,inode)
          end do
       end do
       ee(ielem) = abs(gplap)
       
    end do
    
  end subroutine error_estime_hessian
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Convert error
  !> @details Convert the error into a mesh size
  !> 
  !-----------------------------------------------------------------------

  subroutine error_estime_e_to_h(ee,hh,hh_opt,FACTOR,NELEM_OPT)
    
    real(rp),              pointer, intent(in)    :: ee(:)
    real(rp),              pointer, intent(in)    :: hh(:)
    real(rp),              pointer, intent(inout) :: hh_opt(:)
    real(rp),    optional,          intent(in)    :: FACTOR
    integer(ip), optional,          intent(in)    :: NELEM_OPT
    real(rp)                                      :: rr,kk,nn,dd,xfact
    real(rp)                                      :: beta,gamma,alpha
    integer(ip)                                   :: ielem
    real(rp)                                      :: xx
    
    if( .not. associated(hh_opt) ) & 
         call memory_alloca(memor_dom,'HH_OPT','mod_error_estimate',hh_opt,nelem)
    !
    ! Target number of elements
    !
    xx = optional_argument(1.0_rp,FACTOR)
    if( present(NELEM_OPT) ) then
       nn = real(nelem_opt,rp)/real(max(1_ip,npart),rp)*xx
    else
       nn = real(nelem,rp)*xx
       call PAR_SUM(nn)
    end if
    
    dd    = real(ndime,rp)
    kk    = 2.0_rp
    alpha = 2.0_rp * kk / dd
    xfact = 1.0_rp/(dd*(1.0_rp+alpha))

    beta  = 0.0_rp
    gamma = 2.0_rp/(1.0_rp+alpha)
    !
    ! beta = Sum_{i=1}^N ei^{2/(1+alpha)}
    !
    do ielem = 1,nelem
       beta = beta + abs(ee(ielem)) ** gamma
    end do
    call PAR_SUM(beta)
     
    beta = beta * ( alpha ** ((2.0_rp+alpha)/(1.0_rp+alpha)) + alpha ** (1.0_rp/(1.0_rp+alpha)) )
    beta = ( (1.0_rp + alpha ) * nn / beta ) ** (1.0_rp/dd)
    beta = beta * ( alpha ** xfact )
   
    do ielem = 1,nelem      
       rr            = ( abs(ee(ielem)) ** (2.0_rp*xfact) ) * beta
       hh_opt(ielem) = hh(ielem) / (rr+zeror)
    end do

  end subroutine error_estime_e_to_h

end module mod_error_estimate
!> @}
 
