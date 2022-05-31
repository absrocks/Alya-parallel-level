!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_element_integration.f90
!> @author  houzeaux
!> @date    2018-11-30
!> @brief   Element integration
!> @details Subroutines to deal with element integration, like
!>          computing Cartesian derivatives of shape functions, etc.
!-----------------------------------------------------------------------

module mod_element_integration

  use def_kintyp, only : ip,rp
#ifndef VECTOR_SIZE
  use def_master, only : VECTOR_SIZE
#endif
  use def_domain, only : ndime,ntens,mnode
  use def_domain, only : kfl_savda,elmda
  use def_domain, only : memor_dom
  use def_domain, only : elm
  use mod_elmgeo, only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo, only : elmgeo_hessian
  use mod_elmgeo, only : element_type
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none

  private


  interface element_shape_function_derivatives_jacobian
     module procedure element_shape_function_derivatives_jacobian_scalar,&
          &           element_shape_function_derivatives_jacobian_vector,&
          &           element_shape_function_derivatives_jacobian_vector_2
  end interface element_shape_function_derivatives_jacobian

  public :: element_shape_function_derivatives_jacobian
  public :: element_shape_function_derivatives_jacobian_vector_2
  public :: element_integration
  public :: element_integration_initialization

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-30
  !> @brief   Compute shape functions, derivatives and Hessian
  !> @details This routine calculates:
  !>          GPCAR: Cartesian derivatives
  !>          GPVOL: Unit volume
  !>          GPHES: Hessian matrix
  !>          Scalar and vector subroutines are available
  !> 
  !-----------------------------------------------------------------------

  subroutine element_shape_function_derivatives_jacobian_vector_2(&
       pnode,pgaus,plapl,weigp,shapf,deriv,heslo,&
       elcod,gpvol,gpsha,gpder,gpcar,gphes,&
       list_elements)
    
    integer(ip), intent(in)             :: pnode
    integer(ip), intent(in)             :: plapl
    integer(ip), intent(in)             :: pgaus

    real(rp),    intent(in)             :: weigp                          !< Weight of integration paints
    real(rp),    intent(in)             :: shapf(pnode)                   !< Shape functions at integration points
    real(rp),    intent(in)             :: deriv(ndime,pnode)             !< Shape function derivatives at integration points
    real(rp),    intent(in)             :: heslo(ntens,pnode)             !< Shape function Hessian at integration points

    real(rp),    intent(in)             :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)            :: gpvol(VECTOR_SIZE,*)
    real(rp),    intent(out)            :: gpsha(VECTOR_SIZE,pnode,*)
    real(rp),    intent(out)            :: gpder(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)            :: gpcar(VECTOR_SIZE,ndime,mnode,*)
    real(rp),    intent(out), optional  :: gphes(VECTOR_SIZE,ntens,mnode,*)
    integer(ip), intent(in)             :: list_elements(VECTOR_SIZE)
    real(rp)                            :: dummr(1)

    if( pgaus /= 1 ) call runend('THIS SUBROUTINE IS INTENDED FOR 1 INTEGRATION POINT')
    dummr = weigp
    call element_shape_function_derivatives_jacobian_vector(&
       pnode,pgaus,plapl,dummr,shapf,deriv,heslo,&
       elcod,gpvol,gpsha,gpder,gpcar,gphes,&
       list_elements)
    
  end subroutine element_shape_function_derivatives_jacobian_vector_2

  subroutine element_shape_function_derivatives_jacobian_vector(&
       pnode,pgaus,plapl,weigp,shapf,deriv,heslo,&
       elcod,gpvol,gpsha,gpder,gpcar,gphes,&
       list_elements)

    !-----------------------------------------------------------------------
    !****f* domain/elmca2
    ! NAME
    !    elmca2
    ! DESCRIPTION
    !    This routine calculates:
    !    GPCAR: Cartesian derivatives
    !    GPVOL: Unit volume
    !    GPHES: Hessian matrix
    ! USES
    !    invmtx
    ! USED BY
    !    ***_elmope
    !    extnor
    ! SOURCE
    !-----------------------------------------------------------------------

    integer(ip), intent(in)             :: pnode
    integer(ip), intent(in)             :: plapl
    integer(ip), intent(in)             :: pgaus

    real(rp),    intent(in)             :: weigp(*)                       !< Weight of integration paints
    real(rp),    intent(in)             :: shapf(pnode,*)                 !< Shape functions at integration points
    real(rp),    intent(in)             :: deriv(ndime,pnode,*)           !< Shape function derivatives at integration points
    real(rp),    intent(in)             :: heslo(ntens,pnode,*)           !< Shape function Hessian at integration points

    real(rp),    intent(in)             :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)            :: gpvol(VECTOR_SIZE,*)
    real(rp),    intent(out)            :: gpsha(VECTOR_SIZE,pnode,*)
    real(rp),    intent(out)            :: gpder(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)            :: gpcar(VECTOR_SIZE,ndime,mnode,*)
    real(rp),    intent(out), optional  :: gphes(VECTOR_SIZE,ntens,mnode,*)
    integer(ip), intent(in)             :: list_elements(VECTOR_SIZE)
    integer(ip)                         :: igaus,inode,idime
    integer(ip)                         :: itens,ivect,ielem
    real(rp)                            :: gpcod,gpdet(VECTOR_SIZE,pgaus),d2sdx(27)
    real(rp)                            :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)

    !if( lelch(ielem) == ELCUT ) then
    !   !
    !   ! Cut elements
    !   !
    !   !call elmcar_cut(&
    !   !     ltype(ielem),pnode,plapl,ielem,elcod,pgaus,gpvol,gpsha,gpder,gpcar,gphes)
    !   call runend('CUT ELEMENT WITH VECTORIZATION NOT READY YET')
    !
    !else
    !
    ! GPSHA, GPDER
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          gpsha(1:VECTOR_SIZE,inode,igaus) = shapf(inode,igaus)
          do idime = 1,ndime
             gpder(1:VECTOR_SIZE,idime,inode,igaus) = deriv(idime,inode,igaus) 
          end do
       end do
    end do

    if( kfl_savda == 2 ) then

       !--------------------------------------------------------------
       !
       ! Take values from element database: GPVOL, GPCAR, GPHES
       !
       !--------------------------------------------------------------

       do ivect = 1,VECTOR_SIZE             
          ielem = list_elements(ivect)             
          if( ielem > 0 ) then              
             gpvol(ivect,1:pgaus)                 = elmda(ielem) % gpvol(1:pgaus) 
             gpcar(ivect,1:ndime,1:pnode,1:pgaus) = elmda(ielem) % gpcar(1:ndime,1:pnode,1:pgaus)
             if( present(gphes) ) then
                if( plapl == 1 ) then
                   gphes(ivect,1:ntens,1:pnode,1:pgaus) = elmda(ielem) % gphes(1:ntens,1:pnode,1:pgaus) 
                else
                   gphes(ivect,1:ntens,1:pnode,1:pgaus) = 0.0_rp
                end if
             end if
          end if
       end do

    else

       !--------------------------------------------------------------
       !
       ! Compute sobre la marcha
       !
       !--------------------------------------------------------------

       !
       ! GPCAR, and GPVOL
       !
       call elmgeo_cartesian_derivatives_jacobian(ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

       do igaus = 1,pgaus
          gpvol(1:VECTOR_SIZE,igaus) = weigp(igaus) * gpdet(1:VECTOR_SIZE,igaus)
       end do

       if( present(gphes) ) then
          if( ( ndime == 2 .and. pnode == 3 .and. pgaus <= 3 ) .or. ( ndime == 3 .and. pnode == 4 .and. pgaus <= 4 ) ) then
             !
             ! No Hessian
             !
             gphes(1:VECTOR_SIZE,1:ntens,1:pnode,1:pgaus) = 0.0_rp
             
          else if( plapl == 0 ) then
             !
             ! No Hessian
             !         
             gphes(1:VECTOR_SIZE,1:ntens,1:pnode,1:pgaus) = 0.0_rp        
             
          else
             !
             ! GPHES
             !
             call elmgeo_hessian(ndime,mnode,ntens,pnode,pgaus,heslo,gphes,xjaci,deriv,elcod)
             
          end if
       end if
       
    end if

    !end if

  end subroutine element_shape_function_derivatives_jacobian_vector

  subroutine element_shape_function_derivatives_jacobian_scalar(&
       pnode,pgaus,plapl,weigp,shapf,deriv,heslo,&
       elcod,gpvol,gpsha,gpder,gpcar,gphes,&
       ielem)

    integer(ip), intent(in)             :: pnode
    integer(ip), intent(in)             :: plapl
    integer(ip), intent(in)             :: pgaus

    real(rp),    intent(in)             :: weigp(*)                       !< Weight of integration paints
    real(rp),    intent(in)             :: shapf(pnode,*)                 !< Shape functions at integration points
    real(rp),    intent(in)             :: deriv(ndime,pnode,*)           !< Shape function derivatives at integration points
    real(rp),    intent(in)             :: heslo(ntens,pnode,*)           !< Shape function Hessian at integration points

    real(rp),    intent(in)             :: elcod(ndime,pnode)
    real(rp),    intent(out)            :: gpvol(*)
    real(rp),    intent(out)            :: gpsha(pnode,*)
    real(rp),    intent(out)            :: gpder(ndime,pnode,*)
    real(rp),    intent(out)            :: gpcar(ndime,mnode,*)
    real(rp),    intent(out), optional  :: gphes(ntens,mnode,*)
    integer(ip), intent(in),  optional  :: ielem
    integer(ip)                         :: igaus,inode,idime
    integer(ip)                         :: itens,ivect
    real(rp)                            :: gpcod,gpdet(pgaus),d2sdx(27)
    real(rp)                            :: xjaci(ndime,ndime,pgaus)

    !if( lelch(ielem) == ELCUT ) then
    !   !
    !   ! Cut elements
    !   !
    !   !call elmcar_cut(&
    !   !     ltype(ielem),pnode,plapl,ielem,elcod,pgaus,gpvol,gpsha,gpder,gpcar,gphes)
    !   call runend('CUT ELEMENT WITH VECTORIZATION NOT READY YET')
    !
    !else
    !
    ! GPSHA, GPDER
    !
    gpsha(1:pnode,1:pgaus)         = shapf(1:pnode,1:pgaus)
    gpder(1:ndime,1:pnode,1:pgaus) = deriv(1:ndime,1:pnode,1:pgaus)

    if( kfl_savda == 2 .and. present(ielem) ) then

       !--------------------------------------------------------------
       !
       ! Take values from element database: GPVOL, GPCAR, GPHES
       !
       !--------------------------------------------------------------
       
       gpvol(1:pgaus)                 = elmda(ielem) % gpvol(1:pgaus) 
       gpcar(1:ndime,1:pnode,1:pgaus) = elmda(ielem) % gpcar(1:ndime,1:pnode,1:pgaus)
       if( present(gphes) ) then
          if( plapl == 1 ) then
             gphes(1:ntens,1:pnode,1:pgaus) = elmda(ielem) % gphes(1:ntens,1:pnode,1:pgaus) 
          else
             gphes(1:ntens,1:pnode,1:pgaus) = 0.0_rp
          end if
       end if

    else

       !--------------------------------------------------------------
       !
       ! Compute sobre la marcha
       !
       !--------------------------------------------------------------

       !
       ! GPCAR, and GPVOL
       !
       call elmgeo_cartesian_derivatives_jacobian(ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

       do igaus = 1,pgaus
          gpvol(igaus) = weigp(igaus) * gpdet(igaus)
       end do

       if( present(gphes) ) then
          if( ( ndime == 2 .and. pnode == 3 .and. pgaus <= 3 ) .or. ( ndime == 3 .and. pnode == 4 .and. pgaus <= 4 ) ) then
             !
             ! No Hessian
             !
             gphes(1:ntens,1:pnode,1:pgaus) = 0.0_rp

          else if( plapl == 0 ) then
             !
             ! No Hessian
             !         
             gphes(1:ntens,1:pnode,1:pgaus) = 0.0_rp        

          else
             !
             ! GPHES
             !
             call elmgeo_hessian(ndime,mnode,ntens,pnode,pgaus,heslo,gphes,xjaci,deriv,elcod)

          end if
       end if

    end if

    !end if

  end subroutine element_shape_function_derivatives_jacobian_scalar

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-30
  !> @brief   Compute an integration rule
  !> @details Compute integration rule, POSGP and WEIGP as well as
  !>          shape functions, derivatives and Hessian at the integration
  !>          points
  !> 
  !-----------------------------------------------------------------------

  subroutine element_integration(pelty,elmtyp,RULE,INTEGRATION_POINTS)
    
    integer(ip),            intent(in)    :: pelty               !< Element type
    type(elm),              intent(inout) :: elmtyp              !< Element integration
    character(*),           intent(in)    :: RULE                !< Integration rule name
    integer(ip),  optional, intent(in)    :: INTEGRATION_POINTS  !< Element type
    integer(ip)                           :: pdime,pnode,pgaus
    integer(ip)                           :: pquad,ptopo,mdime
    integer(ip)                           :: ierro,prule

    call memory_deallo(memor_dom,'SHAPE','cshder',elmtyp % shape)
    call memory_deallo(memor_dom,'DERIV','cshder',elmtyp % deriv)
    call memory_deallo(memor_dom,'HESLO','cshder',elmtyp % heslo)
    call memory_deallo(memor_dom,'POSGP','cshder',elmtyp % posgp)
    call memory_deallo(memor_dom,'WEIGP','cshder',elmtyp % weigp)

    pnode = element_type(pelty) % number_nodes
    pdime = element_type(pelty) % dimensions
    ptopo = element_type(pelty) % topology

    if( RULE(1:5) == 'CLOSE' ) then
       pquad = 1
       pgaus = pnode
    else
       pquad = 0
       if( present(INTEGRATION_POINTS) ) then
          pgaus = INTEGRATION_POINTS
       else
          pgaus = pnode
       end if
    end if
    prule = 2*ptopo+pquad+1     
    mdime = max(1_ip,pdime)
    
    call memory_alloca(memor_dom,'SHAPE','cshder',elmtyp % shape,pnode,pgaus)
    call memory_alloca(memor_dom,'DERIV','cshder',elmtyp % deriv,ndime,pnode,pgaus)
    call memory_alloca(memor_dom,'HESLO','cshder',elmtyp % heslo,ntens,pnode,pgaus)
    call memory_alloca(memor_dom,'POSGP','cshder',elmtyp % posgp,ndime,pgaus)
    call memory_alloca(memor_dom,'WEIGP','cshder',elmtyp % weigp,pgaus)

    call rulepw(pdime,pgaus,prule,elmtyp % posgp,elmtyp % weigp,ierro)
    call shafal(&
         elmtyp % posgp,pdime,pnode,pgaus,ntens,elmtyp % shape,elmtyp % deriv,&
         elmtyp % heslo,ierro)
    
  end subroutine element_integration
         
  subroutine element_integration_initialization(elmtyp)

    type(elm), intent(inout) :: elmtyp !< Element integration
    
    nullify(elmtyp % shape)
    nullify(elmtyp % deriv)
    nullify(elmtyp % heslo)
    nullify(elmtyp % posgp)
    nullify(elmtyp % weigp)
    
  end subroutine element_integration_initialization
  
end module mod_element_integration
!> @}
