!------------------------------------------------------------------------
!> @addtogroup Matrix_Toolbox
!> @{
!> @name    ToolBox for matrix operations
!> @file    mod_matrix.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for matrix operations
!> @details ToolBox for matrix operations: fill in, etc.
!------------------------------------------------------------------------

module mod_matrices
  use def_kintyp,              only : ip,rp,lg,i1p,r1p
  use mod_memory,              only : memory_alloca,memory_deallo
  use def_master,              only : INOTMASTER,kfl_paral
  use def_kermod,              only : kfl_element_to_csr
  use def_domain,              only : elmar
  use def_domain,              only : mnode
  use def_domain,              only : mgaus
  use def_domain,              only : ndime
  use def_domain,              only : nelem
  use def_domain,              only : nzdom
  use def_domain,              only : r_dom
  use def_domain,              only : c_dom
  use def_domain,              only : ltype
  use def_domain,              only : lgaus
  use def_domain,              only : lezdo
  use def_domain,              only : ntens
  use def_domain,              only : lnods
  use def_domain,              only : lnnod
  use def_domain,              only : coord
  use def_domain,              only : memor_dom
  use mod_element_integration, only : element_shape_function_derivatives_jacobian 
  implicit none

  private

  public :: matrices_gradient_divergence
  public :: matrices_laplacian
  
contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/05/2017
  !> @brief   Compute the gradient matrix
  !> @details Compute the gradient matrix
  !>
  !----------------------------------------------------------------------

  subroutine matrices_gradient_divergence(Grad,Div,kdime,DEALLOCATE_MATRICES) 

    real(rp),              pointer, intent(inout) :: Grad(:,:)
    real(rp),    optional, pointer, intent(inout) :: Div(:,:)
    integer(ip), optional,          intent(in)    :: kdime
    logical(lg), optional,          intent(in)    :: DEALLOCATE_MATRICES
    integer(ip)                                   :: ielem
    integer(ip)                                   :: pnode
    integer(ip)                                   :: plapl
    integer(ip)                                   :: pgaus
    integer(ip)                                   :: pelty
    integer(ip)                                   :: inode

    real(rp)                                      :: elcod(ndime,mnode)
    real(rp)                                      :: gpvol(mgaus)
    real(rp)                                      :: gpsha(mnode,mgaus)
    real(rp)                                      :: gpder(ndime,mnode,mgaus)
    real(rp)                                      :: gpcar(ndime,mnode,mgaus)
    real(rp)                                      :: gphes(ntens,mnode,mgaus)
    logical(lg)                                   :: if_deallocate_matrices

    
    if( INOTMASTER ) then

       if( present(DEALLOCATE_MATRICES) ) then
          if_deallocate_matrices = DEALLOCATE_MATRICES
       else
          if_deallocate_matrices = .false.
       end if

       if( if_deallocate_matrices ) then
          call memory_deallo(memor_dom,'GRAD','matrices_gradient',Grad)
          if( present(Div) ) then
             call memory_deallo(memor_dom,'Div','matrices_gradient',Div)
          end if
       end if 
       
       if( .not. associated(Grad) ) then
          call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,ndime,nzdom)
       else
          Grad = 0.0_rp
       end if
       if( present(Div) ) then
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'Div','matrices_gradient',Div,ndime,nzdom)
          else
             Div = 0.0_rp
          end if
       end if
       
       do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pgaus = lgaus(ielem)
             pnode = lnnod(ielem)
             plapl = 0
             do inode = 1,pnode
                elcod(1:ndime,inode) = coord(1:ndime,lnods(inode,ielem))
             end do
             call element_shape_function_derivatives_jacobian(&
                  pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                  elmar(pelty) % deriv,elmar(pelty) % heslo,&
                  elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
             if( present(Div) ) then
                call matrices_element_gradient(&
                     pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,&
                     Grad,Div)
             else
                call matrices_element_gradient(&
                     pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpsha,gpcar,&
                     Grad)
             end if
          end if
       end do


    else

       !
       ! Allocate of size 1 otherwise I got : Attempt to use pointer GRAD when it is not associated with a target
       !
       if( .not. associated(Grad) ) then
          call memory_alloca(memor_dom,'GRAD','matrices_gradient',Grad,ndime,1_ip)
       else
          Grad = 0.0_rp
       end if
       if( present(Div) ) then
          if( .not. associated(Div) ) then
             call memory_alloca(memor_dom,'Div','matrices_gradient',Div,ndime,1_ip)
          else
             Div = 0.0_rp
          end if
       end if
       
    end if
    
  end subroutine matrices_gradient_divergence

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/05/2017
  !> @brief   Element and global gradient matrices
  !> @details Element and global gradient matrices
  !>
  !----------------------------------------------------------------------

  subroutine matrices_element_gradient(&
       pgaus,pnode,ielem,lnods,gpvol,gpsha,gpcar,Grad,Div)

    integer(ip),                    intent(in)    :: pgaus
    integer(ip),                    intent(in)    :: pnode
    integer(ip),                    intent(in)    :: ielem
    integer(ip),                    intent(in)    :: lnods(pnode)
    real(rp),                       intent(in)    :: gpvol(pgaus)
    real(rp),                       intent(in)    :: gpsha(pnode,pgaus)
    real(rp),                       intent(in)    :: gpcar(ndime,mnode,pgaus)
    real(rp),    pointer,           intent(inout) :: Grad(:,:)
    real(rp),    pointer, optional, intent(inout) :: Div(:,:)
    real(rp)                                      :: elgra(ndime*pnode,pnode) ! Aup
    real(rp)                                      :: eldiv(pnode,ndime*pnode) ! Apu
    integer(ip)                                   :: igaus,inode,idime,idofn,iz
    integer(ip)                                   :: jcolu,jnode,ipoin,jpoin
    integer(ip)                                   :: jdofn
    !
    ! Element matrix
    !
    elgra = 0.0_rp
    eldiv = 0.0_rp
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             idofn = (inode-1)*ndime+idime
             do jnode = 1,pnode
                elgra(idofn,jnode) = elgra(idofn,jnode) &
                     - gpvol(igaus) * gpsha(jnode,igaus) * gpcar(idime,inode,igaus)
                eldiv(jnode,idofn) = eldiv(jnode,idofn) &
                     + gpvol(igaus) * gpsha(jnode,igaus) * gpcar(idime,inode,igaus)
             end do
          end do
       end do
    end do

    !do inode = 1,pnode*ndime
    !   write(*,'(10(1x,e13.6))') elgra(idime,inode,1:pnode)
    !end do
    !call runend('O.K.!')
    !
    ! Assemble global system
    !
    if( kfl_element_to_csr == 1 ) then
      
       do inode = 1,pnode
          do jnode = 1,pnode
             iz = lezdo(inode,jnode,ielem)
             do idime = 1,ndime
                idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                Grad(idime,iz) = Grad(idime,iz) + elgra(idofn,jnode) ! G
             end do
             if( present(Div) ) then
                do idime = 1,ndime
                   jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                   Div(idime,iz) = Div(idime,iz) + eldiv(inode,jdofn) ! D                      
                end do
             end if
          end do
       end do
       
    else

          do inode = 1,pnode
             ipoin = lnods(inode)
             do jnode = 1,pnode
                jpoin = lnods(jnode)
                iz    = r_dom(ipoin)
                jcolu = c_dom(iz)
                do while( jcolu /= jpoin .and. iz < r_dom(ipoin+1)-1 )
                   iz    = iz + 1
                   jcolu = c_dom(iz)
                end do
                if( jcolu == jpoin ) then
                   do idime = 1,ndime
                      idofn = (inode-1) * ndime + idime
#ifdef NO_COLORING
                      !$OMP ATOMIC
#endif
                      Grad(idime,iz) = Grad(idime,iz) + elgra(idofn,jnode)  ! G
                   end do
                   if( present(Div) ) then
                      do idime = 1,ndime
                         jdofn = (jnode-1) * ndime + idime
#ifdef NO_COLORING
                         !$OMP ATOMIC
#endif
                         Div(idime,iz) = Div(idime,iz) + eldiv(inode,jdofn) ! D
                      end do
                   end if
                end if
             end do
          end do
       
    end if

  end subroutine matrices_element_gradient


  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    27/07/2017
  !> @brief   Compute the laplacian matrix
  !> @details Compute the laplacian matrix
  !>
  !----------------------------------------------------------------------

  subroutine matrices_laplacian(Lapl,DEALLOCATE_MATRIX) 

    real(rp),              pointer, intent(inout) :: Lapl(:)
    logical(lg), optional, intent(in)             :: DEALLOCATE_MATRIX

    integer(ip)                                 :: ielem
    integer(ip)                                 :: pnode
    integer(ip)                                 :: plapl
    integer(ip)                                 :: pgaus
    integer(ip)                                 :: pelty
    integer(ip)                                 :: inode

    real(rp)                                    :: elcod(ndime,mnode)
    real(rp)                                    :: gpvol(mgaus)
    real(rp)                                    :: gpsha(mnode,mgaus)
    real(rp)                                    :: gpder(ndime,mnode,mgaus)
    real(rp)                                    :: gpcar(ndime,mnode,mgaus)
    real(rp)                                    :: gphes(ntens,mnode,mgaus)
    logical(lg)                                 :: if_deallocate_matrices

    if( INOTMASTER ) then

       if( present(DEALLOCATE_MATRIX) ) then
          if_deallocate_matrices = DEALLOCATE_MATRIX
       else
          if_deallocate_matrices = .false.
       end if

       if( if_deallocate_matrices ) then
          call memory_deallo(memor_dom,'LAPL','matrices_gradient',Lapl)
       end if
       
       if( .not. associated(Lapl) ) then
          call memory_alloca(memor_dom,'LAPL','matrices_gradient',Lapl,nzdom)
       else
          Lapl = 0.0_rp
       end if

       
       do ielem = 1,nelem
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pgaus = lgaus(ielem)
             pnode = lnnod(ielem)
             plapl = 0
             do inode = 1,pnode
                elcod(1:ndime,inode) = coord(1:ndime,lnods(inode,ielem))
             end do
             call element_shape_function_derivatives_jacobian(&
                  pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                  elmar(pelty) % deriv,elmar(pelty) % heslo,&
                  elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)
             
             call matrices_element_laplacian(&
                  pgaus,pnode,ielem,lnods(:,ielem),gpvol,gpcar,Lapl)
            
          end if
       end do

    else

       !
       ! Allocate of size 1 otherwise I got : Attempt to use pointer L when it is not associated with a target
       !
       if( .not. associated(Lapl) ) then
          call memory_alloca(memor_dom,'LAPL','matrices_gradient',Lapl,1_ip)
       else
          Lapl = 0.0_rp
       end if
       
    end if
    
  end subroutine matrices_laplacian


  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    27/07/2017
  !> @brief   Element and global Laplacian matrix ( grad p , grad q )
  !> @details Element and global Laplacian matrix ( grad p , grad q )
  !>
  !----------------------------------------------------------------------
  
  subroutine matrices_element_laplacian(&
       pgaus,pnode,ielem,lnods,gpvol,gpcar,Lapl)

    integer(ip),                    intent(in)    :: pgaus
    integer(ip),                    intent(in)    :: pnode
    integer(ip),                    intent(in)    :: ielem
    integer(ip),                    intent(in)    :: lnods(pnode)
    real(rp),                       intent(in)    :: gpvol(pgaus)
    real(rp),                       intent(in)    :: gpcar(ndime,mnode,pgaus)
    real(rp),    pointer,           intent(inout) :: Lapl(:)
    real(rp)                                      :: ellap(pnode,pnode) ! L
    integer(ip)                                   :: igaus,inode,idime,iz
    integer(ip)                                   :: jcolu,jnode,ipoin,jpoin
    !
    ! Element matrix
    !
    ellap = 0.0_rp
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                ellap(inode,jnode) = ellap(inode,jnode) &
                     + gpvol(igaus) * gpcar(idime,inode,igaus) * gpcar(idime,jnode,igaus)  ! We use a + sign as in nsi_element_schur
                                ! In JCP 2001 - Codina - Pressuer Stab..   it is defined with a - sign.

             end do
          end do
       end do
    end do

    !do inode = 1,pnode*ndime
    !   write(*,'(10(1x,e13.6))') ellap(inode,1:pnode)
    !end do
    !call runend('O.K.!')
    !
    ! Assemble global system
    !
    if( kfl_element_to_csr == 1 ) then

       do inode = 1,pnode
          do jnode = 1,pnode
             iz = lezdo(inode,jnode,ielem)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             Lapl(iz) = Lapl(iz) + ellap(inode,jnode) ! L
          end do
       end do

    else

       do inode = 1,pnode
          ipoin = lnods(inode)
          do jnode = 1,pnode
             jpoin = lnods(jnode)
             iz    = r_dom(ipoin)
             jcolu = c_dom(iz)
             do while( jcolu /= jpoin .and. iz < r_dom(ipoin+1)-1 )
                iz    = iz + 1
                jcolu = c_dom(iz)
             end do
             if( jcolu == jpoin ) then
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                Lapl(iz) = Lapl(iz) + ellap(inode,jnode)  ! L
             end if
          end do
       end do

    end if

  end subroutine matrices_element_laplacian

  !
  ! Ojo hay una parecida en nsi_element_laplacian
  !  
  
end module mod_matrices
!> @}
