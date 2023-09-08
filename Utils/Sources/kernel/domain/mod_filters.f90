!------------------------------------------------------------------------
!>
!> @addtogroup Projection_Toolbox
!> @{
!> @name    ToolBox for L2 projections 
!> @file    mod_projec.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for matrix operations
!> @details ToolBox for matrix operations: fill in, etc.
!>
!------------------------------------------------------------------------

module mod_filters

  use def_kintyp
  use def_master
  use def_parame
  use def_domain
  use mod_memory,               only : memory_alloca
  use mod_memory,               only : memory_deallo
  use mod_gradie,               only : gradie
  use mod_element_integration,  only : element_shape_function_derivatives_jacobian
  implicit none
  private
  integer(ip) :: knode
  integer(ip) :: ipoin,idime,inode,ielem,igaus,jnode
  integer(ip) :: pnode,pelty,pgaus,jdime,itens
  real(rp)    :: gpdet,gpvol
  real(rp)    :: xjaci(9),xjacm(9),xfact
  integer(8)  :: memor(2)
 
  public :: filters_nodal

contains

  subroutine filters_nodal(unkno,unkno_filtered)

    real(rp), intent(in)  :: unkno(npoin)
    real(rp), intent(out) :: unkno_filtered(npoin)
    real(rp), pointer     :: unkno_grad(:,:)
    integer(ip)           :: ipoin,ielpo,igaus,plapl
    integer(ip)           :: izdom,pnode,ielem,jpoin
    real(rp)              :: gpvol(mgaus)
    real(rp)              :: gpsha(mnode,mgaus)
    real(rp)              :: gpder(ndime,mnode,mgaus)
    real(rp)              :: gpcar(ndime,mnode,mgaus)
    real(rp)              :: gphes(ntens,mnode,mgaus)
    real(rp)              :: elcod(ndime,mnode)
    real(rp)              :: ellap(mnode)
    real(rp)              :: elgra(ndime,mnode)
    real(rp)              :: xfact,h,rdime,xmin,xmax

    real(rp), pointer     :: unkno_old(:)
    integer(ip)           :: ipass
    real(rp)              :: elunk(mnode)
    real(rp)              :: elave

    if( INOTMASTER ) then

       nullify(unkno_grad)
       nullify(unkno_old)
       call memory_alloca(memor_dom,'UNKNO_GRAD','memgeo',unkno_grad,ndime,npoin)
       call memory_alloca(memor_dom,'UNKNO_OLD','memgeo',unkno_old,npoin)
       call gradie(unkno,unkno_grad)
       plapl = 0
       rdime = 1.0_rp / real(ndime,rp)
       
       do ipass = 1,3

          if( ipass == 1 ) then
             unkno_old(1:npoin)      = unkno(1:npoin) 
          else
             unkno_old(1:npoin)      = unkno_filtered(1:npoin) 
             unkno_filtered(1:npoin) = 0.0_rp
          end if

          do ielem = 1,nelem

             pelty = ltype(ielem)

             if( pelty > 0 ) then

                pgaus = ngaus(pelty)
                pnode = nnode(pelty)
                elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
                elgra(1:ndime,1:pnode) = unkno_grad(1:ndime,lnods(1:pnode,ielem))
                elunk(1:pnode)         = unkno_old(lnods(1:pnode,ielem))

                call element_shape_function_derivatives_jacobian(&
                     pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                     elmar(pelty) % deriv,elmar(pelty) % heslo,&
                     elcod,gpvol,gpsha,gpder,gpcar,gphes,ielem)

                call filters_laplacian(pnode,pgaus,elgra,gpvol,gpsha,gpcar,ellap)

                !h = (sum(gpvol(1:pgaus))) ** rdime 
                !do inode = 1,pnode     
                !   ipoin = lnods(inode,ielem)
                !   unkno_filtered(ipoin) = unkno_filtered(ipoin) + ellap(inode) !* h * h
                !end do

                elave = sum(elunk(1:pnode))/real(pnode,rp)
                do inode = 1,pnode     
                   ipoin = lnods(inode,ielem)
                   do igaus = 1,pgaus
                      unkno_filtered(ipoin) = unkno_filtered(ipoin) + elave * gpsha(inode,igaus) * gpvol(igaus)
                   end do
                end do

             end if

          end do

          call rhsmod(1_ip,unkno_filtered)

          xmax = maxval(unkno(1:npoin))
          xmin = minval(unkno(1:npoin))

          do ipoin = 1,npoin
             unkno_filtered(ipoin) = unkno_filtered(ipoin) / vmass(ipoin)
             unkno_filtered(ipoin) = max(xmin,unkno_filtered(ipoin))
             unkno_filtered(ipoin) = min(xmax,unkno_filtered(ipoin))
          end do

       end do
       return

       !xfact = 1.0_rp / ( real(ndime,rp) + 2.0_rp )
       xfact = 1.0_rp / 24.0_rp
       do ipoin = 1,npoin
          h = vmass(ipoin) ** rdime
          unkno_filtered(ipoin) = unkno(ipoin) - xfact / vmass(ipoin) * unkno_filtered(ipoin) * h * h
          xmax = -huge(1.0_rp)
          xmin =  huge(1.0_rp)
          do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
             jpoin = c_dom(ipoin)
             xmax  = max(xmax,unkno(jpoin))
             xmin  = min(xmin,unkno(jpoin))
          end do
          unkno_filtered(ipoin) = max(xmin,unkno_filtered(ipoin))
          unkno_filtered(ipoin) = min(xmax,unkno_filtered(ipoin))
       end do

       call memory_deallo(memor_dom,'UNKNO_GRAD','memgeo',unkno_grad)

    end if

  end subroutine filters_nodal

  subroutine filters_laplacian(pnode,pgaus,elgra,gpvol,gpsha,gpcar,ellap)

    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: pgaus
    real(rp),    intent(in)  :: elgra(ndime,pnode)
    real(rp),    intent(in)  :: gpvol(pgaus)
    real(rp),    intent(in)  :: gpsha(pnode,pgaus)
    real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
    real(rp),    intent(out) :: ellap(pnode)
   integer(ip)               :: igaus,inode,idime
    real(rp)                 :: gpgra(ndime)

    ellap = 0.0_rp

    do igaus = 1,pgaus
       !
       ! Interpolate gradient at Gauss point
       !
       gpgra = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             gpgra(idime) = gpgra(idime) + gpsha(inode,igaus) * elgra(idime,inode)
          end do                                                                                                                                                                  
       end do
       do inode = 1,pnode
          do idime = 1,ndime
             ellap(inode) = ellap(inode) + gpcar(idime,inode,igaus) * gpgra(idime) * gpvol(igaus)
          end do
       end do

    end do

  end subroutine filters_laplacian

end module mod_filters
!> @}
