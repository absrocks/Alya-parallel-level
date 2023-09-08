!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmop0.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix and RHS for the pressure inicialization - very similar por_solite but only lapl & grav.
!> @details ORDER=1:\n
!!            Pressure equation, elemental operations:\n
!!              1. Compute elemental matrix and RHS \n
!!              2. Impose Dirichlet boundary conditions\n
!!              3. Assemble them\n
!!          ORDER=2:\n
!!            Update the subgrid scale - for the momemnt not ready\n
!> @} 
!------------------------------------------------------------------------
subroutine por_elmop0(order)
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper 
  use def_porous

  implicit none

  integer(ip), intent(in) :: order                     !> 2: compute SGS only

  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  integer(ip) :: ielem,igaus                           ! Indices and dimensions
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo

  real(rp)    :: elpre(mnode,ncomp_por)                ! Gather 
  real(rp)    :: elswa(mnode,ncomp_por)                ! Gather 
  real(rp)    :: elcod(ndime,mnode)

  real(rp)    :: tragl(ndime,ndime)                    ! Stabilization
  real(rp)    :: hleng(ndime)

  real(rp)    :: gpcod(ndime,mgaus)
  real(rp)    :: gpvol(mgaus)                          ! |J|*w
  real(rp)    :: gpdif(ndime,mgaus)                          
  real(rp)    :: gpden(mgaus) 
  real(rp)    :: gpgra(ndime,mgaus)                          ! gravity terms                      
  real(rp)    :: gppre(mgaus,ncomp_por)                ! Pressure
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gpsta(mgaus)                          ! tau
  real(rp)    :: dtmin
  real(rp)    :: dummr(1)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  !
  ! Initialization
  !
  do igaus = 1,mgaus
     gpdif(1:ndime,igaus) = 0.0_rp
     gpden(igaus) = 0.0_rp
     gpsta(igaus) = 0.0_rp
  end do
  !
  ! Loop over elements
  !  
  dtmin = 1.0e6
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)

     if( pelty > 0 ) then
        pnode = lnnod(ielem)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        !
        ! Gather operations
        !        
        call por_elmgat(&
             pnode,lnods(1,ielem),elpre,elswa,elcod,dummr,0_ip)
        !
        ! hleng and tragl at center of gravity
        !
        call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Local time step DTINV_POR
        !
        if(kfl_timco==2) then
            call runend('por_elmop0: local time step not ready')
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        gphes=0.0_rp
        call elmcar(&
             pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
             elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
             gphes,ielem)
        !
        ! Equation coefficients in GP - also gauss point values for press. & water sat.
        !
        call por_elmpre(&
             pnode,pgaus,&
             elmar(pelty) % shape,elpre,elswa,elcod,&
             gpcod,gpden,gpdif,gpgra,gppre,ielem)
       
        call por_elmpm0(&
             pnode,pgaus,gpdif,gpgra,&
             elmar(pelty) % shape,gpcar,gpvol,&
             elmat,elrhs)

        if( order == 1 ) then
           !
           ! Prescribe Dirichlet boundary conditions
           !
!           if( solve(kprsa_por) % kfl_iffix == 0 ) &   ! Fixity not controled by solver
!             call por_elmdir(    ! for the moment no dirich bcs for pressure
! for the initialization perhaps I would need a dirichlet BC  an other option would be to correct by a constant once I 
! have obtainded the result. For the moment I will test it withoun any of the the previous to see what happens.

           if( solve_sol(1)  %  kfl_algso == 9 ) then
                call runend('por_elmop0:   solve_sol(1)  %  kfl_algso == 9  not ready')

           else
              !
              ! Assembly
              !
              call assrhs(&
                   solve_sol(1) % ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
              call assmat(&
                   solve_sol(1) % ndofn,pnode,pnode,npoin,solve_sol(1) % kfl_algso,&
                   ielem,lnods(1,ielem),elmat,amatr)   
           end if

        end if

     end if

  end do elements

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine por_elmop0
