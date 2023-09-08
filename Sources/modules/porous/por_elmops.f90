!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmops.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix and RHS for the water saturation
!> @details ORDER=1:
!!            Satuartion equation, elemental operations:
!!              1. Compute elemental matrix and RHS (non transient part only - K)
!!              2. multiply by elemental values from n-1 and asemble in rhs
!!              3. add transient terms using lumped mass matrix
!!          ORDER=2:
!!            Update the subgrid scale - for the momemnt not ready
!> @} 
!------------------------------------------------------------------------
subroutine por_elmops(order)
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper 
  use def_porous
  use def_solver
  implicit none

  integer(ip), intent(in) :: order                     !> 2: compute SGS only

  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  integer(ip) :: ielem,igaus,inode,jnode                           ! Indices and dimensions
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo

  real(rp)    :: elpre(mnode,ncomp_por)                ! Gather 
  real(rp)    :: elswa(mnode,ncomp_por)                ! Gather 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)

  real(rp)    :: tragl(ndime,ndime)     ! Stabilization
  real(rp)    :: chale(2)
  real(rp)    :: hleng(ndime)

  real(rp)    :: gpcod(ndime,mgaus)
  real(rp)    :: gpvol(mgaus)                          ! |J|*w
  real(rp)    :: gprhs(mgaus)
  real(rp)    :: gpdif(ndime,mgaus)                          
  real(rp)    :: gpadv(ndime,mgaus)                          
  real(rp)    :: gpden(mgaus) 
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gpsta(mgaus)                          ! tau
  real(rp)    :: gplap(mnode,mgaus)                    ! Laplacian - as in tem it is calculated by elmadr but is not used ??
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection - only used for int by part of convection
                                                       ! Something that we will not use for porous
  real(rp)    :: gpsgs(mgaus,2)                        ! Sw'  - not used becaus subgrid tracking not ready

  real(rp)    :: gpdif_aux(mgaus),gpgrd_aux(ndime,mgaus) 
  real(rp)    :: gprea_aux(mgaus)
  real(rp)    :: gpden_aux(mgaus)
  real(rp)    :: gppro_aux(mgaus)
  real(rp)    :: dtinv_aux
  real(rp)    :: dtmin,dummr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  if (order /= 1) call runend('por_elmops: for the moment only order == 1 is ready')
  gpdif_aux = 0.0_rp
  gprea_aux = 0.0_rp
  gpden_aux = 1.0_rp
  gppro_aux = 0.0_rp   ! for OSS - for the moment we leave it to 0
  gpgrd_aux = 0.0_rp   ! diffusion gradient - it is zero because the diffusion is zero
  dtinv_aux = 0.0_rp   ! set to zero so that the transient terms are not assembled and we can do an explict with elmadr
  if ( bemol_por /= 0 ) call runend ('por_elmops:  bemol_por /= 0 not ready')
  !
  ! Initialization
  !
  do igaus = 1,mgaus
     gpdif(1:ndime,igaus) = 0.0_rp
     gpden(igaus) = 0.0_rp
     gpsta(igaus) = 0.0_rp   ! elmadr calculates it
  end do
  !
  ! Obtains aditional terms to lumpped mass matrix using the nodal values of porosity, Bw and dT - stored in pmatr
  !
  if(  solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON .and. &
       solve_sol(1) % kfl_preco == SOL_LOCAL_DIAGONAL) call por_massma()  
  !
  ! Loop over elements
  !  
  dtmin = 1.e6
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
             pnode,lnods(1,ielem),elpre,elswa,elcod,elvel,1_ip)
        !
        ! hleng and tragl at center of gravity
        !
        call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE
        !
        if ( (kfl_ellen_por < 0_ip) .or. (kfl_ellen_por > 2_ip) )   &
             call runend('por_elmops:for the moment only 0 min, 1 max or 2 ave are valid')
        call elmchl(tragl,hleng,dummr,dummr,dummr,chale,pnode,&
             porde,hnatu(pelty),0_ip,kfl_ellen_por)
        !
        ! Local time step DTINV_POR
        !
        if(kfl_timco==2) then
            call runend('por_elmops: locat time step not ready')
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
        call por_elmprs(&
             pnode,pgaus,&
             elmar(pelty) % shape,gpcar,elpre,elswa,elvel,elcod,&
             gpcod,gpadv,gprhs,ielem)
        !
        ! Add diffusion at the wells 
        !
        if ( kfl_difwe_por /=0 ) call por_difwel(&
             pnode,pgaus,lnods(1,ielem),gpdif_aux)
        !
        ! Equation coefficients in GP
        !
        kfl_ortho_por(2)=-1 ! OJO HERBERT QUE HE PUESTO SU AND NOT ASGS
        
        if( solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
           dtinv_aux = 0.0_rp ! Useless because it cancels in por_elswar
        else
           dtinv_aux = dtinv_por 
        end if
        !
        ! Assembly
        !
        call elmadr(& 
             one,ielem,pnode,pgaus,ptopo,plapl,pelty,1_ip,1_ip,lnods(1,ielem),kfl_shock_por(2),&
             kfl_taust_por(2),kfl_ortho_por(2),0_ip,kfl_sgsti_por(2),kfl_limit_por(2),staco_por,gpdif_aux,&
             gprea_aux,gpden_aux,gpadv,gppro_aux,gpvol,gpgrd_aux,gprhs,elmar(pelty) % shape,gpcar,gphes,&
             elswa,elcod,chale,gpsta,dtinv_aux,shock_por,bemol_por,gpdiv,gplap,gpsgs,&
             elmat,elrhs)
!            elswa I understand must not be used it dtinv=0 but I send it anyway  - it should do nothing with it
        
        if( order == 1 ) then

           if( solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
              !
              ! Send elmat * elswa(2) to rhs
              !
              call por_elswar(pnode,ncomp_por,elrhs,elmat,elswa)
              !
              ! obtains aditional terms to lumpped mass matrix using the nodal values of porosity, Bw and dT - stored in pmatr
              !
              call assrhs(&
                   solve(kprsa_por) % ndofn,pnode,lnods(1,ielem),elrhs,rhsid)   

           else
              !
              ! obtains aditional terms to lumpped mass matrix using the nodal values of porosity, Bw and dT - stored in pmatr
              !
              call assrhs(&
                   solve(kprsa_por) % ndofn,pnode,lnods(1,ielem),elrhs,rhsid)      
              call assmat(&
                   solve(kprsa_por) % ndofn,pnode,pnode,npoin,solve(kprsa_por) % kfl_algso,&
                   ielem,lnods(1,ielem),elmat,amatr)

           end if

        end if

     end if

  end do elements
#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine por_elmops
