recursive subroutine rad_elmope(order)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_elmope
  ! NAME 
  !    rad_elmope
  ! DESCRIPTION
  !    ORDER=1:
  !      Radiation heat transfer equation, elemental operations:
  !      1. Compute elemental matrix and RHS 
  !      2. Impose Robin boundary conditions
  !      3. Assemble them
  !    ORDER=2:
  !      Update the subgrid scale
  ! USES
  ! USED BY
  !    rad_matrix
  !------------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_radiat
  implicit none

  integer(ip), intent(in) :: order                     ! =2: compute SGS only

  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  integer(ip) :: ielem,igaus,idime                     ! Indices and dimensions
  integer(ip) :: pelty,pmate,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo

  real(rp)    :: elrad(mnode,ncomp_rad)                ! Gather 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)  !!F Advection in the P1 model is zero

  real(rp)    :: tragl(ndime,ndime),chave(ndime,2)     ! Stabilization
  real(rp)    :: chale(2),hleng(ndime)

  real(rp)    :: gpdtc                                 ! Values at Gauss points
  real(rp)    :: gpvol(mgaus),gprea(mgaus)             ! |J|*w
  real(rp)    :: gpabs(mgaus),gpvel(ndime,mgaus)       ! r, a
  real(rp)    :: gpcon(mgaus),gpcod(ndime,mgaus)       ! k
  real(rp)    :: gpbbr(mgaus),gptem(mgaus)             ! Black body radiation term coming from particles
  real(rp)    :: gpdif(mgaus),gpgrd(ndime,mgaus)       ! k+kt, grad(k+kt)
  real(rp)    :: gppro(mgaus)                          ! Weighted residual L2-projection
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection
  real(rp)    :: gprhs(mgaus)                          ! f (all terms)
  real(rp)    :: gpden(mgaus)                          ! rho and then rho*cp
  real(rp)    :: gpsph(mgaus)                          ! cp
  real(rp)    :: gprad(mgaus,ncomp_rad)                ! T
  real(rp)    :: gpsou(mgaus)                          ! Q
  real(rp)    :: gpsgs(mgaus,2)                        ! G'
  real(rp)    :: gpsgv(ndime,mgaus)                    ! u'
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gplap(mnode,mgaus)                    ! Laplacian
  real(rp)    :: gpsta(mgaus)                          ! tau
  integer(ip) :: dummi,i,j,k
  real(rp)    :: dummr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  !
  ! Initialization
  !
  do igaus = 1,mgaus
     gpdif(igaus) = 0.0_rp
     gpden(igaus) = 1.0_rp
     gprea(igaus) = 0.0_rp
     gpsta(igaus) = 0.0_rp
     gppro(igaus) = 0.0_rp
     gpdiv(igaus) = 0.0_rp
     gprhs(igaus) = 0.0_rp
     do idime = 1,ndime
        gpgrd(idime,igaus) = 0.0_rp
        gpvel(idime,igaus) = 0.0_rp
     end do
  end do
  !
  ! Loop over elements
  !
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
        ! Subgrid scale. radiation: GPSGS 
        !
        call rad_sgsope(&
             1_ip,ielem,pgaus,gpsgs,dummr)   
        !
        ! Gather operations
        !
        call rad_elmgat(&
             pnode,lnods(1,ielem),elrad,elcod) 
        !
        ! hleng and tragl at center of gravity   
        !
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&  !!F   WHAT IS RELEVANCE OF THESE FOR RADIAT
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE   !!F   WHAT IS RELEVANCE OF THESE FOR RADIAT SINCE ELVEL=0
        !
        call elmchl(tragl,hleng,elcod,elvel,chave,chale,pnode,&
             porde,hnatu(pelty),0_ip,kfl_ellen_rad)
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
             gphes,ielem)
        !
        ! Radiation: Gather from element to Gauss points GPRAD+GPSGS.
        !
        call gather(&
             2_ip,pgaus,pnode,1_ip,dummi,&
             elmar(pelty)%shape,elrad,gprad)
        call rad_sgsope(&
             2_ip,ielem,pgaus,gpsgs,gprad)
        !
        !  Properties: Prepare absorption and scattering coefficients. GPDIF (Gamma) and GPREA
        !
        call rad_elmpro(&
             ielem,pmate,pnode,pgaus,1_ip,pgaus,&
             elmar(pelty)%shape,gpcar,gpdif,gpabs,gpbbr,gptem,gpgrd)   
        !
        ! Prepare Equation coefficients
        !
        call rad_elmpr2(&
             pnode,pgaus,pmate,gpsgs,&
             elmar(pelty)%shape,elrad,&
             elcod,gppro,gpcod,&
             lnods(1,ielem))
        !
        ! Exact solution: GPRHS
        !
        call rad_elmexa(&
             pgaus,gpcod,gpden,gpdif,gprea,gpgrd,&
             gpvel,gpsou,gprhs)
        !
        ! Without exact solution, RHS is just the black body radiation
        !
        do igaus = 1,mgaus
           gprhs(igaus) =  4.0_rp* steph_rad * gpabs(igaus) * gptem(igaus)**4 + gpbbr(igaus)  
        enddo
        !
        ! Reaction term is the absorption coefficient
        !
        do igaus = 1,mgaus
           gprea(igaus) = gpabs(igaus)
        enddo        
        !
        ! Subgrid scale residual and update TESGS
        !
        !call rad_sgsope(&
        !     3_ip,ielem,pgaus,gpsgs,dummr)
        call elmadr(& 
             order,ielem,pnode,pgaus,ptopo,plapl,pelty,1_ip,1_ip,lnods(1,ielem),0_ip,& ! 0 is kfl_shock_rad
             kfl_taust_rad,kfl_ortho_rad,0_ip,0_ip,kfl_limit_rad,staco_rad,gpdif,& ! 0 is kfl_sgsti_rad
             gprea,gpden,gpvel,gppro,gpvol,gpgrd,gprhs,elmar(pelty)%shape,gpcar,gphes,&
             elrad,elcod,chale,gpsta,dtinv_rad,0.0_rp,0.0_rp,gpdiv,gplap,gpsgs,& !0.0 is shock parameter and bemol
             elmat,elrhs)

        if( order == 1 ) then
           !
           ! Assembly
           !
           call assrhs(&
                solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
           call assmat(&
                solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
                ielem,lnods(1,ielem),elmat,amatr)   

        end if

     end if

  end do elements

!!$ !Output for Mathematica
!!$print *,'BULK Matrix:'
!!$
!!$print *,'{'
!!$k=1_ip
!!$do i=1,(size(r_sol)-1)
!!$   do j=r_sol(i),(r_sol(i+1)-1)
!!$   print *,'{',i,',',c_sol(j),'} ->',amatr(k),','
!!$   k=k+1_ip
!!$end do
!!$end do
!!$print *,'}'
!!$
!!$print *,'RHS :'
!!$print *,'{'
!!$do j=1,size(rhsid,1)
!!$      print *,rhsid(j),','
!!$enddo
!!$print *,'}'

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine rad_elmope


