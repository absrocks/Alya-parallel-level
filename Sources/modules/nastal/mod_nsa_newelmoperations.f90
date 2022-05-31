module mod_nsa_newelmoperations
  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    04/06/2015
  !> @brief   Module containing elementary operations
  !> @details Module containing elementary operations
  !> @} 
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_kermod
  use      mod_ker_proper

  implicit none

  private
  real(rp) :: &
       chale_elm(2),qufac_elm,hleng_elm(3),&
       hunkn_elm(5,3,3),hunkn_newton_elm(5,3,3),xjaci_elm(3,3),&
       hmini_elm,hmaxi_elm,shmet_elm(3,3,5)

  real(rp) :: &
       epsilon_newton(5) = 0.000001_rp

  real(rp) , pointer      :: &  
       elrhs(:)                 ! nevat_nsa 
  real(rp) , pointer      :: &  
       elmat(:,:)               ! nevat_nsa , nevat_nsa 

  type elm_onnode_gather_nsa
     real(rp)               ::      dtinv_eqs(5,2)=0.0_rp             ! ndofn_nsa, 2                  , mnode 
     real(rp)               ::      eldtt(5,2)=0.0_rp             ! ndofn_nsa, 2                  , mnode 
     real(rp)               ::      elunk(5,10)=0.0_rp             ! ndofn_nsa , mnode , ncomp_nsa (10 max) 
     real(rp)               ::      elsub(5)=0.0_rp               ! ndofn_nsa , mnode 
     real(rp)               ::      elbve(5) =0.0_rp              ! ndofn_nsa , mnode 
     real(rp)               ::      elphy(5)=0.0_rp               ! ndofn_nsa , mnode 
     real(rp)               ::      elcod(3)=0.0_rp               ! ndime , mnode
     real(rp)               ::      elvel(3)=0.0_rp              ! ndime , mnode
     real(rp)               ::      elmsh(3)=0.0_rp               ! ndime , mnode
     real(rp)               ::      elpre=0.0_rp                 ! mnode
     real(rp)               ::      eltem =0.0_rp                 ! mnode
     real(rp)               ::      elvis=0.0_rp                  ! mnode
     real(rp)               ::      elhcp=0.0_rp                ! mnode
     real(rp)               ::      elwme=0.0_rp              ! mnode     
  end type elm_onnode_gather_nsa

  type elm_onnode_comput_nsa
     real(rp)               ::      eldif(5,5,3,3)=0.0_rp         !(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode)
     real(rp)               ::      elcon(5,5,3)=0.0_rp          !(ndofn_nsa,ndofn_nsa,ndime,      mnode)
  end type elm_onnode_comput_nsa

  type elm_ongaus_interp_nsa
     real(rp)               :: xconv_der(5,5,5,3) = 0.0_rp       !(ndofn_nsa,ndofn_nsa,ndofn_nsa,ndime,mgaus),&
     real(rp)               :: xdiff(5,5,3,3) = 0.0_rp          !(ndofn_nsa,ndofn_nsa,ndime,ndime,mgaus),
     real(rp)               :: ddiff(5,5,3,2) = 0.0_rp           !(ndofn_nsa,ndofn_nsa,ndime,    2,mgaus),
     real(rp)               :: xconv(5,5,3) = 0.0_rp             !(ndofn_nsa,ndofn_nsa,ndime,mgaus),&     
     real(rp)               :: xvofo(5,5) = 0.0_rp              !(ndofn_nsa,ndofn_nsa,mgaus), &          
     real(rp)               :: dconv(5,5) = 0.0_rp              !(ndofn_nsa,ndofn_nsa,mgaus), &          
     real(rp)               :: gunkn(5,3) = 0.0_rp               !(ndofn_nsa,ndime,mgaus),&               
     real(rp)               :: gtunk(5,3) = 0.0_rp               !(ndofn_nsa,ndime,mgaus),&               
     real(rp)               :: gsube(5,3) = 0.0_rp              !(ndofn_nsa,ndime,mgaus), &              
     real(rp)               :: xsube(5,3) = 0.0_rp               !(ndofn_nsa,3, mgaus),&               
     real(rp)               :: xdtix(5,3) = 0.0_rp              !(ndofn_nsa,2, mgaus), &              
     real(rp)               :: xtime(5) = 0.0_rp                !(ndofn_nsa,mgaus), &                    
     real(rp)               :: xresi(5) = 0.0_rp               ! (ndofn_nsa,mgaus), &     
     real(rp)               :: xunkn(5,3) = 0.0_rp             ! (ndofn_nsa,3,mgaus), &   
     real(rp)               :: xtunk(5) = 0.0_rp               ! (ndofn_nsa,mgaus), &   
     real(rp)               :: xtide(5) = 0.0_rp               ! (ndofn_nsa,mgaus), &     
     real(rp)               :: taudi(5) = 0.0_rp               ! (ndofn_nsa,mgaus), &     
     real(rp)               :: shocktau_local(5) = 0.0_rp      ! (ndofn_nsa,mgaus), &     
     real(rp)               :: xlopr_conservative(5,5) = 0.0_rp ! (ndofn_nsa,ndofn_nsa,mgaus)
     real(rp)               :: gvelo(3,3) = 0.0_rp              ! (ndime,ndime,mgaus),&    
     real(rp)               :: xldve(3) = 0.0_rp                    ! (ndime,mgaus), &         
     real(rp)               :: xvelo(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: xvmsh(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: gpres(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: gtemp(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: gvisc(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: htrad(3) = 0.0_rp                    ! (ndime,mgaus),&          
     real(rp)               :: dvolu = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: dvelo = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xsoun = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xpres = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xtemp = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xvisc = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xdith = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xlade = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: velmo = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: dhtra = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: sgsdi = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: heats = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xmowe = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xadgam = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xrgacv = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xrgasc = 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xheatcp= 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xheatcv= 0.0_rp                          ! (mgaus),&                
     real(rp)               :: xnutu  = 0.0_rp                          !(mgaus)                                   
  end type elm_ongaus_interp_nsa

  type elm_onnode_onnode_matrix_nsa
     real(rp)               ::         advec_matrix(5,5)= 0.0_rp ! Advective part
     real(rp)               ::         diffu_matrix(5,5)= 0.0_rp ! Diffusive part
     real(rp)               ::         shote_matrix(5,5)= 0.0_rp ! Shock capturing part
     real(rp)               ::         stabi_matrix(5,5)= 0.0_rp ! Stabilization part = Adjoint matrix * Diagonal subscale part
     real(rp)               ::         subdi_matrix(5,5)= 0.0_rp ! Diagonal subscale part
     real(rp)               ::         timas_matrix(5,5,2)= 0.0_rp ! Temporal part for time and pseudotime 
  end type elm_onnode_onnode_matrix_nsa


  !
  ! Public stuff
  ! 
  public nsa_newelmoperations
  

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Perform elementary operations
  !> @details Perform elementary operations
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmoperations
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    type(elm_onnode_gather_nsa)        :: gath(mnode,2)
    type(elm_onnode_comput_nsa)        :: compu(mnode)
    type(elm_ongaus_interp_nsa)        :: gaus(mgaus),gaus_init
    type(elm_onnode_onnode_matrix_nsa) :: matri(mnode,mnode)

    integer(ip) :: ielem,inode,igaus,&
         pelty,pnode,pgaus,plapl,pface,pevat,ipoin,kpoin,icount_newton,ncount_newton

    real(rp)    :: elrhs(nevat_nsa),elmat(nevat_nsa,nevat_nsa)
    ! recall that nevat_nsa = ndofn_nsa * mnode

    ! FEM-data computed from parent domain
    real(rp) :: &
         cartd(ndime,mnode,mgaus),hessi(ntens,mnode,mgaus),xjacm(ndime,ndime),&
         xjaci(ndime,ndime),tragl(ndime,ndime),elcod_local(ndime,mnode),detjm,&
         d2sdx(ndime,ndime,ndime)

    !
    ! Initialise elementary loop
    !
    call nsa_newelmstart

    !
    ! Elementary loop
    !
    elements_loop: do ielem = 1,nelem

       ! Element properties and dimensions
       pelty=ltype(ielem)
       pnode=nnode(pelty)
       pgaus=ngaus(pelty)
       plapl=llapl(pelty)
       pface=nface(pelty)
       pevat = ndofn_nsa*pnode


       qufac_elm = 1.0_rp
       if((ndime.eq.2).and.(pnode.ge.4)) then
          qufac_elm = 0.5_rp 
       end if
       if((ndime.eq.3).and.(pnode.ge.5))then
          qufac_elm = 0.5_rp 
       end if
       
       !       
       ! Compute hleng and tragl at center of gravity
       !       

       call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod_local,hnatu(pelty),hleng_elm)

       ! default values for chale
       chale_elm(1) = hleng_elm(ndime)      ! smallest
       chale_elm(2) = hleng_elm(1)          ! largest

       gaus(1:mgaus) = gaus_init   ! Initialize gaus database to gaus_init, which is zero

       hessi(:,1:mnode,1:mgaus) = 0.0_rp

       
       ! 
       ! Gather values from global to local vectors
       ! 

       call nsa_newelmgather(ielem,pnode,gath(1,1),elcod_local)

       if (kfl_linea_nsa == 2) then
          gath(:,2)= gath(:,1)  ! set the initial values
          call nsa_newelmgather_newton(ielem,pnode,gath(1,2))
       end if

       ! ojo que de momento esto lo voy a hacer solo para la primera pasadita, pero si la malla cambia
       ! puede que tenga que hacerlo para ls dos

       elemental_gauss_points_geometry: do igaus=1,pgaus        

          !
          ! Compute test functions' stuff from the parametric space
          !
          call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod_local,cartd(1,1,igaus),detjm,xjacm,xjaci)
          gaus(igaus)%dvolu=elmar(pelty)%weigp(igaus)*detjm                          
          if(plapl==1) call elmhes(&
               elmar(pelty)%heslo(1,1,igaus),hessi(1,1,igaus),ndime,pnode,ntens,&
               xjaci,d2sdx,elmar(pelty)%deriv(1,1,igaus),elcod_local)     
       end do elemental_gauss_points_geometry

       ncount_newton=1
       if (kfl_linea_nsa == 2) ncount_newton=2 ! to compute the jacobian difference

       do icount_newton=1,ncount_newton
                 
          call nsa_newelmnodevalues(pnode,gath(1,icount_newton),compu)
          
          
          elemental_gauss_points_values: do igaus=1,pgaus        
             
             !
             ! Compute values at gauss (integration) points
             !
             call nsa_newelmgausspointvalues(ielem,igaus,pnode,pgaus,&
                  gath(1,icount_newton),gaus,compu,cartd(1,1,igaus),elmar(pelty)%shape(1:pnode,igaus),hessi(1,1,igaus))
             
             !
             ! Compute local preconditioner at gauss (integration) points
             !          
             call nsa_newelmlocalpreconditioner(igaus,gaus)
             
             !
             ! Compute residuals and sources
             !          
             
             !          if (ielem == 1 .and. igaus==1) then
             !             !       write(6,*) gaus(igaus)%xtide(1:ndofn_nsa)
             !             write(6,*) gaus(igaus)%xunkn(1:ndofn_nsa,ITER_K)
             !             write(6,*) gaus(igaus)%xunkn(1:ndofn_nsa,TIME_N)
             !          end if
             
             
             call nsa_newelmresidualsandsources(igaus,gaus)
             
             
          end do elemental_gauss_points_values

          elemental_gauss_points_stabilization: do igaus=1,pgaus
             
             !
             ! Stabilization by diagonal tau
             !          
             call nsa_newelmvmsdiagonal(ielem,igaus,gaus)
             
             !
             ! Shock capturing: compute shock capturing metrics (cartd is a dummy argument)
             !
             
             call nsa_newelmshocap(igaus,gaus)
             
          end do elemental_gauss_points_stabilization

       end do


!!!! acaaaaa esto esta mal
       if (icount_newton == 1) then
          elrhs= 0.0_rp
          elmat= 0.0_rp
       end if

       elemental_gauss_points_scatter: do igaus=1,pgaus

          do inode=1,pnode

             !
             ! Compute the local rhs and matrix and assemble it to the global rhs and matrix (scatter)
             !

             ! recall that, as defined above, pevat = ndofn_nsa*pnode

             call nsa_newelmatrixcompute(&
                  inode,igaus,ielem,pnode,pevat,&
                  gath(1,1),gaus,cartd(1,1,igaus),elmar(pelty)%shape(1:pnode,igaus),hessi(1,1,igaus),matri)

!!!! acaaaaa ver como restructurar esto


             !
             ! Copy subscale and shocktau to the corresponding global vectors
             !
             call nsa_newelmscattersubscale(inode,igaus,ielem,pnode,gaus,elmar(pelty)%shape(1:pnode,igaus))

          end do

       end do elemental_gauss_points_scatter
       !
       ! Matrix and rhs boundary conditions correction and assembly (for both implicit and explicit cases)
       !
       
       call nsa_newelmatrixsetboundary(ielem,pnode,pevat,lnods(1:pnode,ielem),gath(1,1),elrhs,elmat,kfl_timet_nsa)          


       !
       ! Matrix assembly
       !
       if (kfl_timet_nsa == 2) call assmat(&
            solve(1)%ndofn,pnode,pevat,solve(1)%nunkn,&
            solve(1)%kfl_algso,ielem,lnods(1,ielem),elmat,amatr)
       
       !
       ! RHS assembly
       !
       call assrhs(&
            ndofn_nsa,pnode,lnods(1,ielem),elrhs,rhsid)

    end do elements_loop

    !
    ! Distribute global subscale fields in parallel runs
    !
    call nsa_parall(7_ip) 
    do ipoin = 1,npoin
       umoss_nsa(    1,ipoin,1) = umoss_nsa(    1,ipoin,2)
       umoss_nsa(    2,ipoin,1) = umoss_nsa(    2,ipoin,2)
       if (ndime == 3) umoss_nsa(ndime,ipoin,1) = umoss_nsa(ndime,ipoin,2)
       denss_nsa(      ipoin,1) = denss_nsa(      ipoin,2)
       eness_nsa(      ipoin,1) = eness_nsa(      ipoin,2)
    end do


  end subroutine nsa_newelmoperations

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    12/06/2015
  !> @brief   Start the element loop
  !> @details Start the element loop
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmstart
    use      def_master
    use      def_nastal
    use      def_domain

    implicit none

    integer(ip) :: ipoin,ibopo,idofn

!    dtpseud= 0.0_rp
!    DT_LOCAL_AUXI = DT_PHYSICAL
!    if (kfl_pseud_nsa == 1) then
!       dtpseud= 1.0_rp
!       DT_LOCAL_AUXI = DT_PSEUDO
!    end if

    vdiag_nsa = 1.0_rp
    umoss_nsa(1:ndime,1:npoin,2) = 0.0_rp
    denss_nsa(        1:npoin,2) = 0.0_rp
    eness_nsa(        1:npoin,2) = 0.0_rp
    frequ_nsa(        1:npoin  ) = 0.0_rp
    if (kfl_cotur_nsa <= 0_ip ) turmu(1:npoin) = 0.0_rp

    do ipoin=1,npoin

       ibopo = lpoty(ipoin)
       if (ibopo > 0) then
          !
          ! Initialize rotation and base-change matrices
          !
          do idofn= 1,ndofn_nsa
             jacrot_du_dq_nsa(1:ndofn_nsa,idofn,ibopo)= 0.0_rp
             jacrot_du_dq_nsa(      idofn,idofn,ibopo)= 1.0_rp
             jacrot_dq_du_nsa(1:ndofn_nsa,idofn,ibopo)= 0.0_rp
             jacrot_dq_du_nsa(      idofn,idofn,ibopo)= 1.0_rp
          end do
       end if
       

!!!!! creo que esto no lo necesito mas

!!$       !
!!$       ! PHYSICAl time step:
!!$       ! kfl_dttyp_nsa defines if it is local or not 
!!$       ! it is initialized with dtinv_nsa
!!$       !
!!$       dtinv_eqs(1:5,DT_PHYSICAL)    = dtinv_nsa
!!$       if (kfl_dttyp_nsa(1) > 0 ) then    ! momentum, local time step
!!$          do idime=1,ndime
!!$             itott= (ipoin-1) * ndofn_nsa + idime
!!$             dtinv_eqs(idime,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(1,ipoin,DT_PHYSICAL) 
!!$          end do
!!$       end if
!!$       if (kfl_dttyp_nsa(2) > 0 ) then    ! continuity, local time step
!!$          itott= (ipoin-1) * ndofn_nsa + ndime + 1
!!$          dtinv_eqs(ndime+1,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(2,ipoin,DT_PHYSICAL) 
!!$       end if
!!$       if (kfl_dttyp_nsa(3) > 0 ) then    ! energy, local time step
!!$          itott= (ipoin-1) * ndofn_nsa + ndime + 2
!!$          dtinv_eqs(ndime+2,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(3,ipoin,DT_PHYSICAL) 
!!$       end if
!!$       
!!$       !
!!$       ! PSEUDO time step:
!!$       ! being non-physical, it is always local
!!$       !
!!$       dtinv_eqs(1:ndime,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(1,ipoin,DT_PSEUDO)
!!$       dtinv_eqs(ndime+1,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(2,ipoin,DT_PSEUDO)
!!$       dtinv_eqs(ndime+2,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(3,ipoin,DT_PSEUDO)
!!$       
!!$       ! vdiag is used by explicit schemes
!!$       do idime=1,ndime
!!$          itott= (ipoin-1) * ndofn_nsa + idime           
!!$          !           dtaux= 1.0_rp/dtinv_nsa
!!$          !           if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(1,ipoin,1) 
!!$          dtaux= dtinv_eqs(idime,DT_PHYSICAL) + dtpseud*dtinv_eqs(idime,DT_PSEUDO) 
!!$!          vdiag_nsa(itott) = 1.0_rp / dtaux
!!$          vdiag_nsa(itott) = 1.0_rp 
!!$       end do
!!$       !           dtaux= 1.0_rp/dtinv_nsa
!!$       !        if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(2,ipoin,1) 
!!$       dtaux= dtinv_eqs(ndime+1,DT_PHYSICAL) + dtpseud*dtinv_eqs(ndime+1,DT_PSEUDO) 
!!$!       vdiag_nsa(itott+1) = 1.0_rp / dtaux
!!$       vdiag_nsa(itott+1) = 1.0_rp 
!!$       !           dtaux= 1.0_rp/dtinv_nsa
!!$       !        if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(3,ipoin,1) 
!!$       dtaux= dtinv_eqs(ndime+2,DT_PHYSICAL) + dtpseud*dtinv_eqs(ndime+2,DT_PSEUDO) 
!!$!       vdiag_nsa(itott+2) = 1.0_rp / dtaux
!!$       vdiag_nsa(itott+2) = 1.0_rp 
       
    end do


  end subroutine nsa_newelmstart

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Gather from global to local
  !> @details Gather from global to local
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmgather(ielem,pnode,gath,elcod_local)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip) :: ielem,pnode
    type(elm_onnode_gather_nsa)  :: gath(mnode)

    integer(ip) :: inode,idime,ipoin,dummi,itime_scheme    
    real(rp)    :: &
         prope_tmp(pnode),dummy(ndime,ndime),auxvi(ndofn_nsa),&
         elcod_local(ndime,mnode)


    ! Initialization of dummy variables for viscosity
    dummy = 0.0_rp
    auxvi = 0.0_rp


    !
    ! Properties: viscosity mu, c_p 
    !
    if (kfl_prope /= 0 ) then
       call ker_proper('VISCO','PNODE',dummi,ielem,prope_tmp,pnode)
       gath(1:pnode)%elvis = prope_tmp(1:pnode)
       call ker_proper('SPHEA','PNODE',dummi,ielem,prope_tmp,pnode)
       gath(1:pnode)%elhcp = prope_tmp(1:pnode)
       gath(1:pnode)%elwme  = mowei_nsa
    else
       gath(1:pnode)%elhcp = cpcoe_nsa
       gath(1:pnode)%elwme = mowei_nsa
    endif

    itime_scheme= ITER_K  ! this is required to avoid forbidden memory acces of the global vectors
    if (kfl_tisch_nsa == 2) itime_scheme= TIME_N_MINUS_1

    do inode= 1,pnode

       ipoin= lnods(inode,ielem)

       gath(inode)%elunk(ndime+1,ITER_K) = densi(ipoin,ITER_K)
       gath(inode)%elunk(ndime+2,ITER_K) = energ(ipoin,ITER_K)
       gath(inode)%elunk(ndime+1,TIME_N) = densi(ipoin,TIME_N)
       gath(inode)%elunk(ndime+2,TIME_N) = energ(ipoin,TIME_N)
       gath(inode)%elunk(ndime+1,ITER_AUX) = densi(ipoin,ITER_AUX)
       gath(inode)%elunk(ndime+2,ITER_AUX) = energ(ipoin,ITER_AUX)
       do idime=1,ndime
          gath(inode)%elunk(idime,ITER_K)   =  &
               umome(idime,ipoin,ITER_K)
          gath(inode)%elunk(idime,ITER_AUX) =  &
               umome(idime,ipoin,ITER_AUX)
          gath(inode)%elunk(idime,TIME_N)   =  &
               umome(idime,ipoin,TIME_N)
          gath(inode)%elmsh(idime) = 0.0_rp       ! initialize, because when no mesh motion it must be zero
       end do

       gath(inode)%elsub(ndime+1) = denss_nsa(ipoin,1)
       gath(inode)%elsub(ndime+2) = eness_nsa(ipoin,1)

       !       do idofn = 1,ndofn_nsa
       !          ievat = (inode-1) * ndofn_nsa + idofn
       !          itott = (ipoin-1) * ndofn_nsa + idofn
       !          !!        elsax(ievat) = rhsou_nsa(itott)
       !       end do

       if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then       ! Molecular weight of the mixture
          gath(inode)%elwme  = wmean(ipoin,1)
       endif

       gath(inode)%eldtt(ndime+1,1) = dtieq_nsa(2,ipoin,1)
       gath(inode)%eldtt(ndime+1,2) = dtieq_nsa(2,ipoin,2)
       gath(inode)%eldtt(ndime+2,1) = dtieq_nsa(3,ipoin,1)
       gath(inode)%eldtt(ndime+2,2) = dtieq_nsa(3,ipoin,2)

       gath(inode)%elpre            = press(ipoin,1)
       gath(inode)%eltem            = tempe(ipoin,1)

       do idime= 1,ndime

          gath(inode)%elcod(idime) = coord(idime,ipoin  )          
          elcod_local(idime,inode) = coord(idime,ipoin  )          

          gath(inode)%elsub(idime) = umoss_nsa(idime,ipoin,1)

          gath(inode)%elvel(idime) = gath(inode)%elunk(idime,ITER_K) / gath(inode)%elunk(ndime+1,ITER_K)
          gath(inode)%eldtt(idime,1) = dtieq_nsa(1,ipoin,1)
          gath(inode)%eldtt(idime,2) = dtieq_nsa(1,ipoin,2)
       end do

       !
       ! PHYSICAl time step:
       ! kfl_dttyp_nsa defines if it is local or not 
       ! it is initialized with dtinv_nsa
       !
       gath(inode)%dtinv_eqs(1:5,DT_PHYSICAL)    = dtinv_nsa
       if (kfl_dttyp_nsa(1) > 0 ) then    ! momentum, local time step
          do idime=1,ndime
             gath(inode)%dtinv_eqs(idime,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(1,ipoin,DT_PHYSICAL) 
          end do
       end if
       if (kfl_dttyp_nsa(2) > 0 ) then    ! continuity, local time step
          gath(inode)%dtinv_eqs(ndime+1,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(2,ipoin,DT_PHYSICAL) 
       end if
       if (kfl_dttyp_nsa(3) > 0 ) then    ! energy, local time step
          gath(inode)%dtinv_eqs(ndime+2,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(3,ipoin,DT_PHYSICAL) 
       end if
       
       !
       ! PSEUDO time step:
       ! being non-physical, it is always local
       !
       if (kfl_pseud_nsa == 1) then
          gath(inode)%dtinv_eqs(1:ndime,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(1,ipoin,DT_PSEUDO)
          gath(inode)%dtinv_eqs(ndime+1,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(2,ipoin,DT_PSEUDO)
          gath(inode)%dtinv_eqs(ndime+2,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(3,ipoin,DT_PSEUDO)
       else
          gath(inode)%dtinv_eqs(1:ndime,DT_PSEUDO)  = 0.0_rp
          gath(inode)%dtinv_eqs(ndime+1,DT_PSEUDO)  = 0.0_rp
          gath(inode)%dtinv_eqs(ndime+2,DT_PSEUDO)  = 0.0_rp
       end if       

    end do

    !
    ! Mesh velocity
    !     
    if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          do idime = 1,ndime
             gath(inode)%elmsh(idime) = velom(idime,ipoin)
          end do
       end do
    end if


  end subroutine nsa_newelmgather
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Gather from global to local
  !> @details Gather from global to local
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmgather_newton(ielem,pnode,gath)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip) :: ielem,pnode
    type(elm_onnode_gather_nsa)  :: gath(mnode)

    integer(ip) :: inode,idime,ipoin,itime_scheme,imode,ievat
    real(rp)    :: xhecv,xvelo,xpres,xdens,xtemp,xener,velsq,rdumy,& 
         epsilon_element(nevat_nsa)

    !
    ! Assign initial values to epsilon
    !

    epsilon_element = epsilon_newton

    !
    ! Correct epsilon according to boundary conditions
    !

    ! ESTO DE MOMENTO NO, VAMOS A VER COMO DA SIN CORREGIR...



    itime_scheme= ITER_K  ! this is required to avoid forbidden memory acces of the global vectors
    if (kfl_tisch_nsa == 2) itime_scheme= TIME_N_MINUS_1

    do inode= 1,pnode

       ievat= (inode-1)*ndofn_nsa
       
       ipoin= lnods(inode,ielem)

       xdens = densi(ipoin,ITER_K) + epsilon_element(ievat + ndime + 1)
       gath(inode)%elunk(ndime+1,ITER_K) = xdens
       xener = energ(ipoin,ITER_K) + epsilon_element(ievat + ndime + 2)
       gath(inode)%elunk(ndime+2,ITER_K) = xener
       velsq = 0.0_rp
       do idime=1,ndime
          gath(inode)%elunk(idime,ITER_K)   =  umome(idime,ipoin,ITER_K) + epsilon_element(ievat + idime)
          xvelo = gath(inode)%elunk(idime,ITER_K) / gath(inode)%elunk(ndime+1,ITER_K)
          velsq = velsq + xvelo * xvelo
          gath(inode)%elvel(idime) = xvelo
       end do
       
       ! Recompute derived variables press and tempe
       xhecv = shecp_nsa(ipoin) - runiv_nsa / wmean(ipoin,1)
       
       xtemp= (xener / xdens - 0.5_rp * velsq) / xhecv 

       imode = 0
       if (kfl_prope /= 0 ) imode = 1
       call nsa_stalaw(2,imode,xdens,xpres,xtemp,rdumy,rdumy,wmean(ipoin,1),shecp_nsa(ipoin))

       gath(inode)%elpre            = xpres
       gath(inode)%eltem            = xtemp

    end do

    !
    ! Mesh velocity
    !     
    if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          do idime = 1,ndime
             gath(inode)%elmsh(idime) = velom(idime,ipoin)
          end do
       end do
    end if



  end subroutine nsa_newelmgather_newton

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Compute matrices on nodes 
  !> @details Compute matrices on nodes 
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmnodevalues(pnode,gath,compu)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip) :: pnode
    type(elm_onnode_gather_nsa)  :: gath(mnode)
    type(elm_onnode_comput_nsa)  :: compu(mnode)

    integer(ip) :: inode,idime,kdime
    real(rp)    :: &
         elhcv_aux(pnode),rgacv,velsq,enepe,visci,dicod,dvite,&
         elthe_aux(mnode),eltun_aux(ndofn_nsa,mnode),velno(3)



    do inode=1,pnode

       eltun_aux(ndime+1,inode  ) = gath(inode)%elunk(ndime+1,1)
       eltun_aux(ndime+2,inode  ) = gath(inode)%elunk(ndime+2,1)
       if (kfl_track_nsa == 1) then
          eltun_aux(ndime+1,inode) = gath(inode)%elunk(ndime+1,ITER_K) + gath(inode)%elsub(ndime+1)
          eltun_aux(ndime+2,inode) = gath(inode)%elunk(ndime+2,ITER_K) + gath(inode)%elsub(ndime+2)
       end if

       !
       ! Properties: viscosity mu, c_p 
       !
       if (kfl_prope /= 0 ) then
          elthe_aux(inode) = gath(inode)%elvis * gath(inode)%elhcp / prand_nsa
       else 
          call nsa_lawvis( -1 , 1 ,gath(inode)%elvis,gath(inode)%eltem,dvite) 
          elthe_aux(inode) = gath(inode)%elvis * gath(inode)%elhcp / prand_nsa
       endif

       elhcv_aux(inode) = gath(inode)%elhcp - runiv_nsa / gath(inode)%elwme
       rgacv = runiv_nsa / gath(inode)%elwme / elhcv_aux(inode)
       velsq = 0.0_rp
       do idime=1,ndime
          velsq = velsq + velno(idime)*velno(idime)
          velno(idime) = eltun_aux(idime,inode) / eltun_aux(ndime+1,inode)
          eltun_aux(idime,inode  ) = gath(inode)%elunk(idime,ITER_K) 
          if (kfl_track_nsa == 1) then
             eltun_aux(idime,inode) = gath(inode)%elunk(idime,ITER_K) + gath(inode)%elsub(idime)
          end if
       end do

       enepe = eltun_aux(ndime+2,inode) / eltun_aux(ndime+1,inode) 
       visci = gath(inode)%elvis        / eltun_aux(ndime+1,inode)
       dicod = elthe_aux(inode) / elhcv_aux(inode) / eltun_aux(ndime+1,inode)

       do kdime=1,ndime        
          compu(inode)%elcon(ndime+1,kdime  ,kdime)= 1.0_rp
          compu(inode)%elcon(kdime  ,ndime+2,kdime)= rgacv
          compu(inode)%elcon(kdime  ,ndime+1,kdime)= rgacv * 0.5_rp * velsq
          compu(inode)%elcon(ndime+2,kdime  ,kdime)= ((1.0_rp + rgacv) * enepe - rgacv * 0.5_rp * velsq)
          compu(inode)%elcon(ndime+2,ndime+1,kdime)= &
               - velno(kdime) * ((1.0_rp + rgacv) * enepe - rgacv * velsq)
          compu(inode)%elcon(ndime+2,ndime+2,kdime)= ((1.0_rp + rgacv) * velno(kdime) )
          compu(inode)%eldif(ndime+2,ndime+1,kdime,kdime)= (dicod-visci) * velsq - dicod * enepe
          compu(inode)%eldif(ndime+2,ndime+2,kdime,kdime)= dicod
          do idime=1,ndime
             compu(inode)%elcon(idime,idime  ,kdime)= compu(inode)%elcon(idime,idime  ,kdime) + velno(kdime) 
             compu(inode)%elcon(kdime,idime  ,kdime)= &
                  compu(inode)%elcon(kdime,idime  ,kdime) - rgacv * velno(idime)
             compu(inode)%elcon(idime,kdime  ,kdime)= compu(inode)%elcon(idime,kdime  ,kdime) + velno(idime)
             compu(inode)%elcon(idime,ndime+1,kdime)= &
                  compu(inode)%elcon(idime,ndime+1,kdime) - velno(idime) * velno(kdime)
             compu(inode)%elcon(ndime+2,idime,kdime)= compu(inode)%elcon(ndime+2,idime,kdime) - &
                  rgacv * velno(idime) * velno(kdime)
             compu(inode)%eldif(kdime,kdime,idime,idime  )= visci
             compu(inode)%eldif(kdime,idime,idime,kdime  )= &
                  compu(inode)%eldif(kdime,idime,idime,kdime) + visci
             compu(inode)%eldif(kdime,idime,kdime,idime  )= &
                  compu(inode)%eldif(kdime,idime,kdime,idime) - 2.0_rp * visci / 3.0_rp
             compu(inode)%eldif(kdime,ndime+1,idime,idime)= - visci * velno(kdime)
             compu(inode)%eldif(kdime,ndime+1,idime,kdime)= &
                  compu(inode)%eldif(kdime,ndime+1,idime,kdime) - visci * velno(idime)
             compu(inode)%eldif(kdime,ndime+1,kdime,idime)= &
                  compu(inode)%eldif(kdime,ndime+1,kdime,idime) &
                  + 2.0_rp * visci * velno(idime) / 3.0_rp
             compu(inode)%eldif(ndime+2,idime,kdime,kdime)= (visci-dicod) * velno(idime)
             compu(inode)%eldif(ndime+2,kdime,kdime,idime)= &
                  compu(inode)%eldif(ndime+2,kdime,kdime,idime) &
                  + visci * velno(idime)
             compu(inode)%eldif(ndime+2,kdime,idime,kdime)= &
                  compu(inode)%eldif(ndime+2,kdime,idime,kdime) - &
                  2.0_rp * visci * velno(idime) / 3.0_rp
             compu(inode)%eldif(ndime+2,ndime+1,kdime,idime)= &
                  compu(inode)%eldif(ndime+2,ndime+1,kdime,idime) &
                  + 0.5_rp * visci * velno(kdime) * velno(idime)
          end do
       end do
    end do


  end subroutine nsa_newelmnodevalues
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    09/06/2015
  !> @brief   Interpolate values to the gauss points
  !> @details Interpolate values to the gauss points
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmgausspointvalues(ielem,igaus,pnode,pgaus,gath,gaus,compu,cartigaus,shapigaus,hessigaus)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip), intent(in):: ielem
    integer(ip), intent(in):: igaus
    integer(ip), intent(in):: pnode
    integer(ip), intent(in):: pgaus

    type(elm_onnode_gather_nsa)        :: gath(mnode)
    type(elm_onnode_comput_nsa)        :: compu(mnode)
    type(elm_ongaus_interp_nsa)        :: gaus(mgaus)
    real(rp)                   :: &
         shapigaus(mnode),cartigaus(ndime,mnode),hessigaus(ntens,mnode)

    integer(ip)  ::  idime,idofn,inode,jdofn,jdime,ipoin,ievat

    real(rp)                :: &
         prope_tmp(1),xshai,dwall, &
         gpreo(3),xpreo,dvite,velsq,enepe,dicod,visci,rdumy, &
         hesma(3,3), &
         xmile,seci4

    dwall       = 0.0_rp
    xpreo       = 0.0_rp       ! old pressure and pressure gradient: only 
    rdumy       = 0.0_rp

    hmini_elm = hleng_elm(ndime)  ! hleng(ndime) is the smallest
    hmaxi_elm = hleng_elm(1)      ! hleng(1) is the largest 

    !
    ! gaus(...) is set to zero outside
    !
    
    if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
       gaus(igaus)%htrad(1:ndime) = enthalpy_transport(ielem)%a(1:ndime,igaus,1)
       gaus(igaus)%dhtra          = div_enthalpy_transport(ielem)%a(igaus,1,1)
       gaus(igaus)%heats          = chemical_heat(ielem)%a(igaus,1,1)
    endif
    

    do idime=1,ndime
       gpreo(idime) = 0.0_rp               !     used as limiter when current pressure goes below zero
       do idofn=1,ndofn_nsa
          do jdime=1,ndime
             hunkn_elm(idofn,idime,jdime) = 0.0_rp
          end do
       end do
    end do

    do inode= 1,pnode
       xshai= shapigaus(inode)       
       xpreo= xpreo + xshai * gath(inode)%elpre
       do idime=1,ndime
          do jdime=1,ndime
             hesma(idime,jdime) = hessigaus(nindx_nsa(idime,jdime),inode)
          end do
       end do

       do idofn= 1,ndofn_nsa
          !!        xsube(idofn,igaus,2) = xsube(idofn,igaus,2) + xshai * elsub(idofn,inode)
          gaus(igaus)%xunkn(idofn,ITER_K) = gaus(igaus)%xunkn(idofn,ITER_K) &
               + xshai*gath(inode)%elunk(idofn,ITER_K)
          gaus(igaus)%xunkn(idofn,ITER_AUX) = gaus(igaus)%xunkn(idofn,ITER_AUX) &
               + xshai*gath(inode)%elunk(idofn,ITER_AUX)
          gaus(igaus)%xunkn(idofn,TIME_N) = gaus(igaus)%xunkn(idofn,TIME_N) &
               + xshai*gath(inode)%elunk(idofn,TIME_N)

          gaus(igaus)%xdtix(idofn,1) = gaus(igaus)%xdtix(idofn,1) + xshai*gath(inode)%eldtt(idofn,1)     
          gaus(igaus)%xdtix(idofn,2) = gaus(igaus)%xdtix(idofn,2) + xshai*gath(inode)%eldtt(idofn,2)

          ievat = (inode-1) * ndofn_nsa + idofn

          do jdofn=1,ndofn_nsa
             do idime=1,ndime
                gaus(igaus)%dconv(idofn,jdofn) = gaus(igaus)%dconv(idofn,jdofn) &
                     + cartigaus(idime,inode) * compu(inode)%elcon(idofn,jdofn,idime)
                do jdime=1,ndime
                   gaus(igaus)%ddiff(idofn,jdofn,idime,1) = gaus(igaus)%ddiff(idofn,jdofn,idime,1) &
                        + cartigaus(jdime,inode) * compu(inode)%eldif(idofn,jdofn,jdime,idime)
                   gaus(igaus)%ddiff(idofn,jdofn,idime,2) = gaus(igaus)%ddiff(idofn,jdofn,idime,2) &
                        + cartigaus(jdime,inode) * compu(inode)%eldif(idofn,jdofn,idime,jdime)
                end do
             end do
          end do
       end do


       do idime=1,ndime
          gpreo(idime) = gpreo(idime) + gath(inode)%elpre * cartigaus(idime,inode)
          do idofn = 1,ndofn_nsa
             gaus(igaus)%gunkn(idofn,idime) = gaus(igaus)%gunkn(idofn,idime) &
                  + cartigaus(idime,inode)*gath(inode)%elunk(idofn,ITER_K)
             gaus(igaus)%gsube(idofn,idime) = gaus(igaus)%gsube(idofn,idime) &
                  + cartigaus(idime,inode)*gath(inode)%elsub(idofn)
             do jdime=1,ndime
                hunkn_elm(idofn,idime,jdime) = hunkn_elm(idofn,idime,jdime) &
                     + hesma(idime,jdime) * gath(inode)%elunk(idofn,ITER_K)
             end do
          end do
          ! only the first ndime values of hessi_nsa are used now, i.e. the diagonal
          gaus(igaus)%xlade        = gaus(igaus)%xlade        &
               + hessigaus(idime,inode) * gath(inode)%elunk(ndime+1,ITER_K) 
          gaus(igaus)%xldve(idime) = gaus(igaus)%xldve(idime) &
               + hessigaus(idime,inode) * gath(inode)%elunk(ndime+1,ITER_K)
          ! xvmsh is only different than zero when coupled to alefor
          gaus(igaus)%xvmsh(idime) = gaus(igaus)%xvmsh(idime) + xshai*gath(inode)%elmsh(idime)  
       end do
       !
       ! Turbulent viscosity at gauss point when coupled with TURBUL
       !
       if (kfl_cotur_nsa /= 0 ) then 
          ipoin = lnods(inode,ielem)
          dwall = dwall + xshai * walld(ipoin)    ! Interpolation of the wall distance at gauss point
          if (kfl_cotur_nsa == 1) then
             gaus(igaus)%xnutu = gaus(igaus)%xnutu + xshai * turmu(ipoin)
          endif
       endif
       !
       ! Molecular weight at gauss point when coupled with CHEMIC
       !
       if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
          ipoin = lnods(inode,ielem)
          gaus(igaus)%xmowe = gaus(igaus)%xmowe + xshai * wmean(ipoin,1)
       endif

    end do
            
    !
    ! Mesh velocity
    !     
    !    do inode = 1,pnode
    !       xshai= xshap_nsa(inode,igaus)
    !       ipoin = lnods(inode,ielem)
    !       do idime = 1,ndime
    !       end do
    !    end do

    do idofn=1,ndofn_nsa
       gaus(igaus)%xtunk(idofn) = gaus(igaus)%xunkn(idofn,ITER_K)
       do idime=1,ndime           
          gaus(igaus)%gtunk(idofn,idime) = gaus(igaus)%gunkn(idofn,idime)
       end do
    end do

    if (kfl_track_nsa == 1) then
       do idofn=1,ndofn_nsa
          gaus(igaus)%xtunk(idofn) = gaus(igaus)%xunkn(idofn,ITER_K) + gaus(igaus)%xsube(idofn,2)
          do idime=1,ndime           
             gaus(igaus)%gtunk(idofn,idime) = &
                  gaus(igaus)%gunkn(idofn,idime) + gaus(igaus)%gsube(idofn,idime)
          end do
       end do
    end if


    !
    ! If local preconditioners are used, xtide is overwritten later
    !
    
!    if (kfl_pseud_nsa == 1) then
       do idofn=1,ndofn_nsa
          gaus(igaus)%xtide(idofn) = gaus(igaus)%xunkn(idofn,ITER_K)-gaus(igaus)%xunkn(idofn,TIME_N)
          gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
       end do
!    end if

    ! no shock capturing to x-momentum means no shock capturing at all
    if (kfl_shock_nsa(1) > 0) then
       gaus(igaus)%shocktau_local(1:ndime)=shocktau_nsa(ielem)%a(1,igaus,1) 
       gaus(igaus)%shocktau_local(ndime+1)=shocktau_nsa(ielem)%a(2,igaus,1) 
       gaus(igaus)%shocktau_local(ndime+2)=shocktau_nsa(ielem)%a(3,igaus,1) 
    end if

    
    !
    ! Derived values
    !
   
    velsq= 0.0_rp
    gaus(igaus)%xpres= 0.0_rp
    gaus(igaus)%velmo= 0.0_rp
    do idime=1,ndime
       gaus(igaus)%xvelo(idime) =  gaus(igaus)%xtunk(idime) / gaus(igaus)%xtunk(ndime+1)
       velsq = velsq + gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(idime) 
       gaus(igaus)%xpres = gaus(igaus)%xpres + gaus(igaus)%xtunk(idime) * gaus(igaus)%xtunk(idime)
       do jdime= 1,ndime
          gaus(igaus)%gvelo(idime,jdime) = (gaus(igaus)%gtunk(idime,jdime) &
               - gaus(igaus)%xtunk(idime) * gaus(igaus)%gtunk(ndime+1,jdime) / gaus(igaus)%xtunk(ndime+1)) &
               / gaus(igaus)%xtunk(ndime+1)
       end do
    end do
    gaus(igaus)%velmo = sqrt(velsq)

    if (kfl_prope /= 0 ) then
       call ker_proper('SPHEA','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,shapigaus,cartigaus) 
       gaus(igaus)%xheatcp = prope_tmp(1)
       if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
          gaus(igaus)%xmowe = mowei_nsa
       endif
       gaus(igaus)%xrgasc = runiv_nsa / gaus(igaus)%xmowe  
       gaus(igaus)%xheatcv = gaus(igaus)%xheatcp - gaus(igaus)%xrgasc  !Cv is computed from R & Cp
    else
       gaus(igaus)%xheatcp = cpcoe_nsa 
       gaus(igaus)%xmowe = mowei_nsa
       gaus(igaus)%xrgasc = runiv_nsa / gaus(igaus)%xmowe
       gaus(igaus)%xheatcv = gaus(igaus)%xheatcp - gaus(igaus)%xrgasc  !Cv is computed from R & Cp
    endif

    do idime=1,ndime
       do jdime=1,ndime
          gaus(igaus)%gpres(idime) = gaus(igaus)%gpres(idime) &
               + gaus(igaus)%gvelo(jdime,idime) * gaus(igaus)%xtunk(jdime) &
               + gaus(igaus)%xvelo(jdime) * gaus(igaus)%gtunk(jdime,idime)
       end do
       gaus(igaus)%gpres(idime) = gaus(igaus)%xrgasc * (gaus(igaus)%gtunk(ndime+2,idime) &
            - 0.5_rp * gaus(igaus)%gpres(idime)) / gaus(igaus)%xheatcv 
    end do

    gaus(igaus)%xpres  = gaus(igaus)%xrgasc * (gaus(igaus)%xtunk(ndime+2) - 0.5_rp * gaus(igaus)%xpres &
         / gaus(igaus)%xtunk(ndime+1))  /   gaus(igaus)%xheatcv

    if (gaus(igaus)%xpres .lt. zensa) then
       gaus(igaus)%xpres = xpreo
       do idime=1,ndime
          gaus(igaus)%gpres(idime)= gpreo(idime)
       end do
    end if

    gaus(igaus)%xtemp = gaus(igaus)%xpres / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xrgasc

    gaus(igaus)%xvisc = 0.0_rp
    gaus(igaus)%xdith = 0.0_rp
    dvite = 0.0_rp
    gaus(igaus)%sgsdi = 0.0_rp

    if (kfl_visco_nsa > 0) then 
       !
       ! SGS viscous dissipation for LES
       !
       if (kfl_cotur_nsa < 0) then     
          seci4 = 0.0_rp
          do idime = 1,ndime                     ! 2 S_ij : S_ij
             do jdime = 1,ndime         
                seci4 = seci4 + &
                     gaus(igaus)%gvelo(idime,jdime) &
                     * (gaus(igaus)%gvelo(idime,jdime) + gaus(igaus)%gvelo(jdime,idime))
             end do
          end do
          xmile =  gaus(igaus)%dvolu**0.3333333_rp
          !
          ! SGS_DISSIPATION = C_eps * (K^sgs)**3/2 / V**1/3, K^sgs = sqrt(3/4) * nut * |S|
          !
          gaus(igaus)%sgsdi = 0.916_rp * 0.866025_rp * (sqrt(seci4)**1.5_rp) / xmile 
       endif
       !
       ! Properties from the kernel: viscosity mu, c_p, K
       !
       if (kfl_prope /= 0 ) then
          call ker_proper('VISCO','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,shapigaus,cartigaus) 
          gaus(igaus)%xvisc = prope_tmp(1)
          call ker_proper('CONDU','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,shapigaus,cartigaus) 
          gaus(igaus)%xdith = prope_tmp(1)
          !
          ! Computation viscosity derivative dmu/dT
          !
          if (lawvi_nsa == 1) then               ! power law
             dvite = vispa_nsa(1)*vispa_nsa(2)*gaus(igaus)%xtemp**(vispa_nsa(2)-1.0_rp)        
          else if (lawvi_nsa==2) then            ! sutherland law
             dvite = gaus(igaus)%xvisc &
                  *( 1.5_rp/gaus(igaus)%xtemp - 1.0_rp/(gaus(igaus)%xtemp+vispa_nsa(2)))
          end if

          if (kfl_cotur_nsa /= 0) then           ! If turbulence model ON 
             if (kfl_cotur_nsa < 0) then        ! Turbulent viscosity 
                call nsa_turbul(&
                     dwall,gaus(igaus)%velmo,gaus(igaus)%xvisc,gaus(igaus)%gvelo(1,1),&
                     gaus(igaus)%dvolu,gaus(igaus)%xnutu,&
                     gaus(igaus)%xunkn(1,ITER_K))
                gaus(igaus)%sgsdi = gaus(igaus)%sgsdi*gaus(igaus)%xnutu
                gaus(igaus)%xnutu = gaus(igaus)%xnutu*gaus(igaus)%xtunk(ndime+1)
             endif
             ! If turbulence model ON, xnutu /= 0, otherwise = 0 
             gaus(igaus)%xvisc = gaus(igaus)%xvisc + gaus(igaus)%xnutu  
             gaus(igaus)%xdith = gaus(igaus)%xdith + gaus(igaus)%xnutu  * gaus(igaus)%xheatcp / prand_nsa
          endif

       else
          call nsa_lawvis(-1,1,gaus(igaus)%xvisc,gaus(igaus)%xtemp,dvite)! Dynamic viscosity 

          if (kfl_cotur_nsa < 0) then            ! LES Turbulent viscosity 
             call nsa_turbul(&
                  dwall,gaus(igaus)%velmo,gaus(igaus)%xvisc,gaus(igaus)%gvelo(1,1),&
                  gaus(igaus)%dvolu,gaus(igaus)%xnutu,&
                  gaus(igaus)%xunkn(1,ITER_K))
             gaus(igaus)%sgsdi = gaus(igaus)%sgsdi*gaus(igaus)%xnutu
             gaus(igaus)%xnutu = gaus(igaus)%xnutu*gaus(igaus)%xtunk(ndime+1)
          endif
          gaus(igaus)%xvisc = gaus(igaus)%xvisc + gaus(igaus)%xnutu                  ! If turbulence model ON, xnutu /= 0, otherwise = 0 
          gaus(igaus)%xdith = gaus(igaus)%xvisc * gaus(igaus)%xheatcp / prand_nsa
       endif
    endif

    gaus(igaus)%dvelo = 0.0_rp

    do idime=1,ndime
       gaus(igaus)%gtemp(idime) = &
            ( gaus(igaus)%gpres(idime) -  gaus(igaus)%xpres * gaus(igaus)%gtunk(ndime+1,idime) &
            / gaus(igaus)%xtunk(ndime+1)) &
            / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xrgasc
       gaus(igaus)%gvisc(idime) = dvite *  gaus(igaus)%gtemp(idime)
       gaus(igaus)%dvelo =  gaus(igaus)%dvelo +  gaus(igaus)%gvelo(idime,idime)
    end do
    enepe = gaus(igaus)%xtunk(ndime+2) / gaus(igaus)%xtunk(ndime+1)
    visci = gaus(igaus)%xvisc / gaus(igaus)%xtunk(ndime+1)
    dicod = gaus(igaus)%xdith / gaus(igaus)%xheatcv / gaus(igaus)%xtunk(ndime+1)

    gaus(igaus)%xrgacv = gaus(igaus)%xrgasc / gaus(igaus)%xheatcv
    gaus(igaus)%xadgam = gaus(igaus)%xheatcp / gaus(igaus)%xheatcv

    gaus(igaus)%xsoun= sqrt(gaus(igaus)%xadgam * gaus(igaus)%xrgasc *  gaus(igaus)%xtemp)

    ! xvmsh is the mesh velocity, only different than zero when coupled to alefor


    !
    ! Compute jacobian matrices 
    !

    do jdime=1,ndime        
       gaus(igaus)%xconv(ndime+1,jdime  ,jdime)= 1.0_rp
       gaus(igaus)%xconv(jdime  ,ndime+2,jdime)= gaus(igaus)%xrgacv
       gaus(igaus)%xconv(jdime  ,ndime+1,jdime)= gaus(igaus)%xrgacv * 0.5_rp * velsq
       gaus(igaus)%xconv(ndime+2,jdime  ,jdime)= &
            ((1.0_rp + gaus(igaus)%xrgacv) * enepe - gaus(igaus)%xrgacv * 0.5_rp * velsq)
       gaus(igaus)%xconv(ndime+2,ndime+1,jdime)= &
            - gaus(igaus)%xvelo(jdime) &
            * ((1.0_rp + gaus(igaus)%xrgacv) * enepe - gaus(igaus)%xrgacv * velsq) + gaus(igaus)%htrad(jdime)
       gaus(igaus)%xconv(ndime+2,ndime+2,jdime)= &
            ((1.0_rp + gaus(igaus)%xrgacv) * ( gaus(igaus)%xvelo(jdime) - gaus(igaus)%xvmsh(jdime) ))
       gaus(igaus)%xdiff(ndime+2,ndime+1,jdime,jdime)= (dicod-visci) * velsq - dicod * enepe
       gaus(igaus)%xdiff(ndime+2,ndime+2,jdime,jdime)= dicod
       do idime=1,ndime
          gaus(igaus)%xconv(idime,idime  ,jdime)=  gaus(igaus)%xconv(idime,idime  ,jdime) &
               + ( gaus(igaus)%xvelo(jdime) - gaus(igaus)%xvmsh(jdime) ) 
          gaus(igaus)%xconv(jdime,idime  ,jdime)=  gaus(igaus)%xconv(jdime,idime  ,jdime) &
               - gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idime)
          gaus(igaus)%xconv(idime,jdime  ,jdime)=  gaus(igaus)%xconv(idime,jdime  ,jdime) &
               + gaus(igaus)%xvelo(idime)
          gaus(igaus)%xconv(idime,ndime+1,jdime)=  gaus(igaus)%xconv(idime,ndime+1,jdime) &
               - gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(jdime)
          gaus(igaus)%xconv(ndime+2,idime,jdime)=  gaus(igaus)%xconv(ndime+2,idime,jdime) &
               - gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(jdime)
          gaus(igaus)%xdiff(jdime,jdime,idime,idime)= visci
          gaus(igaus)%xdiff(jdime,idime,idime,jdime  )=  gaus(igaus)%xdiff(jdime,idime,idime,jdime) + visci
          gaus(igaus)%xdiff(jdime,idime,jdime,idime  )=  gaus(igaus)%xdiff(jdime,idime,jdime,idime) &
               - 2.0_rp * visci / 3.0_rp
          gaus(igaus)%xdiff(jdime,ndime+1,idime,idime)= - visci * gaus(igaus)%xvelo(jdime)
          gaus(igaus)%xdiff(jdime,ndime+1,idime,jdime)=  gaus(igaus)%xdiff(jdime,ndime+1,idime,jdime) &
               - visci * gaus(igaus)%xvelo(idime)
          gaus(igaus)%xdiff(jdime,ndime+1,jdime,idime)=  gaus(igaus)%xdiff(jdime,ndime+1,jdime,idime) &
               + 2.0_rp * visci * gaus(igaus)%xvelo(idime) / 3.0_rp
          gaus(igaus)%xdiff(ndime+2,idime,jdime,jdime)= (visci-dicod) * gaus(igaus)%xvelo(idime)
          gaus(igaus)%xdiff(ndime+2,jdime,jdime,idime)=  gaus(igaus)%xdiff(ndime+2,jdime,jdime,idime) &
               + visci * gaus(igaus)%xvelo(idime)
          gaus(igaus)%xdiff(ndime+2,jdime,idime,jdime)=  gaus(igaus)%xdiff(ndime+2,jdime,idime,jdime) &
               - 2.0_rp * visci * gaus(igaus)%xvelo(idime) / 3.0_rp
          gaus(igaus)%xdiff(ndime+2,ndime+1,jdime,idime)=  gaus(igaus)%xdiff(ndime+2,ndime+1,jdime,idime) &
               + 0.5_rp * visci * gaus(igaus)%xvelo(jdime) * gaus(igaus)%xvelo(idime)
       end do
    end do



  end subroutine nsa_newelmgausspointvalues


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    09/06/2015
  !> @brief   Interpolate values to the gauss points
  !> @details Interpolate values to the gauss points
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmlocalpreconditioner(igaus,gaus)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip), intent(in):: igaus

    type(elm_ongaus_interp_nsa)        :: gaus(mgaus)

    integer(ip)  ::  idime,kdime,idofn,jdofn,jdime,kdofn,ldofn,mdofn,odofn,pdofn

    real(rp)     :: &
         xconv_ori(5,5,3), & ! Original convective Jacobian without preconditioning
         xconv_lopre_diff(5,5,3), & ! dP * K for CM local preconditioner
         xdiff_ori(5,5,3,3), & ! Original diffusion matrix without preconditioning
         xlopr_symmet(5,5),&
         xttra(5,5), & ! Symmetrizing transformation of the N-S equations
         xttra_inv(5,5), & ! Inverse of the Symmetrizing transformation of the N-S equations
         xqtra(5,5), & ! Streamline transformation of the N-S equations
         xqtra_inv(5,5), & ! Inverse of the Streamline transformation of the N-S equations
         xftra(5,5), & ! Symmetrizing streamline transformation of the N-S equations
         xftra_inv(5,5), & ! Inverse of the Symmetrizing streamline transformation of the N-S equations
         velmo_xy, & ! sqrt(u_1^2 + u_2^2)
         xmach, & ! Mach number at Gauss points
         xmach_ref, & ! Reference Mach number at Guass points
         aturkel, & ! Turkel's factor on reference mach for preconditioners 
         xbeta, &  ! local preconditioning parameter
         xbeta_visc, &  ! local preconditioning parameter
         xtaup, & ! local preconditioning parameter
         cfl_sound, & ! local preconditioning parameter
         reyno_cell, reyno_cell_inv, & ! cell Reynolds number and its inverse
         gmach(3), gsoun(3), gmach_ref(3), &
         gvemo(3), gbeta(3), gbeta_visc(3), greyn_cell_inv(3), &
         glopr(5,5,3), & ! gradient of CM preconditioner
         termA,termB,termC,termD,termA_der,termB_der,termC_der,termD_der, &
         termE,termF,termG,termH,termI,term1,term2,term3,term4,term5,zensa_lopre,velsq



    if (kfl_lopre_nsa == 0) return 

    velsq= gaus(igaus)%velmo * gaus(igaus)%velmo


    zensa_lopre = 0.1_rp
    xmach = 1.0e-10           ! very small mach number when u=0 to avoid ill definition of xmach
    if (gaus(igaus)%velmo > 1.0e-8_rp) then
       xmach = gaus(igaus)%velmo / gaus(igaus)%xsoun
    end if


    do idofn=1,ndofn_nsa
       do jdofn=1,ndofn_nsa
          xlopr_symmet(idofn,jdofn) = 0.0_rp
          xttra(idofn,jdofn) = 0.0_rp
          xttra_inv(idofn,jdofn) = 0.0_rp
          xqtra(idofn,jdofn) = 0.0_rp
          xqtra_inv(idofn,jdofn) = 0.0_rp

          xftra(idofn,jdofn) = 0.0_rp
          xftra_inv(idofn,jdofn) = 0.0_rp
          gaus(igaus)%xlopr_conservative(idofn,jdofn) = 0.0_rp
          do idime=1,ndime
             xconv_ori(idofn,jdofn,idime) = 0.0_rp           
             xconv_lopre_diff(idofn,jdofn,idime) = 0.0_rp
             glopr(idofn,jdofn,idime) = 0.0_rp
             xconv_ori(idofn,jdofn,idime) = gaus(igaus)%xconv(idofn,jdofn,idime)
             gaus(igaus)%xconv(idofn,jdofn,idime) = 0.0_rp
             do jdime=1,ndime
                xdiff_ori(idofn,jdofn,idime,jdime) = gaus(igaus)%xdiff(idofn,jdofn,idime,jdime)
                gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) = 0.0_rp
             end do
          end do
       end do
    end do

    if (kfl_lopre_nsa == 1) then ! IDENTITY LOCAL PRECONDITIONER (JUST FOR TESTING PURPOUSES)

       do idofn=1,ndofn_nsa
          gaus(igaus)%xlopr_conservative(idofn,idofn) = 1.0_rp           
       end do

       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
!             if (kfl_pseud_nsa == 1) &
                  gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
                  + gaus(igaus)%xlopr_conservative(idofn,jdofn) &
                  * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
          end do
!          if (kfl_pseud_nsa == 1) &
               gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
       end do

       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
             do idime=1,ndime
                gaus(igaus)%xconv(idofn,jdofn,idime) = xconv_ori(idofn,jdofn,idime)
                do jdime=1,ndime
                   gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) = xdiff_ori(idofn,jdofn,idime,jdime)
                end do
             end do
          end do
       end do

    else if (kfl_lopre_nsa == 2) then  ! VAN LEER-LEE-ROE STEADY INVISCID LOCAL PRECONDITIONER IS APPLIED

        if (velsq == 0.0_rp) then

           ! when velsq is strictly zero and VLR is used, do not precondition

           gaus(igaus)%xlopr_conservative= 1.0_rp           
           
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 if (kfl_pseud_nsa == 1) &
                      gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
                      + gaus(igaus)%xlopr_conservative(idofn,jdofn) &
                      * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
              end do
              if (kfl_pseud_nsa == 1) gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
           end do
           
           gaus(igaus)%xconv = xconv_ori
           gaus(igaus)%xdiff = xdiff_ori
           
        else
           
           
           if (xmach < 1.0_rp-zensa_lopre) then
              xbeta = sqrt(1.0_rp-xmach*xmach)
              xtaup = xbeta
           else if (xmach < 1.0_rp) then
              xbeta = sqrt((2.0_rp-zensa_lopre)*zensa_lopre)
              xtaup = xbeta
           else if (xmach < 1.0_rp+zensa_lopre) then
              xbeta = sqrt((2.0_rp+zensa_lopre)*zensa_lopre)
              xtaup = xbeta / xmach
           else
              xbeta = sqrt(xmach*xmach-1.0_rp)
              xtaup = xbeta / xmach
           end if
           
!!! BEGIN: VLR preconditioner with conservative variables (any change of variables is applied)
           termA = xtaup / xbeta / xbeta
           termB = 1.0_rp - xmach * xmach
           termC = 0.5_rp * gaus(igaus)%xrgacv * xmach * xmach
           termD = 1.0_rp / gaus(igaus)%xsoun / gaus(igaus)%xsoun
           termE = termA * termD
           termF = gaus(igaus)%xrgacv * termD
           termG = termA * xmach * xmach
           termH = gaus(igaus)%xrgacv * termD * (1.0_rp - termG)
           termI = termC * (termG - 1.0_rp)
           
           term1 = (1.0_rp+termA-xtaup) &
                / velsq + (gaus(igaus)%xrgacv*termA*termB-termA+gaus(igaus)%xrgacv) * termD
           term2 = - termA * termB * (1.0_rp+termC) - termC
           term3 = - gaus(igaus)%xrgacv * (termA*termB+1.0_rp) * termD
           term4 = - termE + termH
           term5 = 1.0_rp + termA - termA/gaus(igaus)%xrgacv &
                - 3.0_rp*termG/2.0_rp + gaus(igaus)%xrgacv*termG + termC*(1.0_rp-termG)
           
           do idofn=1,ndime
              gaus(igaus)%xlopr_conservative(idofn,ndime+1) = term2 * gaus(igaus)%xvelo(idofn)
              gaus(igaus)%xlopr_conservative(idofn,ndime+2) = term3 * gaus(igaus)%xvelo(idofn)
              gaus(igaus)%xlopr_conservative(ndime+1,idofn) = term4 * gaus(igaus)%xvelo(idofn)
              gaus(igaus)%xlopr_conservative(ndime+2,idofn) = term5 * gaus(igaus)%xvelo(idofn)
              do jdofn=idofn,ndime
                 gaus(igaus)%xlopr_conservative(idofn,jdofn) = &
                      term1 * gaus(igaus)%xvelo(idofn) * gaus(igaus)%xvelo(jdofn)
                 gaus(igaus)%xlopr_conservative(jdofn,idofn) = gaus(igaus)%xlopr_conservative(idofn,jdofn)
              end do
              gaus(igaus)%xlopr_conservative(idofn,idofn) = gaus(igaus)%xlopr_conservative(idofn,idofn) + xtaup 
           end do
           gaus(igaus)%xlopr_conservative(ndime+1,ndime+1) = 1.0_rp + termG + termC * (termG-1.0_rp)
           gaus(igaus)%xlopr_conservative(ndime+1,ndime+2) = - termH
           gaus(igaus)%xlopr_conservative(ndime+2,ndime+1) = &
                (1.0_rp/gaus(igaus)%xrgacv-termB) * termA * velsq - 0.5_rp * velsq &
                + 0.5_rp * termI * velsq - termC * velsq * termA
           gaus(igaus)%xlopr_conservative(ndime+2,ndime+2) = (1.0_rp-gaus(igaus)%xrgacv) * termG + termI
!!! END: VLR preconditioner with conservative variables (any change of variables is applied)
           
           ! Compute P_conservative * A
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 !             if (kfl_pseud_nsa == 1) &
                 gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
                      + gaus(igaus)%xlopr_conservative(idofn,jdofn) &
                      * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
                 do idime=1,ndime
                    do pdofn=1,ndofn_nsa
                       gaus(igaus)%xconv(idofn,jdofn,idime) = &
                            gaus(igaus)%xconv(idofn,jdofn,idime) &
                            + gaus(igaus)%xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                    end do
                 end do
              end do
              !          if (kfl_pseud_nsa == 1) &
              gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
           end do
           
        end if

    else if (kfl_lopre_nsa == 3) then  ! CHOI & MERKLE STEADY/UNSTEADY VISCOUS/INVISCID PRECONDITIONER IS APPLIED

       do idime=1,ndime
          gvemo(idime) = 0.0_rp
          greyn_cell_inv(idime) = 0.0_rp
          if (gaus(igaus)%velmo > 0.0_rp) then
             gvemo(idime) = gvemo(idime) / gaus(igaus)%velmo
             do jdime=1,ndime
                gvemo(idime) = gvemo(idime) + gaus(igaus)%xvelo(jdime) * gaus(igaus)%gvelo(jdime,idime)
             end do
             !              greyn_cell_inv(idime) = - xvisc * (gaus(igaus)%gtunk(ndime+1,idime)*gaus(igaus)%velmo+gaus(igaus)%xtunk(ndime+1)*gvemo(idime)) &
             !                   / hmaxi_elm / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xtunk(ndime+1) / velsq
             greyn_cell_inv(idime) = - gaus(igaus)%xvisc &
                  * (gaus(igaus)%gtunk(ndime+1,idime)*gaus(igaus)%velmo+gaus(igaus)%xtunk(ndime+1)&
                  *gvemo(idime)) &
                  / hmini_elm / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xtunk(ndime+1) / velsq
          end if

          gsoun(idime) = 0.5_rp * gaus(igaus)%xadgam * gaus(igaus)%xrgasc * gaus(igaus)%gtemp(idime) &
               / gaus(igaus)%xsoun
          gmach(idime) = (gvemo(idime) - xmach * gsoun(idime)) / gaus(igaus)%xsoun 
          gmach_ref(idime) = gmach(idime)

       end do

       !
       !  Compute epsilon (xmach_ref):
       !  
       !  epsilon = min(1 , max (M_limit, M_computed, a*M_inflow))
       !  
       !  the third case is suggested by Turkel for more robust behaviour deep down boundary layers
       !  see paper from colin 2001. a=1 is fine.
       !

       xmach_ref = xmach  ! C&M STEADY CASE
       if (kfl_pseud_nsa == 1) then  ! C&M UNSTEADY CASE        
          cfl_sound = gaus(igaus)%xsoun / dtinv / hmini_elm
          xmach_ref = sqrt(xmach * xmach + 1.0_rp / cfl_sound / cfl_sound)
       end if

       !
       ! M_inflow is rmach_nsa (Turkel's third term )
       !
       aturkel= 2.0_rp
       if (xmach_ref < aturkel * rmach_nsa) xmach_ref= rmach_nsa


       if (xmach_ref < 0.00001_rp) then 
          xmach_ref = 0.00001_rp
          gmach_ref = 0.0_rp
       else if (xmach_ref > 1.0_rp) then
          xmach_ref = 1.0_rp
          gmach_ref = 0.0_rp
       end if

       if (kfl_visco_nsa == 0) then  ! C&M INVISCID CASE
          xbeta_visc = 1.0_rp
          gbeta_visc = 0.0_rp
       else  ! C&M VISCOUS CASE
          !           reyno_cell     = gaus(igaus)%xtunk(ndime+1) * gaus(igaus)%velmo * hmaxi_elm / gaus(igaus)%xvisc
          reyno_cell     = gaus(igaus)%xtunk(ndime+1) * gaus(igaus)%velmo * hmini_elm / gaus(igaus)%xvisc
          reyno_cell_inv = 1.0e10   ! very large 1/Re when u=0
          if (reyno_cell > 0.0_rp) reyno_cell_inv = 1.0_rp / reyno_cell
          termA = xmach_ref * xmach_ref * (reyno_cell_inv - 1.0_rp + 1.0_rp / xmach / xmach)
          xbeta_visc = reyno_cell_inv * (reyno_cell_inv - 1.0_rp) / termA
          do idime=1,ndime
             termB = 2.0_rp * xmach_ref * gmach_ref(idime) &
                  * (reyno_cell_inv - 1.0_rp + 1.0_rp / xmach / xmach) &
                  + xmach_ref * xmach_ref * (greyn_cell_inv(idime) &
                  - 2.0_rp * gmach(idime) / xmach / xmach / xmach)
             gbeta_visc(idime) = (greyn_cell_inv(idime) * (2.0_rp * reyno_cell_inv - 1.0_rp) &
                  - reyno_cell_inv * (reyno_cell_inv - 1.0_rp) * termB / termA) / termA
          end do
          if (xbeta_visc < 1.0_rp) then
             xbeta_visc = 1.0_rp
             gbeta_visc = 0.0_rp
          end if
       end if

       xbeta = xbeta_visc * gaus(igaus)%xsoun * gaus(igaus)%xsoun
       xtaup = 1.0_rp
       do idime=1,ndime
          gbeta(idime) = gbeta_visc(idime) * gaus(igaus)%xsoun * gaus(igaus)%xsoun &
               + 2.0_rp * gaus(igaus)%xsoun * xbeta_visc * gsoun(idime)
       end do

!!! BEGIN: CHOI & MERKLE preconditioner with primitive variables (p, u T)
       termA = xbeta * xmach_ref * xmach_ref
       termB = 0.5_rp * velsq - gaus(igaus)%xheatcp * gaus(igaus)%xtemp + termA * xtaup

       do idime=1,ndime
          xlopr_symmet(idime,idime) = 1.0_rp / gaus(igaus)%xtunk(ndime+1)
          xlopr_symmet(idime,ndime+1) = - gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1)
          xlopr_symmet(ndime+2,idime) = &
               - gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xheatcp
       end do
       xlopr_symmet(ndime+1,ndime+1) = termA
       xlopr_symmet(ndime+2,ndime+1) = termB / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xheatcp
       xlopr_symmet(ndime+2,ndime+2) = 1.0_rp / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xheatcp
!!! END: CHOI & MERKLE preconditioner with primitive variables (p, u T)

!!! BEGIN: CHOI & MERKLE preconditioner with conservative variables (we apply a change of variables)
       do idime=1,ndime
          xftra(idime,idime) = gaus(igaus)%xtunk(ndime+1)
          xftra(idime,ndime+1) = gaus(igaus)%xtunk(idime) / gaus(igaus)%xpres
          xftra(idime,ndime+2) = - gaus(igaus)%xtunk(idime) / gaus(igaus)%xtemp
          xftra(ndime+2,idime) = gaus(igaus)%xtunk(idime)
       end do
       xftra(ndime+1,ndime+1) = gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xpres
       xftra(ndime+1,ndime+2) = - gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xtemp
       xftra(ndime+2,ndime+1) = &
            0.5_rp * velsq / gaus(igaus)%xrgasc / gaus(igaus)%xtemp + 1.0_rp / gaus(igaus)%xrgacv
       xftra(ndime+2,ndime+2) = - 0.5_rp * gaus(igaus)%xtunk(ndime+1) * velsq / gaus(igaus)%xtemp

       do idime=1,ndime
          xftra_inv(idime,idime) = 1.0_rp / gaus(igaus)%xtunk(ndime+1)
          xftra_inv(idime,ndime+1) = - gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1)
          xftra_inv(ndime+1,idime) = - gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idime)
          xftra_inv(ndime+2,idime) = &
               - gaus(igaus)%xvelo(idime) / gaus(igaus)%xheatcv / gaus(igaus)%xtunk(ndime+1)
       end do
       xftra_inv(ndime+1,ndime+1) = 0.5_rp * gaus(igaus)%xrgacv * velsq
       xftra_inv(ndime+1,ndime+2) = gaus(igaus)%xrgacv
       xftra_inv(ndime+2,ndime+1) = (velsq - gaus(igaus)%xtunk(ndime+2) / gaus(igaus)%xtunk(ndime+1)) &
            / gaus(igaus)%xheatcv / gaus(igaus)%xtunk(ndime+1)
       xftra_inv(ndime+2,ndime+2) = 1.0_rp / gaus(igaus)%xheatcv / gaus(igaus)%xtunk(ndime+1)

!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              do ldofn=1,ndofn_nsa
!!$                 do mdofn=1,ndofn_nsa
!!$                    gaus(igaus)%xlopr_conservative(idofn,jdofn) = gaus(igaus)%xlopr_conservative(idofn,jdofn) + xftra_inv(idofn,ldofn) &
!!$                         *xlopr_symmet(ldofn,mdofn)*xftra(mdofn,jdofn)
!!$                 end do
!!$              end do
!!$           end do
!!$        end do


       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
             do ldofn=1,ndofn_nsa
                gaus(igaus)%xlopr_conservative(idofn,jdofn) = &
                     gaus(igaus)%xlopr_conservative(idofn,jdofn) + xftra(idofn,ldofn) &
                     *xlopr_symmet(ldofn,jdofn)
             end do
          end do
       end do

!!! END: CHOI & MERKLE preconditioner with conservative variables (we apply a change of variables) 

!!! BEGIN: dP_conservative computation
       termB = 1.0_rp + 0.5_rp * velsq / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
       termC = xbeta * xmach_ref * xmach_ref / gaus(igaus)%xrgacv - 0.5_rp * velsq
       termD = 1.0_rp / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
       do kdime=1,ndime
          term1 = 0.0_rp
          do idime=1,ndime
             term1 = term1 + gaus(igaus)%xvelo(idime) * gaus(igaus)%gvelo(idime,kdime)
          end do
          termB_der = (term1 - 0.5_rp * velsq * gaus(igaus)%gtemp(kdime) / gaus(igaus)%xtemp ) &
               / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
          termC_der = (gbeta(kdime)*xmach_ref*xmach_ref + 2.0_rp*xbeta*xmach_ref * gmach_ref(kdime)) & 
               / gaus(igaus)%xrgacv - term1
          termD_der = - gaus(igaus)%gtemp(kdime)/gaus(igaus)%xtemp/gaus(igaus)%xtemp/gaus(igaus)%xheatcp
          do idime=1,ndime
             termA = gaus(igaus)%xvelo(idime) / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
             termA_der = &
                  (gaus(igaus)%gvelo(idime,kdime)-gaus(igaus)%xvelo(idime)*gaus(igaus)%gtemp(kdime)/gaus(igaus)%xtemp) / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
             do jdime=1,ndime
                glopr(idime,jdime,kdime) = &
                     ( ( gaus(igaus)%xvelo(idime)*gaus(igaus)%gvelo(jdime,kdime)+gaus(igaus)%xvelo(jdime)*gaus(igaus)%gvelo(idime,kdime) ) &
                     - gaus(igaus)%xvelo(idime)*gaus(igaus)%xvelo(jdime)*gaus(igaus)%gtemp(kdime)/gaus(igaus)%xtemp ) / gaus(igaus)%xtemp / gaus(igaus)%xheatcp
             end do
             glopr(ndime+1,idime,kdime) = termA_der
             glopr(ndime+2,idime,kdime) = gaus(igaus)%gvelo(idime,kdime) * termB + gaus(igaus)%xvelo(idime) * termB_der
             glopr(idime,ndime+1,kdime) = termA_der * termC + termA * termC_der
             glopr(idime,ndime+2,kdime) = - termA_der
          end do
          glopr(ndime+1,ndime+1,kdime) = termD_der * termC + termD * termC_der
          glopr(ndime+1,ndime+2,kdime) = - termD_der
          glopr(ndime+2,ndime+1,kdime) = termC_der * termB + termC * termB_der
          glopr(ndime+2,ndime+2,kdime) = - termB_der
       end do
!!! END: dP_conservative computation

       ! Compute: P_conservative * A + dP_conservative * K
       !          P_conservative * K

       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
             do jdime=1,ndime
                do pdofn=1,ndofn_nsa
                   do idime=1,ndime
                      xconv_lopre_diff(idofn,jdofn,jdime) = xconv_lopre_diff(idofn,jdofn,jdime) + &
                           glopr(idofn,pdofn,idime) * xdiff_ori(pdofn,jdofn,idime,jdime)
                   end do
                end do
             end do
          end do
       end do

       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
!             if (kfl_pseud_nsa == 1) &
                  gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
                  + gaus(igaus)%xlopr_conservative(idofn,jdofn) * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
             do idime=1,ndime
                do pdofn=1,ndofn_nsa
                   gaus(igaus)%xconv(idofn,jdofn,idime) = &
                        gaus(igaus)%xconv(idofn,jdofn,idime) + gaus(igaus)%xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                   do jdime=1,ndime
                      gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) = gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) + &
                           gaus(igaus)%xlopr_conservative(idofn,pdofn)*xdiff_ori(pdofn,jdofn,idime,jdime)
                   end do
                end do
                gaus(igaus)%xconv(idofn,jdofn,idime) = gaus(igaus)%xconv(idofn,jdofn,idime) + xconv_lopre_diff(idofn,jdofn,idime)
             end do
          end do
!          if (kfl_pseud_nsa == 1) &
               gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
       end do


!!$        ! Compute: P_conservative * A
!!$        !          P_conservative * K
!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              if (kfl_pseud_nsa == 1) &
!!$                   gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
!!$                   + gaus(igaus)%xlopr_conservative(idofn,jdofn) * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
!!$              do idime=1,ndime
!!$                 do pdofn=1,ndofn_nsa
!!$                    gaus(igaus)%xconv(idofn,jdofn,idime) = &
!!$                         gaus(igaus)%xconv(idofn,jdofn,idime) + gaus(igaus)%xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
!!$                    do jdime=1,ndime
!!$                       gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) = gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) + &
!!$                            gaus(igaus)%xlopr_conservative(idofn,pdofn)*xdiff_ori(pdofn,jdofn,idime,jdime)
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$           if (kfl_pseud_nsa == 1) gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
!!$        end do

    else if (kfl_lopre_nsa == 4) then  ! CHOI & MERKLE INVISCID UNSTEADY PRECONDITIONER IS APPLIED

       xmach_ref = xmach
       if (xmach_ref < 0.0001_rp) then 
          xmach_ref = 0.0001_rp
       else if (xmach_ref > 1.0_rp) then
          xmach_ref = 1.0_rp
       end if

       xbeta = gaus(igaus)%xsoun * gaus(igaus)%xsoun
       xtaup = 1.0_rp

!!! BEGIN: CHOI & MERKLE INVISCID UNSTEADY preconditioner with conservative variables (any change of variables is applied)
       termA = 1.0_rp / gaus(igaus)%xheatcp / gaus(igaus)%xtemp
       termB = xbeta * xmach_ref * xmach_ref
       termC = 0.5_rp * velsq - gaus(igaus)%xheatcp * gaus(igaus)%xtemp + termB * xtaup
       termD = 1.0_rp / gaus(igaus)%xrgasc / gaus(igaus)%xtemp
       termE = 1.0_rp + 0.5_rp * termA * velsq

       do idime=1,ndofn_nsa
          do jdime=1,ndime
             gaus(igaus)%xlopr_conservative(idime,jdime) = gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(jdime) * termA
          end do
          gaus(igaus)%xlopr_conservative(idime,idime) = gaus(igaus)%xlopr_conservative(idime,idime) + 1.0_rp
          gaus(igaus)%xlopr_conservative(idime,ndime+1) = gaus(igaus)%xvelo(idime) * (termB*termD-1.0_rp-termA*termC)
          gaus(igaus)%xlopr_conservative(idime,ndime+2) = - gaus(igaus)%xvelo(idofn) * termA
          gaus(igaus)%xlopr_conservative(ndime+1,idime) = gaus(igaus)%xvelo(idime) * termA
          gaus(igaus)%xlopr_conservative(ndime+2,idime) = gaus(igaus)%xvelo(idime) * termE
       end do
       gaus(igaus)%xlopr_conservative(ndime+1,ndime+1) = termB * termD - termA * termC
       gaus(igaus)%xlopr_conservative(ndime+1,ndime+2) = - termA
       gaus(igaus)%xlopr_conservative(ndime+2,ndime+1) = 1.0_rp / gaus(igaus)%xrgacv + 0.5_rp * gaus(igaus)%xadgam * xmach * xmach &
            - velsq * (1.0_rp + 0.5_rp * termA * termC)
       gaus(igaus)%xlopr_conservative(ndime+2,ndime+2) = 0.5_rp * velsq * termA
!!! END: CHOI & MERKLE INVISCID UNSTEADY preconditioner with conservative variables (any change of variables is applied)

       ! Compute P_conservative * A
       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
!             if (kfl_pseud_nsa == 1) &
                  gaus(igaus)%xtide(idofn) = &
                  gaus(igaus)%xtide(idofn) + gaus(igaus)%xlopr_conservative(idofn,jdofn) * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
             do idime=1,ndime
                do pdofn=1,ndofn_nsa
                   gaus(igaus)%xconv(idofn,jdofn,idime) = gaus(igaus)%xconv(idofn,jdofn,idime) + gaus(igaus)%xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                end do
             end do

          end do
!          if (kfl_pseud_nsa == 1) &
               gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
       end do



    else ! OTHER LOCAL PRECONDITIONERS ARE APPLIED

       ! to program!! other local preconditioners!!! (MM)

       velmo_xy = sqrt(gaus(igaus)%xvelo(1)*gaus(igaus)%xvelo(1)+gaus(igaus)%xvelo(2)*gaus(igaus)%xvelo(2))
!!$     if (gaus(igaus)%velmo < zensa) then           
!!$        gaus(igaus)%velmo = zensa 
!!$     end if
!!$     if (velmo_xy < zensa) then           
!!$        velmo_xy = zensa 
!!$     end if

       if (xmach < 1.0_rp-zensa_lopre) then
          xbeta = sqrt(1.0_rp-xmach*xmach)
          xtaup = xbeta
       else if (xmach < 1.0_rp) then
          xbeta = sqrt((2.0_rp-zensa_lopre)*zensa_lopre)
          xtaup = xbeta
       else if (xmach < 1.0_rp+zensa_lopre) then
          xbeta = sqrt((2.0_rp+zensa_lopre)*zensa_lopre)
          xtaup = xbeta / xmach
       else
          xbeta = sqrt(xmach*xmach-1.0_rp)
          xtaup = xbeta / xmach
       end if

       ! program another local preconditioner HERE!!!!
!!! BEGIN: VLR preconditioner with symmetrizing variables and streamline coordinates
       xlopr_symmet(1,1) = 1.0_rp + xtaup / xbeta / xbeta
       do idime=2,ndime
          xlopr_symmet(idime,idime) = xtaup
       end do
       xlopr_symmet(1,ndime+1) = - xtaup * xmach / xbeta / xbeta
       xlopr_symmet(ndime+1,1) = xlopr_symmet(1,ndime+1)
       xlopr_symmet(ndime+1,ndime+1) = xtaup * xmach * xmach / xbeta / xbeta
       xlopr_symmet(ndime+2,ndime+2) = 1.0_rp
!!! END: VLR preconditioner with symmetrizing variables and streamline coordinates

!!! BEGIN: VLR preconditioner with conservative variables (we apply a change of variables)
!!$        do idofn=1,ndime
!!$           xftra(1,idofn) = gaus(igaus)%xvelo(idofn) / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo
!!$           xftra(ndime+1,idofn) = - gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idofn) / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
!!$           xftra(ndime+2,idofn) = - gaus(igaus)%xrgasc * gaus(igaus)%xvelo(idofn) / gaus(igaus)%xpres
!!$        end do
!!$        xftra(1,ndime+1) = - gaus(igaus)%velmo / gaus(igaus)%xtunk(ndime+1)
!!$        xftra(2,1) = - gaus(igaus)%xvelo(2) / gaus(igaus)%xtunk(ndime+1) / velmo_xy
!!$        xftra(2,2) = gaus(igaus)%xvelo(1) / gaus(igaus)%xtunk(ndime+1) / velmo_xy
!!$        if (ndime==3) then
!!$           xftra(3,1) = - gaus(igaus)%xvelo(1) * gaus(igaus)%xvelo(3) /  gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo / velmo_xy
!!$           xftra(3,2) = - gaus(igaus)%xvelo(2) * gaus(igaus)%xvelo(3) /  gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo / velmo_xy
!!$           xftra(3,3) = velmo_xy /  gaus(igaus)%xtunk(ndime+1) / velmo
!!$        end if
!!$        xftra(ndime+1,ndime+1) = 0.5_rp * gaus(igaus)%xrgacv * gaus(igaus)%xsoun * xmach * xmach / gaus(igaus)%xtunk(ndime+1)
!!$        xftra(ndime+2,ndime+1) = - gaus(igaus)%xheatcv * gaus(igaus)%xsoun * gaus(igaus)%xsoun / gaus(igaus)%xpres + 0.5_rp * gaus(igaus)%xrgasc * gaus(igaus)%xsoun * gaus(igaus)%xsoun * xmach * xmach / gaus(igaus)%xpres
!!$        xftra(ndime+1,ndime+2) = gaus(igaus)%xrgacv / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun 
!!$        xftra(ndime+2,ndime+2) = gaus(igaus)%xrgasc / gaus(igaus)%xpres
!!$        
!!$        do idofn=1,ndime
!!$           xftra_inv(idofn,1) = gaus(igaus)%xvelo(idofn) * gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo
!!$           xftra_inv(idofn,ndime+1) = gaus(igaus)%xvelo(idofn) * gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
!!$           xftra_inv(idofn,ndime+2) = - gaus(igaus)%xvelo(idofn) * gaus(igaus)%xpres / gaus(igaus)%xheatcv / gaus(igaus)%xsoun / gaus(igaus)%xsoun
!!$        end do
!!$        xftra_inv(ndime+2,1) = gaus(igaus)%velmo * gaus(igaus)%xtunk(ndime+1)
!!$        xftra_inv(1,2) = - gaus(igaus)%xvelo(2) * gaus(igaus)%xtunk(ndime+1) / velmo_xy
!!$        xftra_inv(2,2) = gaus(igaus)%xvelo(1) * gaus(igaus)%xtunk(ndime+1) / velmo_xy
!!$        if (ndime==3) then
!!$           xftra_inv(1,3) = - gaus(igaus)%xvelo(1) * gaus(igaus)%xvelo(3) *  gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo / velmo_xy
!!$           xftra_inv(2,3) = - gaus(igaus)%xvelo(2) * gaus(igaus)%xvelo(3) *  gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo / velmo_xy
!!$           xftra_inv(3,3) = velmo_xy *  gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%velmo
!!$        end if
!!$        xftra_inv(ndime+1,ndime+1) = gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
!!$        xftra_inv(ndime+2,ndime+1) = gaus(igaus)%xtunk(ndime+1) * gaus(igaus)%xsoun * (gaus(igaus)%xheatcv/gaus(igaus)%xrgasc+0.5_rp*xmach*xmach)
!!$        xftra_inv(ndime+1,ndime+2) = - gaus(igaus)%xpres / gaus(igaus)%xheatcv / gaus(igaus)%xsoun / gaus(igaus)%xsoun 
!!$        xftra_inv(ndime+2,ndime+2) = - 0.5_rp * gaus(igaus)%xpres * xmach * xmach / gaus(igaus)%xheatcv
!!$        
!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              do ldofn=1,ndofn_nsa
!!$                 do mdofn=1,ndofn_nsa
!!$                    gaus(igaus)%xlopr_conservative(idofn,jdofn) = gaus(igaus)%xlopr_conservative(idofn,jdofn) + xftra_inv(idofn,ldofn) &
!!$                         *xlopr_symmet(ldofn,mdofn)*xftra(mdofn,jdofn)
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!! END: VLR preconditioner with conservative variables (we apply a change of variables) 


!!! BEGIN: VLR preconditioner with conservative variables (we apply two change of variables)
       xttra(ndime+1,ndime+2) = gaus(igaus)%xrgacv / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
       xttra(ndime+2,ndime+1) = - gaus(igaus)%xheatcv * gaus(igaus)%xsoun * gaus(igaus)%xsoun / gaus(igaus)%xpres
       xttra(ndime+2,ndime+2) = gaus(igaus)%xrgasc / gaus(igaus)%xpres

       xttra_inv(ndime+1,ndime+1) = gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
       xttra_inv(ndime+1,ndime+2) = - gaus(igaus)%xpres / gaus(igaus)%xsoun / gaus(igaus)%xsoun / gaus(igaus)%xheatcv
       xttra_inv(ndime+2,ndime+1) = gaus(igaus)%xtunk(ndime+1) * gaus(igaus)%xsoun / gaus(igaus)%xrgacv

       xqtra(ndime+1,ndime+1) = 1.0_rp
       xqtra(ndime+2,ndime+2) = 1.0_rp
       xqtra(2,1) = - gaus(igaus)%xvelo(2) / velmo_xy
       xqtra(2,2) = gaus(igaus)%xvelo(1) / velmo_xy

       if (velmo_xy < zensa) then
          xqtra(2,1) = -1.0_rp !!(MM) Veure signe!!!
          xqtra(2,2) = 1.0_rp !!(MM) Veure signe!!!
       end if

       xqtra_inv(ndime+1,ndime+1) = 1.0_rp
       xqtra_inv(ndime+2,ndime+2) = 1.0_rp
       xqtra_inv(1,2) = xqtra(2,1)
       xqtra_inv(2,2) = xqtra(2,2)

       do idime=1,ndime
          xttra(idime,idime) = 1.0_rp / gaus(igaus)%xtunk(ndime+1)
          xttra(idime,ndime+1) = - gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1)
          xttra(ndime+1,idime) = - gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
          xttra(ndime+1,ndime+1) = xttra(ndime+1,ndime+1) + 0.5_rp * gaus(igaus)%xrgacv * gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(idime) / gaus(igaus)%xtunk(ndime+1) / gaus(igaus)%xsoun
          xttra(ndime+2,idime) = - gaus(igaus)%xrgasc * gaus(igaus)%xvelo(idime) / gaus(igaus)%xpres
          xttra(ndime+2,ndime+1) = xttra(ndime+2,ndime+1) + 0.5_rp * gaus(igaus)%xrgasc * gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(idime) / gaus(igaus)%xpres

          xttra_inv(idime,idime) = gaus(igaus)%xtunk(ndime+1)
          xttra_inv(idime,ndime+1) = gaus(igaus)%xtunk(idime) / gaus(igaus)%xsoun
          xttra_inv(idime,ndime+2) = - gaus(igaus)%xvelo(idime) * gaus(igaus)%xpres / gaus(igaus)%xsoun / gaus(igaus)%xsoun / gaus(igaus)%xheatcv
          xttra_inv(ndime+2,idime) = gaus(igaus)%xtunk(idime)
          xttra_inv(ndime+2,ndime+1) = xttra_inv(ndime+2,ndime+1) + 0.5_rp * gaus(igaus)%xtunk(ndime+1) * gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(idime) / gaus(igaus)%xsoun
          xttra_inv(ndime+2,ndime+2) = xttra_inv(ndime+2,ndime+2) - 0.5_rp * gaus(igaus)%xpres * gaus(igaus)%xvelo(idime) * gaus(igaus)%xvelo(idime) / gaus(igaus)%xheatcv / gaus(igaus)%xsoun / gaus(igaus)%xsoun

          xqtra(1,idime) = gaus(igaus)%xvelo(idime) / gaus(igaus)%velmo           
          if (gaus(igaus)%velmo < zensa) then
             xqtra(1,idime) = 1.0_rp !!(MM) Veure signe!!!
          end if
          xqtra_inv(idime,1) = xqtra(1,idime)           
       end do

       if (ndime == 3) then
          xqtra(3,1) = - gaus(igaus)%xvelo(1) * gaus(igaus)%xvelo(3) / gaus(igaus)%velmo / velmo_xy
          xqtra(3,2) = - gaus(igaus)%xvelo(2) * gaus(igaus)%xvelo(3) / gaus(igaus)%velmo / velmo_xy
          xqtra(3,3) = velmo_xy / gaus(igaus)%velmo

          if (velmo_xy < zensa) then
             xqtra(3,1) = -1.0_rp  !!(MM) Veure signe!!!
             xqtra(3,2) = -1.0_rp  !!(MM) Veure signe!!!
          end if
          if (gaus(igaus)%velmo < zensa) then
             xqtra(3,3) = 1.0_rp
          end if

          xqtra_inv(1,3) = xqtra(3,1)
          xqtra_inv(2,3) = xqtra(3,2)
          xqtra_inv(3,3) = xqtra(3,3)
       end if

       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
             do kdofn=1,ndofn_nsa
                do ldofn=1,ndofn_nsa
                   do mdofn=1,ndofn_nsa
                      do odofn=1,ndofn_nsa
                         gaus(igaus)%xlopr_conservative(idofn,jdofn) &
                              = gaus(igaus)%xlopr_conservative(idofn,jdofn) + xttra_inv(idofn,kdofn)*xqtra_inv(kdofn,ldofn) &
                              *xlopr_symmet(ldofn,mdofn)*xqtra(mdofn,odofn)*xttra(odofn,jdofn)
                      end do
                   end do
                end do
             end do
          end do
       end do
!!! END: VLR preconditioner with conservative variables (we apply two change of variables) 

       ! Compute P_conservative * A
       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa
!             if (kfl_pseud_nsa == 1) &
                  gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) &
                  + gaus(igaus)%xlopr_conservative(idofn,jdofn) * (gaus(igaus)%xunkn(jdofn,ITER_K)-gaus(igaus)%xunkn(jdofn,TIME_N))
             do idime=1,ndime
                do pdofn=1,ndofn_nsa
                   gaus(igaus)%xconv(idofn,jdofn,idime) = &
                        gaus(igaus)%xconv(idofn,jdofn,idime) + gaus(igaus)%xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                end do
             end do
          end do
!          if (kfl_pseud_nsa == 1) &
               gaus(igaus)%xtide(idofn) = gaus(igaus)%xtide(idofn) * dtinv
       end do

!!$     do idofn=1,ndofn_nsa
!!$        do jdofn=1,ndofn_nsa
!!$           do idime=1,ndime
!!$              do kdofn=1,ndofn_nsa
!!$                 do ldofn=1,ndofn_nsa
!!$                    do mdofn=1,ndofn_nsa
!!$                       do odofn=1,ndofn_nsa
!!$                          do pdofn=1,ndofn_nsa
!!$                             gaus(igaus)%xconv(idofn,jdofn,idime) = gaus(igaus)%xconv(idofn,jdofn,idime) + xttra_inv(idofn,kdofn)*xqtra_inv(kdofn,ldofn) &
!!$                                  *xlopr_symmet(ldofn,mdofn)*xqtra(mdofn,odofn)*xttra(odofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
!!$                          end do
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end do

    end if


  end subroutine nsa_newelmlocalpreconditioner


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    09/06/2015
  !> @brief   Compute residuals and sources at the gauss points
  !> @details Compute residuals and sources at the gauss points
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmresidualsandsources(igaus,gaus)
    use      def_master
    use      def_domain
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip), intent(in):: igaus

    type(elm_ongaus_interp_nsa)        :: gaus(mgaus)

    integer(ip)  ::  idofn,jdofn,idime,jdime

    !
    ! volume forces: gravity, dissipation (if LES) and diffusion of enthalpy (if species)
    !

    gaus(igaus)%xresi(ndime+2) = gaus(igaus)%heats  ! Source term of chemical reactions
    gaus(igaus)%xvofo(ndime +2,ndime + 1) = gaus(igaus)%sgsdi + gaus(igaus)%dhtra

    do idofn=1,ndofn_nsa
!       if (kfl_pseud_nsa == 1) &
            gaus(igaus)%xresi(idofn) = gaus(igaus)%xresi(idofn) - gaus(igaus)%xtide(idofn)
       do jdofn=1,ndofn_nsa
          gaus(igaus)%xvofo(idofn,jdofn) = gaus(igaus)%xvofo(idofn,jdofn) + gravm_nsa(idofn,jdofn)

          gaus(igaus)%xresi(idofn)= gaus(igaus)%xresi(idofn) + gaus(igaus)%xvofo(idofn,jdofn)*gaus(igaus)%xtunk(jdofn)
          do idime=1,ndime
             gaus(igaus)%xresi(idofn)= gaus(igaus)%xresi(idofn) - gaus(igaus)%xconv(idofn,jdofn,idime) * gaus(igaus)%gunkn(jdofn,idime)
             do jdime=1,ndime
                if (kfl_lopre_nsa < 2) then
                   gaus(igaus)%xresi(idofn) = gaus(igaus)%xresi(idofn) + gaus(igaus)%ddiff(idofn,jdofn,idime,1) * gaus(igaus)%gunkn(jdofn,jdime) + &
                        gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) * hunkn_elm(jdofn,idime,jdime)
                else
                   ! (MM) This is an approximation. gaus(igaus)%ddiff(idofn,jdofn,idime,1) * gaus(igaus)%gunkn(jdofn,jdime) tendría que estar
                   ! Pero ddiff no tiene el precondicionador dentro !!!! hay que programarlo !!!!!
!!!                   gaus(igaus)%xresi(idofn) = gaus(igaus)%xresi(idofn) + gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) * hunkn_elm(jdofn,idime,jdime)
                end if
                if  (kfl_taudi_nsa >= 5) then   ! non-diagonal tau
                   call runend('NSA_ELMOPERATIONS: NON-DIAGONAL TAU ONLY POSSIBLE IN THE OLD ELCONS VERSION')
                end if
             end do
          end do

       end do
    end do


!    if (ielem == 1 .and. igaus==1) then
!!       write(6,*) gaus(igaus)%xtide(1:ndofn_nsa)
!       write(6,*) gaus(igaus)%xunkn(1:ndofn_nsa,ITER_K)
!       write(6,*) gaus(igaus)%xunkn(1:ndofn_nsa,TIME_N)
!    end if

  end subroutine nsa_newelmresidualsandsources



  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Compute VMS stabilisation terms (la monyos, the diagonal tau)
  !> @details Compute VMS stabilisation terms (la monyos, the diagonal tau)
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmvmsdiagonal(ielem,igaus,gaus)
    use      def_master
    use      def_domain
    use      def_parame
    use      def_nastal
    use      def_kermod
    use      mod_ker_proper

    implicit none

    integer(ip), intent(in):: ielem
    integer(ip), intent(in):: igaus

    type(elm_ongaus_interp_nsa)        :: gaus(mgaus)

    integer(ip)  :: ndofc,ndofe,idime,idofn

    real(rp)    ::  &
         hconv,hleti,hsoun,reate, &
         tauen,taupa,tauco,taush, & 
         tauxi(ndofn_nsa,5),xauxi,xsute,xunte,&
         gdens(ndime),tiint,tiine,tiinc,xdens,&
         cfl_veloc, &
         zensa_lopre_stabi,velmo_lopre_stabi


    if (kfl_stabi_nsa == 0) return

    reate = (1.0_rp + gaus(igaus)%xrgasc / gaus(igaus)%xheatcv)*gaus(igaus)%dvelo

    reate = reate + gaus(igaus)%heats  ! add heat source term from chemical reactions
    ! heats = 0.0 if chemistry off
    if (gaus(igaus)%dvelo < 0.0_rp) then
       reate= -reate
       !     reate= 0.0_rp
    end if

    ! chales are defined in nsa_elconsxy

    hconv = hmini_elm
    hsoun = hmini_elm
    hleti = hmini_elm*hmaxi_elm
    if (kfl_higha_nsa == 1) hleti= hmini_elm*hmini_elm

    ndofc= ndime+1
    ndofe= ndime+2

    xdens        = gaus(igaus)%xunkn(ndofc,    1)
    gdens(    1) = gaus(igaus)%gunkn(ndofc,    1)
    gdens(    2) = gaus(igaus)%gunkn(ndofc,    2)
    gdens(ndime) = gaus(igaus)%gunkn(ndofc,ndime)

    tauxi= 0.0_rp


    xauxi= 0.5_rp + gaus(igaus)%xheatcv/gaus(igaus)%xadgam / gaus(igaus)%xrgasc
    xauxi= 1.0_rp / xauxi

    ! c = sqrt(gamma R T) = sqrt(gamma p / rho) => p = c*c rho / gamma


    zensa_lopre_stabi= 1.0e-6
    velmo_lopre_stabi= gaus(igaus)%velmo

    if(kfl_lopre_nsa > 1 .and. gaus(igaus)%velmo < zensa_lopre_stabi) then           
       !       velmo_lopre_stabi = zensa_lopre_stabi
       !       it looks like it is much better and robust to put velmo_lopre_stabi to zero
       velmo_lopre_stabi = 0.0_rp
    end if

    if (kfl_taufa_nsa(1,2) == 1) then
       tauxi(      1,1) = velmo_lopre_stabi/hconv
       tauxi(ndime+1,1) = velmo_lopre_stabi/hconv
       tauxi(ndime+2,1) = velmo_lopre_stabi/hconv
    end if

    if (kfl_reate_nsa == 0) reate=0.0_rp

    if (kfl_taufa_nsa(2,2) == 1) then
       tauxi(      1,2) = reate
       tauxi(ndime+1,2) = reate
       tauxi(ndime+2,2) = reate
    end if

    if (kfl_taufa_nsa(3,2) == 1) then
       tauxi(      1,3) = gaus(igaus)%xsoun / hsoun
       tauxi(ndime+1,3) = gaus(igaus)%xsoun / hsoun
       tauxi(ndime+2,3) = gaus(igaus)%xsoun / hsoun
    end if
    if (kfl_taufa_nsa(4,2) == 1) then     
       tauxi(      1,4) = 4.0_rp * gaus(igaus)%xvisc / hleti / xdens
       tauxi(ndime+2,4) = 4.0_rp * gaus(igaus)%xdith / hleti / xdens / gaus(igaus)%xheatcp
    end if

    ! esto es para usar la velocidad del choque

    !  tauxi(      1,5) = sspee / hleng(ndime)
    !  tauxi(ndime+1,5) = sspee / hleng(ndime)
    !  tauxi(ndime+2,5) = sspee / hleng(ndime)

    if (kfl_lopre_nsa > 1) then
       !!tauxi(      1,4) = 0.0_rp
       !!tauxi(ndime+2,4) = 0.0_rp 
       if (kfl_lopre_nsa == 3) then ! CHOI & MERKLE PRECONDITIONER IS APPLIED
          if (kfl_pseud_nsa == 1 .and. kfl_taufa_nsa(1,2) == 1) then !  PSEUDO TIME IS APPLIED
             
             cfl_veloc = velmo_lopre_stabi / dtinv / hconv

! esto de aca abajo va mucho peor:
!             dtinv_phys_local=  (tauxi(1,1) + tauxi(1,3) + tauxi(1,4))/qufac_elm
 !            cfl_veloc = velmo_lopre_stabi / dtinv_phys_local / hconv 

             if (cfl_veloc < zensa) cfl_veloc = zensa
!           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc) / 2.0_rp
!           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc) / 4.0_rp  ! too diffusive for the implicit
           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc)           ! very low diffusivity, but not enough perhaps for sharp gradients


             tauxi(      1,1) = tauxi(      1,1) * cfl_veloc 
             tauxi(ndime+1,1) = tauxi(ndime+1,1) * cfl_veloc 
             tauxi(ndime+2,1) = tauxi(ndime+2,1) * cfl_veloc 
          end if
       end if
       tauxi(      1,3) = 0.0_rp
       tauxi(ndime+1,3) = 0.0_rp
       tauxi(ndime+2,3) = 0.0_rp
    end if

    !  tiint= tauxi(      1,1) + tauxi(      1,2) + tauxi(      1,3) + tauxi(      1,4) + tauxi(      1,5) 
    !  tiinc= tauxi(ndime+1,1) + tauxi(ndime+1,2) + tauxi(ndime+1,3)                    + tauxi(ndime+1,5)  
    !  tiine= tauxi(ndime+2,1) + tauxi(ndime+2,2) + tauxi(ndime+2,3) + tauxi(ndime+2,4) + tauxi(ndime+2,5)

    tiint= tauxi(      1,1) + tauxi(      1,3) + tauxi(      1,4)
    tiinc= tauxi(ndime+1,1) + tauxi(ndime+1,3)                    
    tiine= tauxi(ndime+2,1) + tauxi(ndime+2,3) + tauxi(ndime+2,4) 


    taush= 0.0_rp
    !    if (sspee > sound(igaus)) then
    !       taush= qufac_elm * hleng(ndime) / sspee      
    !    end if

    tauen= 0.0_rp
    taupa= 0.0_rp
    tauco= 0.0_rp

    if (tiine > zensa) tauen= qufac_elm / tiine
    if (tiint > zensa) taupa= qufac_elm / tiint
    if (tiinc > zensa) tauco= qufac_elm / tiinc

!!$  ! Tau parameter as computed in López&Nigro's paper
!!$  tauen= 0.0_rp
!!$  taupa= 0.0_rp
!!$  tauco= 0.0_rp
!!$  xmach = velmo(igaus) / sound(igaus)
!!$  cfl_sound = sound(igaus) / dtinv / hsoun
!!$  xmach_ref = sqrt(xmach * xmach + 1.0_rp / cfl_sound / cfl_sound)
!!$  a = 1.0_rp + xmach_ref * xmach_ref
!!$  sound_pseud = velmo(igaus) * velmo(igaus) * a * a + &
!!$       4.0_rp * xmach_ref * xmach_ref * (sound(igaus) * sound(igaus) - velmo(igaus) * velmo(igaus))
!!$  sound_pseud = sqrt(sound_pseud)
!!$  b = (velmo(igaus) * a + sound_pseud) / hleng(ndime)
!!$  tauen = b * b + 4.0_rp * dtinv * dtinv
!!$  tauen = 1.0_rp / sqrt(tauen)
!!$  taupa = tauen
!!$  tauco = tauen

    if (kfl_stabi_nsa == 2) then
       call runend('NSA_MONYOS: CG DEPRECATED - USE MULTISCALE STABILIZATION!!')
    end if

    gaus(igaus)%taudi(1:ndime) = taupa + taush
    gaus(igaus)%taudi(ndofc)   = tauco + taush
    gaus(igaus)%taudi(ndofe)   = tauen + taush


    if (kfl_isent_nsa == 1) gaus(igaus)%taudi(ndofe)=0.0_rp

    xsute= 0.0_rp
    xunte= 0.0_rp

    !
    ! Update local subscale values
    !
    do idofn=1,ndofn_nsa
       !     xsute = xsube(idofn,igaus,2) / xdtix(idofn,igaus,2)  ! it does not work with this term     
       gaus(igaus)%xsube(idofn,1) = gaus(igaus)%taudi(idofn) * ( xsute + gaus(igaus)%xresi(idofn) - xunte) 
    end do

    !
    ! Update global subscale values
    !
    ! xsube(..,..,1) = new subscale, coming from nsa_monyos
    ! xsube(..,..,2) = subscale of the last time step, coming from nsa_gauvalxy
    
    do idime=1,ndime
       umosg_nsa(idime,ielem,igaus,2) = umosg_nsa(idime,ielem,igaus,1) 
       umosg_nsa(idime,ielem,igaus,1) = gaus(igaus)%xsube(idime,1)
    end do

    densg_nsa(ielem,igaus,2) = densg_nsa(ielem,igaus,1)
    enesg_nsa(ielem,igaus,2) = enesg_nsa(ielem,igaus,1)
    densg_nsa(ielem,igaus,1) = gaus(igaus)%xsube(ndime+1,1)
    enesg_nsa(ielem,igaus,1) = gaus(igaus)%xsube(ndime+2,1)


    !if (ielem == 1 .and. igaus == 1) then
!!$print*
!!$print*,'dtinv_nsa',itinn(modul),1.0_rp/dtinv_nsa,1.0_rp/dtinv
!!$print*,'xsube',xsube(ndime+1,igaus,1)
    !   write(6,*) 'tololo',dtinv,velmo(igaus),sound(igaus)
    !end if

  end subroutine nsa_newelmvmsdiagonal

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Compute shock capturing operations
  !> @details Compute shock capturing operations
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmshocap(igaus,gaus)
    use      def_kintyp
    use      def_domain
    use      def_nastal

    implicit none

    integer(ip), intent(in):: igaus

    type(elm_ongaus_interp_nsa)        :: gaus(mgaus)

    integer(ip) :: if_shocap(5),napsh
    integer(ip)           :: idime,idofn

    real(rp)              :: &
         grau2(ndofn_nsa),grmod(ndofn_nsa),remod(ndofn_nsa),difeq(ndofn_nsa),&
         ficve(ndofn_nsa),shpec,shfau,zesh2,hconv,xmach

    real(rp)              :: epsdc,epssu,epssl,fiso(ndofn_nsa),faniso(ndofn_nsa),uvel,vvel,wvel,xvel2,zesho


    zesho= 1000.0_rp * zensa

    !
    ! itask = 1 -> compute shock capturing difusion tensor for igaus 
    !
    ! initialize
    shmet_elm = 0.0_rp

    ! no shock capturing to x-momentum means no shock capturing at all, then return
    if (kfl_shock_nsa(1) == 0) return

    hconv = chale_elm(1)

    if_shocap = 1             ! default: apply shock capturing
    xvel2 = gaus(igaus)%velmo*gaus(igaus)%velmo
    zesh2 = zesho*zesho
    !     if (xvel2 < zesho) if_shocap = 0            ! velocity is the first threshold  
    if (shock_nsa < zesho) then
       if_shocap = 0
       return
    end if

    napsh = 0
    do idofn = 1,ndofn_nsa
       grau2(idofn)= 0.0_rp
       do idime = 1,ndime                                     
          grau2(idofn) = grau2(idofn) &
               + gaus(igaus)%gunkn(idofn,idime)*gaus(igaus)%gunkn(idofn,idime)            !evaluate grad modul
       end do
       grmod(idofn) = sqrt(grau2(idofn))
       remod(idofn) =  abs(gaus(igaus)%xresi(idofn))
       ficve(idofn) = 0.0_rp
       if (grau2(idofn) < zesho) then
          if_shocap(idofn) = 0
       else
          ficve(idofn) = remod(idofn) / grmod(idofn)
       end if
       xmach= gaus(igaus)%velmo / gaus(igaus)%xsoun

       if (ficve(idofn) < zesho) if_shocap(idofn)= 0

       !          if (ficve(idofn) > zesho) then
       !           write (6,200) idofn,ficve(idofn),sound,velmo
       !200          format (i4,10(2x,f8.4))
       !          end if
       napsh = napsh + if_shocap(idofn)
    end do

    if (napsh == 0) return       ! no one needs shock capturing, return

    ! momentum equation viscosity
    difeq(1:ndime) = gaus(igaus)%xvisc / gaus(igaus)%xunkn(ndime+1,ITER_K)   
    ! continuity (no viscosity terms)
    difeq(ndime+1) = 0.0_rp       
    ! energy (adding thermal diffusion)
    difeq(ndime+2) = gaus(igaus)%xdith / gaus(igaus)%xunkn(ndime+1,ITER_K) / gaus(igaus)%xheatcp

    uvel= gaus(igaus)%xvelo(    1)
    vvel= gaus(igaus)%xvelo(    2)
    wvel= gaus(igaus)%xvelo(ndime)

    do idofn=1,ndofn_nsa
       shfau = shock_nsa
       if (difeq(idofn) > zesh2) then
          shpec = ficve(idofn) * hconv / 2.0_rp / difeq(idofn)           
          if (shpec > 0.0_rp) shfau = shock_nsa - 1.0 / shpec
          if (shfau < 0.0_rp) shfau = 0.0_rp
       end if
       if (idofn == ndime+1) shfau = shock_nsa  ! no diffusion for continuity
       if (resid_nsa(4) < shtol_nsa) then
          epsdc= gaus(igaus)%shocktau_local(idofn)
       else
          epsdc= 0.5_rp * shfau * hconv * ficve(idofn)      
          gaus(igaus)%shocktau_local(idofn) = epsdc
       end if
       fiso(idofn)    = epsdc

       if (kfl_shock_nsa(idofn) == 1) then
          !
          ! Anisotropic shock capturing (DEFAULT): compares epsdc with supg-like difusion and 
          ! compute an anisotropic difusion tensor
          !
          epssu = gaus(igaus)%taudi(idofn) * xvel2
          epssl = epsdc - epssu
          if (epssl < 0.0_rp) then
             epssl= 0.0_rp
          end if
          if (xvel2 > 0.0_rp) then
             faniso(idofn)= (epssl - epsdc)/xvel2     
          else
             faniso(idofn)= 0.0_rp
          end if
          
          shmet_elm(1,1,idofn) =  faniso(idofn)*uvel*uvel + fiso(idofn)
          shmet_elm(2,2,idofn) =  faniso(idofn)*vvel*vvel + fiso(idofn)        
          shmet_elm(1,2,idofn) =  faniso(idofn)*uvel*vvel
          shmet_elm(2,1,idofn) =  shmet_elm(1,2,idofn)
          
          if (ndime == 3) then           
             shmet_elm(3,3,idofn)= faniso(idofn)*wvel*wvel + fiso(idofn)
             shmet_elm(1,3,idofn)= faniso(idofn)*uvel*wvel
             shmet_elm(2,3,idofn)= faniso(idofn)*vvel*wvel
             shmet_elm(3,1,idofn)= shmet_elm(1,3,idofn)
             shmet_elm(3,2,idofn)= shmet_elm(2,3,idofn)
          end if

       else if (kfl_shock_nsa(idofn) == 2) then
          !
          ! Isotropic shock capturing: isotropic epdsc, don't compare
          !
          shmet_elm(1,1,idofn) =  fiso(idofn)
          shmet_elm(2,2,idofn) =  fiso(idofn)        
          shmet_elm(1,2,idofn) =  0.0_rp
          shmet_elm(2,1,idofn) =  0.0_rp
          
          if (ndime == 3) then           
             shmet_elm(3,3,idofn)= fiso(idofn)
             shmet_elm(1,3,idofn)= 0.0_rp
             shmet_elm(2,3,idofn)= 0.0_rp
             shmet_elm(3,1,idofn)= 0.0_rp
             shmet_elm(3,2,idofn)= 0.0_rp
          end if



       end if
          

    end do



  end subroutine nsa_newelmshocap


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmatrixcompute.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Compute jacobians and mass matrices
  !> @details Compute jacobians and mass matrices
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmatrixcompute(&
       inode,igaus,ielem,pnode,pevat,gath,gaus,cartigaus,shapigaus,hessigaus,matri)
    use      def_master
    use      def_nastal
    use      def_domain

    implicit none

    integer(ip), intent(in):: inode
    integer(ip), intent(in):: igaus
    integer(ip), intent(in):: ielem
    integer(ip), intent(in):: pnode
    integer(ip), intent(in):: pevat

    type(elm_ongaus_interp_nsa), intent(in)        :: gaus(mgaus)
    type(elm_onnode_gather_nsa), intent(in)        :: gath(mnode)

    real(rp), intent(in)                   :: &
         cartigaus(ndime,mnode), shapigaus(mnode),hessigaus(ntens,mnode)

    type(elm_onnode_onnode_matrix_nsa), intent(out) :: matri(mnode,mnode)
    type(elm_onnode_onnode_matrix_nsa)              :: matri_init

    integer(ip) :: &
         ITER_NEWTON,ipoin,idime,jdime,itott,idofn,jdofn,kdofn,jnode,ievat,jevat,&
         if_assemble_mat,if_assemble_rhs

!!$    real(rp) :: &
!!$         advec_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Advective part
!!$         diffu_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Diffusive part
!!$         shote_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Shock capturing part
!!$         stabi_matrix(ndofn_nsa,ndofn_nsa,mnode)  ,& ! Stabilization part = Adjoint matrix * Diagonal subscale part
!!$         timas_matrix(ndofn_nsa,ndofn_nsa,mnode,2),& ! Temporal part for time and pseudotime 
!!$         subdi_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Diagonal subscale part


    real(rp) :: &
         adjoi_matrix(ndofn_nsa,ndofn_nsa),&       ! Adjoint matrix
         subes(ndofn_nsa),&                        ! Subscale
         state(ndofn_nsa)  ,galte(ndofn_nsa)  ,rhsmat(ndofn_nsa),&
         state_n(ndofn_nsa),galte_n(ndofn_nsa),shote_n(ndofn_nsa),shote(ndofn_nsa),&
         aumat,autim(ndofn_nsa),autim_n(ndofn_nsa),ausax(ndofn_nsa),&         
         xconv_newton(ndofn_nsa,ndofn_nsa,ndime),&
         resid_n, elunk_value,diff_fact,supg_fact,&
         ximpl_visc, ximpl_conv, ximpl_stab, ximpl_shot,&
         xx_taudi,cn_left,cn_pseudo,explicit_pseudo(5)

    cn_left= 1.0_rp  
    if (kfl_tisch_nsa == 3) cn_left= 0.5_rp  ! crank-nicolson 0.5
    cn_pseudo= 1.0_rp - cn_left
    explicit_pseudo= 1.0_rp

    ximpl_visc= 0.0_rp
    if (kfl_ximpl_nsa(1)==1) ximpl_visc= 0.5_rp
    ximpl_conv= 0.0_rp
    if (kfl_ximpl_nsa(2)==1) ximpl_conv= 0.5_rp
    ximpl_stab= 0.0_rp
    if (kfl_ximpl_nsa(3)==1) ximpl_stab= 0.5_rp
    ximpl_shot= 0.0_rp
    if (kfl_ximpl_nsa(4)==1) ximpl_shot= 0.5_rp

    supg_fact= 1.0_rp
    if (kfl_stabi_nsa == 6) supg_fact=0.0_rp  !generalized supg

    ITER_NEWTON= ITER_K

    adjoi_matrix(1:ndofn_nsa,idofn) = 0.0_rp ! Adjoint matrix
    xconv_newton(1:ndofn_nsa,idofn,1:ndime) = 0.0_rp 

    matri(inode,1:pnode)= matri_init 

!!$    do idofn=1,ndofn_nsa
!!$       do jdofn=1,ndofn_nsa          
!!$          do jnode=1,pnode
!!$             matri(inode,1:pnode)%advec_matrix= 0.0_rp ! Advective part
!!$             matri(inode,1:pnode)%diffu_matrix= 0.0_rp ! Diffusive part
!!$             matri(inode,1:pnode)%shote_matrix= 0.0_rp ! Shock capturing part
!!$             matri(inode,1:pnode)%subdi_matrix= 0.0_rp ! Diagonal subscale part
!!$             matri(inode,1:pnode)%stabi_matrix= 0.0_rp ! Stabi part = Adjoint matrix * Diagonal subscale part
!!$             matri(inode,1:pnode)%timas_matrix= 0.0_rp ! Temporal part
!!$          end do
!!$       enddo
!!$    enddo

    if (kfl_linea_nsa == 2) then
       !
       ! Compute NR contribution
       !

    end if

    do idime=1,ndime
       do idofn=1,ndofn_nsa
          do jdofn=1,ndofn_nsa  
             xx_taudi = gaus(igaus)%taudi(jdofn)
             do jnode=1,pnode
                ! Advective part
                matri(inode,jnode)%advec_matrix(idofn,jdofn) = matri(inode,jnode)%advec_matrix(idofn,jdofn) &
                     + (gaus(igaus)%xconv(idofn,jdofn,idime) + xconv_newton(idofn,jdofn,idime)) &
                     * cartigaus(idime,jnode) * shapigaus(inode)
                ! Diagonal subscale part
                diff_fact= supg_fact * gaus(igaus)%ddiff(idofn,jdofn,idime,1)
                matri(inode,jnode)%subdi_matrix(idofn,jdofn) = matri(inode,jnode)%subdi_matrix(idofn,jdofn) &
                     + gaus(igaus)%xconv(idofn,jdofn,idime) * cartigaus(idime,jnode) * xx_taudi &
                     - diff_fact * cartigaus(idime,jnode) * xx_taudi

                do jdime = 1,ndime
                   diff_fact= supg_fact *gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) 
                   matri(inode,jnode)%subdi_matrix(idofn,jdofn) = matri(inode,jnode)%subdi_matrix(idofn,jdofn) &
                        - diff_fact * hessigaus(nindx_nsa(idime,jdime),jnode) * xx_taudi
                   ! Diffusive part                   
                   matri(inode,jnode)%diffu_matrix(idofn,jdofn) = matri(inode,jnode)%diffu_matrix(idofn,jdofn) &
                        + gaus(igaus)%xdiff(idofn,jdofn,idime,jdime) &
                        * cartigaus(jdime,jnode) * cartigaus(idime,inode)
                end do

                if (idofn == jdofn) then
                   do jdime = 1,ndime
                      ! Shock capturing
                      matri(inode,jnode)%shote_matrix(idofn,jdofn) = matri(inode,jnode)%shote_matrix(idofn,jdofn) &
                           + shmet_elm(jdime,idime,idofn) * cartigaus(jdime,jnode) * cartigaus(idime,inode)
                   end do
                end if

             end do
          end do
       end do
    end do

    do idofn=1,ndofn_nsa        
       do jdofn=1,ndofn_nsa        
          ! Adjoint matrix        
          if (kfl_lopre_nsa == 0) then ! (MM) This 'if' should disappear in the future
             adjoi_matrix(jdofn,idofn) = gaus(igaus)%dconv(jdofn,idofn) * shapigaus(inode)
          else
             adjoi_matrix(jdofn,idofn) = 0.0_rp           
          end if
          do idime=1,ndime
             diff_fact= supg_fact *gaus(igaus)%ddiff(jdofn,idofn,idime,2)
             adjoi_matrix(jdofn,idofn) = adjoi_matrix(jdofn,idofn) &
                  + (gaus(igaus)%xconv(jdofn,idofn,idime) + diff_fact) * cartigaus(idime,inode)
             do jdime = 1,ndime
                diff_fact= supg_fact * gaus(igaus)%xdiff(idofn,jdofn,idime,jdime)
                adjoi_matrix(idofn,jdofn) = adjoi_matrix(idofn,jdofn) &
                     + diff_fact * hessigaus(nindx_nsa(idime,jdime),inode)
             end do

          end do
       end do
    end do

!!$
    do jnode=1,pnode
       do idofn= 1,ndofn_nsa
          ! Temporal part
          matri(inode,jnode)%timas_matrix(idofn,idofn,DT_PSEUDO)    = &
               shapigaus(inode) * shapigaus(jnode) * gath(inode)%dtinv_eqs(idofn,DT_PSEUDO)

          if (kfl_lopre_nsa > 0) then  
             !
             ! If preconditioning and pseudo-time, then transform physical time mass matrix. 
             !  System matrix comes transformed from nsa_gauvalxy.
             !
             do jdofn= 1,ndofn_nsa
                matri(inode,jnode)%timas_matrix(idofn,jdofn,DT_PHYSICAL)  = &
                     shapigaus(inode) * shapigaus(jnode) &
                     * gath(inode)%dtinv_eqs(idofn,DT_PHYSICAL) * gaus(igaus)%xlopr_conservative(idofn,jdofn)
             end do
          else
             matri(inode,jnode)%timas_matrix(idofn,idofn,DT_PHYSICAL)  = &
                  shapigaus(inode) * shapigaus(jnode) * gath(inode)%dtinv_eqs(idofn,DT_PHYSICAL)
          end if

          if (kfl_stabi_nsa >= 1) then
             do jdofn= 1,ndofn_nsa
                do kdofn= 1,ndofn_nsa
                   ! Stabilization part = Adjoint matrix * Diagonal subscale part
                   matri(inode,jnode)%stabi_matrix(idofn,jdofn)= matri(inode,jnode)%stabi_matrix(idofn,jdofn) &
                        + adjoi_matrix(idofn,kdofn) * matri(inode,jnode)%subdi_matrix(kdofn,jdofn)
                end do
             end do
          end if
       end do
    end do

  end subroutine nsa_newelmatrixcompute

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmatrixassemble.f90
  !> @author  Mariano Vazquez
  !> @date    08/06/2015
  !> @brief   Assemble elementary matrices from jacobians and mass matrix
  !> @details Assemble elementary matrices from jacobians and mass matrix
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmatrixassemble(&
       inode,igaus,pnode,pevat,gath,gaus,matri,elrhs,elmat)
    use      def_master
    use      def_nastal
    use      def_domain

    implicit none

    integer(ip), intent(in):: inode
    integer(ip), intent(in):: igaus
    integer(ip), intent(in):: pnode
    integer(ip), intent(in):: pevat

    type(elm_ongaus_interp_nsa), intent(in)        :: gaus(mgaus)
    type(elm_onnode_gather_nsa), intent(in)        :: gath(mnode)
    type(elm_onnode_onnode_matrix_nsa), intent(in) :: matri(mnode,mnode)

    real(rp),    intent(out):: elrhs(pevat),elmat(pevat,pevat)

    ! IMPORTANT: IF ELRHS AND ELMAT ARE HERE DEFINED USING PEVAT, KEEP IT ALL THE WAY UNTIL ASSEMBLY

    integer(ip) :: &
         ITER_NEWTON,idime,jdime,itott,idofn,jdofn,kdofn,jnode,ievat,jevat,&
         if_assemble_mat,if_assemble_rhs

    real(rp) :: &
         subes(ndofn_nsa),&                        ! Subscale
         state(ndofn_nsa)  ,galte(ndofn_nsa)  ,rhsmat(ndofn_nsa),&
         state_n(ndofn_nsa),galte_n(ndofn_nsa),shote_n(ndofn_nsa),shote(ndofn_nsa),&
         aumat,autim(ndofn_nsa),autim_n(ndofn_nsa),ausax(ndofn_nsa),&         
         xconv_newton(ndofn_nsa,ndofn_nsa,ndime),&
         resid_n, elunk_value,diff_fact,supg_fact,&
         ximpl_visc, ximpl_conv, ximpl_stab, ximpl_shot,&
         xx_taudi,cn_left,cn_pseudo,explicit_pseudo(5)


    cn_left= 1.0_rp  
    if (kfl_tisch_nsa == 3) cn_left= 0.5_rp  ! crank-nicolson 0.5
    cn_pseudo= 1.0_rp - cn_left
    explicit_pseudo= 1.0_rp

    ximpl_visc= 0.0_rp
    if (kfl_ximpl_nsa(1)==1) ximpl_visc= 0.5_rp
    ximpl_conv= 0.0_rp
    if (kfl_ximpl_nsa(2)==1) ximpl_conv= 0.5_rp
    ximpl_stab= 0.0_rp
    if (kfl_ximpl_nsa(3)==1) ximpl_stab= 0.5_rp
    ximpl_shot= 0.0_rp
    if (kfl_ximpl_nsa(4)==1) ximpl_shot= 0.5_rp

    supg_fact= 1.0_rp
    if (kfl_stabi_nsa == 6) supg_fact=0.0_rp  !generalized supg

    ITER_NEWTON= ITER_K


    state= 0.0_rp
    shote= 0.0_rp    ! shote is computed here from shmet_nsa matrix
    galte= 0.0_rp
    autim= 0.0_rp     
    rhsmat= 0.0_rp     
    state_n= 0.0_rp
    galte_n= 0.0_rp
    shote_n= 0.0_rp
    autim_n= 0.0_rp     
    subes= 0.0_rp
    ausax= 0.0_rp     

    do idofn=1,ndofn_nsa
!       if (itinn(modul) > 1) then
          !
          ! if pseudo is NOT active, dtinv_eqs(... , DT_PSEUDO) = 0.0, as defined above
          !
          explicit_pseudo(idofn) = gath(inode)%dtinv_eqs(idofn,DT_PSEUDO) + gath(inode)%dtinv_eqs(idofn,DT_PHYSICAL) 
          if (explicit_pseudo(idofn) > 0.0_rp) then
             explicit_pseudo(idofn)= 1.0_rp / explicit_pseudo(idofn)
          end if
!       end if

    enddo

    do jnode=1,pnode
       do jdofn=1,ndofn_nsa
          do idofn= 1,ndofn_nsa
             ievat= (inode-1) * ndofn_nsa + idofn
             jevat= (jnode-1) * ndofn_nsa + jdofn

             if_assemble_rhs= 1
             if (kfl_linea_nsa == 3) then
                if_assemble_rhs = 0
                if (itinn(modul) > 1) then
                   if (jnode == inode) then
                      if (jdofn == idofn) if_assemble_rhs = 1
                   end if
                end if
             end if

             ! explicit formulation terms
             galte(idofn)= galte(idofn) &
                  + (matri(inode,jnode)%advec_matrix(idofn,jdofn) &
                  +  matri(inode,jnode)%diffu_matrix(idofn,jdofn)) * gath(jnode)%elunk(jdofn,ITER_NEWTON)

             state(idofn)= state(idofn) &
                  + matri(inode,jnode)%stabi_matrix(idofn,jdofn) * gath(jnode)%elunk(jdofn,ITER_NEWTON) 

             shote(idofn)= shote(idofn) &
                  + matri(inode,jnode)%shote_matrix(idofn,jdofn) * gath(jnode)%elunk(jdofn,ITER_NEWTON) 

             subes(jdofn)= subes(jdofn) &
                  - matri(inode,jnode)%subdi_matrix(jdofn,idofn) * gath(jnode)%elunk(idofn,ITER_NEWTON)

             if_assemble_mat = 1
             if (kfl_linea_nsa == 3) then
                if (if_assemble_rhs == 1) if_assemble_mat = 0
             end if
             aumat= 0.0_rp

             if (if_assemble_mat == 1) then
                aumat = &
                     matri(inode,jnode)%advec_matrix(idofn,jdofn) &
                     + matri(inode,jnode)%diffu_matrix(idofn,jdofn) &
                     + matri(inode,jnode)%stabi_matrix(idofn,jdofn) &
                     + matri(inode,jnode)%shote_matrix(idofn,jdofn) 

                ! adding elmat because this line computes the contribution of each gauss point
                elmat(ievat,jevat) = elmat(ievat,jevat) &
                     + gaus(igaus)%dvolu * (cn_left * aumat + matri(inode,jnode)%timas_matrix(idofn,jdofn,DT_PHYSICAL))


                !        esto lo hace mas parecido al explicito antiguo
                !              elmat(ievat,jevat) = elmat(ievat,jevat) &
                !                   + dvolu_nsa(igaus) * (cn_left * aumat)

                !
                ! In the case of newton or jacobi, timas_matrix(...., DT_PSEUDO) is for the current time step
                ! When tau-newton or tau-jacobi, it is for the real current iter and we have to add the DT_PHYSICAL.
                !
                elmat(ievat,jevat) = elmat(ievat,jevat) &
                     + gaus(igaus)%dvolu * matri(inode,jnode)%timas_matrix(idofn,jdofn,DT_PSEUDO)


             end if


             elunk_value = gath(jnode)%elunk(jdofn,TIME_N)           
             if (kfl_pseud_nsa == 1) then
                elunk_value = &
                     cn_left * gath(jnode)%elunk(jdofn,ITER_NEWTON) &
                     + cn_pseudo * gath(jnode)%elunk(jdofn,TIME_N)
             end if

             if (kfl_linea_nsa==1) then
                ! this is the right hand side contribution for the delta form Jacobi 
                galte_n(idofn)= galte_n(idofn) &
                     + (matri(inode,jnode)%advec_matrix(idofn,jdofn) &
                     +  matri(inode,jnode)%diffu_matrix(idofn,jdofn)) * elunk_value
                state_n(idofn)= state_n(idofn) &
                     + (matri(inode,jnode)%stabi_matrix(idofn,jdofn)) * elunk_value
                shote_n(idofn)= shote_n(idofn) &
                     + (matri(inode,jnode)%shote_matrix(idofn,jdofn)) * elunk_value
                if (kfl_pseud_nsa == 1) then
                   autim_n(idofn) = autim_n(idofn) &
                        + matri(inode,jnode)%timas_matrix(idofn,jdofn,DT_PHYSICAL) &
                        * (gath(jnode)%elunk(jdofn,ITER_NEWTON) - gath(jnode)%elunk(jdofn,TIME_N))
                end if
             else if (kfl_linea_nsa==2) then
                ! this is the right hand side contribution for the delta form Newton-Raphson 
                autim_n(idofn) = autim_n(idofn) &
                     + matri(inode,jnode)%timas_matrix(idofn,jdofn,DT_PSEUDO) &
                     * (gath(jnode)%elunk(jdofn,ITER_NEWTON) - gath(jnode)%elunk(jdofn,TIME_N))
                galte_n(idofn)= galte_n(idofn) &
                     + (matri(inode,jnode)%advec_matrix(idofn,jdofn) &
                     +  matri(inode,jnode)%diffu_matrix(idofn,jdofn)) * gath(jnode)%elunk(jdofn,ITER_NEWTON)
                state_n(idofn)= state_n(idofn) &
                     + (matri(inode,jnode)%stabi_matrix(idofn,jdofn)) * gath(jnode)%elunk(jdofn,ITER_NEWTON)
                shote_n(idofn)= shote_n(idofn) &
                     + (matri(inode,jnode)%shote_matrix(idofn,jdofn)) * gath(jnode)%elunk(jdofn,ITER_NEWTON)
             end if
          end do
       end do

       ! When no pseudo time step, time contribution to RHS is zero

       if (kfl_pseud_nsa == 0) then
          autim_n(idofn) = 0.0_rp
       end if

       !  esto seria usando xtide_nsa
!!!     autim_n(idofn)= xtide_nsa(idofn,igaus) * xshap_nsa(inode,igaus)

    end do

    if (kfl_stabi_nsa == 0) then
       state   = 0.0_rp     
       state_n = 0.0_rp     
       shote   = 0.0_rp     
       shote_n = 0.0_rp     
    end if

    if (kfl_timet_nsa == 1) then        
       ! explicit: 
       ! no time contribution to elrhs 
       ! only galte, state and shote explicitly computed terms
       ! elmat is not used
       ! shote comes from outside, it is a positive term
       do idofn = 1,ndofn_nsa  
          ievat = (inode-1) * ndofn_nsa + idofn
          ! adding elrhs because this line computes the contribution of each gauss point

          resid_n = gaus(igaus)%dvolu * (galte_n(idofn) + state_n(idofn) + shote_n(idofn) + autim_n(idofn))

!          elrhs(ievat) = elrhs(ievat) &
!               + gaus(igaus)%dvolu * explicit_pseudo(idofn) * &
!               (- shapigaus(inode)*gaus(igaus)%xtide(idofn) - galte(idofn) - state(idofn) - shote(idofn))

          elrhs(ievat) = elrhs(ievat) - resid_n * explicit_pseudo(idofn) 

       end do


    else if (kfl_timet_nsa == 2) then   
       ! implicit: 
       ! time contribution in elrhs (time n) and elmat (time n+1) 
       ! no galte, state and shote explicitly computed terms
       ! elmat used

!!!     do idofn = 1,ndofn_nsa  
!!!        ievat = (inode-1) * ndofn_nsa + idofn
!!!        ! adding elrhs because this line computes the contribution of each gauss point
!!!        elrhs(ievat) = elrhs(ievat) + dvolu_nsa(igaus) * autim(idofn)
!!!     end do

       if (kfl_linea_nsa >= 1) then ! jacobi or inexact newton

          !        iauxi= 0
          !        if (itinn(modul) == 1) then
          !           iauxi= 1
          !        end if
          !        if (iauxi== 1) then
          do idofn = 1,ndofn_nsa  
             ievat = (inode-1) * ndofn_nsa + idofn
             resid_n   = 0.0_rp
             if (kfl_linea_nsa == 1) then
                ! jacobi, delta form
                resid_n   = gaus(igaus)%dvolu * (galte_n(idofn) + state_n(idofn) + shote_n(idofn) + autim_n(idofn))
             else if (kfl_linea_nsa == 2) then
                ! newton raphson
                resid_n   = gaus(igaus)%dvolu * (galte_n(idofn) + state_n(idofn) + shote_n(idofn) + autim_n(idofn))
             end if
             elrhs(ievat) = elrhs(ievat) - resid_n
          end do
          !        end if

       end if

    end if

!!$
!!$
!!$

!!!!!!!!!!!!!!!!
!!!! OJOOOOO QUE ESTO ES UNA PRUEBA PARA HACER EXPLICITOS!!! NO TIENE QUE IR ASI!!
!!$  auele= 0.0_rp
!!$  do idofn = 1,ndofn_nsa  
!!$     ievat = (inode-1) * ndofn_nsa + idofn
!!$     do jdofn = 1,ndofn_nsa  
!!$        do jnode=1,pnode
!!$           jevat = (jnode-1) * ndofn_nsa + jdofn           
!!$           elunk_jevat  = gath(jnode)%elunk(jdofn,ITER_K)
!!$           auele(ievat) = auele(ievat) + elmat(ievat,jevat) * elunk_jevat
!!$        end do
!!$     end do
!!$  end do
!!$  do idofn = 1,ndofn_nsa  
!!$     ievat = (inode-1) * ndofn_nsa + idofn
!!$     auele_rhs = elrhs(ievat) + dvolu_nsa(igaus) * ( shote(idofn) - galte(idofn) - state(idofn))
!!$     auele(ievat) = auele(ievat) - auele_rhs
!!$  end do
!!$
!!$  if ((ielem == 2).and.(igaus==4) .and. (inode==4)) then
!!$     write(6,*)
!!$     do ievat=1,nevat_nsa
!!$        write(6,*) auele(ievat)
!!$     end do
!!$     stop
!!$
!!$  end if

!!!!!!!!!!!!!!!!

    ! para minicucu
    !  if (inode == 4 .and. igaus == 4 .and. itinn(modul)== 1) then
    !     write(6677,100) itinn(modul),elunk(1,2,ITER_NEWTON), elunk(1,2,TIME_N),autim_n(1),galte_n(1),&
    !          cn_left * elunk(1,2,ITER_NEWTON) + cn_pseudo * elunk(1,2,TIME_N)
    !     write(6677,100) kfl_pseud_nsa,dtinv_eqs(1,DT_PSEUDO),dtinv_eqs(1,DT_PHYSICAL)
    !  end if
    !100 format(i,10(2x,e))



  end subroutine nsa_newelmatrixassemble


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    11/06/2015
  !> @brief   Boundary conditions assembly for implicit scheme
  !> @details Boundary conditions assembly for implicit schemes \n
  !!          Conditions are always set on physical variables.
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmatrixsetboundary(&
       ielem,pnode,pevat,lnode,gath,elrhs,elmat,kfl_matvec)
    use      def_master
    use      def_nastal, only: cvcoe_nsa,ncomp_nsa,ndofn_nsa,kfl_fixrs_nsa,jacrot_du_dq_nsa,&
         jacrot_dq_du_nsa,runiv_nsa,kfl_linea_nsa
    use      def_domain, only: ndime,mnode,skcos,exnor,lpoty,lnods

    implicit none

    integer(ip), intent(in):: ielem
    integer(ip), intent(in):: pnode
    integer(ip), intent(in):: pevat
    integer(ip), intent(in):: lnode(pnode) 
    integer(ip), intent(in):: kfl_matvec   !> matrix-vector or only vector?
    real(rp),    intent(out):: elrhs(pevat),elmat(pevat,pevat)

    ! IMPORTANT: IF ELRHS AND ELMAT ARE HERE DEFINED USING PEVAT, KEEP IT ALL THE WAY UNTIL ASSEMBLY

    type(elm_onnode_gather_nsa), intent(in)        :: gath(mnode)

    integer(ip)    :: inode,jnode,ipoin,ibopo,idime,jdime,ievat,jevat,if_jacobians
    integer(ip)    :: kfl_elfix(ndofn_nsa,mnode),kinfl,iroty,idofn,jdofn,nofix,nofix_vel,iffix

    real(rp)  :: &
         elbve(ndofn_nsa,mnode),elphy(ndofn_nsa,mnode),velsq,xveln(ndime),xadve(ndime),&
         xpres,xdens,xtemp,xener,xmome(ndime), &
         rotgl(ndime,ndime),rotlg(ndime,ndime),rotqu_aux(ndofn_nsa,ndofn_nsa),&
         rotuq_aux(ndofn_nsa,ndofn_nsa),xmach_chkinf, adiag,&
         rgacv,rgasc,xhecv


    !
    ! Boundary conditions: Dirichlet (essential B.Cs)
    !
    ! PROGRAMMED BOUNDARY DIRICHLET BOUNDARY CONDITIOS:
    !
    ! 00000 : free, no boundary, (U,rho,E)
    ! 11111 : all prescribed, (U,rho,E)
    ! 50000 : checking inflow/outflow 
    ! 20000 : velocity local base prescribed, (U,rho,E)
    ! 20100 : velocity local base prescribed, (U,rho,E) "the famous condition 20100"
    ! 11101 : no slip NS, (U,rho,T) or subsonic inflow 
    ! 11103 : no slip NS, (U,rho,T) 
    ! 11100 : no slip NS, (U,rho,T) but no prescription on T 
    ! 00010 : subsonic outflow
    ! 11110 : dani 1
    ! 11102 : subsonic inflow, prescribing p instead of T (u,rho,p)
    ! 00002 : subsonic outflow, prescribing p instead of rho (u,rho,p)
    !

    iffix= 0
    do inode=1,pnode
       ipoin= lnods(inode,ielem)
       do idofn=1,ndofn_nsa
          kfl_elfix(idofn,inode)= kfl_fixno_nsa(idofn,ipoin)
          iffix= iffix + kfl_elfix(idofn,inode)
          ! in elmoperations, both explicit and implicit use the delta-form 
          elbve(idofn,inode)= 0.0_rp  
       end do
    end do
    
    if (iffix == 0) return

    do inode = 1,pnode

       ipoin = lnode(inode)
       ibopo = lpoty(ipoin)

       elphy(1:ndime,inode) = gath(inode)%elvel(1:ndime)
       elphy(ndime+1,inode) = gath(inode)%elpre
       elphy(ndime+2,inode) = gath(inode)%eltem

       !
       ! Check if this node has some kind of boundary condition
       !
       nofix = 0
       do idofn= 1,ndofn_nsa
          if (kfl_elfix(idofn,inode) .gt. 0) then
             nofix= nofix + 1
          end if
       end do
       nofix_vel= 0
       do idime= 1,ndime
          if (kfl_elfix(idime,inode) .gt. 0) then
             nofix_vel= nofix_vel + 1
          end if
       end do

       if (nofix > 0 .and. ibopo > 0) then 

          ! this is a boundary node with a condition

          rgasc = runiv_nsa / gath(inode)%elwme
          xhecv = gath(inode)%elhcp - rgasc                      !Cv is computed from R & Cp
          rgacv = runiv_nsa / gath(inode)%elwme / xhecv

          xdens= gath(inode)%elunk(ndime+1,ITER_K)
          xener= gath(inode)%elunk(ndime+2,ITER_K)
          xpres= elphy(ndime+1,inode)
          xtemp= elphy(ndime+2,inode)
          velsq= 0.0_rp
          do idime= 1,ndime
             xmome(idime) = gath(inode)%elunk(idime,ITER_K)
             xadve(idime) = elphy(idime,inode)
             xveln(idime) = elphy(idime,inode)
             velsq = velsq + xveln(idime)*xveln(idime)
          end do
          ! When coupled with alefor, substract the mesh velocity velom to the advection velocity xadve
          if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
             do idime=1,ndime
                xadve(idime) = xadve(idime) - velom(idime,ipoin)
             end do
          end if


          if(kfl_elfix(1,inode)==5) then         
             !
             ! 1. Check wether inflow or outflow, sub or supersonic
             !

             if( solve(1) % kfl_iffix .ne. 0 ) &
                  call runend("NSA_ASSBOU: BC 50000 NOT PREPARED FOR ZEROFIX, TO BE PROGRAMMED.")

             call nsa_chkinf(&
                  kinfl,xmach_chkinf,ipoin,xadve(1:ndime),xveln(1:ndime),xpres,xdens)

             if (kinfl==1) then
                ! inflow
                do idime=1,ndime
                   kfl_elfix(idime,inode)= 1
                end do
                kfl_elfix(ndime+2,inode)= 1
                if (xmach_chkinf >= 1.0_rp) then 
                   ! supersonic, so fix also the density
                   kfl_elfix(ndime+1,inode)  = 1
                end if
             else
                ! outflow
                if (xmach_chkinf < 1.0_rp) then 
                   ! subsonic, so fix only the density
                   kfl_elfix(ndime+1,inode)  = 1
                end if
             end if

          end if

          !
          ! Change variables from conservative to physical 
          !
          ! 1. Compute base change jacobians dU/dQ and dQ/dU
          !
          ! These are the boundary conditions options so far:
          !
          ! (a) u,rho,T (so we use U,rho,E as unknowns) : no jacobians computed -> nofix_vel  > 0 and nofix == ndofn_nsa
          ! (b) u, T    : jacobians computed    -> nofix_vel  > 0 and nofix  < ndofn_nsa
          ! (c) rho     : no jacobians computed -> nofix_vel == 0 and nofix  < ndofn_nsa
          ! (d) u       : jacobians computed, but only for the velocity local frame -> nofix_vel > 0 and nofix < ndofn_nsa,
          !                                                                            but kfl_elfix(ndime+1,inode) == 0
          ! (e) u, rho  : no jacobians computed    -> nofix_vel  > 0 and nofix  < ndofn_nsa
          ! (f) u, p    : jacobians computed       -> nofix_vel  > 0 and nofix  < ndofn_nsa
          ! (g) p       : jacobians computed       -> nofix_vel == 0 and nofix  < ndofn_nsa
          
          if_jacobians= 0 ! (a), (c), (e)

          if (nofix < ndofn_nsa .and. kfl_elfix(ndime+1,inode) == 0) if_jacobians= 1  !(b) and (d)
          if (kfl_elfix(ndime+2,inode) == 2) if_jacobians= 1                          !(f)

          
          if (if_jacobians == 1) then
             ! (b), (d) and (f)
             ! Initialize rotation matrices. As this is done over ibopos but within an element loop,
             ! the rotation matrices are computed and stored more than once, redundantly. 
             ! Pero bueno, la vida es asín. 
             do idofn= 1,ndofn_nsa
                do jdofn= 1,ndofn_nsa
                   jacrot_dq_du_nsa(idofn,jdofn,ibopo)=0.0_rp
                   jacrot_du_dq_nsa(idofn,jdofn,ibopo)=0.0_rp
                end do
                jacrot_dq_du_nsa(idofn,idofn,ibopo)=1.0_rp
                jacrot_du_dq_nsa(idofn,idofn,ibopo)=1.0_rp              
             end do

             ! Compute rotqu and rotuq only when one or more dof are not prescribed, unless only density is prescribed
             !           if (nofix_vel > 0) then

             if (kfl_elfix(ndime+2,inode) == 2) then   ! U=(U, rho, E) and Q=(u, rho, p)
                
                do idime=1,ndime
                   jacrot_du_dq_nsa(  idime,  idime,ibopo) =   xdens         ! dU_i / du_i 
                   jacrot_du_dq_nsa(  idime,ndime+1,ibopo) =   xveln(idime)  ! dU_i / drho
                   jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =   xmome(idime)  ! dE / du_i   
                   
                   jacrot_dq_du_nsa(  idime,  idime,ibopo) =   1.0_rp / xdens            ! du / dU
                   jacrot_dq_du_nsa(  idime,ndime+1,ibopo) = - xveln(idime) / xdens      ! du / drho           
                   jacrot_dq_du_nsa(ndime+2,  idime,ibopo) = - xveln(idime) * rgacv      ! dp / dU_i                 
                end do
                
                jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = 0.5_rp * velsq                 ! dE / drho
                jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = 1.0 / rgacv                    ! dE / dp
                jacrot_du_dq_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                         ! drho / drho
                
                jacrot_dq_du_nsa(ndime+2,ndime+1,ibopo) = rgacv * 0.5_rp * velsq         ! dp / drho
                jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = rgacv                          ! dp / dE
                jacrot_dq_du_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                         ! drho / drho 
                
             else if (kfl_elfix(ndime+1,inode) == 0) then    ! U=(U, rho, E) and Q=(u, rho, T)

                ! it is ELSE because density and pressure cannot be prescribed simultaneously
   
                do idime=1,ndime
                   jacrot_du_dq_nsa(  idime,  idime,ibopo) =   xdens                       ! dU_i / du_i = d (rho u_i) / du_i 
                   jacrot_du_dq_nsa(  idime,ndime+1,ibopo) =   xveln(idime)                ! dU_i / drho = d (rho u_i) / drho 
                   jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =   xmome(idime)                ! dE / du_i   = d (rho c_v T + rho 0.5 u^2) / du_i     
                   
                   jacrot_dq_du_nsa(  idime,  idime,ibopo) =   1.0_rp / xdens
                   jacrot_dq_du_nsa(  idime,ndime+1,ibopo) = - xveln(idime)/ xdens                   
                   jacrot_dq_du_nsa(ndime+2,  idime,ibopo) = - xveln(idime) / xhecv / xdens 
                   
                end do
                jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = xtemp * xhecv + 0.5_rp * velsq    ! dE / drho = d (rho c_v T + rho 0.5 u^2) / drho
                jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = xdens * xhecv                     ! dE / dT   = d (rho c_v T + rho 0.5 u^2) / dT
                jacrot_du_dq_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                            ! drho / drho 
                
                jacrot_dq_du_nsa(ndime+2,ndime+1,ibopo) = (velsq - xener / xdens) / xhecv / xdens
                jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = 1.0_rp / xhecv / xdens 
                jacrot_dq_du_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                            ! drho / drho 
                
                !           end if
                
                
                !           if( kfl_elfix(ndime+2,inode) > 0 ) then                       ! (b), temperature prescribed             
                !              do idime=1,ndime
                !                 jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =    xmome(idime) / xdens
                !                 jacrot_dq_du_nsa(  idime,ndime+2,ibopo) =   -xmome(idime) / xhecv / xdens / xdens
                !              end do
                !              jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = xtemp * xhecv - 0.5_rp * velsq
                !              jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = xdens * xhecv
                !              jacrot_dq_du_nsa(ndime+1,ndime+2,ibopo) = (velsq - xener / xdens) / xhecv / xdens 
                !              jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = xener /  xhecv / xdens 
                !           end if

             end if
             !
             ! 2. Check the specific prescription
             !
             if( kfl_elfix(1,inode) == 2 ) then                             ! compute local frame
                !
                ! 1. First, check local base velocity prescriptions
                !
                !   rotau(1:ndime,1:ndime,1)  is  GL , goes in rotuq
                !   rotau(1:ndime,1:ndime,2)  is  LG , goes in rotqu
                !
                !   Choose the proper local basis
                !
                iroty=kfl_fixrs_nsa(ibopo)
                if( iroty == -1 ) then                                    ! Tangent system
                   do idime=1,ndime
                      do jdime= 1,ndime
                         rotlg(idime,jdime)= exnor(idime,jdime,ibopo)
                         rotgl(idime,jdime)= exnor(jdime,idime,ibopo)
                      end do
                   end do

                else if( iroty >= 1 ) then                                ! Given system
                   do idime=1,ndime
                      do jdime= 1,ndime
                         rotlg(idime,jdime)= skcos(idime,jdime,iroty)
                         rotgl(idime,jdime)= skcos(jdime,idime,iroty)
                      end do
                   end do
                   !           else if( iroty == -2 ) then                               ! Given system
                else if( iroty == -3 ) then                               ! Geometrical normal
                   do idime=1,ndime
                      do jdime= 1,ndime
                         rotlg(idime,jdime)= skcos(idime,jdime,ibopo)
                         rotgl(idime,jdime)= skcos(jdime,idime,ibopo)
                      end do
                   end do
                end if


                !
                ! 2. Correct jacrot_du_dq_nsa and jacrot_dq_du_nsa accordingly: 
                !    F   = dU/dQ *   LG   in jacrot_du_dq_nsa(...)  
                !    F^⁻1=  GL  * dQ/dU   in jacrot_dq_du_nsa(...)

                rotqu_aux= 0.0_rp
                rotuq_aux= 0.0_rp
                do idime=1,ndime
                   do idofn= 1,ndofn_nsa
                      do jdime=1,ndime
                         rotqu_aux(idofn,idime)= &
                              rotqu_aux(idofn,idime) + jacrot_du_dq_nsa(idofn,jdime,ibopo) * rotlg(jdime,idime) 
                         rotuq_aux(idime,idofn)= &
                              rotuq_aux(idime,idofn) + rotgl(idime,jdime) * jacrot_dq_du_nsa(jdime,idofn,ibopo)
                      end do
                   end do
                end do

                do idime=1,ndime
                   do idofn= 1,ndofn_nsa
                      jacrot_du_dq_nsa(idofn,idime,ibopo)=rotqu_aux(idofn,idime) 
                      jacrot_dq_du_nsa(idime,idofn,ibopo)=rotuq_aux(idime,idofn) 
                   end do
                end do

             end if

             !
             ! Rotate elmat and elrhs
             !
             
             call nsa_rotsys(1_ip,&
                  inode,pnode,ndofn_nsa,pevat,elmat,elrhs,&
                  jacrot_du_dq_nsa(1,1,ibopo),jacrot_dq_du_nsa(1,1,ibopo),kfl_linea_nsa,kfl_matvec)

             !           if( kfl_elfix(1,inode) == 2) then
             !              matri_debu=0.0_rp
             !              matri_debu_ndime=0.0_rp
             !              matri_debu=matmul(jacrot_du_dq_nsa(:,:,ibopo),jacrot_dq_du_nsa(:,:,ibopo))
             !              matri_debu_ndime=matmul(rotgl(:,:),rotlg(:,:))
             !              matri_debu_ndime=0.0_rp
             !           end if

          end if


          if( solve(1) % kfl_iffix == 0 ) then


             !
             ! Correct boundary conditions
             !
             !
             !matrix and rhsid modified now to account for the boundary conditions
             !

             ! Momentum dof's
             do idime = 1,ndime
                if( kfl_elfix(idime,inode) == 1 .or. kfl_elfix(idime,inode) == 2) then
                   ievat = (inode-1)*ndofn_nsa + idime
                   
                   if (kfl_matvec == 2) then   ! both matrix and vector (implicit schemes)
                      ! either local frame or cartesian frame, both are equally treated
                      adiag = elmat(ievat,ievat)
                      do jnode = 1,pnode 
                         do jdofn = 1,ndofn_nsa
                            jevat = (jnode-1)*ndofn_nsa + jdofn
                            elmat(ievat,jevat) = 0.0_rp
                            elmat(jevat,ievat) = 0.0_rp
                         end do
                      end do
                      elmat(ievat,ievat)        = adiag
                   end if

                   elrhs(ievat)              = 0.0_rp
                   !              if (kfl_linea_nsa == 2) elsou(ievat) = adiag * xvalu
                end if
             end do


             ! Continuity dof
             if( kfl_elfix(ndime+1,inode) == 1 ) then
                ievat = (inode-1)*ndofn_nsa + ndime + 1
                
                if (kfl_matvec == 2) then   ! both matrix and vector (implicit schemes)
                   adiag = elmat(ievat,ievat)
                   do jnode = 1,pnode 
                      do jdofn = 1,ndofn_nsa
                         jevat = (jnode-1)*ndofn_nsa + jdofn
                         elmat(ievat,jevat) = 0.0_rp
                         elmat(jevat,ievat) = 0.0_rp
                      end do
                   end do
                   elmat(ievat,ievat)        = adiag
                end if
                elrhs(ievat)              = 0.0_rp
                !           if (kfl_linea_nsa == 2) elsou(ievat) = adiag * xvalu
             end if

             ! Energy or Temperature
             if( kfl_elfix(ndime+2,inode) > 0) then
                ievat = (inode-1)*ndofn_nsa + ndime + 2

                if (kfl_matvec == 2) then   ! both matrix and vector (implicit schemes)
                   adiag = elmat(ievat,ievat)
                   do jnode = 1,pnode 
                      do jdofn = 1,ndofn_nsa
                         jevat = (jnode-1)*ndofn_nsa + jdofn
                         elmat(ievat,jevat) = 0.0_rp
                         elmat(jevat,ievat) = 0.0_rp
                      end do
                   end do
                   elmat(ievat,ievat)        = adiag
                end if

                elrhs(ievat)              = 0.0_rp
             end if

          end if

       end if   ! nofix

    end do

  end subroutine nsa_newelmatrixsetboundary


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> @addtogroup NastalElmoperations
  !> @{
  !> @file    mod_nsa_newelmoperations.f90
  !> @author  Mariano Vazquez
  !> @date    11/06/2015
  !> @brief   Scatter subscale to the elements
  !> @details Scatter subscale to the elements
  !> @} 
  !-----------------------------------------------------------------------
  subroutine nsa_newelmscattersubscale(inode,igaus,ielem,pnode,gaus,shapigaus)
    use      def_master
    use      def_nastal
    use      def_domain

    implicit none

    integer(ip) :: inode
    integer(ip) :: igaus
    integer(ip) :: ielem
    integer(ip) :: pnode
    type(elm_ongaus_interp_nsa)  :: gaus(mgaus)
    real(rp)    :: shapigaus(pnode)

    integer(ip) :: ipoin
    real(rp)    :: asfac
    
    ipoin= lnods(inode,ielem)
    asfac= gaus(igaus)%dvolu * shapigaus(inode) / vmass(ipoin)
    umoss_nsa(    1,ipoin,2) = &
         umoss_nsa(    1,ipoin,2) + asfac * umosg_nsa(    1,ielem,igaus,1)
    umoss_nsa(    2,ipoin,2) = &
         umoss_nsa(    2,ipoin,2) + asfac * umosg_nsa(    2,ielem,igaus,1)
    if (ndime == 3) &
         umoss_nsa(ndime,ipoin,2) =&
         umoss_nsa(ndime,ipoin,2) + asfac * umosg_nsa(ndime,ielem,igaus,1)
    denss_nsa(      ipoin,2) = &
         denss_nsa(      ipoin,2) + asfac * densg_nsa(      ielem,igaus,1)
    eness_nsa(      ipoin,2) = &
         eness_nsa(      ipoin,2) + asfac * enesg_nsa(      ielem,igaus,1)

    ! no shock capturing to x-momentum means no shock capturing at all
    if (kfl_shock_nsa(1) > 0) then
       shocktau_nsa(ielem)%a(1,igaus,1) = gaus(igaus)%shocktau_local(1)
       shocktau_nsa(ielem)%a(2,igaus,1) = gaus(igaus)%shocktau_local(ndime+1)
       shocktau_nsa(ielem)%a(3,igaus,1) = gaus(igaus)%shocktau_local(ndime+2)
    end if


  end subroutine nsa_newelmscattersubscale

end module mod_nsa_newelmoperations
