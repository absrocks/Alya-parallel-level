!==============================================================================!
!
!< Dic08, 2014. (0) created  
!< Dic09, 2014. (1) convective matrix -> characteristics 
!
! Dont put all your eggs in one basket !!
!
! Patm = 103400.0            !> [Nm2/s2/m3]
! Eint = 5.0/2.0*287.0*300.0 !> [J/kg] = [J/kg/K][K], Eint = 215250.0 (gas ideal)
!
!==============================================================================!
module mod_nsa_euler
  use def_kermod, only: kfl_prope
  use def_parame, only: ip, rp  
  use def_domain, only: ltype, lnods, coord, hnatu
  use def_domain, only: nnode, ngaus, mnode, mgaus, nelem, ndime, npoin
  use def_master, only: inotmaster
  use def_domain, only: vmasc, vmass, elmar
  implicit none
  !
  logical(ip)           :: cons
  integer(ip)           :: ndofn, nevat 
  integer(ip)           :: dt_inv 
  real(rp), pointer     :: gamma_nsa(:)
  !
  type LODI
!   integer(ip)           :: idofn           =  1 
    integer(ip)           :: ndofn           =  0 
    real(rp)              :: gamme           =  1.2_rp
    !
    real(rp), pointer     :: xdvol(:,:)      => null()
    real(rp), pointer     :: xchrc(:,:,:,:)  => null()
    real(rp), pointer     ::  chrc(:,:,:)    => null()
    !
    real(rp), pointer     ::  densi(:)       => null()
    real(rp), pointer     ::  veloc(:,:)     => null()
    real(rp), pointer     ::  press(:)       => null()
  end type 
  type(LODI), save :: EULER
  !
  private
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
    public :: EULER
    public :: nsa_euler_allocate
    public :: nsa_euler_deallocate
    public :: nsa_euler_calc
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
contains 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  ! +
  ! |_Xxxxxx
  !   |_xxx_turnon
  !
  subroutine nsa_euler_allocate( STRUCT )
  use def_master, only: veloc, densi, press
  use def_domain, only: mnode
  use def_nastal, only: adgam_nsa
  implicit none 
  type(LODI),  intent(inout) :: STRUCT
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
      !
      cons         = .false. 
      ndofn        =  ndime+2
      nevat        =  ndofn * mnode
      !
      STRUCT%ndofn =  ndofn
      STRUCT%gamme =  adgam_nsa
      !
      STRUCT%densi => densi(        1:npoin,1)
      STRUCT%veloc => veloc(1:ndime,1:npoin,1)
      STRUCT%press => press(        1:npoin,1)
      !
      allocate( STRUCT%xchrc(ndofn, nelem, mgaus, ndime) )
      !call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xchrc_nsa)
      !
      allocate( STRUCT%chrc(ndofn, npoin, ndime) )
      !call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', chrc_nsa)
      !
      allocate( STRUCT%xdvol(nelem, mgaus) )

      STRUCT%xchrc = 0.0_rp
      STRUCT% chrc = 0.0_rp
      STRUCT%xdvol = 0.0_rp
  endif 
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsa_euler_deallocate( STRUCT )
  implicit none 
  type(LODI),  intent(inout) :: STRUCT
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    deallocate( STRUCT%xchrc )
    deallocate( STRUCT% chrc )
    deallocate( STRUCT%xdvol )
  endif 
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine nsa_euler_calc( STRUCT )
  use def_master, only: rhsid
  use def_nastal, only: dtinv_nsa
  implicit none 
  type(LODI), intent(inout) :: STRUCT
  !-----------------------------------------------------------------------||---!
  real(rp)    :: hmini
  real(rp)    :: elcod(ndime, mnode)
  real(rp)    :: detjm, xshap(mgaus) 
  real(rp)    :: xjaci(ndime, ndime), xjacm(ndime, ndime), xcartd(ndime, mnode, mgaus)
  real(rp)    :: tragl(ndime, ndime), hleng(ndime) 
  !
  real(rp)    :: xdphit(ndofn, ndime, mgaus)  
  real(rp)    ::   phit(ndofn, mnode), xphit(ndofn, mgaus)
  real(rp)    ::  xvelmo(mgaus), xsound(mgaus) 
  real(rp)    :: Sxinv(1:ndofn,1:ndofn)
  real(rp)    :: elgamme(mnode)
  real(rp)    ::  xgamme(mgaus)
  !
  integer(ip) :: ielem, inode, pelty, pgaus
  integer(ip) :: pnode, igaus, idofn, idime, ipoin
  integer(ip) :: itime, jtime 
  logical(ip) :: physc 
  !
  !< subscale 
  real(rp)    ::   hinv
  real(rp)    ::   xres(ndofn, mgaus)
  real(rp)    :: ssxphi(ndofn, mgaus)
  real(rp)    ::   xtau(ndofn, mgaus)
  !
  !< shock
  real(rp)    :: xdifeq(ndofn) 
  real(rp)    :: xkapsh(ndofn,mgaus)
  real(rp)    :: xshmet(ndime, ndime, ndofn, mgaus) 
  real(rp)    :: xshote(ndofn, mgaus)
  !
  !< convection
  real(rp)    ::  elconv(ndofn, ndofn, ndime, mnode)
  real(rp)    ::   xconv(ndofn, ndofn, ndime, mgaus)
  real(rp)    ::  xdconv(ndofn, ndofn, mgaus)
  real(rp)    ::   xchrc(ndofn, ndime, mgaus)
  real(rp)    ::   xSkL(ndofn, ndime, mgaus)
  !
  !< solver 
  real(rp)    ::  elmat(nevat,nevat)
  real(rp)    ::  elrhs(nevat)
  real(rp)    :: elunks(nevat)
  integer(ip) :: idofg, idofl, ievat

  !-----------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------------| inisol |---!
  itime  = 1_ip
  jtime  = 2_ip
  physc  = .not.cons 
  dt_inv =  dtinv_nsa
  !
  call inisol() 
  !
  !-----------------------------------------------------------------| main |---!
  master:& 
  if(INOTMASTER) then
  !-----------------------------------------------------------------------||---!
  elementary_loop: & 
  do ielem = 1,nelem 
    !---------------------------------------------------------------------||---!
    pelty = ltype(ielem)
    pnode = nnode(pelty)
    pgaus = ngaus(pelty)
    !----------------------------------------------| PHIt -> [rho, vel, p] |---!
    node_loop01:&
    do inode = 1,pnode
      !-------------------------------------------------------------------||---!
      ipoin                  = lnods(inode, ielem) 
      elcod(1:ndime,inode)   = coord(1:ndime,ipoin)
      !-------------------------------------------------------------------||---!
      elgamme(inode)         = STRUCT%gamme
      phit(1,         inode) = STRUCT%densi(         ipoin)  
      phit(2:ndime+1, inode) = STRUCT%veloc(1:ndime, ipoin)
      phit(  ndime+2, inode) = STRUCT%press(         ipoin) 
      !---------------------------------------------------------| tracking |---!
      call set_subscale(physc, jtime, ipoin, elgamme(inode), phit(1:ndofn,inode) )
      !-----------------------------------------------------------| elconv |---!
      call convection_jacobian(physc, elgamme(inode), phit(1:ndofn,inode), &
                               elconv(1:ndofn,1:ndofn,1:ndime,inode) ) 
      !-------------------------------------------------------------------||---!
    enddo node_loop01 
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    call elmlen(ndime, pnode, elmar(pelty)%dercg, tragl, elcod, hnatu(pelty), hleng)
    hmini = hleng(ndime) !> the smallest
    hinv  = 1.0_rp/hleng(ndime)
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    gauss_loop01: &
    do igaus = 1,pgaus
      !-------------------------------------------------------------------||---!
      call elmder(pnode, ndime, elmar(pelty)%deriv(1,1,igaus), elcod, xcartd(1,1,igaus), detjm, xjacm, xjaci)
      STRUCT%xdvol(ielem,igaus) = elmar(pelty)%weigp(igaus) * detjm 
      !-------------------------------------------------------------------||---!
      !
      !------------------------------------------------------------| dPHIt |---!
      xdphit(1:ndofn,1:ndime,igaus) = 0.0_rp
      do idofn = 1,ndofn
        do idime = 1,ndime
        xdphit(idofn,idime,igaus) = dot_product( phit(idofn,1:pnode), xcartd(idime,1:pnode,igaus)       )
        enddo
      enddo
      !-------------------------------------------------------------| PHIt |---!
      do idofn = 1,ndofn
        xphit(idofn,       igaus) = dot_product( phit(idofn,1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      enddo
      !-----------------------------------------------------------| xgamme |---!
      xgamme(igaus) = dot_product( elgamme(1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      !---------------------------------------------------------| tracking |---! 
      call set_xsubscale(physc, itime, ielem, igaus, xgamme(igaus), xphit(1:ndofn,igaus))
      !---------------------------------------------------| xsound, xvelmo |---!
      call psound(xphit(1:ndofn,igaus), xgamme(igaus), xsound(igaus), xvelmo(igaus))
      !-------------------------------------------------------------------||---!
      !
      !-----------------------------| L = (\Lambda^k * S_k^{-1}) * dU/dr_k |---!
      !do idime = 1,ndime
      ! call SINVMTRX( Sxinv(1:ndofn,1:ndofn), xphit(1:ndofn,igaus), idime, xsound(igaus) )
      ! STRUCT%xchrc(1:ndofn,ielem,igaus,idime) = matmul( Sxinv(1:ndofn,1:ndofn), xdphit(1:ndofn,idime,igaus) )
      !enddo
      !-------------------------------------------------------------------||---!
      !
      !----------------------------------------------------| xconv, xdconv |---! 
      call xconvection_jacobian(physc, pnode,                            & 
                                xgamme(                        igaus),   &
                                xsound(                        igaus),   &
                                xdphit(1:ndofn,1:ndime,        igaus),   &
                                 xphit(1:ndofn,                igaus),   &
                                xcartd(1:ndime,1:mnode,        igaus),   & 
                                elconv(1:ndofn,1:ndofn,1:ndime,1:mnode), &
                                 xconv(1:ndofn,1:ndofn,1:ndime,igaus),   & 
                                xdconv(1:ndofn,1:ndofn,        igaus),   & 
                                 xchrc(1:ndofn,1:ndime,        igaus),   &
                                  xSkL(1:ndofn,1:ndime,        igaus)    )
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------------------| xres |---! 
      xres(1:ndofn,igaus) = 0.0  
      do idime = 1,ndime
        xres(1:ndofn,igaus) = xres(1:ndofn,igaus) + &
          matmul( -xconv(1:ndofn,1:ndofn,idime,igaus), xdphit(1:ndofn,idime,igaus) ) 
      enddo
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------| Shock capturing  |---!
      xdifeq = 0.0
      xshmet(1:ndime,1:ndime,1:ndofn,igaus) = 0.0
      xkapsh(                1:ndofn,igaus) = 0.0
      call init_shock( xkapsh(1:ndofn,        igaus),  &
                       xdphit(1:ndofn,1:ndime,igaus),  &
                         xres(1:ndofn,        igaus),  & 
                                              xdifeq,  &
                                       xvelmo(igaus),  &
                                               hmini,  &  
               xshmet(1:ndime,1:ndime,1:ndofn,igaus))
      !
      xshote = 0.0 
      call calc_shock(         xkapsh(1:ndofn,igaus), & 
                       xdphit(1:ndofn,1:ndime,igaus), & 
                         xcartd(1:ndime,inode,igaus), & 
               xshmet(1:ndime,1:ndime,1:ndofn,igaus), &
                               xshote(1:ndofn,igaus) )
      !-------------------------------------------------------------------||---!
      !
      !---------------------------------------------------------| tracking |---! 
      call get_xsubscale(physc, itime, ielem, igaus,  &
                                       xgamme(igaus), & 
                                 xres(1:ndofn,igaus), &
                                       xvelmo(igaus), & 
                                       xsound(igaus), &
                                         hinv, pnode, &
                               ssxphi(1:ndofn,igaus), &
                                 xtau(1:ndofn,igaus) )
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------------------------||---!
    enddo gauss_loop01
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    call euler_calc_matrices(physc, ielem, & 
                                              xtau(1:ndofn,1:mgaus), &
                                            xshote(1:ndofn,1:mgaus), &
                                        STRUCT%xdvol(ielem,1:mgaus), &
                                    xcartd(1:ndime,1:mnode,1:mgaus), &
                             xconv(1:ndofn,1:ndofn,1:ndime,1:mgaus), &  
                                      xSkL(1:ndofn,1:ndime,1:mgaus), &
                                             elmat(1:nevat,1:nevat), & 
                                             elrhs(1:nevat        )  ) 
    !---------------------------------------------------------------------||---!
    !
    !-----------------------------------------------------------| Assembly |---!
    if(.true.) then
      do inode = 1,pnode
        do idofn = 1,ndofn
          ievat = (inode-1) * ndofn + idofn
          elunks(ievat) = phit(idofn,inode) ! -> [RHO, VELi, PRE] 
        enddo
      enddo
      !
      elrhs = matmul(elmat,elunks) - elrhs 
    endif 
    !-----------------------------------------------| kernel/mathru/assrhs |---!
    do inode = 1,pnode
      ipoin = lnods(inode,ielem)
      idofg = (ipoin-1)*ndofn
      idofl = (inode-1)*ndofn
      rhsid(idofg+1:idofg+ndofn) = rhsid(idofg+1:idofg+ndofn) + elrhs(idofl+1:idofl+ndofn)
    enddo
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end do elementary_loop
  !-----------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------| gauss2nodes |---<! 
  !do idofn = 1,ndofn
  !  do idime = 1,ndime
  !  call gauss2nodes( STRUCT%xchrc(idofn,1:nelem,1:mgaus,idime), STRUCT%xdvol(1:nelem,1:mgaus), STRUCT%chrc(idofn,1:npoin,idime) ) 
  !  enddo
  !enddo
  !-----------------------------------------------------------------------||---!
  !
  !-------------------------------------------------------------| tracking |---!
  call get_subscale(itime, jtime, STRUCT%xdvol(1:nelem,1:mgaus))
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  call phit2unkno(physc, jtime, STRUCT%gamme)
  call    solving(physc, STRUCT%gamme)
  !call nsa_upcons
  !if(.not.kfl_track_nsa==0) call runend("ERROR lodi: unavailable tracking stabilization (set STABILIZATION:MULTISCALE)!")
  !call inisol()


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  endif master
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------------| upcons |---!
  subroutine phit2unkno(phys, i_time, gamme)
  use def_master, only: unkno
  use def_master, only: veloc, umome, densi, energ, densi, press
  implicit none
  logical(ip), intent(in)  :: phys
  integer(ip), intent(in)  :: i_time 
  real(rp), intent(in)     :: gamme
  !
  integer(ip) :: i_poin, i_evat, i_dime !, i_time 
  !i_time = 2_ip
  !---------------------------------------------------------------------||---!
  !
  !---------------------------------------------------------------------||---!
  if(INOTMASTER) then
    do i_poin = 1,npoin
      i_evat  = (i_poin-1)*ndofn
      do i_dime = 1,ndime
      unkno(i_evat+i_dime ) = umome(i_dime,i_poin,i_time) 
      end do
      unkno(i_evat+ndime+1) = densi(       i_poin,i_time)
      unkno(i_evat+ndime+2) = energ(       i_poin,i_time) 
      if(phys) call cons2phys( unkno(i_evat+1:i_evat+ndofn), gamme, unkno(i_evat+1:i_evat+ndofn) )
    enddo
    !print *, "01.rho:", minval(densi(1:npoin,1)), maxval(densi(1:npoin,1)), sum(densi(1:npoin,1))/npoin   
    !print *, "01.pre:", minval(press(1:npoin,1)), maxval(press(1:npoin,1)), sum(press(1:npoin,1))/npoin    
    !print *, ""
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------------------|  |---!
  subroutine solving(phys, gamme)
  use def_master, only: unkno, rhsid,  solve, solve_sol, npoin_type
  implicit none
  logical(ip), intent(in)  :: phys
  real(rp), intent(in)     :: gamme
  !
  integer(ip) :: n_dofn, i_totn, i_dofn, i_evat
  integer(ip) :: i_poin
  real(rp)    :: diag, fact1
  !-----------------------------------------------------------------------||---!
  !
  !------------------------------| Solving DU/Dt -> U(n+1)-U(n) = RHS * Dt |---!
  if(INOTMASTER) then
    solve(1) % xdiag = 1.0_rp/dt_inv
    if(solve_sol(1) % kfl_recov /= 2) call pararr('SLX', NPOIN_TYPE, solve_sol(1) % ndofn * npoin, rhsid)
    !---------------------------------------------------------------------||---!
    n_dofn = solve_sol(1) % ndofn
    diag   = solve_sol(1) % xdiag
    do i_poin = 1,npoin
      i_totn = (i_poin-1) *n_dofn
      fact1  = diag/vmasc(i_poin)
      do i_dofn = 1,n_dofn
        unkno(i_totn+i_dofn) = unkno(i_totn+i_dofn) + rhsid(i_totn+i_dofn) * fact1
      end do
    enddo
    !---------------------------------------------------------------------||---!
    if(phys) then 
      do i_poin = 1,npoin
        i_evat  = (i_poin-1)*ndofn
        call phys2cons( unkno(i_evat+1:i_evat+ndofn), gamme, unkno(i_evat+1:i_evat+ndofn)) 
      enddo 
    endif 
    !---------------------------------------------------------------------||---!
    call nsa_setvar(3_ip,1_ip)
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !---------------------------------------------------------| SEMI-PRIVATE |---!
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine euler_calc_matrices(phys, ielem, xtau, xshote, xdvol, xcartd, xconv, xSkL, mtx, rhs)
  use def_nastal, only: nevat_nsa
  implicit none
  logical(ip), intent(in)   :: phys
  integer(ip), intent(in)   :: ielem 
  !
  real(rp),    intent(in )  ::   xtau(ndofn, mgaus)
  real(rp),    intent(in )  :: xshote(ndofn, mgaus)
  real(rp),    intent(in )  ::  xdvol(mgaus)
  real(rp),    intent(in )  :: xcartd(ndime,mnode,mgaus)
 !real(rp),    intent(in )  ::
 !real(rp),    intent(in )  ::
 !real(rp),    intent(in )  ::
  real(rp),    intent(in )  ::  xconv(ndofn,ndofn,ndime,mgaus)
  real(rp),    intent(in )  ::   xSkL(ndofn,      ndime,mgaus)
  real(rp),    intent(out)  :: mtx(nevat_nsa, nevat_nsa)
  real(rp),    intent(out)  :: rhs(nevat_nsa)
  !
  real(rp)    ::  xmtx01(ndofn,ndofn)
  real(rp)    ::  xmtx02(ndofn,ndofn)
  !
  real(rp)    ::  xadvec_matrix(ndofn, ndofn, mnode)
  real(rp)    ::  xsubdi_matrix(ndofn, ndofn, mnode)
  real(rp)    ::  xstabi_matrix(ndofn, ndofn, mnode)
  real(rp)    ::  xtimas_matrix(ndofn, ndofn, mnode)
  real(rp)    ::  xadjoi_matrix(ndofn, ndofn)
  real(rp)    :: xsource(ndofn)
  !
  integer(ip) :: pelty, pnode, pgaus
  integer(ip) :: ipoin, idime, igaus, inode, jnode, jdofn
  integer(ip) :: icols(ndofn), irows(ndofn) 
  integer(ip) :: kk
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  mtx(1:nevat_nsa,1:nevat_nsa) = 0.0_rp
  rhs(1:nevat_nsa            ) = 0.0_rp
  !
  pelty = ltype(ielem)
  pnode = nnode(pelty)
  pgaus = ngaus(pelty)
  do igaus = 1,pgaus
    do inode = 1,pnode
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------------------------||---!
      xtimas_matrix(1:ndofn,1:ndofn,1:pnode) = 0.0_rp 
      xadvec_matrix(1:ndofn,1:ndofn,1:pnode) = 0.0_rp
      xstabi_matrix(1:ndofn,1:ndofn,1:pnode) = 0.0_rp ! stabi = subdi * adjoi 
      xsubdi_matrix(1:ndofn,1:ndofn,1:pnode) = 0.0_rp
      xadjoi_matrix(1:ndofn,1:ndofn        ) = 0.0_rp 
      xsource(1:ndofn) = 0.0_rp
      !
      do idime = 1,ndime 
        xmtx01(1:ndofn,1:ndofn)  = xconv(1:ndofn,1:ndofn,idime,igaus) * elmar(pelty) % shape(inode,igaus)
        !
        do jdofn = 1,ndofn
        xmtx02(1:ndofn,  jdofn)  = xconv(1:ndofn,  jdofn,idime,igaus) * xtau(jdofn, igaus)
        enddo
        !
        do jnode = 1,pnode
        xtimas_matrix(:,:,jnode) = elmar(pelty)%shape(inode,igaus) * elmar(pelty)%shape(jnode,igaus) * dt_inv 
       !xadvec_matrix(:,:,jnode) = xadvec_matrix(:, :,jnode) + xcartd(idime,jnode,igaus) * xmtx01(1:ndofn,1:ndofn) !**
        xsubdi_matrix(:,:,jnode) = xsubdi_matrix(:, :,jnode) - xcartd(idime,jnode,igaus) * xmtx02(1:ndofn,1:ndofn)
        enddo
        !
        xadjoi_matrix(:,:)       = xadjoi_matrix(:,:       ) + xcartd(idime,inode,igaus) * xconv(1:ndofn,1:ndofn,idime,igaus)
        !
        xsource(1:ndofn) = xsource(1:ndofn) + xSkL(1:ndofn,idime,igaus) * elmar(pelty)%shape(inode,igaus) !< (1)
        !
      enddo
      !
      do jnode = 1,pnode
        xstabi_matrix(1:ndofn,1:ndofn,jnode) = matmul( xadjoi_matrix(1:ndofn,1:ndofn), xsubdi_matrix(1:ndofn,1:ndofn,jnode) ) !**
      enddo
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------------------------||---!
      if(.true.) then 
        irows      = (inode-1) * ndofn + (/ (kk, kk=1,ndofn) /)
        rhs(irows) = rhs(irows) + xdvol(igaus) *  xshote(1:ndofn,igaus) !< work ?? 
        rhs(irows) = rhs(irows) + xdvol(igaus) * xsource(1:ndofn)       !< (1) 
        !
        do jnode = 1,pnode
        icols                   =  (jnode-1) * ndofn + (/ (kk, kk=1,ndofn) /)
        xmtx01(1:ndofn,1:ndofn) = -xadvec_matrix(1:ndofn,1:ndofn,jnode) + xstabi_matrix(1:ndofn,1:ndofn,jnode) 
        mtx(irows,icols)        =  mtx(irows,icols) + xmtx01(1:ndofn,1:ndofn) * xdvol(igaus)
        enddo
        ! 
      endif 
      !-------------------------------------------------------------------||---!
      !
      !-------------------------------------------------------------------||---!
    enddo
    ! 
  enddo 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine convection_jacobian(phys, gamme, el_phit, el_conv)
    implicit none
    logical(ip), intent(in)  :: phys
    real(rp), intent( in)    :: gamme
    real(rp), intent( in)    :: el_phit(ndofn)
    real(rp), intent(out)    :: el_conv(ndofn, ndofn, ndime)
    !
    integer(ip) :: i_dime
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    if(phys) then
        do i_dime = 1,ndime
          call PHYMTRX( el_conv(1:ndofn,1:ndofn,i_dime), el_phit(1:ndofn), i_dime, gamme)
        enddo
    else
        do i_dime = 1,ndime
          call   CMTRX( el_conv(1:ndofn,1:ndofn,i_dime), el_phit(1:ndofn), i_dime, gamme)
        enddo
    endif 
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine xconvection_jacobian(phys, p_node, gamme, sound, x_dphit, x_phit, x_cartd, el_conv, x_conv, x_dconv, x_chrc, x_SkL)
    implicit none
    logical(ip), intent(in)  :: phys
    integer(ip), intent(in)  :: p_node
    real(rp), intent( in)    :: gamme
    real(rp), intent( in)    :: sound 
    real(rp), intent( in)    :: x_dphit(ndofn,ndime) 
    real(rp), intent( in)    ::  x_phit(ndofn)
    real(rp), intent( in)    :: x_cartd(ndime, mnode)
    real(rp), intent( in)    :: el_conv(ndofn, ndofn, ndime, mnode)
    real(rp), intent(out)    ::  x_conv(ndofn, ndofn, ndime)
    real(rp), intent(out)    :: x_dconv(ndofn, ndofn)
    real(rp), intent(out)    ::  x_chrc(ndofn, ndime) 
    real(rp), intent(out)    ::   x_SkL(ndofn, ndime)
    !
    real(rp)    :: Sx_inv(ndofn,ndofn      )
    real(rp)    ::     Sx(ndofn,ndofn,ndime) 
    !
    integer(ip) :: i_dofn, j_dofn, i_dime 
    !
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    if(phys) then
      do i_dime = 1,ndime
        !<   L = \Lambda^k * S_k^{-1} * dU/dr_k
        call SINVMTRX( Sx_inv(1:ndofn,1:ndofn       ), x_phit(1:ndofn), i_dime, sound )
        x_chrc(1:ndofn,i_dime) = matmul( Sx_inv(1:ndofn,1:ndofn       ), x_dphit(1:ndofn,i_dime) )
        !< d_k = S_k * L  
        call    SMTRX(     Sx(1:ndofn,1:ndofn,i_dime), x_phit(1:ndofn), i_dime, sound )
!x_chrc(ndofn,i_dime) = -x_chrc(1,i_dime) 
!x_chrc(1,i_dime) = -x_chrc(ndofn,i_dime) 
!x_chrc(ndofn,3) = -x_chrc(1,3) 
!x_chrc(1,3) = -x_chrc(ndofn,3)           
        x_SkL( 1:ndofn,i_dime) = matmul(     Sx(1:ndofn,1:ndofn,i_dime),  x_chrc(1:ndofn,i_dime) ) !< (1)
      enddo  
      do i_dime = 1,ndime
        call PHYMTRX( x_conv(1:ndofn,1:ndofn,i_dime), x_phit(1:ndofn), i_dime, gamme )
      enddo 
    else  
      do i_dime = 1,ndime
        call CMTRX( x_conv(1:ndofn,1:ndofn,i_dime), x_phit(1:ndofn), i_dime, gamme )
      enddo
    endif 
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------| dA_{\alpha}/dr_i |---!
    x_dconv(1:ndofn,1:ndofn) = 0.0
    do i_dofn = 1,ndofn
        do j_dofn = 1,ndofn
          do i_dime = 1,ndime
            x_dconv(i_dofn,j_dofn) = x_dconv(i_dofn,j_dofn) + &
                                     dot_product( el_conv(i_dofn,j_dofn,i_dime,1:p_node), x_cartd(i_dime,1:p_node) )
          enddo
        enddo
    enddo
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine set_subscale( phys, i_time, i_poin, gamme, phi )
    use def_nastal, only: kfl_track_nsa, umoss_nsa, denss_nsa, eness_nsa
    implicit none
    logical(ip), intent(in) :: phys
    integer(ip), intent(in) :: i_poin, i_time
    real(rp), intent(in   ) :: gamme
    real(rp), intent(inout) :: phi(ndofn) 
    !
    real(rp), pointer :: umoss(:,:) => null()
    real(rp), pointer :: denss(  :) => null()
    real(rp), pointer :: eness(  :) => null()
    !
    real(rp) ::  sub(ndofn)
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    sub(1:ndofn) = 0.0
    if(kfl_track_nsa == 1) then
        umoss => umoss_nsa(1:ndime,1:npoin,i_time) 
        denss => denss_nsa(        1:npoin,i_time) 
        eness => eness_nsa(        1:npoin,i_time) 
        !
        sub(1:ndime) = umoss(1:ndime,i_poin)
        sub(ndime+1) = denss(        i_poin)
        sub(ndime+2) = eness(        i_poin)
        if(phys) call cons2phys( sub(1:ndofn), gamme, sub(1:ndofn)) 

        umoss(1:ndime,i_poin) = 0.0
        denss(        i_poin) = 0.0
        eness(        i_poin) = 0.0
    end if
    phi(1:ndofn) = phi(1:ndofn) + sub(1:ndofn)    
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine get_subscale(i_time, j_time, xdvol)
    use def_nastal, only: umoss_nsa, denss_nsa, eness_nsa
    use def_nastal, only: umosg_nsa, densg_nsa, enesg_nsa, kfl_goite_nsa, miinn_nsa
    implicit none
    integer(ip), intent(in) :: i_time, j_time 
    real(rp),    intent(in) :: xdvol(nelem,mgaus)
    !
    integer(ip) :: i_dime
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    do i_dime = 1,ndime 
      call gauss2nodes( umosg_nsa(i_dime,1:nelem,1:mgaus,i_time), xdvol(1:nelem,1:mgaus), umoss_nsa(i_dime,1:npoin,j_time) )
    enddo 

    call gauss2nodes( densg_nsa(1:nelem,1:mgaus,i_time), xdvol(1:nelem,1:mgaus), denss_nsa(1:npoin,j_time) )
    call gauss2nodes( enesg_nsa(1:nelem,1:mgaus,i_time), xdvol(1:nelem,1:mgaus), eness_nsa(1:npoin,j_time) )
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine set_xsubscale( phys, i_time, i_elem, i_gaus, gamme, phi )
    use def_nastal, only: kfl_track_nsa, umosg_nsa, densg_nsa, enesg_nsa
    implicit none
    logical(ip), intent(in) :: phys
    integer(ip), intent(in) :: i_time, i_elem, i_gaus
    real(rp), intent(in   ) :: gamme
    real(rp), intent(inout) :: phi(ndofn) 
    !
    real(rp), pointer :: umoss(:,:,:) => null()
    real(rp), pointer :: denss(  :,:) => null()
    real(rp), pointer :: eness(  :,:) => null()
    !
    real(rp) ::  sub(ndofn)
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    sub(1:ndofn) = 0.0
    if(kfl_track_nsa == 1) then
        umoss => umosg_nsa(1:ndime,1:nelem,1:mgaus,i_time)
        denss => densg_nsa(        1:nelem,1:mgaus,i_time) 
        eness => enesg_nsa(        1:nelem,1:mgaus,i_time) 
        !
        sub(1:ndime) = umoss(1:ndime,i_elem,i_gaus)
        sub(ndime+1) = denss(        i_elem,i_gaus)
        sub(ndime+2) = eness(        i_elem,i_gaus)
        if(phys) call cons2phys( sub(1:ndofn), gamme, sub(1:ndofn)) 
    end if
    phi(1:ndofn) = phi(1:ndofn) + sub(1:ndofn)    
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine get_xsubscale(phys, i_time, i_elem, i_gaus, gamme, &
                           x_res, x_velmo, x_sound, h_inv, p_node, x_sub, x_tau)
    use def_nastal, only: kfl_taufa_nsa, zensa
    use def_nastal, only: umosg_nsa, densg_nsa, enesg_nsa
    implicit none
    logical(ip), intent(in)  ::  phys
    integer(ip), intent(in)  ::  i_time, i_elem, i_gaus
    real(rp), intent( in)    ::  gamme
    real(rp), intent( in)    ::  x_res(ndofn)
    real(rp), intent( in)    ::  x_velmo
    real(rp), intent( in)    ::  x_sound
    real(rp), intent( in)    ::  h_inv
    integer(ip), intent( in) ::  p_node
    real(rp), intent(out)    ::  x_sub(ndofn)
    real(rp), intent(out)    ::  x_tau(ndofn)
    !
    real(rp) :: xtau(ndofn), aux(ndofn)
    real(rp) :: qufac 
    integer(ip) :: idofn 
    !---------------------------------------------------------------------||---!
    !
    !-------------------------------------------------------| ssxphi, xtau |---! 
    qufac = 1.0_rp
    if((ndime==2).and.(p_node>=4)) qufac = 0.5_rp
    if((ndime==3).and.(p_node>=5)) qufac = 0.5_rp
    !
    xtau = 0.0_rp
    if(kfl_taufa_nsa(1,2) == 1) xtau = xtau + x_velmo * h_inv
    if(kfl_taufa_nsa(3,2) == 1) xtau = xtau + x_sound * h_inv
    !
    x_sub(1:ndofn) = 0.0
    do idofn = 1,ndofn
      if(xtau(idofn) > zensa) x_sub(idofn) = qufac/xtau(idofn) * x_res(idofn)
    end do
    x_tau(1:ndofn) = qufac/xtau(1:ndofn) !-> xsubdi_matrix
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
    !> PHIsg_nsa(:,1) is used by nsa_updtss.
    if(phys) call phys2cons( x_sub(1:ndofn), gamme, aux(1:ndofn) ) 
    umosg_nsa(1:ndime,i_elem,i_gaus,i_time) = aux(1:ndime)
    densg_nsa(        i_elem,i_gaus,i_time) = aux(ndime+1)
    enesg_nsa(        i_elem,i_gaus,i_time) = aux(ndime+2)
    !---------------------------------------------------------------------||---!
    !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !--------------------------------------------------------------| PRIVATE |---!
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine SMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(inout) :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: Wk(ndofn), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: inv_c, inv_c2, inv_c_rho 
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    inv_c     = 1.0_rp/cvel
    inv_c2    = inv_c*inv_c
    inv_c_rho = 1.0_rp/Wk(1)*inv_c
    Ak(1:ndofn,1:ndofn)   = 0.0_rp
    !> |
    Ak(1,(/1_ip,ndime+2/)) = inv_c2*0.5_rp 
    Ak(1,         kk+1) = inv_c2
    Ak(ndime+2, (/1_ip,ndime+2/)) = 0.5_rp
    !> _
    Ak(kk+1,      1) = -0.5_rp*inv_c_rho
    Ak(kk+1,ndime+2) =  0.5_rp*inv_c_rho 
    !> |_ 
    do ii = 2,ndime+1
      Ak(ii,ii) = 1
    enddo    
    Ak(kk+1,kk+1) = 0.0_rp 
  end subroutine 
  !!
  subroutine SINVMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(out)   :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: Wk(ndofn), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: c2, c_rho, rho, vk 
    real(rp)    :: lambda(ndofn)  
    integer(ip) :: ii
    ! W  = [rho, u1, u2, u3, p]^T 
    ! S^{-1}_k * F^{k} S_k = \Lambda_k
    ! F^k = [rho*u_i*u_k+p, rho*u_k, (rho*e+p)*u_k]^T  
    rho   = Wk(1)
    vk    = Wk(kk+1)
    c2    = cvel*cvel
    c_rho = cvel*rho
    Ak(:,:) = 0.0_rp
    ! |
    Ak(1,ndime+2) =  1.0_rp 
    Ak(1,   kk+1) = -c_rho 
    Ak(ndime+2,   kk+1) = c_rho 
    Ak(ndime+2,ndime+2) = 1.0_rp 
    ! _
    Ak(kk+1,      1) =  c2
    Ak(kk+1,ndime+2) = -1.0_rp
    ! |_ 
    do ii = 2,ndofn
      Ak(ii,ii) = 1
    enddo    
    Ak(kk+1,kk+1) = 0.0_rp
    !
    lambda(2:ndime+1) = vk
    lambda(        1) = vk - cvel 
    lambda(ndofn) = vk + cvel
    !
    do ii = 1,ndofn
      Ak(1:ndofn,ii) = Ak(1:ndofn,ii)*lambda(1:ndofn)
    enddo
  end subroutine SINVMTRX
  !!
  subroutine PHYMTRX(Ak, Wk, kk, gamme)
    implicit none
    real(rp), intent(inout) :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: Wk(ndofn), gamme
    integer(ip), intent(in) :: kk
    real(rp)    :: vkk, rho, pre
    integer(ip) :: ii
    !
    !! Ak = dF/dUk, Ak = A(row, col, k) 
    !! W  = [rho, u1, u2, u3, p]^T 
    !! e  = cp*T - p/rho =  cv*T -> rho*e = p/(gamma-1)  
    !! E  = e + 0.5*vi*vi
    !! p  = rho*R*T 
    !! cp = cv + R, gamma = cp/cv
    !
    vkk = Wk(kk+1)
    rho = Wk(1)
    pre = Wk(ndime+2)
    !
    Ak(:,:)         = 0.0
    Ak(1,kk+1)      = rho
    Ak(ndofn,kk+1)  = gamme*pre
    Ak(kk+1,ndofn)  = 1.0/rho
    !
    do ii = 1,ndofn
      Ak(ii,ii) = vkk
    enddo
  end subroutine PHYMTRX  
  !!
  subroutine CMTRX(Ak, U, kk, gamme) 
    implicit none
    real(rp), intent(out)   :: Ak(ndofn, ndofn)
    real(rp), intent(in)    :: U(ndofn)
    real(rp), intent(in)    :: gamme 
    integer(ip), intent(in) :: kk 
    integer(ip) :: jj
    real(rp)    :: rho, ene, uk(ndime), ukuk, Rcv  
    !!
    !! U  = [rho*u1, rho*u2, rho*u3, rho, rho*E]^T  
    !! F  = [rho*ui*uk+p, rho*uk, (rho*E+p)*uk]^T, E =  e + vivi/2 
    !! Ak = dF/dUk, Ak = A(row, col, k) 
    !! Rcv = R/Cv = gamma - 1
    !!
    Ak(1:ndofn,1:ndofn) = 0.0
    rho  = U(ndime+1)
    ene  = U(ndime+2)/rho
    uk   = U(1:ndime)/rho 
    ukuk = sum(uk(1:ndime)*uk(1:ndime))
    Rcv  = gamme - 1.0  
    !!  
    Ak(     kk, ndime+1) = Rcv * 0.5_rp * ukuk                                 ! A(1, 4, 1) = 1/2 * R/Cv * Ui*Ui 
    Ak(     kk, ndime+2) = Rcv                                                 ! A(1, 5, 1) = R/Cv 
    Ak(ndime+2, ndime+1) = -uk(kk) * ((1.0_rp + Rcv) * ene - Rcv * ukuk)       ! A(5, 4, 1) = -Ui * ((1+R/Cv)*E/rho - R/Cv*Ui*Ui)
    Ak(ndime+2, ndime+2) = (1.0_rp + Rcv) * uk(kk)                             ! A(4, 4, 1) = (1+R/Cv) * Ui 
    Ak(ndime+1,      kk) = 1.0_rp                                              ! A(4, 1, 1) = 1.0 
    Ak(ndime+2,      kk) = (1.0_rp + Rcv) * ene - Rcv * 0.5_rp * ukuk          ! A(5, 1, 1) = (1+R/Cv) * E/rho - 1/2 * R/Cv Ui*Ui 
    do jj = 1,ndime
      Ak(     kk,      jj) = Ak(     kk,      jj) - Rcv * uk(jj)           ! A(1, 1:3, 1)   += -R/Cv * Vj  
      Ak(     jj,      jj) = Ak(     jj,      jj) + uk(kk)                 ! A(1:3, 1:3, 1) += Vi 
      Ak(     jj,      kk) = Ak(     jj,      kk) + uk(jj)                 ! A(1:3, 1, 1)   += Vj 
      Ak(     jj, ndime+1) = Ak(     jj, ndime+1) - uk(jj) * uk(kk)        ! A(1:3, 4, 1)   += -Vi*Vj  
      Ak(ndime+2,      jj) = Ak(ndime+2,      jj) - Rcv * uk(jj) * uk(kk)  ! A(5, 1:3, 1)   += -R/Cv * Vi*Vj  
    end do
  end subroutine CMTRX
  !!
  subroutine PMTRX(Ak, Wk, gamme, vkvk) !! ??
    implicit none
    real(rp), intent(out) :: Ak(ndofn, ndofn)
    real(rp), intent(in)  :: Wk(ndofn), gamme, vkvk
    real(rp)    :: rho 
    integer(ip) :: ii

    rho = Wk(1)

    Ak(:,:) = 0.0_rp
    !> |_ 
    do ii = 2,ndime+1 
      Ak(     ii,ii) = 0.5_rp*rho
      Ak(     ii, 1) = 0.5_rp*Wk(ii)
      Ak(ndime+2,ii) = Wk(ii)*rho 
    enddo
    Ak(      1,      1) =  1.0_rp
    Ak(ndime+2,ndime+2) =  1.0_rp/(gamme-1.0_rp)
    Ak(ndime+2,      1) = -0.5_rp*vkvk*vkvk 
  end subroutine 
  !!
  subroutine cons2phys(Uk, gamme, Wk)
    implicit none
    real(rp), intent(in)  :: Uk(ndofn), gamme  
    real(rp), intent(out) :: Wk(ndofn) 
    real(rp) :: aux(ndofn) 
    !> P = (gamma-1.0) * rho * (Et - 1/2 |V|) 
    !>   = (gamme-1.0) * (ENE - 0.5*sum(MOM(1:ndime)*MOM(1:ndime))/RHO)
    !Wk(1:ndime) = Uk(1:ndime)/Uk(ndime+1) 
    !Wk(ndime+1) = Uk(ndime+1) 
    !Wk(ndime+2) = (gamme-1.0)*(Uk(ndime+2) - 0.5/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime)))  
    !> [mom, rho, et] -> [rho, vel, p] 
    aux(      1)   = Uk(ndime+1) 
    aux(2:ndime+1) = Uk(1:ndime)/Uk(ndime+1) 
    aux(ndime+2)   = (gamme-1.0_rp)*(Uk(ndime+2) - 0.5_rp/Uk(ndime+1)*sum(Uk(1:ndime)*Uk(1:ndime))) 
    Wk = aux 
  end subroutine cons2phys
  !!
  subroutine phys2cons(Wk, gamme, Uk)
    implicit none
    real(rp), intent(in)  :: Wk(ndofn), gamme
    real(rp), intent(out) :: Uk(ndofn)
    real(rp) :: aux(ndofn)
    !> rho*Et = P/(gamma-1) + 1/2*rho*|V| 
    !>        =   rho*Cv*T  + 1/2*rho*|V|
    !Uk(1:ndime) = Wk(1:ndime)*Wk(ndime+1) 
    !Uk(ndime+1) = Wk(ndime+1)
    !Uk(ndime+2) = Wk(ndime+2)/(gamme-1.0) + 0.5*Wk(ndime+1)*sum(Wk(1:ndime)*Wk(1:ndime))
    !> [rho, vel, p] -> [mom, rho, et]
    aux(1:ndime) = Wk(2:ndime+1)*Wk(1)
    aux(ndime+1) = Wk(1)
    aux(ndime+2) = Wk(ndime+2)/(gamme-1.0_rp) + 0.5_rp*Wk(1)*sum(Wk(2:ndime+1)*Wk(2:ndime+1))
    Uk = aux 
  end subroutine phys2cons  
  !!
  subroutine csound(Uk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Uk(ndofn), gamme
    real(rp), intent(out) :: c, vmod  
    !> c2 = gamma*P/rho = gamma(gamma-1)(Et-1/2|V|) 
    vmod = sum(Uk(1:ndime)/Uk(ndime+1)*Uk(1:ndime)/Uk(ndime+1)) 
    c    = Uk(ndime+2)/Uk(ndime+1) - 0.5_rp*vmod
    vmod = sqrt(vmod) 
    c    = sqrt(gamme*(gamme-1.0_rp)*c)
  end subroutine csound
  !!
  subroutine psound(Wk, gamme, c, vmod)
    implicit none
    real(rp), intent(in)  :: Wk(ndofn), gamme
    real(rp), intent(out) :: c, vmod 
    !> c^2 = gamma*p/rho = gamma*(gamma-1)*cv*T = gamma*(gamma-1)*e
    !> W = [rho, vel, p] 
    c    = Wk(ndime+2)/Wk(1) 
    c    = sqrt( gamme*c )
    vmod = sqrt( dot_product(Wk(2:ndime+1), Wk(2:ndime+1)) )
  end subroutine psound
  !!
  subroutine calc_shock(kapsh, grad, cartd, shmet, vecsk) 
    use def_nastal, only: zensa, shock_nsa
    implicit none
    real(rp), intent(inout) :: kapsh(ndofn)
    real(rp), intent(in)    :: grad(ndofn, ndime), cartd(ndime)
    real(rp), intent(in)    :: shmet(ndime, ndime, ndofn) 
    real(rp), intent(out)   :: vecsk(ndofn)
    integer(ip) :: ii, jj, kk 

    if(shock_nsa < 1000.0_rp*zensa) then
      kapsh = 0
      return
    end if
    ! 
    do kk = 1,ndofn
      if(kapsh(kk) > 0) then
        do jj = 1,ndime
          do ii = 1,ndime
            vecsk(kk) = vecsk(kk) + cartd(ii)*shmet(jj,ii,kk)*grad(kk,jj)
          end do
        end do
        vecsk(kk) = -vecsk(kk)
      end if
    end do
  end subroutine  
  !!
  subroutine init_shock(kapsh, grad, resid, difeq, xvel2, hmin, shmet)
    use def_nastal, only: zensa, shock_nsa
    implicit none
    real(rp), intent(inout) :: kapsh(ndofn)
    real(rp), intent(in)    :: grad(ndofn,ndime), resid(ndofn)
    real(rp), intent(in)    :: difeq(ndofn), xvel2, hmin
    real(rp), intent(out)   :: shmet(ndime, ndime, ndofn)

    real(rp)    :: ficve(ndofn) 
    real(rp)    :: grad2(ndofn), grmod(ndofn), remod(ndofn)
    real(rp)    :: shpec, shfau, zesho, zesho2, aux01 
    integer(ip) :: napsh
    integer(ip) :: idofn 
    
     zesho = zensa*1000.0_rp*zensa 

     kapsh  = 1
     shmet  = 0.0_rp

     if(xvel2 < zesho) kapsh = 0

     if(shock_nsa < zesho) then
        kapsh = 0
        return
     end if

     napsh = 0
     do idofn = 1,ndofn
        grad2(idofn) = sum(grad(idofn,1:ndime)*grad(idofn,1:ndime))
        grmod(idofn) = sqrt(grad2(idofn))
        remod(idofn) =  abs(resid(idofn))
        !
        if(grad2(idofn) < zesho) then
           kapsh(idofn) = 0
           ficve(idofn) = 0.0_rp
        else
           ficve(idofn) = remod(idofn) / grmod(idofn)
        end if
        !
        if(ficve(idofn) < zesho) kapsh(idofn) = 0
        napsh = napsh + kapsh(idofn)
     end do

     if(napsh == 0) return

     zesho2 = zesho*zesho  
     do idofn = 1,ndofn
        !
        shfau = shock_nsa
        if(difeq(idofn) > zesho2) then
           shpec = ficve(idofn) * 0.5_rp * hmin/difeq(idofn)
           if(shpec > 0.0_rp) shfau = shock_nsa - 1.0/shpec
           if(shfau < 0.0_rp) shfau = 0.0_rp
        end if
        if(idofn == ndime+1)  shfau = shock_nsa
        !
        aux01 = 0.5_rp * shfau * hmin * ficve(idofn)

        if(ndime == 2) then
          shmet(1,1,idofn) = aux01
          shmet(1,2,idofn) = 0.0
          shmet(2,1,idofn) = shmet(1,2,idofn)
          shmet(2,2,idofn) = aux01
        endif

       if(ndime == 3) then
          shmet(1,1,idofn) = aux01
          shmet(1,2,idofn) = 0.0
          shmet(1,3,idofn) = 0.0

          shmet(2,1,idofn) = shmet(1,2,idofn)
          shmet(2,2,idofn) = aux01
          shmet(2,3,idofn) = 0.0

          shmet(3,1,idofn) = shmet(1,3,idofn)
          shmet(3,2,idofn) = shmet(2,3,idofn)
          shmet(3,3,idofn) = aux01
       end if
     enddo 
  end subroutine
  !!
  subroutine gauss2nodes(gprop, dvol, nprop) 
    implicit none
    real(rp), intent(in ) :: gprop(nelem, mgaus) 
    real(rp), intent(in ) ::  dvol(nelem, mgaus)
    real(rp), intent(out) :: nprop(npoin)
    real(rp)    :: ficve(ndofn)
    real(rp)    :: gshape(mgaus)
    real(rp)    ::   gfac(mgaus)
    integer(ip) :: col, dimsize
    integer(ip) :: inode, ielem, ipoin, pelty, pgaus

    col     = 1
    dimsize = 1 

    nprop(1:npoin) = 0.0_rp

    do ielem = 1,nelem
      pelty = ltype(ielem)

      pgaus = ngaus(pelty)
      gfac(1:pgaus) = dvol(ielem,1:pgaus) * gprop(ielem,1:pgaus)

      do inode = 1,nnode(pelty)
        gshape(1:pgaus) = elmar(pelty)%shape(inode,1:pgaus)

        ipoin = lnods(inode,ielem)
        nprop(ipoin) = nprop(ipoin) + dot_product( gshape(1:pgaus), gfac(1:pgaus) )
      enddo
    enddo

    call rhsmod( dimsize, nprop(1:npoin) )
    nprop(1:npoin) = nprop(1:npoin)/vmass(1:npoin)
  end subroutine 
  !!
  subroutine getgammaker(Runiv, shecp, ggamme)
    use def_master,     only: ID_NASTAL, ID_CHEMIC, kfl_coupl
    use def_master,     only: wmean
    use mod_ker_proper, only: ker_proper
    implicit none
    real(rp), intent(in   ) :: Runiv
    real(rp), intent(inout) :: shecp(npoin)
    real(rp), intent(out)   :: ggamme(npoin)
    real(rp)    :: dummy(ndime,ndime) 
    integer(ip) :: dummi, rdim 

    !> Specificic heat Cp
    if(kfl_prope/=0) then
      call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp(1:npoin))
      !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp(1:npoin),dummi,dummi,dummy(:,1),dummy(:,1))
      !> Mean molecular weight (wmean) calcuted from chemic
      if(kfl_coupl(ID_NASTAL,ID_CHEMIC)==1) ggamme(1:npoin) = 1.0_rp/(1.0_rp - runiv/wmean(1:npoin,1)/shecp(1:npoin))

      !> thermal conductivity, kappa, used by diffussive matrix 
      !call ker_proper('CONDU','IGAUS',1_ip,ielem,xdith,pnode,1_ip)
    endif

    !> Visco
    !if(kfl_prope/=0) call ker_proper('VISCO','NPOIN',dummi,dummi,visco(:,icomp),dummi,dummi,xrano,xrano)
  end subroutine
!!
  subroutine getmaxmach(mom, rho, ene, gamme, mach)
    implicit none
    real(rp), intent(in ) :: mom(ndime, npoin), rho(npoin), ene(npoin), gamme(npoin)
    real(rp), intent(out) :: mach(npoin)
    real(rp)    :: Uk(ndofn), c, vmod
    integer(ip) :: ipoin

    !if(INOTMASTER) then
      mach(1:npoin) = -666.666_rp

      do ipoin = 1,npoin
        Uk(1:ndofn) = (/ mom(1:ndime,ipoin), rho(ipoin), ene(ipoin) /)
        call csound(Uk, gamme(ipoin), c, vmod)
        if(c/=0.0_rp) mach(ipoin) = vmod/c
      end do
    !endif
  end subroutine
!!
  subroutine getvorticity(dvel, vort)
    implicit none
    real(rp), intent(in ) :: dvel(ndime,ndime)
    real(rp), intent(out) :: vort 
    real(rp) :: omega(3)

    omega(1:3) = 0.0_rp 
    omega(3) = dvel(2,1) - dvel(1,2) !> dv_y/r_x - dv_x/r_y 
    if(ndime==3) then 
    omega(2) = dvel(1,3) - dvel(3,1) !> dv_x/r_z - dv_z/r_x
    omega(1) = dvel(3,2) - dvel(2,3) !> dv_z/r_y - dv_y/r_z
    endif
 
    vort = sqrt( dot_product(omega(1:3), omega(1:3)) )
  end subroutine
!!
!==============================================================================! 
end module  mod_nsa_euler
