!==============================================================================!
!!
!! V0.1: conservative and primitive forms  
!! V0.2: conservative and primitive forms. improved  
!!
!! Dont put all your eggs in one basket !!
!!
subroutine nsa_lodi 
  use def_nastal
  use def_master, only: inotmaster, veloc, umome, densi, energ, densi, press
  use def_master, only: rhsid, unkno, solve, solve_sol, npoin_type   
  use def_domain, only: ltype, nnode, ngaus, lnods, coord, mnode, mgaus, nelem, ndime, npoin
  use def_domain, only: elmar, hnatu, vmass, vmasc, lboel, LNODB, NBOUN, LPOTY, exnor 
  use def_parame, only: ip, rp  
  use def_kermod, only: kfl_prope
  use def_master, only: ID_NASTAL, ID_CHEMIC, kfl_coupl, kfl_paral

 !use mod_memchk, only: memchk 
 !use def_master, only: enthalpy_transport, div_enthalpy_transport, chemical_heat
 !use mod_ker_proper, only: ker_proper
  implicit none 

  ! starting letters for variables
  ! integer: i-o
  ! reales: otherwise

  integer(ip)      :: ielem, inode, pelty, pgaus, pnode, igaus, idofn, idime
  real(rp)         :: hmini, qufac
  real(rp)         :: elcod(ndime, mnode)
  real(rp)         :: detjm, xshap(mgaus) 
  real(rp)         :: xjaci(ndime, ndime), xjacm(ndime, ndime), cartd(ndime, mnode, mgaus)
  real(rp)         :: tragl(ndime, ndime), hleng(ndime) 

  integer(ip) ::  tprev, tcurr, tnext, jdofn, ievat, ipoin 
  real(rp)    ::  elconv(ndofn_nsa, ndofn_nsa, ndime, mnode)
  real(rp)    ::   xconv(ndofn_nsa, ndofn_nsa, ndime, mgaus)
  real(rp)    ::  xdconv(ndofn_nsa, ndofn_nsa, mgaus)
  real(rp)    :: xdphit(ndofn_nsa, ndime, mgaus)  
  real(rp)    ::  phit(ndofn_nsa, mnode), xphit(ndofn_nsa, mgaus)
  real(rp)    :: ssphi(ndofn_nsa, mnode), ssxphi(ndofn_nsa, mgaus)
  real(rp)    :: elvelmo, elsound
  real(rp)    ::  xvelmo(mgaus), xsound(mgaus) 
  real(rp)    ::  xres(ndofn_nsa, mgaus)
  real(rp)    ::  xtau(ndofn_nsa)
  !real(rp)    :: xdvol(nelem, mgaus)
  real(rp)    :: xshape(mgaus)   
  real(rp)    :: hinv, ggamma
  real(rp)    :: Patm, Eint 

  !> shock
 !real(rp) :: difeq(ndofn_nsa) 
 !real(rp) :: xkapsh(ndofn_nsa,mgaus)
 !real(rp) :: xshmet(ndime, ndime, ndofn_nsa, mgaus) 
 !real(rp) :: shote(ndofn_nsa)

  !> solver 
  real(rp)    :: state(ndofn_nsa), galte(ndofn_nsa), elrhs(nevat_nsa)
  integer(ip) :: ndofn, itotn
  integer(ip) :: phys_ok, cons_ok
  real(rp)    :: xdiag
  real(rp)    :: fact1, resf, idofg 
  real(rp)    :: xmat(ndofn_nsa, ndofn_nsa), xrhs(ndofn_nsa)

  !> riemman/characteristics 
  integer(ip) :: rmnn_ok, ipoity         !> actived
  integer(ip) :: chrc_id, chrc_idime, ichrc  
  real(rp)    ::  xaux(ndofn_nsa)
  real(rp)    :: Maux01(ndofn_nsa,ndofn_nsa), Maux02(ndofn_nsa,ndofn_nsa), phitaux01(ndofn_nsa)
  real(rp)    :: xvort 
  real(rp)    :: Sxinv(1:ndofn_nsa,1:ndofn_nsa)

  !> mach 
  real(rp)    :: maxmach !, mach(npoin)
!  real(rp)    :: sigma, maxmach, Kfact

  !> chemical reactions
  !real(rp)    :: gamme(npoin)
  real(rp)    :: elgamme(mnode)
  real(rp)    ::  xgamme(mgaus)
  !real(rp)    :: xWk(mgaus), xdHk(mgaus), xHk(ndime)
  !real(rp)    :: xSource(mgaus)

  integer(ip) :: dummi, rdim 
  real(rp)    :: dummy(ndime,ndime) 

  !> -----------------------------------------------------------------------------| PRE-INIT |---<!
  phys_ok   = 0 
  cons_ok   = 1 

!  Patm  = 103400.0            !> [Nm2/s2/m3]
!  Eint  = 5.0/2.0*287.0*300.0 !> [J/kg] = [J/kg/K][K], Eint = 215250.0 (gas ideal)
  !> -----------------------------------------------------------------------------| PRE-INIT |---<!

!  !> ------------------------------------------------------------------------------| Riemman |---<!
!  xgamme = adgam_nsa
!   gamme = adgam_nsa
!
!  if(INOTMASTER) then 
!    if(chem_ok) call getgammaker(gamme(1:npoin))
!    call getmaxmach(umome(1:ndime,1:npoin,1), densi(1:npoin,1), energ(1:npoin,1), gamme(1:npoin), mach(1:npoin))
!  endif 
!  call getmaxval(mach(1:npoin), maxmach)
!
!  sigma = 0.25/0.01   !> sigma/L
!  Kfact = sigma*(1.0-maxmach*maxmach)
!  !> ------------------------------------------------------------------------------| Riemman |---<!

  !================================================================================================================!
  if_master: &
  if(INOTMASTER) then
  !
  tprev = 1
  tcurr = 2
  tnext = 2
  !
  elcod  = -666.666_rp 
  elconv = -666.666_rp 
  xconv  = -666.666_rp 
  cartd  = -666.666_rp 
  xsound = -666.666_rp
  !xdvol  = -666.666_rp 
  xvelmo = -666.666_rp
  ! 
  ssphi  = 0.0_rp
  ssxphi = 0.0_rp

  xgamme = adgam_nsa
!   gamme = adgam_nsa
  if(kfl_coupl(ID_NASTAL,ID_CHEMIC) /= 0) call getgammaker( gamma_nsa(1:npoin) )

  !> Loop over elements
  elementary_loop: & 
  do ielem = 1,nelem 
    pelty = ltype(ielem)
    pnode = nnode(pelty)
    pgaus = ngaus(pelty)

    !> -------------------------------------------------------| NODE LOOP (SEE nsa_gaconsxy) |---<! 
    node_loop01: &
    do inode = 1,pnode
      ipoin = lnods(inode, ielem) 

      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
      elgamme(inode)       = gamma_nsa(ipoin) 

      phit(1:ndime, inode) = umome(1:ndime, ipoin, tprev) 
      phit(ndime+1, inode) = densi(         ipoin, tprev) 
      phit(ndime+2, inode) = energ(         ipoin, tprev) 

      !> -------------------------------------------------------------------------| Tracking |---<!
      !if(kfl_track_nsa==1) then
      !  ssphi(1:ndime, inode) = umoss_nsa(1:ndime, ipoin, tcurr)
      !  ssphi(ndime+1, inode) = denss_nsa(         ipoin, tcurr)
      !  ssphi(ndime+2, inode) = eness_nsa(         ipoin, tcurr)
      !
      !  umoss_nsa(1:ndime, ipoin, tcurr) = 0.0
      !  denss_nsa(         ipoin, tcurr) = 0.0
      !  eness_nsa(         ipoin, tcurr) = 0.0
      !
      !  phit(1:ndofn_nsa,inode) = phit(1:ndofn_nsa,inode) + ssphi(1:ndofn_nsa,inode)
      !end if
      !> ----------------------------------------------------------------------------------------<!
      call cons2phys( phit(1:ndofn_nsa,inode), elgamme(inode), phit(1:ndofn_nsa,inode) )
      call    psound( phit(1:ndofn_nsa,inode), elgamme(inode), elsound, elvelmo )
      if(elsound /= 0.0) then 
        vmach_nsa(ipoin) = elvelmo/elsound
        sound_nsa(ipoin) = elsound
      endif 
    enddo node_loop01
    !> ------------------------------------------------------------------------------------------<! 


    !> -------------------------------------------------------------------------| GAUSS LOOP |---<! 
    !> hleng and tragl at center of gravity
    call elmlen(ndime, pnode, elmar(pelty)%dercg, tragl, elcod, hnatu(pelty), hleng)
    hmini = hleng(ndime) !> hleng(ndime) is the smallest
    hinv  = 1.0_rp/hleng(ndime)

    elrhs = 0.0_rp  
    gauss_loop01: &
    do igaus = 1,pgaus
      call elmder(pnode, ndime, elmar(pelty)%deriv(1,1,igaus), elcod, cartd(1,1,igaus), detjm, xjacm, xjaci)
      xdvol_nsa(ielem,igaus) = elmar(pelty)%weigp(igaus) * detjm 

      !> ----------------------------------------------------------------------------------------<! 
      xdphit(1:ndofn_nsa,1:ndime,igaus) = 0.0_rp
      do idofn = 1,ndofn_nsa
        do idime = 1,ndime
          xdphit(idofn,idime,igaus) = dot_product( phit(idofn,1:pnode), cartd(idime,1:pnode,igaus) )
        enddo
      enddo

      !> chlieren/vorticity 
      xschlrn_nsa(ielem,igaus) = sqrt( dot_product(xdphit(1,1:ndime,igaus), xdphit(1,1:ndime,igaus)) )
      !call getvorticity( xdphit(2:ndime+1,1:ndime,igaus), xschlrn_nsa(ielem,igaus) )

      do idofn = 1,ndofn_nsa
        xphit(idofn,igaus) = dot_product( phit(idofn,1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      enddo
      !
      if (kfl_coupl(ID_NASTAL,ID_CHEMIC) /= 0) &
           xgamme(igaus) = dot_product( elgamme(1:pnode), elmar(pelty)%shape(1:pnode,igaus) )
      !
      !> ----------------------------------------------------------------------------------------<! 

      !> -------------------------------------------------------------------------| Tracking |---<!
      !if(kfl_track_nsa == 1) then
      !  ssxphi(1:ndime, igaus) = umosg_nsa(1:ndime,ielem, igaus, tprev)
      !  ssxphi(ndime+1, igaus) = densg_nsa(        ielem, igaus, tprev)
      !  ssxphi(ndime+2, igaus) = enesg_nsa(        ielem, igaus, tprev)
      !
      !  xphit(1:ndofn_nsa,igaus) = xphit(1:ndofn_nsa,igaus) + ssxphi(1:ndofn_nsa,igaus)
      !end if
      !ssxphi(1:ndofn_nsa, igaus) = 0.0 !> xphi = tau*residue 
      !!call cons2phys( xphit(1:ndofn_nsa,igaus), xgamme(igaus), xphit(1:ndofn_nsa,igaus) )
      !> ----------------------------------------------------------------------------------------<!
      call psound(xphit(1:ndofn_nsa,igaus), xgamme(igaus), xsound(igaus), xvelmo(igaus))
      !> ----------------------------------------------------------------------------------------<! 

      !> -------------------------------------------------| L = \Lambda^k * S_k^{-1} * dU/dr_k |---<! 
      if(cons_ok == 1) then
        !> primitives (V) to conservative (U). P = dU/dV  
        call PMTRX(Maux01(1:ndofn_nsa,1:ndofn_nsa), xphit(1:ndofn_nsa,igaus), xgamme(igaus), xvelmo(igaus))
        do idime = 1,ndime
!call SINVMTRX(Sxinv_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xphit(1:ndofn_nsa,igaus), idime, xsound(igaus))
!xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) = matmul( Sxinv_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xdphit(1:ndofn_nsa,idime,igaus) )

call SINVMTRX(Sxinv(1:ndofn_nsa,1:ndofn_nsa), xphit(1:ndofn_nsa,igaus), idime, xsound(igaus))
xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) = matmul( Sxinv(1:ndofn_nsa,1:ndofn_nsa), xdphit(1:ndofn_nsa,idime,igaus) )

          call SMTRX(Maux02(1:ndofn_nsa,1:ndofn_nsa), xphit(1:ndofn_nsa,igaus), idime, xsound(igaus))
          Sx_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem) = matmul(Maux01(1:ndofn_nsa,1:ndofn_nsa), Maux02(1:ndofn_nsa,1:ndofn_nsa))
        enddo 
      endif

      !if(phys_ok) then
      !  do idime = 1,ndime
      !    call SINVMTRX(Sxinv_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xphit(1:ndofn_nsa,igaus), idime, xsound(igaus))
      !    xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) = matmul( Sxinv_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xdphit(1:ndofn_nsa,idime,igaus) )
      !
      !    call SMTRX(Sx_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xphit(1:ndofn_nsa,igaus), idime, xsound(igaus))
      !  enddo
      !endif 
      !> ------------------------------------------------------------------------------------------<! 
    enddo gauss_loop01
    !> --------------------------------------------------------------------------------------------<! 
  end do elementary_loop

  !> -------------------------------------------------------------------------| gauss2nodes |---<! 
  do idofn = 1,ndofn_nsa
    do idime = 1,ndime
    !  call gauss2nodes(xchrc_nsa(idofn,1:nelem,1:mgaus,idime), xdvol_nsa(1:nelem,1:mgaus), chrc_nsa(idofn,1:npoin,idime)) 
    enddo
  enddo

  !call gauss2nodes(xschlrn_nsa(1:nelem,1:mgaus), xdvol_nsa(1:nelem,1:mgaus), schlrn_nsa(1:npoin))
  !> --------------------------------------------------------------------------------------------<! 
  endif if_master

  !> -----------------------------------------------------------------------------| max mach |---<!

  !call getmaxval( vmach_nsa(:), maxmach )
  call getmaxval( maxmach )

  !> K = sigma/L * c * (1 - M*M); c = ?? 
  kfact_lodi_nsa = sigma_lodi_nsa * (1.0_rp - maxmach * maxmach)
  !> ------------------------------------------------------------------------------| Riemman |---<!

  999 format(5f11.2)
!================================================================================================================!
!================================================================================================================!
contains
!!
  subroutine SMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(inout) :: Ak(ndofn_nsa, ndofn_nsa)
    real(rp), intent(in)    :: Wk(ndofn_nsa), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: inv_c, inv_c2, inv_c_rho 
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    inv_c     = 1.0_rp/cvel
    inv_c2    = inv_c*inv_c
    inv_c_rho = 1.0_rp/Wk(1)*inv_c
    Ak(1:ndofn_nsa,1:ndofn_nsa)   = 0.0_rp
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
 !subroutine SINVMTRX(Wk, kk, cvel)
 subroutine SINVMTRX(Ak, Wk, kk, cvel)
    implicit none
    real(rp), intent(out)   :: Ak(ndofn_nsa, ndofn_nsa)
    real(rp), intent(in)    :: Wk(ndofn_nsa), cvel
    integer(ip), intent(in) :: kk
    real(rp)    :: c2, c_rho, rho, vk 
    real(rp)    :: lambda(ndofn_nsa)  
    integer(ip) :: ii
    !> W  = [rho, u1, u2, u3, p]^T 
    !> S^{-1}_k * F^{k} S_k = \Lambda_k
    !> F^k = [rho*u_i*u_k+p, rho*u_k, (rho*e+p)*u_k]^T  
    rho   = Wk(1)
    vk    = Wk(kk+1)
    c2    = cvel*cvel
    c_rho = cvel*rho
    Ak(:,:) = 0.0_rp
    !> |
    Ak(1,ndime+2) =  1.0_rp 
    Ak(1,   kk+1) = -c_rho 
    Ak(ndime+2,   kk+1) = c_rho 
    Ak(ndime+2,ndime+2) = 1.0_rp 
    !> _
    Ak(kk+1,      1) =  c2
    Ak(kk+1,ndime+2) = -1.0_rp
    !> |_ 
    do ii = 2,ndofn_nsa
      Ak(ii,ii) = 1
    enddo    
    Ak(kk+1,kk+1) = 0.0_rp

    lambda(2:ndime+1) = vk
    lambda(        1) = vk - cvel 
    lambda(ndofn_nsa) = vk + cvel

    do ii = 1,ndofn_nsa
      Ak(1:ndofn_nsa,ii) = Ak(1:ndofn_nsa,ii)*lambda(1:ndofn_nsa)
    enddo
  end subroutine SINVMTRX
!!
  subroutine PMTRX(Ak, Wk, gamme, vkvk) 
    implicit none
    real(rp), intent(out) :: Ak(ndofn_nsa, ndofn_nsa)
    real(rp), intent(in)  :: Wk(ndofn_nsa), gamme, vkvk
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
    real(rp), intent(in)  :: Uk(ndofn_nsa), gamme  
    real(rp), intent(out) :: Wk(ndofn_nsa) 
    real(rp) :: aux(ndofn_nsa) 
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
    real(rp), intent(in)  :: Wk(ndofn_nsa), gamme
    real(rp), intent(out) :: Uk(ndofn_nsa)
    real(rp) :: aux(ndofn_nsa)
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
    real(rp), intent(in)  :: Uk(ndofn_nsa), gamme
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
    real(rp), intent(in)  :: Wk(ndofn_nsa), gamme
    real(rp), intent(out) :: c, vmod 
    !> c^2 = gamma*p/rho = gamma*(gamma-1)*cv*T = gamma*(gamma-1)*e
    !> W = [rho, vel, p] 
    c    = Wk(ndime+2)/Wk(1) 
    c    = sqrt( gamme*c )
    vmod = sqrt( dot_product(Wk(2:ndime+1), Wk(2:ndime+1)) )
  end subroutine psound
!!
  subroutine gauss2nodes(gprop, dvol, nprop) 
!    use def_domain, only: mgaus  
    implicit none

    real(rp), intent(in ) :: gprop(nelem, mgaus) 
    real(rp), intent(in ) ::  dvol(nelem, mgaus)
    real(rp), intent(out) :: nprop(npoin)
    real(rp)    :: ficve(ndofn_nsa)
    real(rp)    :: gshape(mgaus)
    real(rp)    ::   gfac(mgaus)
    integer(ip) :: col, dimsize  

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
        nprop(ipoin) = nprop(ipoin) + dot_product(gshape(1:pgaus),gfac(1:pgaus))
      enddo
    enddo

    call rhsmod(dimsize, nprop(1:npoin))
    nprop(1:npoin) = nprop(1:npoin)/vmass(1:npoin)
  end subroutine 
!!
  subroutine getgammaker(ggamme)
    use def_master,     only: ID_NASTAL, ID_CHEMIC, kfl_coupl
    use def_master,     only: wmean
    use mod_ker_proper, only: ker_proper
    implicit none
    real(rp), intent(out) :: ggamme(npoin)

    !> Specificic heat Cp
    if(kfl_prope/=0) then
      !call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa(1:npoin),dummi,dummi,dummy(:,1),dummy(:,1))
      call ker_proper('SPHEA','NPOIN',dummi,dummi,shecp_nsa(1:npoin))
      !> Mean molecular weight (wmean) calcuted from chemic
      if(kfl_coupl(ID_NASTAL,ID_CHEMIC)==1) ggamme(1:npoin) = 1.0_rp/(1.0_rp - runiv_nsa/wmean(1:npoin,1)/shecp_nsa(1:npoin))

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
    real(rp) :: Uk(ndofn_nsa), c, vmod

    !if(INOTMASTER) then
      mach(1:npoin) = -666.666_rp

      do ipoin = 1,npoin
        Uk(1:ndofn_nsa) = (/ mom(1:ndime,ipoin), rho(ipoin), ene(ipoin) /)
        call csound(Uk, gamme(ipoin), c, vmod)
        if(c/=0.0_rp) mach(ipoin) = vmod/c
      end do
    !endif
  end subroutine
!!
!  subroutine getmaxval(array, maxv)
  subroutine getmaxval(maxv)
    use def_master, only: IPARALL, ISLAVE ! I am a parallel, master or slave
    use def_master, only: nparr           ! Number of real model parameters, parre buffer size
    use def_master, only: parre           ! Working real array
    implicit none
!    real(rp), intent(in ) :: array(:)
    real(rp), intent(out) :: maxv
    real(rp), target      :: vmax(1)

    vmax(1) = 0.0_rp 

    if(INOTMASTER) vmax(1) = maxval(vmach_nsa(1:npoin))

    if(IPARALL) then
      nparr =  1
      parre => vmax
      call par_operat(2_ip) ! seek minimul values in parallel
      maxv = vmax(1)
    end if
    !print *, maxv
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
end subroutine 
!==============================================================================!


!==============================================================================!
subroutine get_normaltype_lodi()
  use def_nastal
  use def_domain, only: ltype, nnode, ngaus, lnods, coord
  use def_domain, only: mnode, mgaus, nelem, ndime, npoin
  use def_domain, only: lboel, LNODB, NBOUN, LPOTY, exnor
  use def_parame, only: ip, rp
  implicit none
  integer(ip) :: ipoin, iboun, idime 
  integer(ip) :: ipoity, ielem, pelty, pnode, pgaus, inode, igaus 
  real(rp)    ::  xaux(ndofn_nsa)
!
!  type lodi_type
!    integer(ip) :: id, idime, ichrc
!  endtype lodi_type
!  type(lodi_type) normal(npoin)
!
!
!  !> Loop over elements
!  !> ----------------------------------------------------------------------------------------------<!
!  elementary_loop: &
!  do ielem = 1,nelem
!    pelty = ltype(ielem)
!    pnode = nnode(pelty)
!    pgaus = ngaus(pelty)
!    elem_lodi_nsa(ielem)%id = -666
!
!    !> --------------------------------------------------------------------------------------------<!
!    node_loop01: &
!    do inode = 1,pnode
!      ipoin  = lnods(inode,ielem)
!      iboun  = lpoty(ipoin)       
!
!      if(normal(ipoin)%id>0) then 
!        elem_lodi_nsa(ielem)%id = ielem 
!
!        do igaus = 1,pgaus
!          elem_lodi_nsa(ielem)%normal(igaus)%id    = normal(ipoin)%id  
!          elem_lodi_nsa(ielem)%normal(igaus)%idime = normal(ipoin)%idime
!          elem_lodi_nsa(ielem)%normal(igaus)%ichrc = normal(ipoin)%ichrc 
!        enddo 
!      endif 
!    enddo node_loop01
!    !> --------------------------------------------------------------------------------------------<!
!  enddo elementary_loop 
!  !> ----------------------------------------------------------------------------------------------<!
!
end subroutine
!==============================================================================!

!==============================================================================!
subroutine nsa_lodi_normaltype()
    !> nsa_elconsxy
    !>
    !> NOTA: Subsonic outflow boundary
    !>       The incoming Riemann invariant is overwritten,
    !>       All other characteristic variables retain their
    !>       boundary values computed at the predictor step.
    !>
    use def_nastal
    use def_domain, only: ltype, nnode, ngaus, lnods, coord, mnode, mgaus, nelem, ndime, npoin
    use def_domain, only: lboel, LNODB, NBOUN, LPOTY, exnor
    use def_master, only: veloc, umome, densi, energ, densi, press, inotmaster
    use def_parame, only: ip, rp
    implicit none
    real(rp)    :: xaux(ndofn_nsa)
    real(rp)    :: inor(ndime)
    real(rp)    :: ivel(ndime), irot_vel(ndime), vel_n
    integer(ip) :: ipoin, zpoin, ipoity, idime  
    !> 
    !> type lodi_type
    !>  integer(ip) :: id, idime, ichrc
    !> endtype lodi_type
    !> type(lodi_type) normal_nsa(npoin)
    !>
    !do ipoin = 1,npoin
  if(INOTMASTER) then
     do ipoin = 1,npoin
      ipoity = lpoty(ipoin)                                              !> lpoty(ipoin)>0 means this is boundary node 
      normal_nsa(ipoin)%id    = -ipoin
      normal_nsa(ipoin)%idime = -69
      normal_nsa(ipoin)%ichrc = -69

      idime = -69
      boundary: &
      if(ipoity /= 0) then                                                !> Is it a boundary?
        !
        if(      kfl_fixno_nsa(1,ipoin)==6 ) then                         !> x's
           idime = 1
        else if( kfl_fixno_nsa(1,ipoin)==7 ) then                         !> y's
           idime = 2
        else if( kfl_fixno_nsa(1,ipoin)==8 ) then                         !> z's 
           idime = 3
        endif
        !
        if(idime>0) then                                                  !> Is the boundary actived?
           ivel(1:ndime) = veloc(1:ndime,ipoin,1)

           call nsa_rotvec(1_rp, ivel(1:ndime), exnor(:,1,ipoity), irot_vel(1:ndime), ndime)
           vel_n = irot_vel(1) !dot_product(ivel(1:ndime), inor(1:ndime))

           if(abs(vel_n)>1e-6 .and. abs(vel_n-ivel(idime))<1.e-6 ) then   !> Is ALL of flux go out?
             inor(1:ndime) = exnor(:,1,ipoity)
             if(vel_n>=0.0_rp) then !                                        <=== ??
               if(inor(idime)>=0.0_rp) then
                 normal_nsa(ipoin)%id    = ipoin
                 normal_nsa(ipoin)%idime = idime
                 normal_nsa(ipoin)%ichrc = 1                              !> L1 characteristic
               else if(inor(idime)<0.0_rp) then
                 normal_nsa(ipoin)%id    = ipoin
                 normal_nsa(ipoin)%idime = idime
                 normal_nsa(ipoin)%ichrc = ndofn_nsa                      !> L5 characteristic
               endif
             else
!DMM               call runend("Lodi: Recirculation found!! ")
!DMM               write(*,*)'Recirculation found at the boundary'
             endif
           endif
        endif
        !
      endif boundary
    enddo
    !  
  endif 
!> 
!> +Nastal 
!> |_nsa_turnon
!>  |_nsa_inibcs 
!>    kfl_fixrs => kfl_fixrs_nsa
!>    kfl_fixno => kfl_fixno_nsa
!>    bvess     => bvess_nsa(:,:,1)
!> +Nastal 
!> |_nsa_iniunk
!>   if(kfl_fixno_nsa(idofn,ipoin) <= 0) bvess_nsa(idofn,ipoin,1) = xfiel(-nfiel_nsa(XXX))%a(1,ipoin) 
!>
!> +Nastal 
!> |_nsa_turnon
!>   |_nsa_reanut 
!>     do while(words(1)/='ENDNU') 
!>       if(words(1) == 'STABI') then
!>         if(exists('MULTI')) then   !> /kernel/defmod/def_inpout/exists()
!>           if(exists('TRACK')) kfl_track_nsa = 1   !Subscales tracking
!>         endif 
!>       end if
!>     enddo 
!>
end subroutine
!==============================================================================!

!> def_nastal.f90
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
!    real(rp), pointer :: xschlrn_nsa(:,:)      !> schlieren 
!    real(rp), pointer ::  schlrn_nsa(:)        !> schlieren
!    real(rp), pointer ::   xsrc_nsa(:,:,:)     !> gauss point source 
!    real(rp), pointer ::  xchrc_nsa(:,:,:,:)   !> gauss point characteristic 
!    real(rp), pointer ::   chrc_nsa(:,:,:)     !> nodal point characteristic 
!    real(rp), pointer ::  sxinv_nsa(:,:,:,:,:) !> matrix S**-1
!    real(rp), pointer ::     sx_nsa(:,:,:,:,:) !> matrix S 
!    real(rp), pointer ::  xgalte_nsa(:,:,:,:) !> 1:ndofn_nsa,igaus,inode,ielem  
!
!    integer(ip)       :: ichrc_nsa, rmnn_nsa
!
!    !> nsa_lodi_normaltype   
!    type lodi_type
!      integer(ip) :: id, idime, ichrc
!    endtype lodi_type
!
!    type(lodi_type), pointer :: normal(:)
!
!    type type_lodi
!      integer(ip)               :: id
!       type(lodi_type), pointer :: normal(:)
!    endtype type_lodi
!
!    type(type_lodi), dimension(:), allocatable :: elem_lodi_nsa
    !-----------------------------------------------------------------------! 

!> nsa_membcs
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------!
!    allocate(elem_lodi_nsa(nelem))
!    do ielem = 1,nelem
!      allocate(elem_lodi_nsa(ielem)%normal(mgaus))
!    enddo
!
!    allocate(normal(npoin), stat=istat)
!    normal(1:npoin)%id    = -1
!    normal(1:npoin)%idime = -1
!    normal(1:npoin)%ichrc = -1
!
!    allocate( xchlrn_nsa(nelem,mgaus) ,stat=istat)
!    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xchlrn_nsa)
!
!    allocate( chlrn_nsa(npoin) ,stat=istat)
!    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs',  chlrn_nsa)
!
!    allocate(  xsrc_nsa(ndofn_nsa,mgaus,nelem) ,stat=istat)
!    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xsrc_nsa)
!
!    allocate( xchrc_nsa(ndofn_nsa,nelem,mgaus,ndime), stat=istat)
!    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', xchrc_nsa)
!
!    allocate( chrc_nsa(ndofn_nsa,npoin,ndime), stat=istat)
!    call memchk(zero,istat,mem_modul(1:2,modul),'LODI_NSA','nsa_membcs', chrc_nsa)
!
!    allocate( sxinv_nsa(ndofn_nsa,ndofn_nsa,ndime,mgaus,nelem) ,stat=istat)
!    allocate(    sx_nsa(ndofn_nsa,ndofn_nsa,ndime,mgaus,nelem) ,stat=istat)
!
!    allocate( xgalte_nsa(ndofn_nsa,mgaus,mnode,nelem), stat=istat )
!    xgalte_nsa = -666.666
    !-----------------------------------------------------------------------!

!> nsa_elerhs.f90
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
!    xsour(1:ndofn_nsa) = 0.0_rp

!    if(rmnn_nsa) then
!      if(elem_lodi_nsa(ielem)%id>0) then
!         idime = elem_lodi_nsa(ielem)%normal(igaus)%idime
!
!         if(idime>0) then
!           !xconv(1:ndofn_nsa,1:ndofn_nsa,idime) = 0.0
!
!           !xsour(1:ndofn_nsa) = -matmul(Sx_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) ) 
!           !xsrc_nsa(1:ndofn_nsa,igaus,ielem) = -matmul(Sx_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) )
!           xsour(1:ndofn_nsa) = xsrc_nsa(1:ndofn_nsa,igaus,ielem)
!         endif
!      endif
!    endif
    !
    !> elrhs(ievat) = elrhs(ievat) + dvolu*xshai*xsour(jdofn)
    !> -----------------------------------------------------------------------------------------<! 

!> nsa_cvgunk
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
    !call minmax(one, npoin, zero, chrc_nsa(ichrc_nsa,ndime,:), vamxm_nsa(1,4), vamxm_nsa(2,4)) !__my_change__
    !-----------------------------------------------------------------------!

!> nsa_inivar
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
    !../../kernel/defmod/def_kintyp.f90. nvarp = postprocess variables  
!     postp(1) % wopos( 1,61) = 'CHARA'  !> LODI      
!     postp(1) % wopos( 2,61) = 'VECTO'
!
!     postp(1) % wopos( 1,62) = 'SCHLI'  !> SCHLIEREN SCHLIEREN  
!     postp(1) % wopos( 2,62) = 'SCALA'
    !-----------------------------------------------------------------------! 

!> nsa_output
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------!
!    if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) call nsa_outwit()
    !-----------------------------------------------------------------------!

!> nsa_outvar
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
!  case( 61_ip )
!    if(INOTMASTER) then
!      do idime = 1,ndime
!        gevec(idime,1:npoin) = chrc_nsa(1,1:npoin,idime)
!      enddo
!    endif
!    !gesca => unkno
!
!  case( 62_ip )
!    if(INOTMASTER) unkno(1:npoin) = chlrn_nsa(1:npoin)
!    gesca => unkno

!> nsa_reanut 
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
!        else if(words(1)=='NSCBC') then
!           if(exists('ON   ')) then
!              lodi_nsa       = 1
!              sigma_lodi_nsa = getrea('VALUE', 0.7_rp, 'Shock capturing parameter')
!           end if
    !-----------------------------------------------------------------------!
!

!> nsa_sendat 
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
!    call iexcha( lodi_nsa )
!    call rexcha( sigma_lodi_nsa )
!    call rexcha( press_lodi_nsa )
    !-----------------------------------------------------------------------! 
!
!==============================================================================! 
