subroutine tur_elmmsudiff_spacon(&
     kfl_ortho, kfl_shock,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
     gprea,gptur,gpgrd, gprhs,gpden,gpsha,gpcar,gpvol,h, elunk, gphes, sreac, &
     gppro, gpprr, gppgr,gpmut,eltur,elwal,gpvis,elres_diff,gprhs_diff,gprea_diff)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmmsudiff_spacon
  ! NAME 
  !    tur_elmmsudiff_spacon
  ! DESCRIPTION
  !    
  !      This calculates the partial derivative of the residual w.r.t. constants of turbulent visco (dR/dD)
  ! USES
  ! USED BY
  !    
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur,iunkn_tur,kfl_timei_tur,&
       &                        kfl_weigh_tur,dtinv_tur,kfl_advec_tur, kfl_taust_tur, staco_tur,param_tur,kfl_clipp_tur
  use def_master, only       :  kfl_lumped
  use def_kermod, only       : kfl_adj_prob, kfl_ndvars_opt
  use mod_matrix
  use mod_tauadr, only       :  tauadr
  
  implicit none
  integer(ip), intent(in)    :: pnode,plapl,pgaus, kfl_ortho, kfl_shock
  real(rp),    intent(in)    :: gpvel(ndime,pgaus), gptur(nturb_tur,3,pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus), shock_tur
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), h, elunk(pnode, 2),sreac(pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)             ! dNk/dxidxj
  real(rp),    intent(in)    :: gppro(pgaus), gpprr(pgaus), gppgr(ndime, pgaus)
  real(rp),    intent(in)    :: gpmut(pgaus)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)    :: elwal(pnode)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(out)   :: elres_diff(pnode,ndime,pnode)
  real(rp),    intent(in)    :: gprhs_diff(ndime,pnode,pgaus),gprea_diff(ndime,pnode,pgaus)
  
  integer(ip)                :: inode,jnode,idime,igaus,knode,kdime
  real(rp)                   :: fact1,fact2,fact3, resid(pnode), gpper(pnode), rhsit , uscoe
  real(rp)                   :: xmuit, tau, gpadv(pnode), gpnve, grvgr, resi2(pnode)
  real(rp)                   :: rhnve, gpad1(pnode), rresi, ugrau, grtur(ndime), grnor,  xmui3, switc
  real(rp)                   :: CD, gplap, rhnvv, auxve(ndime), gppe2(pnode)
  real(rp)                   :: grtur_forw(pgaus,ndime)
  
  real(rp)                   :: elmat_diff(pnode,pnode,ndime,pnode),elrhs_diff(pnode,ndime,pnode)
  real(rp)                   :: tau_diff(ndime,pnode)
  real(rp)                   :: resid_diff(pnode,ndime,pnode),gpad1_diff(pnode,ndime,pnode)
  real(rp)                   :: gpadv_diff(pnode,ndime,pnode),gpper_diff(pnode,ndime,pnode)
  real(rp)                   :: rhsit_diff(ndime,pnode),xmuit_diff,fact2_diff
  real(rp)                   :: rhsitwal_diff(ndime,pnode),residwal_diff(pnode,ndime,pnode)
  real(rp)                   :: elmatwal_diff(pnode,pnode,ndime,pnode),elrhswal_diff(pnode,ndime,pnode)
  
  
  elres_diff = 0.0_rp
      
  
  !-------------------------------------------------------------------
  !
  ! Assembly
  !
  !-------------------------------------------------------------------
  do igaus = 1,pgaus
    elrhs_diff = 0.0_rp 
    elmat_diff = 0.0_rp 
    !
    ! Stabilization parameter without reactive (nonlinear) term :
    ! tau = 1.0_rp /( 4.0_rp*gpdif(igaus)/chale(2)/chale(2) + 2.0_rp*rhnve/chale(1) ) 
    !
    call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip)      
    rhnve = gpden(igaus) * gpnve 
    call tauadr(&
         kfl_taust_tur,staco_tur,rhnve,gpdif(igaus),sreac(igaus),&
         h,h,tau)
    !
    ! Calculus of residual resid and perturbation function Pi=gppre
    !
    ! Rj  = rho*Nj/dt + rho*a.grad(Nj) + s*Nj
    ! R2j = -grad(k).grad(Nj) - k*Lap(Nj)
    ! Pi  = Ni*(1-tau*s) + tau*a.grad(Ni)
    !
    do inode = 1,pnode
      gpad1(inode) = 0.0_rp
      do idime = 1,ndime
        gpad1(inode) = gpad1(inode) + gpvel(idime,igaus) * gpcar(idime,inode,igaus)   
      end do
      gpadv(inode) = gpad1(inode) * gpden(igaus)
      gpper(inode) = ( gpsha(inode,igaus)*(1.0_rp-tau*sreac(igaus)) + tau*gpadv(inode)) * gpvol(igaus)        
    end do !inode
        
    do inode = 1,pnode
      resid_diff(inode,1,1) = gprea_diff(1,1,igaus) * gpsha(inode,igaus)
    end do !inode
    rhsit_diff(1,1) = gprhs_diff(1,1,igaus)
    !
    ! Assembly of the matrix and rhs for ASGS, SUPG or Full OSS 
    !
    if( kfl_ortho /= 2 ) then 
      do inode = 1,pnode
        elrhs_diff(inode,1,1) = elrhs_diff(inode,1,1) + rhsit_diff(1,1) * gpper(inode)
        do jnode = 1,pnode
          elmat_diff(inode,jnode,1,1) = elmat_diff(inode,jnode,1,1) + resid_diff(jnode,1,1) * gpper(inode)
        enddo !jnode
      enddo!inode           
    end if !kfl_ortho
        
    do inode = 1,pnode
      do jnode = 1,pnode
        elres_diff(inode,1,1) = elres_diff(inode,1,1) + elmat_diff(inode,jnode,1,1)*eltur(1,jnode)
      enddo
      elres_diff(inode,1,1) = elres_diff(inode,1,1) - elrhs_diff(inode,1,1)
    enddo
        
  end do !igaus
   
  
end subroutine tur_elmmsudiff_spacon
