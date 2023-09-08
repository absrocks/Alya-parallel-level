subroutine got_elmop2()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmop2
  ! NAME 
  !    got_elmop2
  ! DESCRIPTION
  !    Elemental operations:
  !    1. Compute elemental matrix and RHS 
  !    2. Impose Dirichlet boundary conditions
  !    3. Assemble them
  ! USES
  !    got_elmgat
  !    elmlen
  !    elmchl
  !    got_elmtss
  !    got_elmpre
  !    got_elmpro
  !    got_elmsgs
  !    got_elmrhs
  !    got_exacso
  !    got_elmtes
  !    got_elmres
  !    got_elmmat
  !    got_elmdir
  !    got_assrhs
  !    got_assmat
  ! USED BY
  !    got_matrix
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_gotita
  implicit none
  !
  ! Element matrices and vectors (stiffness and preconditioner)
  !
  real(rp)    :: elmat(nevat_got,nevat_got)
  real(rp)    :: elrhs(nevat_got)
  !
  ! Working arrays (created here to avoid automatic arrays)
  !
  real(rp)    :: xjaci(9),xjacm(9) 
  !
  ! Indices and dimensions
  !
  integer(ip) :: ielem,igaus,ndofn,ndofr
  integer(ip) :: pnode,pgaus,pevat,pelty,porde
  !
  ! Gauss point values
  !
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
  real(rp)    :: gpvol(mgaus),gpdet                    ! w*|J|, |J|
  real(rp)    :: gppor(mgaus),gpcod(ndime,mgaus)       ! sig, x
  real(rp)    :: gpdif(mgaus),gpvno(mgaus)             ! k, |u|
  real(rp)    :: gpst1(mgaus),gpst2(mgaus)             ! tau1, tau2
  real(rp)    :: gprhs(ndime,mgaus)                    ! RHS
  real(rp)    :: gpvdr(ndime,mgaus,2)                  ! u
  real(rp)    :: gpvel(ndime,mgaus)                    ! u_a
  real(rp)    :: gpgvd(ndime,ndime,mgaus)              ! grad(u)
  real(rp)    :: gpdiv(mgaus)                          ! div(u)
  real(rp)    :: gpugu(ndime,mgaus)                    ! (u.grad)u)
  real(rp)    :: gpcdr(mgaus)
  !
  ! Element characteristics (to compute tau1 and tau2)
  !
  real(rp)    :: tragl(9),chave(ndime,2)               ! Stabilization
  real(rp)    :: chale(2),hleng(3)
  !
  ! Perturbation and residuals
  !
  real(rp)    :: resim(9,mnode,mgaus)                  ! Residual momentum
  real(rp)    :: tesmv(mnode,mgaus)
  real(rp)    :: dtcri,dummr(ndime*ndime*mgaus*mnode)  ! Others
  !
  ! Size of the residual
  !
  if(kfl_linea_got==1) then
     ndofr=1
  else
     ndofr=ndime*ndime
  end if
  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem
     !
     ! Element properties and dimensions
     !
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     porde=lorde(pelty)
     ndofn=ndime
     pevat=ndofn*pnode 

     call got_elmgat(&
          1_ip,1_ip,pnode,lnods(1,ielem),elvdr_got,elcdr_got,&
          elcod_got,elvel_got,eldif_got)
     !
     ! HLENG and TRAGL at center of gravity 
     !
     call elmlen(&
          ndime,pnode,elmar(pelty)%dercg,tragl,elcod_got,&
          hnatu(pelty),hleng)
     !
     ! Compute the characteristic length: CHALE
     ! 
     call elmchl(&
          tragl,hleng,elcod_got,elvdr_got,chave,chale,pnode,&
          lorde(pelty),hnatu(pelty),1.0_ip,kfl_ellen_got)
     !
     ! Local time step: DTINV_NSI
     !
     if(kfl_timco==2) then 
        call got_elmtss(&
             pnode,porde,elcod_got,elvel_got,elvdr_got,elcdr_got,&
             eldif_got,elmar(pelty)%shacg,elmar(pelty)%dercg,&
             elmar(pelty)%weicg,dtcri)
        dtinv_got=1.0_rp/(dtcri*safet_got)
        if(kfl_stead_got==1) dtinv_got = 0.0_rp
        if(kfl_timei_got==0) dtinv_got = 0.0_rp
     end if
     !
     ! 1st and 2nd order Cartesian derivatives GPCAR, GPVOL=dV=|J|*wg
     !
     do igaus=1,pgaus     
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
             elcod_got,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
        gpvol(igaus)=elmar(pelty)%weigp(igaus)*gpdet            ! |J|*wg
     end do
     !
     ! Preliminaries: GPCOD, GPPRE, GPGPR
     !
     call got_elmpre(&
          pnode,pgaus,elmar(pelty)%shape,gpcar,elcod_got,&
          dummr,elvdr_got,elvel_got,dummr,gpcod,gpvdr,gpvel,&
          gpgvd,dummr,gpdiv,gpugu,gpvno)
     !
     ! Properties: GPDEN, GPVIS, GPPOR, GPGVI and GPCOD
     !
     call got_elmpro(&                                    
          pnode,pgaus,1_ip,pgaus,eldif_got,elmar(pelty)%shape,&
          gpvdr,gpvel,gpvno,gpcdr,chale,gppor,gpdif)
     !
     ! Residual: RESIM
     !
     call got_elmre2(&
          pnode,pgaus,ndofr,elmar(pelty)%shape,gpcar,gppor,gpvdr,&
          gpgvd,resim) 
     !
     ! Exact solution: GPRHS
     !
     call got_elmexa(&
          2_ip,pgaus,gpcod,gpvel,gprhs)
     !
     ! Right-hand side: GPRHS
     !
     call got_elmrhs(&
          pgaus,ndofn,gppor,dummr,gpvdr,dummr,gpdiv,gpugu,gpvel,gprhs)
     !
     ! Stabilization parameters: GPST1, GPST2
     !
     call got_elmsgs(& 
          pgaus,gpcdr,gpvno,gpdiv,gppor,chale,gpst1,gpst2)
     !
     ! Test function: TESMV
     !
     call got_elmte2(&
          pnode,pgaus,elmar(pelty)%shape,gpcar,gppor,&
          gpvdr,gpvel,gpst1,tesmv) 
     !
     ! Assemble elemental matrix
     !
     call got_elmma2(&
          pgaus,pnode,pevat,ndofr,gpvol,resim,tesmv,gprhs,&
          elmat,elrhs)
     !
     ! Shock capturing
     !
     if(kfl_shocm_got/=0) then
        call got_elmsm2(&
             pnode,pgaus,pevat,ndofr,gpgvd,gpvdr,gpst1,gprhs,resim,gpvol,&
             gpcar,elvdr_got,dtinv_got,chale,elmat)
     end if
     !
     ! Prescribe Dirichlet boundary conditions
     !
     call got_elmdir(&
          1_ip,pnode,pevat,ndofn,lnods(1,ielem),elmat,elrhs)
     !
     ! Assemble LHS and RHS
     !
     call assrhs(&
          ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
     call assmat(&
          ndofn,pnode,pevat,nunkn_got(ivari_got),solve(ivari_got)%kfl_algso,&
          ielem,lnods(1,ielem),elmat,amatr)

  end do elements

end subroutine got_elmop2
