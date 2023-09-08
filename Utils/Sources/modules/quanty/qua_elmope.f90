subroutine qua_elmope(itask)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_elmope
  ! NAME 
  !    qua_elmope
  ! DESCRIPTION
  !    Compute elemental matrix and RHS 
  !    itask==0   for shoridenger equation
  !    itask==1   for Poisson equation (DFT - All electron)
  ! USES
  ! USED BY
  !    qua_matrix
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  implicit none

  integer(ip),intent(in) :: itask

  real(rp)    :: elmat(mnode,mnode)
  real(rp)    :: elrhs(mnode),prueba
  integer(ip) :: ielem,igaus                           ! Indices and dimensions
  integer(ip) :: pelty,pnode,porde,ptopo,inode
  integer(ip) :: pgaus,plapl,jnode,knl,jnl,ii,jj
  real(rp)    :: elpro(mnode)                          ! Gather 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpvol(mgaus)                          ! |J|*w
  real(rp)    :: gpelf(mgaus)                          ! r, a
  real(rp)    :: gpion(mgaus)                          ! Pot Ionico 
  real(rp)    :: gphar(mgaus)                          ! Pot Hartree
  real(rp)    :: gppxc(mgaus)                          ! Pot XC
  real(rp)    :: gprhs(mgaus)                          ! f (all terms)
  real(rp)    :: gpcon(mgaus)                          ! k
  real(rp)    :: gposc(mgaus)
  real(rp)    :: gpcod(ndime,mgaus)
  real(rp)    :: gpdif(mgaus)                          ! k
  real(rp)    :: gpden(mgaus)                          ! rho and then rho*cp
  real(rp)    :: gpcou(mgaus)                          ! Coulomb
  real(rp)    :: gpaxi(mgaus)                          ! axisymetric
  real(rp)    :: gpesp(mgaus)                          ! Espheric
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: Anl(mnode,mnode)                      ! nonlocal PP
  !
  ! Initialization
  !
  do igaus = 1,mgaus
     gpaxi(igaus) = 0.0_rp
     gpesp(igaus) = 0.0_rp
     gpcon(igaus) = 0.0_rp
     gposc(igaus) = 0.0_rp
     gppxc(igaus) = 0.0_rp
     gpelf(igaus) = 0.0_rp
  end do
  !
  ! Loop over elements
  !  
  elements: do ielem=1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)
     plapl = 0
     porde = lorde(pelty)
     ptopo = ltopo(pelty)
     !
     ! Gather operations
     !
     call qua_elmgat(&
          pnode,lnods(1,ielem),elpro,elcod)
     !
     ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
     !
     call elmcar(&
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
          gphes,ielem)
     !
     !  Properties:  GPDIF  coeficiente de difusion
     !
     do igaus = 1,pgaus
        gpdif(igaus) = 0.5_rp * (1.0_rp + real(itask)) 
        ! itask==0  schrodinger eq.  => gpdif = 1/2
        ! itask==1  Poisson          => gpdif = 1
     end do
     !
     ! Equation coefficients
     !
     call qua_elmpre(itask,&
          pnode,pgaus,gpden,gpdif,elmar(pelty)%shape,gpcar,elpro,&
          elcod,gprhs,gpcod,gpcou,gpaxi,lnods(1,ielem),gposc,&
          gpelf,gpesp,gppxc,gphar,gpion)
     !
     ! Assemble element matrix and right-hand side ELMAT, ELRHS
     !     
     call qua_elmmat(&
          pnode,pgaus,gpdif,gposc,gpvol,elmar(pelty)%shape,gpcar,&
          gpcod,gprhs,gpcou,gpaxi,elmat,gpelf,gpesp,elrhs,gpion,&
          gphar,gppxc) 


     ! OJO!!!!!!!!!!!!!!!!******
     DO inode = 1,pnode
        DO jnode = 1,pnode
           Anl(inode,jnode) = 0.0_rp
        enddo
     enddo

     if(kfl_dftgs_qua==1) then
        do knl=1,lmaximo+1
           do jnl=1,2*(knl-1)+1
              if(nl_denom(knl,jnl).ne.0.0) then
                 DO inode=1,pnode
                    DO jnode=1,pnode
                       ii=lnods(inode,ielem)
                       jj=lnods(jnode,ielem)               
                       prueba =nonloc(knl,jnl,ii)*CONJG(nonloc(knl,jnl,jj))
                       Anl(inode,jnode)=Anl(inode,jnode)+ (prueba)/nl_denom(knl,jnl)
                    enddo
                 enddo
              endif
           enddo
        enddo
     endif

     do inode = 1,pnode
        do jnode = 1,pnode
           !elmat(inode,jnode) = elmat(inode,jnode) + Anl(inode,jnode)
        enddo
     enddo
     !
     ! Prescribe Dirichlet boundary conditions
     !
     call qua_elmdir(itask,pnode,lnods(1,ielem),elmat,elrhs)
     !
     ! Assembly
     !
     if( itask == 1_ip ) then    
        call assrhs(solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
     endif

!if( itask == 0 ) then
!   if(ielem == 721 ) then
!      print*,elmat
!      stop
!   end if
!end if
     call assmat(&
          solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
          ielem,lnods(1,ielem),elmat,amatr)

  end do elements

end subroutine qua_elmope
!------------------------------------------------------------------------
! NOTES
! 
! Governing equation
!
!***
!------------------------------------------------------------------------
