subroutine hlm_elmope()

  !------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_elmope.f90
  ! NAME 
  !    hlm_elmope
  ! DESCRIPTION
  !      This routine
  !      1. Computes element matrix and element RHS for each element in a mesh;
  !      2. Assembles elemet equations into the system equations.
  ! USES
  !    hlm_assemb
  ! USED BY
  !    hlm_matrix
  !------------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz

  implicit none

  complex(rp)         :: elmat(4*mnode,4*mnode),elrhs(4*mnode),summ,sumr   !Element matrix, element RHS
  complex(rp)         :: pvecpo(ndime,mnode)                     !Primary magnetic vector potential
  complex(rp)         :: pscapo(mnode)                           !Primary electric scalar potential
  integer(ip)         :: ielem,igaus,idime,inode,ii,jj            !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,ipoin,poin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  real(rp)            :: elcod(ndime,mnode),abselmat(4*mnode,4*mnode)
  real(rp)            :: gpvol(mgaus)                            !w * |J|(gaus point)
  complex(rp)         :: gprea(ncond_hlm,mgaus),dgprea(ncond_hlm,mgaus)  !Reaction terms
  complex(rp)         :: gprhs(4,mnode)                          !f (all terms)
  real(rp)            :: gpcar(ndime,mnode,mgaus)                !dNk/dxj ...
  real(rp)            :: gphes(ntens,mnode,mgaus)                !dNk/dxidx ...

  real(rp)            :: dmax,dmin,rcons

  amatx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
  rhsix(:)=cmplx(0.0_rp,0.0_rp,kind=rp)

  !Loop over elements
  elements: do ielem = 1,nelem

	!Element dimensions
	pelty = ltype(ielem)       !Element type	
	pnode = nnode(pelty)       !Number of element nodes 
	pgaus = ngaus(pelty)       !Number of element Gauss points
	plapl = 0_ip               !Existence of Laplasian in formulation	
	porde = lorde(pelty)       !Element order
	ptopo = ltopo(pelty)       !Element topology

	!Check material
	pmate = 1_ip
	if ( nmate > 1_ip ) then
	pmate = lmate(ielem)
	end if


	!Gather
	do inode = 1,pnode
		ipoin = lnods(inode,ielem)
		do idime = 1,ndime
			elcod(idime,inode) = coord(idime,ipoin)       !Assignment of global nodes' coordinates to local element nodes
		end do
	end do
 
	!Calculate Cartesian derivatives of element trial functions, Hessian matrix and volume: GPCAR, GPHES and PGVOL
	call elmcar(pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
		elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
		gphes,ielem)
  
	!Calculate reaction term from material properties
	do igaus = 1,pgaus
		do ii=1,ncond_hlm
			gprea(ii,igaus) = cmplx(0.0_rp,anguf_hlm*perma_hlm(pmate)*sigma_hlm(ii,pmate),kind=rp)       !Reaction term: i * omega * mu * sigma 
			dgprea(ii,igaus) = cmplx(0.0_rp,anguf_hlm*perma_hlm(pmate)*dsigma_hlm(ii,pmate),kind=rp)     !Reaction term: i * omega * mu * dsigma 
		end do
	end do


	!!! This approach has to be re-written to avoid calculating them in each node for each element !!!
	!Create GPRHS vector that presents f (interior, volume, loads) using values of primary potentials in element nodes
	call hlm_prim_vec_pot(pnode,elcod,lnods(1,ielem),pvecpo)     !Calculate values of primary vector potential in element nodes
	call hlm_prim_sca_pot(pnode,elcod,pscapo)                    !Calculate values of primary scalar potential in element nodes

	do inode = 1,pnode
		gprhs(1,inode) = pvecpo(1,inode)
		gprhs(2,inode) = pvecpo(2,inode)
		gprhs(3,inode) = pvecpo(3,inode)
		gprhs(4,inode) = pscapo(inode)  

!!!!For testing
     !ipoin = lnods(inode,ielem) 
     !pmgvp_hlm(1,ipoin) = pvecpo(1,inode)
     !pmgvp_hlm(2,ipoin) = pvecpo(2,inode)
     !pmgvp_hlm(3,ipoin) = pvecpo(3,inode)
     !pelsp_hlm(ipoin)   = pscapo(inode)
!!!!For testing
	end do

	!Creation of element equations: ELMAT and ELRHS
	call hlm_elm_equs(pnode,pgaus,gprea,dgprea,gpvol,gprhs,elmar(pelty)%shape,gpcar,elmat,elrhs)
!if (lnods(1,ielem) == 5750 ) write(*,*) gpvol(1), ' ', elmar(pelty)%shape(1,1), ' ', gpcar(1,1,1), ' ', gpcar(2,1,1), ' ', gpcar(3,1,1)
!if (lnods(2,ielem) == 5750 ) write(*,*) gpvol(1), ' ', elmar(pelty)%shape(2,1), ' ', gpcar(1,2,1), ' ', gpcar(2,2,1), ' ', gpcar(3,2,1)
!if (lnods(3,ielem) == 5750 ) write(*,*) gpvol(1), ' ', elmar(pelty)%shape(3,1), ' ', gpcar(1,3,1), ' ', gpcar(2,3,1), ' ', gpcar(3,3,1)
!if (lnods(4,ielem) == 5750 ) write(*,*) gpvol(1), ' ', elmar(pelty)%shape(4,1), ' ', gpcar(1,4,1), ' ', gpcar(2,4,1), ' ', gpcar(3,4,1)
!if (lnods(4,ielem) == 5750) then
!  sumr = (0.0_rp,0.0_rp)
  summ = (0.0_rp,0.0_rp)
  do jj=1,16
  	do ii=jj,16
!		abselmat(ii, jj) = abs(elmat(ii, jj))
		if (elmat(ii, jj) /= elmat(jj, ii)) summ = 1.0_rp
  	end do
  end do
  if (summ == 1.0_rp) then
	write (*,*) 'NOT symmetric!'
!  else
!	write (*,*) 'ITS OK!!!'
  end if
!  do jj=1,16
!	sumr = sumr + elrhs(jj)
!  	do ii=1,16
!		summ = summ + elmat(ii, jj)
!		write (*,'(f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4,f10.4)') abselmat(jj, 1:16)
!  	end do
!  end do
!  write (*,*) summ, sumr, sqrt(elcod(1,1)*elcod(1,1)+elcod(2,1)*elcod(2,1)+elcod(3,1)*elcod(3,1))
!end if

	!Assembly of element equations into system equations: AMATX and RHSIX
	call hlm_assemb(1_ip,pnode,lnods(1,ielem),elmat,elrhs)

  end do elements

end subroutine hlm_elmope

