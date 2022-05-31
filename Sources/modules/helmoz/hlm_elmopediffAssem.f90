
subroutine hlm_elmopediffAssem(indvars)

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

  integer(ip), intent(in)   :: indvars

  complex(rp)         :: elmat(4*mnode,4*mnode),elrhs(4*mnode)   !Element matrix, element RHS
  complex(rp)         :: pvecpo(ndime,mnode)                     !Primary magnetic vector potential
  complex(rp)         :: pscapo(mnode)                           !Primary electric scalar potential
  integer(ip)         :: ielem,igaus,idime,inode,jnode,ii              !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,tnode,ipoin,poin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  real(rp)            :: elcod(ndime,mnode)
  real(rp)            :: gpvol(mgaus)                            !w * |J|(gaus point)
  complex(rp)         :: gprea(ncond_hlm,mgaus),dgprea(ncond_hlm,mgaus)  !Reaction terms
  complex(rp)         :: gprhs(4,mnode)                          !f (all terms)
  real(rp)            :: gpcar(ndime,mnode,mgaus)                !dNk/dxj ...
  real(rp)            :: gphes(ntens,mnode,mgaus)                !dNk/dxidx ...

  real(rp)            :: dmax,dmin,rcons
  real(rp)            :: dsigmatmp
  real(rp)            :: centerX, centerY, centerZ

  real(rp)            :: countInner

  !Loop over elements
  elements: do ielem = 1,nelem

        centerX=0.0_rp
        centerY=0.0_rp
        centerZ=0.0_rp

        dsigmatmp = 0.0_rp
        countInner=0.0_rp

        !Element dimensions
        pelty = ltype(ielem)       !Element type	
        pnode = nnode(pelty)       !Number of element nodes 
        tnode = 4_ip * pnode       !pnode * 4 unknown parameters in a node (Asx, Asy, Asz, Psis)
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
                !if(lninv_loc(ipoin)==indvars .and. diffj_isInside(lninv_loc(ipoin))>0_ip)then
                if(lninv_loc(ipoin)==indvars)then
                        !dsigmatmp = 1.0_rp
                        !dsigmatmp = 0.25_rp
                        dsigmatmp = 1.0_rp/real(pnode)
                end if
        end do

        if(dsigmatmp/=0.0_rp)then

           !Calculate Cartesian derivatives of element trial functions, Hessian matrix and volume: GPCAR, GPHES and PGVOL
           call elmcar(pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
           elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
           gphes,ielem)
  
           !Calculate reaction term from material properties
           do igaus = 1,pgaus
              do ii=1,ncond_hlm
                 gprea(ii,igaus) = cmplx(0.0_rp,anguf_hlm*perma_hlm(pmate)*dsigmatmp ,kind=rp)       !Reaction term: i * omega * mu * sigma 
                 dgprea(ii,igaus) = cmplx(0.0_rp,anguf_hlm*perma_hlm(pmate)*dsigmatmp,kind=rp)     !Reaction term: i * omega * mu * dsigma 
              end do
           end do

           !Create GPRHS vector that presents f (interior, volume, loads) using values of primary potentials in element nodes
           call hlm_prim_vec_pot(pnode,elcod,lnods(1,ielem),pvecpo)     !Calculate values of primary vector potential in element nodes
           call hlm_prim_sca_pot(pnode,elcod,pscapo)                    !Calculate values of primary scalar potential in element nodes

           do inode = 1,pnode
              gprhs(1,inode) = pvecpo(1,inode)
              gprhs(2,inode) = pvecpo(2,inode)
              gprhs(3,inode) = pvecpo(3,inode)
              gprhs(4,inode) = pscapo(inode)  
           end do

           !Creation of element equations: ELMAT and ELRHS
           !call hlm_elm_equs(pnode,pgaus,gprea,dgprea,gpvol,gprhs,elmar(pelty)%shape,gpcar,elmat,elrhs)
           call hlm_elm_equsdiff(pnode,pgaus,gprea,dgprea,gpvol,gprhs,elmar(pelty)%shape,gpcar,elmat,elrhs)


           !Assembly of element equations into system equations: DAMAT and DRHSI
           call hlm_assembdiff(1_ip,pnode,lnods(1,ielem),elmat,elrhs,indvars)


        else

        ! set elmat and elmrhs to zero
        
           do jnode = 1,tnode,4
              elrhs(jnode)   = (0.0_rp,0.0_rp)
              elrhs(jnode+1) = (0.0_rp,0.0_rp)
              elrhs(jnode+2) = (0.0_rp,0.0_rp)
              elrhs(jnode+3) = (0.0_rp,0.0_rp)
              do inode = 1,tnode,4
                 elmat(inode,   jnode) = (0.0_rp,0.0_rp)
                 elmat(inode+1, jnode) = (0.0_rp,0.0_rp)
                 elmat(inode+2, jnode) = (0.0_rp,0.0_rp)
                 elmat(inode+3, jnode) = (0.0_rp,0.0_rp)
              end do
              do inode = 1,tnode,4
                 elmat(inode,   jnode+1) = (0.0_rp,0.0_rp)
                 elmat(inode+1, jnode+1) = (0.0_rp,0.0_rp)
                 elmat(inode+2, jnode+1) = (0.0_rp,0.0_rp)
                 elmat(inode+3, jnode+1) = (0.0_rp,0.0_rp)
              end do
              do inode = 1,tnode,4
                 elmat(inode,   jnode+2) = (0.0_rp,0.0_rp)
                 elmat(inode+1, jnode+2) = (0.0_rp,0.0_rp)
                 elmat(inode+2, jnode+2) = (0.0_rp,0.0_rp)
                 elmat(inode+3, jnode+2) = (0.0_rp,0.0_rp)
              end do
              do inode = 1,tnode,4
                 elmat(inode,   jnode+3) = (0.0_rp,0.0_rp)
                 elmat(inode+1, jnode+3) = (0.0_rp,0.0_rp)
                 elmat(inode+2, jnode+3) = (0.0_rp,0.0_rp)
                 elmat(inode+3, jnode+3) = (0.0_rp,0.0_rp)
              end do
           end do

        end if

        !!Assembly of element equations into system equations: AMATX and RHSIX
        !call hlm_assembdiff(1_ip,pnode,lnods(1,ielem),elmat,elrhs,indvars)

  end do elements

end subroutine hlm_elmopediffAssem

