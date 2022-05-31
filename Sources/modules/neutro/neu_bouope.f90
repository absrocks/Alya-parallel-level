!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_bouope.f90
!> @author  Guillaume Houzeaux
!> @brief   Boundary assembly
!> @details This assembles the boundary elements into the matrix of the system
!> @} 
!-----------------------------------------------------------------------

subroutine neu_bouope()

  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper 
  use def_domain
  use def_neutro
  use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix, only : matrix_assemble_element_RHS
  implicit none
  real(rp)    :: elmat(mnode,mnode) !> Elementary matrix
  real(rp)    :: elrhs(mnode) !> Elementary rhs
  real(rp)    :: baloc(ndime,ndime,mgaub) !> Local directions, normal to boundary element (baloc(1:ndime, ndime))  and two tangencial (baloc(1:ndime,1) and baloc(1:ndime,2) for 3D
  real(rp)    :: elcod(ndime,mnode) !> Element code
  real(rp)    :: bocod(ndime,mnodb) !> Boundary code
  real(rp)    :: borad(num_energies_neu,num_directions_neu,mnodb) !> Boundary radiation for each energy, direction and boundary node
  integer(ip) :: ielem,inode,ipoin,kfl_gobou,pgaub
  integer(ip) :: igaus,igaub,iboun,inodb,pblty,idime
  integer(ip) :: pnodb,pmate,pnode,pelty,pgaus
  integer(ip) :: dummi,idirection,ienergy
  real(rp)    :: eucta,tmatr,gbsur(mgaub),gpdet,adotn !> Eucta?? Tangent space??, not used, surface, determinant, scalar product between velocity and normal?
  real(rp)    :: gbsph,gbden,gbcon,gbtem,gbvel(3)
  real(rp)    :: gpsph(mgaus),gpden(mgaus),gpdif(mgaus) !> These variables are unused
  real(rp)    :: gprea(mgaus) !> Unused as well
  real(rp)    :: gptur(mgaus)                                 ! Turbulent viscosity (unused)
  real(rp)    :: arobi,trobi,qrobi,twall,acvis !> (unused)
  real(rp)    :: para1,para2,para3,para4 !> (unused)
  real(rp)    :: xmrhs,xmmat !> (unused)
  real(rp)    :: eledd(mnode),gpcar(ndime,mnode,mgaus) !> (unused only the first one)
  real(rp)    :: gpcon(mgaus),dummr(ndime*mnode),gpcod, gpvis(mgaus) !> (unused)
  real(rp)    :: xjaci(9),xjacm(9)
  !
  ! Loop over elements of the boundary
  !
  boundaries: do iboun = 1,nboun

     if( kfl_fixbo_neu(current_energy_neu,current_direction_neu) % l(iboun) == 1 ) then

        pblty = ltypb(iboun) !> We put the boundary type of the current boundary element in pblty
        pnodb = nnode(pblty) !> From the boundary type we get the number of nodes in the boundary
        ielem = lelbo(iboun) !> We determine the identifier for the current element touching the boundary
        pelty = ltype(ielem) !> We determine the type of the current element

        if( pelty > 0 ) then !> If the element is nonzero

           pnode = nnode(pelty) !> Number of nodes in current element
           pgaus = ngaus(pelty) !> Number of Gauss points in current element
           pgaub = ngaus(pblty) !> Number of Gauss points in the boundary
           pmate = 1 !> We consider only one material
           if( nmate > 1 ) pmate = lmate(ielem) !> If there is more than one material we determine the current element's
           !
           ! Inititalize: setting the matrix and the vector column to zero
           !
           elmat = 0.0_rp 
           elrhs = 0.0_rp
           !
           ! Gather operations
           !
           bocod(1:ndime,1:pnodb) = coord(1:ndime,lnodb(1:pnodb,iboun)) !> Boundary code determination
           elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem)) !> Element code determination
           borad(1:num_energies_neu,1:num_directions_neu,1:pnodb) &
                = neutr(1:num_energies_neu,1:num_directions_neu,lnodb(1:pnodb,iboun),1) !> We equate the flux in the energies and directions and nodes of the boundary nodes to the boundary radiation
           !
           ! Cartesian derivatives
           !
           do igaus = 1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty) % deriv(1,1,igaus),&      ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian (elmar is the element data base matrix)
           end do
           do igaub = 1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % deriv(1,1,igaub),&
                   bocod,baloc(1,1,igaub),eucta) !> This routine calculates the baloc (normal and tangencial directions in the element) and eucta
              gbsur(igaub) = elmar(pblty) % weigp(igaub) * eucta !> Here we calculate the surface for the current Gauss point, CHECK THIS
              call chenor(pnode,baloc(1,1,igaub),bocod,elcod)          ! Check normal
           end do
           !
           ! Element matrix and rhs
           !
           call neu_boumat(&
                pnodb,pgaub,pnode,lboel(1,iboun),elmar(pblty) % shape,gbsur,baloc,borad,&
                bvnat_neu(current_energy_neu,current_direction_neu) % a(1,iboun,1),&
                elmat,elrhs) !> We call this subroutine to calculate and implement the physical boundary conditions
           !
           ! Assembly
           !
           call matrix_assemble_element_RHS(&
                1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid) !> We assemble the RHS into the system RHS
           call matrix_assemble_element_matrix_to_CSR(&
                kfl_element_to_csr,1_ip,pnode,pnode,&
                ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo) !> We assemble the elementary matrix into the system matrix

        end if

     end if

  end do boundaries

end subroutine neu_bouope
