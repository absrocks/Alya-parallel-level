subroutine cderda()
  !-----------------------------------------------------------------------
  !****f* Domain/cderda
  ! NAME
  !    cderda
  ! DESCRIPTION
  !    This routine defines the derivated parameters of the element
  ! OUTPUT
  !    NDIMB ... Boundary dimension
  !    NTENS ... # of independent Hessian matrix values
  !    NINER ... # of independent tensor of inertia values
  !    LRULE ... 1 Quad/Hexa:  open   
  !          ... 2 Quad/Hexa:  closed 
  !          ... 3 Tria/Tetra: open
  !          ... 4 Tria/Tetra: closed
  !          ... 5   - /Penta: open
  !          ... 6   - /Penta: closed
  !          ... 7   - /Pyram: open
  !          ... 8   - /Pyram: closed
  !    HNATU ... 2 Quad/Hexa
  !          ... 1 Others
  !    MNODE ... Maximum # of element nodes in the mesh
  !    MGAUS ... Maximum # of element Gauss points in the mesh
  !    MNODB ... Maximum # of boundary nodes in the mesh
  !    MGAUB ... Maximum # of boundary Gaus points in the mesh
  !    MLAPL ... 1 if at least one element needs Laplacian
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memory, only : memory_alloca
  use mod_elmgeo, only : elmgeo_element_type_retrocompatibility
  implicit none
  integer(ip) :: ielty,pface,pgaus,pquad,ndife
  !
  ! Where element and boundary element types start and stop
  ! We put IBSTA_DOM=1 to treat BARD3D elements, which boundary
  ! is the 1D POINT element
  !
  if( ndime == 1 ) then
     iesta_dom =  2
     iesto_dom =  9
     ibsta_dom =  1
     ibsto_dom =  1
  else if( ndime == 2 ) then
     iesta_dom = 10
     iesto_dom = 29
     ibsta_dom =  2
     ibsto_dom =  9
  else
     iesta_dom = 30
     iesto_dom = 60
     ibsta_dom = 10
     ibsto_dom = 29
  end if
  if( lexis(BAR3D) /= 0 ) ibsta_dom = 1
  !
  ! Gobal dimensions
  !
  ndimb = ndime-1
  if( ndime == 1 ) then
     ntens = 1
  else
     ntens = 3 * ndime - 3
  end if
  if( ndime == 3 ) then
     niner = 3
  else
     niner = 1
  end if
  !
  ! Compute LRULE and HNATU
  !
  do ielty = 1,nelty
     lrule(ielty) = 2*ltopo(ielty)+lquad(ielty)+1     
     if( ltopo(ielty) == -1 ) then
        hnatu(ielty) = 2.0_rp                     ! BAR
     else if( ltopo(ielty) == 0 ) then
        hnatu(ielty) = 2.0_rp                     ! QUA/HEX
     else  
        hnatu(ielty) = 1.0_rp                     ! OTHERS
     end if
  end do
  !
  ! Divide mesh
  !
  if( kfl_divid == 1 ) then
     ndife = 0
     do ielty = 1,nelty
        ndife = ndife + lexis(ielty)
     end do
     if( ndife /= 1 )        call runend('READIM: CANNOT DIVIDE MESH')
     if( lexis(HEX08) == 0 ) call runend('READIM: CAN ONLY DIVIDE HEXAHEDRA')
     lexis(TET04) = 1
  end if
  !
  ! If Volume Gauss points have not been assigned, put default option
  !
  do ielty = 1,nelty     
     if( lexis(ielty) /= 0 .and. ngaus(ielty) <= 0 ) then
        if( lquad(ielty) == 0 ) then
           !
           ! Open rule
           !
           if(      ielty == POINT ) then
              ngaus(ielty) = 1
           else if( ielty == BAR04 ) then
              ngaus(ielty) = 11
           else if( ielty == QUA08 ) then
              ngaus(ielty) = 9
           else if( ielty == QUA09 ) then
              ngaus(ielty) = 9
           else if( ielty == HEX20 ) then
              ngaus(ielty) = 27
           else if( ielty == HEX27 ) then
              ngaus(ielty) = 27
           else if( ielty == TET10 ) then
              ngaus(ielty) = 11
           else
              ngaus(ielty) = nnode(ielty) 
           end if
        else
           !
           ! Close rule
           !
           ngaus(ielty) = nnode(ielty)
        end if
     end if
  end do
  !
  ! Treat elements of dimension ndime-1: NGAUS, LQUAD and LEXIS
  !
  do ielty = iesta_dom,iesto_dom 
     if( lexis(ielty) == 1 ) then
        pgaus = ngaus(ielty)
        pquad = lquad(ielty)
        call bouele(nelty,pgaus,pquad,ielty,ngaus,lquad,lexis)
     end if
  end do
  !
  ! Maximum values
  ! MNODE can be optionally given in dimension field
  !
  mnode = max(mnode,-1_ip)
  mgaus = -1
  mnodb = -1
  mgaub = -1
  mlapl = -1
  do ielty = iesta_dom,iesto_dom
     if( lexis(ielty) == 1 ) then
        mnode = max(mnode,nnode(ielty))
        mgaus = max(mgaus,ngaus(ielty))
        mlapl = max(mlapl,llapl(ielty))
     end if
  end do
  do ielty = ibsta_dom,ibsto_dom
     if( lexis(ielty) == 1 ) then
        mnodb = max(mnodb,nnode(ielty))
        mgaub = max(mgaub,ngaus(ielty))
     end if
  end do
  if( lexis(SHELL) == 1 ) mnodb = max(2_ip,mnodb)
  if( lexis(BAR3D) == 1 ) mnodb = max(2_ip,mnodb)
  mnoga = max(mnode,mgaus)
  !
  ! Allocate memory for face list
  !
  mface = maxval(nface)
  do ielty = iesta_dom,iesto_dom
     pface = nface(ielty)
     if( lexis(ielty) == 1 .and. pface /= 0 ) then
        call memory_alloca(memor_dom,'LTYPF','cderda',ltypf(ielty) % l,pface)
        call memory_alloca(memor_dom,'NNODF','cderda',nnodf(ielty) % l,pface)
        call memory_alloca(memor_dom,'LFACE','cderda',lface(ielty) % l,mnodb,pface)
        call domfac(&
             ielty,mnodb,pface,lface(ielty)%l,&
             ltypf(ielty)%l,nnodf(ielty)%l)
     end if
  end do
  !
  ! Element data base, check retro-compatibility
  !
  call elmgeo_element_type_retrocompatibility(&
       iesta_dom,iesto_dom,nface,lexis,needg,leedg,&
       ldime,lorde,ltopo,nnode(1:),nnodf,ltypf,lface)
  !
  ! Check errors
  !  
  if( mnode == -1 ) call runend('DOMAIN: NO ELEMENT TYPE HAS BEEN DECLARED')
  if( mgaus == -1 ) call runend('DOMAIN: NO NUMERICAL INTEGRATION HAS BEEN DEFINED')

end subroutine cderda
