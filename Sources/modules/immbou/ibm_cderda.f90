subroutine ibm_cderda()
  !-----------------------------------------------------------------------
  !****f* Domain/ibm_cderda
  ! NAME
  !    ibm_cderda
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
  use mod_memchk
  use def_immbou
  implicit none
  integer(ip)              :: ielty,iimbo

  if( INOTSLAVE ) then
     !
     ! Default
     !
     do ielty = 1,nelty     
        if( lexib(ielty) /= 0 ) then
           if( ngaib(ielty) <= 0 ) then
              lquib(ielty) = 0               ! open rule
              if( ielty == BAR04 ) then
                 ngaib(ielty) = 11
              else if(ielty == QUA08 ) then
                 ngaib(ielty) = 9
              !else if( ielty == HEX20 ) then
              !   ngaib(ielty) = 27
              else
                 ngaib(ielty) = nnode(ielty) ! ngaus=nnode
              end if
           end if
        end if
     end do
  end if
  !
  ! Immersed boundary method: Maximum values and allocate memory
  ! Existence of Dirichlet and/or force IB
  !
  ielty = 0
  kfl_embed_ibm = 0                         ! Existence of Embedded mesh IB
  kfl_diric_ibm = 0                         ! Existence of Dirichlet IB
  kfl_force_ibm = 0                         ! Existence of Force IB
  if( INOTSLAVE ) then
     do iimbo = 1,nimbo
        if( imbou(iimbo) % kfl_typeb >= 1 .and. imbou(iimbo) % kfl_coupl == 0 ) then
           ielty = 1
        else
           kfl_embed_ibm = 1
        end if
        if( imbou(iimbo) % kfl_coupl == 0 .and. imbou(iimbo) % kfl_typeb <= 0 ) then
           kfl_diric_ibm = kfl_diric_ibm + 1
        else if( imbou(iimbo) % kfl_coupl == 1 ) then
           kfl_force_ibm = kfl_force_ibm + 1
        end if
     end do
  end if
  call parari('BCT',0_ip,1_ip,kfl_embed_ibm)
  call parari('BCT',0_ip,1_ip,kfl_diric_ibm)
  call parari('BCT',0_ip,1_ip,kfl_force_ibm)

  if( INOTSLAVE ) then
     if( ielty == 0 ) then
        mgaib = -1
        !do ielty = ibsta_dom,ibsto_dom
        do ielty = 1,nelty
           if( lexib(ielty) == 1 ) then
              mgaib = max(mgaib,ngaib(ielty))
           end if
        end do
        if( mnoib == -1 ) call runend('DOMAIN: NO IB TYPE HAS BEEN DECLARED')
        if( mgaib == -1 ) call runend('DOMAIN: NO IB NUMERICAL INTEGRATION HAS BEEN DEFINED')
        do ielty = 1,nelty
           lruib(ielty) = 2 * ltopo(ielty) + lquib(ielty) + 1     
        end do
     end if

  end if

end subroutine ibm_cderda
