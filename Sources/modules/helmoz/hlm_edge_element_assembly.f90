!------------------------------------------------------------------------
!> @addtogroup HelmozMatrixAssembly
!> @ingroup    Helmoz
!> @{
!> @file    hlm_edge_element_assembly.f90
!> @date    08/02/2016
!> @author  Guillaume Houzeaux
!> @brief   Assembly
!> @details Assembly
!>
!> @} 
!------------------------------------------------------------------------
subroutine hlm_edge_element_assembly()

  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use def_helmoz
  implicit none

  real(rp)    :: Ke(meshe(ndivi) % medge,meshe(ndivi) % medge) 
  real(rp)    :: Me(meshe(ndivi) % medge,meshe(ndivi) % medge) 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: ellen(meshe(ndivi) % medge)
  real(rp)    :: elsig(meshe(ndivi) % medge)

  integer(ip) :: ielem,inode,pelty,pnode
  integer(ip) :: pedge,ipoin,iedgg,iedge
  !
  ! OpenMP declarations 
  ! 
  ! $OMP PARALLEL DO                                         &
  ! $OMP SCHEDULE      ( DYNAMIC , par_omp_nelem_chunk )     & 
  ! $OMP DEFAULT       ( NONE )                              &
  ! $OMP PRIVATE       (  )                                  &
  ! $OMP SHARED        (  )                                  &
  ! $OMP REDUCTION     ( +: )                                &
  ! $OMP REDUCTION     ( MIN: )                              &  
  ! $OMP REDUCTION     ( MAX: )                   
  ! 
  ! Loop over elements
  !     
  do ielem = 1,nelem

     pelty = ltype(ielem)
     pnode = lnnod(ielem)
     pedge = meshe(ndivi) % lnned(ielem)

     if( pelty > 0 ) then

        Ke = 0.0_rp
        Me = 0.0_rp
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           elcod(1:ndime,inode) = coord(1:ndime,ipoin)
        end do
        do iedge = 1,pedge
           iedgg = meshe(ndivi) % ledgs(iedge,ielem)
           ellen(iedge) = length_edges_hlm(iedgg)
        end do
        elsig(1:pedge) = sign_edges_hlm(1:pedge,ielem)

        if( pelty == TET04 ) then
           !
           ! Tetrahedra
           !
           call elemental_matrices(elcod,ellen,elsig,Ke,Me)
           do iedge = 1,pedge
              write(6,'(6(1x,e12.6))') Ke(iedge,1:pedge)
           end do
           call runend('O.K.!')
        else
           call runend('ELEMENT NOT CODED')
        end if

     end if
  end do  
  ! $OMP END PARALLEL DO


end subroutine hlm_edge_element_assembly

subroutine elemental_matrices(eleNodes,lengthE,signsEle,Ke,Me)
  use def_kintyp, only : ip,rp
  implicit none
  ! Functions
  real(rp)                                               :: m33Det
  ! Parameters
  integer(ip), parameter                                 :: orderEle = 6
  ! Input 
  real(rp),    intent(in),  dimension(3,4)               :: eleNodes
  real(rp),    intent(in),  dimension(orderEle)          :: lengthE, signsEle
  ! Output
  real(rp),    intent(out), dimension(orderEle,orderEle) :: Me, Ke
  ! Variables
  real(rp),    dimension(4)                              :: b, c, d, f1, f2, f3, f4
  integer(ip), dimension(7)                              :: tmp = (/ 1, 2, 3, 4, 1, 2, 3 /)
  real(rp),    dimension(3,3)                            :: tmpM
  integer(ip), dimension(orderEle,3)                     :: Nele
  real(rp)                                               :: A1, A2
  real(rp)                                               :: eleVol
  integer(ip)                                            :: i, j
  !
  ! Compute volume
  !
  eleVol = 1.0_rp

  ! Compute Coefficients
  do i = 1,4
     ! Coefficient b
     tmpM(1,:) = (/1, 1, 1/)
     tmpM(2,:) = (/eleNodes(2,tmp(i+1)), eleNodes(2,tmp(i+2)), eleNodes(2,tmp(i+3))/)
     tmpM(3,:) = (/eleNodes(3,tmp(i+1)), eleNodes(3,tmp(i+2)), eleNodes(3,tmp(i+3))/)
     b(i)      = m33Det(tmpM)
     ! Coefficient c
     tmpM(1,:) = (/1, 1, 1/)
     tmpM(2,:) = (/eleNodes(1,tmp(i+1)), eleNodes(1,tmp(i+2)), eleNodes(1,tmp(i+3))/)
     tmpM(3,:) = (/eleNodes(3,tmp(i+1)), eleNodes(3,tmp(i+2)), eleNodes(3,tmp(i+3))/)
     c(i)      = m33Det(tmpM)
     ! Coefficient d
     tmpM(1,:) = (/1, 1, 1/)
     tmpM(2,:) = (/eleNodes(1,tmp(i+1)), eleNodes(1,tmp(i+2)), eleNodes(1,tmp(i+3))/)
     tmpM(3,:) = (/eleNodes(2,tmp(i+1)), eleNodes(2,tmp(i+2)), eleNodes(2,tmp(i+3))/)
     d(i)      = m33Det(tmpM)
  end do

  ! Add signs
  b(1) = b(1)*(-1.0_rp)
  b(3) = b(3)*(-1.0_rp)
  c(2) = c(2)*(-1.0_rp)
  c(4) = c(4)*(-1.0_rp)
  d(1) = d(1)*(-1.0_rp)
  d(3) = d(3)*(-1.0_rp)

  ! Global notation of edges within each elemental_matrices
  ! [edge node_i node_j]
  Nele(:,1) = (/1, 2, 3, 4, 5, 6/)
  Nele(:,2) = (/1, 1, 1, 2, 4, 3/)
  Nele(:,3) = (/2, 3, 4, 3, 2, 4/)

  do i=1, 4
     f1(i) = b(1)*b(i) + c(1)*c(i) + d(1)*d(i)
     f2(i) = b(2)*b(i) + c(2)*c(i) + d(2)*d(i)
     f3(i) = b(3)*b(i) + c(3)*c(i) + d(3)*d(i)
     f4(i) = b(4)*b(i) + c(4)*c(i) + d(4)*d(i)
  end do

  A1 = (360_rp*eleVol)**(-1)
  A2 = (720_rp*eleVol)**(-1)
  ! Mass matrix computation
  Me(1,1) = lengthE(1)**2 * A1 * (f2(2) - f1(2) + f1(1))
  Me(1,2) = lengthE(1) * lengthE(2) * A2 * (2*f2(3)-f1(2)-f1(3)+f1(1))
  Me(1,3) = lengthE(1) * lengthE(3) * A2 * (2*f2(4)-f1(2)-f1(4)+f1(1))
  Me(1,4) = lengthE(1) * lengthE(4) * A2 * (f2(3)-f2(2)-2*f1(3)+f1(2))
  Me(1,5) = lengthE(1) * lengthE(5) * A2 * (f2(2)-f2(4)-f1(2)+2*f1(4))
  Me(1,6) = lengthE(1) * lengthE(6) * A2 * (f2(4)-f2(3)-f1(4)+f1(3))
  Me(2,2) = lengthE(2)**2 * A1 * (f3(3) - f1(3) + f1(1))
  Me(2,3) = lengthE(2) * lengthE(3) * A2 * (2*f3(4)-f1(3)-f1(4)+f1(1))
  Me(2,4) = lengthE(2) * lengthE(4) * A2 * (f3(3)-f2(3)-f1(3)+2*f1(2))
  Me(2,5) = lengthE(2) * lengthE(5) * A2 * (f2(3)-f3(4)-f1(2)+f1(4))
  Me(2,6) = lengthE(2) * lengthE(6) * A2 * (f1(3)-f3(3)-2*f1(4)+f3(4))
  Me(3,3) = lengthE(3)**2 * A1 * (f4(4) - f1(4) + f1(1))
  Me(3,4) = lengthE(3) * lengthE(4) * A2 * (f3(4)-f2(4)-f1(3)+f1(2))
  Me(3,5) = lengthE(3) * lengthE(5) * A2 * (f2(4)-f4(4)-2*f1(2)+f1(4))
  Me(3,6) = lengthE(3) * lengthE(6) * A2 * (f4(4)-f3(4)-f1(4)+2*f1(3))
  Me(4,4) = lengthE(4)**2 * A1 * (f3(3) - f2(3) + f2(2))
  Me(4,5) = lengthE(4) * lengthE(5) * A2 * (f2(3)-2*f3(4)-f2(2)+f2(4))
  Me(4,6) = lengthE(4) * lengthE(6) * A2 * (f3(4)-f3(3)-2*f2(4)+f2(3))
  Me(5,5) = lengthE(5)**2 * A1 * (f2(2) - f2(4) + f4(4))
  Me(5,6) = lengthE(5) * lengthE(6) * A2 * (f2(4)-2*f2(3)-f4(4)+f3(4))
  Me(6,6) = lengthE(6)**2 * A1 * (f4(4) - f3(4) + f3(3))

  do i=1, orderEle
     do j=1, orderEle
        Me(j,i) = Me(i,j)
     end do
  end do

  ! Add signs of edges
  do i=1, orderEle
     do j=1, orderEle
        Me(i,j) = Me(i,j) * signsEle(i) * signsEle(j)
     end do
  end do

  ! Stiffness matrix computation
  do i=1, orderEle
     do j=1, orderEle
        Ke(i,j) = lengthE(i) * lengthE(j) * &
             ( (c(NEle(i,2))*d(NEle(i,3)) - d(NEle(i,2))*c(NEle(i,3))) * (c(NEle(j,2))*d(NEle(j,3)) - d(NEle(j,2))*c(NEle(j,3))) &
             + (d(NEle(i,2))*b(NEle(i,3)) - b(NEle(i,2))*d(NEle(i,3))) * (d(NEle(j,2))*b(NEle(j,3)) - b(NEle(j,2))*d(NEle(j,3))) &
             + (b(NEle(i,2))*c(NEle(i,3)) - c(NEle(i,2))*b(NEle(i,3))) * (b(NEle(j,2))*c(NEle(j,3)) - c(NEle(j,2))*b(NEle(j,3))) ) &
             * signsEle(i) * signsEle(j)
     end do
  end do

  Ke = Ke * (4.0_rp*eleVol)/((6.0_rp*eleVol)**4)

end subroutine elemental_matrices

! Compute determinant
function m33Det (M) result (det)
  use def_kintyp, only : ip,rp

  implicit none
  
  real(rp), dimension(3,3), intent(in)  :: M
  real(rp) :: det
  
  det =    M(1,1)*M(2,2)*M(3,3)  &
       & - M(1,1)*M(2,3)*M(3,2)  &
       & - M(1,2)*M(2,1)*M(3,3)  &
       & + M(1,2)*M(2,3)*M(3,1)  &
       & + M(1,3)*M(2,1)*M(3,2)  &
       & - M(1,3)*M(2,2)*M(3,1)
  
  return
  
end function m33Det
