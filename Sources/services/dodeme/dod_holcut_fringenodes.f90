!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut_fringenodes.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Identify fringe nodes
!> @details Identify fringe node, only when a hole has been created
!>
!>          \verbatim
!>          Input:
!>          ------
!>                     Free    hole  fringe
!>          Chimera      0   -jsubd  -jsubd
!>          Patch        0      x     jsubd
!>          Prescribed   0      x     jsubd
!> 
!>          Output:
!>          -------
!>                     Free    hole  fringe
!>          Chimera      0   -jsubd   jsubd
!>          Patch        0      x     jsubd
!>          Prescribed   0      x     jsubd
!>          \endverbatim
!>
!> @} 
!-----------------------------------------------------------------------
subroutine dod_holcut_fringenodes()
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use def_kermod
  use mod_kdtree
  use mod_graphs
  implicit none
  integer(ip) :: isubd,ipoin,inodb,pnodb,iboun,npoin_subd
  !
  ! LNSUB_NPOIN(IPOIN) < 0 for holes
  ! Loop over extension boundaries to set it to positive
  !
  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)

     if( ihole_dod(isubd) == 1 ) then
        npoin_subd = current_subdomain % npoin
        do iboun = 1,current_subdomain % nboun
           if( current_subdomain % lboch(iboun) == BOEXT ) then
              pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
              do inodb = 1,pnodb
                 ipoin = current_subdomain % lnodb(inodb,iboun)
                 current_subdomain % lsubd_npoin(ipoin) = abs(current_subdomain % lsubd_npoin(ipoin))
              end do
           end if
        end do

     end if

  end do

end subroutine dod_holcut_fringenodes
 
