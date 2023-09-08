!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_kdtree.f90
!> @author  Guillaume Houzeaux
!> @date    27/02/2013
!> @brief   Compute KD-Tree for patch meshes
!> @details Compute KD-Tree for patch meshes
!> @} 
!-----------------------------------------------------------------------

subroutine dod_kdtree(itask)
  use def_parame
  use def_master
  use def_domain
  use def_dodeme
  use mod_memory
  use mod_kdtree
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: isubd
  integer(ip)             :: ipoin,idime,iboun,inodb 
  integer(ip)             :: pnodb,ipoin_bou
  !
  ! Compute KD-Tree for patch meshes
  !
  if( INOTSLAVE ) then

     call livinf(0_ip,'CONSTRUCT SKD-TREE STRUCTURE',0_ip)  

     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)

        if( itask == 1 .or. ( itask == 2 .and. ihole_dod(isubd) /= 0 ) ) then
           !
           ! Allocate and boundary as if the boundary mesh were indpendent 
           ! LNODB_BOU(INODB,IBOUN) ....... IPOIN_BOU in boundary numbering
           ! COORD_BOU(IDIME,IPOIN_BOU) ... Node coordinates
           !
           call memgen(1_ip,current_subdomain % npoin,0_ip)
           current_subdomain % npoin_bou = 0
           do iboun = 1,current_subdomain % nboun
              pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
              do inodb = 1,pnodb
                 ipoin = current_subdomain % lnodb(inodb,iboun)
                 if( gisca(ipoin) == 0 ) then
                    current_subdomain % npoin_bou = current_subdomain % npoin_bou + 1
                    gisca(ipoin) = current_subdomain % npoin_bou
                 end if
              end do
           end do
           call memory_alloca(mem_servi(1:2,servi),'COORD_BOU','dod_kdtree',current_subdomain % coord_bou,ndime,current_subdomain % npoin_bou,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_servi(1:2,servi),'LNODB_BOU','dod_kdtree',current_subdomain % lnodb_bou,mnodb,current_subdomain % nboun,    'DO_NOT_INITIALIZE')
           do ipoin = 1,current_subdomain % npoin
              ipoin_bou = gisca(ipoin)
              if( ipoin_bou > 0 ) then
                 do idime = 1,ndime
                    current_subdomain % coord_bou(idime,ipoin_bou) = current_subdomain % coord(idime,ipoin)
                 end do
              end if
           end do
           do iboun = 1,current_subdomain % nboun
              pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
              do inodb = 1,pnodb
                 ipoin = current_subdomain % lnodb(inodb,iboun)
                 ipoin_bou = gisca(ipoin)
                 current_subdomain % lnodb_bou(inodb,iboun) = ipoin_bou
              end do
           end do

           call memgen(3_ip,current_subdomain % npoin,0_ip)
           !
           ! Allocte and coompute SKD-Tree structure
           !
           call kdtree(&
                1_ip,mnodb,&
                current_subdomain % npoin_bou , current_subdomain % nboun , current_subdomain % coord_bou , &
                current_subdomain % lnodb_bou , current_subdomain % ltypb , current_subdomain % fabox     , &
                current_subdomain % bobox     , current_subdomain % sabox , current_subdomain % blink     , &
                current_subdomain % stru2     , current_subdomain % ldist , current_subdomain % lnele     )

        else if( itask == 3 ) then
           !
           ! Deallocate everything
           !
           call memory_deallo(mem_servi(1:2,servi),'COORD_BOU','dod_kdtree',current_subdomain % coord_bou)
           call memory_deallo(mem_servi(1:2,servi),'LNODB_BOU','dod_kdtree',current_subdomain % lnodb_bou)
           call kdtree(&
                2_ip,mnodb,&
                current_subdomain % npoin_bou , current_subdomain % nboun , current_subdomain % coord_bou, &
                current_subdomain % lnodb_bou , current_subdomain % ltypb , current_subdomain % fabox,     &
                current_subdomain % bobox     , current_subdomain % sabox , current_subdomain % blink,     &
                current_subdomain % stru2     , current_subdomain % ldist , current_subdomain % lnele      ) 
           
        end if

     end do

  end if

end subroutine dod_kdtree
 
