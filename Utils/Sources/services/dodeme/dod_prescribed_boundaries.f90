 !-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_prescribed_boundaries.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Mark extension nodes for prescribed interfaces
!> @details For prescribed interfaces, read in dod_readat(), 
!>          assign the subdomain number to the interface node 
!>          PRESCRIBED_BOUNDARIES(IBOUN) 
!>          =>
!>          SUBDOMAIN(ISUBD) % LSUBD_NPOIN(IPOIN) = JSUBD
!>          SUBDOMAIN(ISUBD) % LBOCH(IBOUN) = BOEXT
!> @} 
!-----------------------------------------------------------------------
subroutine dod_prescribed_boundaries()
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_dodeme
  implicit none
  integer(ip) :: isubd,ipoin,inodb,pnodb,iboun,jsubd,iboun_global
  integer(ip) :: intij,ksubd

  do isubd = 1,nsubd
     if( ipres_dod(isubd) == 1 ) then
        current_subdomain => subdomain(isubd)
        do iboun = 1,current_subdomain % nboun
           iboun_global = current_subdomain % lbper(iboun)
           if( iboun_global == 0 ) then
              jsubd = current_subdomain % lsubd_nboun(iboun)
           else 
              jsubd = prescribed_boundaries(iboun_global)
           end if
           if( jsubd > 0 ) then
              intij = intyp_dod(isubd,jsubd)
              if( intij == 0 ) then
                 call runend('DOD_PRESCRIBED_BOUNDARIES: A PRESCRIBED BOUDNARIES DOES NOT BELONG TO AN INTERFACE')
              else
                 if( ictop_dod(intij) /= DOD_PRESCRIBED ) then
                    call runend('DOD_PRESCRIBED_BOUNDARIES: A PRESCRIBED BOUDNARIES DOES NOT BELONG TO A PRESCRIBED INTERFACE')
                 else
                    current_subdomain % lboch(iboun) = BOEXT  
                    pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
                    do inodb = 1,pnodb
                       ipoin = current_subdomain % lnodb(inodb,iboun)
                       current_subdomain % lsubd_npoin(ipoin) = jsubd
                    end do
                 end if
              end if
           end if
        end do

     end if
  end do

 
  return
  isubd = 2
  current_subdomain => subdomain(isubd)
  do iboun = 1,current_subdomain % nboun
     current_subdomain % lboch(iboun) = BOFEM
     pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
     do inodb = 1,pnodb
        ipoin = current_subdomain % lnodb(inodb,iboun)
        current_subdomain % lsubd_npoin(ipoin) = 0
     end do
  end do


end subroutine dod_prescribed_boundaries
