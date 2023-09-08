!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut_inversehole.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Finish hole cutting taking into account solid
!> @details Finish hole cutting:
!!
!!          1. Look for solid elements in a recursive way
!!
!!          Elements o are hole. Element x are solid, because their
!!          nodes were identified as free (they fall inside a solid
!!          of the patch). The procedure to declare them as holes 
!!          is to mark free element recursively from the background boundary.
!!
!!           +----------------------------------+
!!           |                                  |
!!           |         o o o o o o o            |
!!           |         o o o o o o o            |
!!           |         o o x x o o o            |
!!           |         o o x x o o o            |
!!           |         o o o o o o o            |
!!           |                                  |
!!           +----------------------------------+
!!
!!                          ||
!!                          \/
!!
!!           +----------------------------------+
!!           | . . . . . . . . . . . . . . . . .|
!!           | . . . . o o o o o o o . . . . . .|
!!           | . . . . o o o o o o o . . . . . .|
!!           | . . . . o o o o o o o . . . . . .|
!!           | . . . . o o o o o o o . . . . . .|
!!           | . . . . o o o o o o o . . . . . .|
!!           | . . . . . . . . . . . . . . . . .|
!!           +----------------------------------+
!!
!!           2. Define fringe nodes
!!
!!           +----------------------------------+
!!           |                                  |
!!           |         f f f f f f f            |
!!           |         f           f            |
!!           |         f           f            |
!!           |         f           f            |
!!           |         f f f f f f f            |
!!           |                                  |
!!           +----------------------------------+
!!
!> @} 
!-----------------------------------------------------------------------

subroutine dod_holcut_inversehole()

  use def_parame
  use def_master
  use def_domain
  use def_dodeme
  use mod_memchk
  use mod_memory
  use mod_postpr
  use mod_messages, only : livinf
  implicit none
  integer(ip)          :: isubd,iboun,nelem_subd,ielem,ielel,jelem
  integer(ip)          :: istack,nstack,pnodb,inode,ipoin,pnode
  integer(ip), pointer :: lstack(:)
  
  nullify(lstack)

  call livinf(0_ip,'MARK INVERSE HOLE TO DETECT SOLIDS',0_ip)

  do isubd = 1,nsubd

     if( ihole_dod(isubd) /= 0 ) then

        !----------------------------------------------------------------
        !
        ! ISUBD has a Chimera interface
        !
        !----------------------------------------------------------------

        current_subdomain  => subdomain(isubd)
        nelem_subd         =  current_subdomain % nelem

        !----------------------------------------------------------------
        !
        ! Mark all free elements by stopping at the border, whenever the
        ! border has no hole node
        ! GISCA(IELEM) = 1 ... Element was marked because it is free
        !
        !----------------------------------------------------------------

        call memory_alloca(mem_servi(1:2,servi),'LSTACK','dod_holcut_inversehole',lstack,nelem_subd)
        call memgen(1_ip,nelem_subd,0_ip) 

        do ielem = 1,nelem_subd
           lstack(ielem) = 0
           gisca(ielem)  = 0 
        end do

        do iboun = 1,current_subdomain % nboun
           pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
           ielem = current_subdomain % lelbo(iboun) 
           if( gisca(ielem) == 0 ) then
              !
              ! Check if element has all its nodes free
              !
              pnode = current_subdomain % lnnod(ielem)
              inode = 0
              do while( inode < pnode )
                 inode = inode + 1
                 ipoin = current_subdomain % lnods(inode,ielem)
                 if( current_subdomain % lsubd_npoin(ipoin) < 0 ) inode = pnode + 1
              end do

              if( inode == pnode ) then
              nstack       = 1
              lstack(1)    = ielem
              gisca(ielem) = 1
              istack       = 0                
              do 
                 if( istack == nstack ) exit
                 istack = istack+1   
                 ielem  = lstack(istack)
                 do ielel = current_subdomain % pelel(ielem), current_subdomain % pelel(ielem+1)-1
                    jelem = current_subdomain % lelel(ielel)
                    if( gisca(jelem) == 0 ) then
                       if( current_subdomain % lsubd_nelem(jelem) == 0 ) then
                          gisca(jelem)   = 1
                          nstack         = nstack + 1
                          lstack(nstack) = jelem
                       end if
                    end if
                 end do
              end do
           end if

           end if
        end do
        
        !----------------------------------------------------------------
        !
        ! Mark all non-free elements as candidate hole elements
        !
        !----------------------------------------------------------------

        do ielem = 1,nelem_subd
           if( gisca(ielem) == 0 ) then
              !
              ! Element was not marked
              !
              if( current_subdomain % lsubd_nelem(ielem) == 0 ) then
                 current_subdomain % lsubd_nelem(ielem) = DOD_SOLID                 
              end if
           else
              !
              ! Element was marked so it's free
              !
              current_subdomain % lsubd_nelem(ielem) = 0
           end if 
        end do

        call memgen(3_ip,nelem_subd,0_ip) 
        call memory_deallo(mem_servi(1:2,servi),'LSTACK','dod_holcut_inversehole',lstack)

     end if
  end do

end subroutine dod_holcut_inversehole




