!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_trabcs.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Correct bc for trailing edges
!> @details Correct bc for trailing edges
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_trabcs(pnode,lnode,elrhs,elmat)
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_memchk

  use def_nastal
  implicit none

  integer(ip)  :: &
       pnode,inode,jnode,jdofn,idofn,ievat,jevat,lnode(pnode),ipoin,icorrect,lcorrect(pnode)
  real(rp)     :: elrhs(nevat_nsa),diago(nevat_nsa),elmat(nevat_nsa,nevat_nsa)

  icorrect= 0
  lcorrect= 0
  do inode= 1,pnode
     ipoin= lnode(inode)
     if (kfl_fixno_nsa(1,ipoin) == 4) then
        icorrect= icorrect + 1
        lcorrect(inode)= ipoin
        do idofn=1,ndime
           ievat= (inode-1) * ndofn_nsa + idofn
           elrhs(ievat) = 0.0_rp
        end do
     end if
  end do

  if (icorrect > 0) then
     do inode= 1,pnode        
        if (lcorrect(inode) > 0) then
           do jnode= 1,pnode
              if (lcorrect(jnode) == 0) then
                 do jdofn=1,ndime 
                    jevat= (jnode-1) * ndofn_nsa + jdofn
                    ievat= (inode-1) * ndofn_nsa + jdofn
                    elrhs(ievat) = elrhs(ievat) + elrhs(jevat) / real(pnode - icorrect)
                 end do
              end if
           end do
        end if
     end do
  end if

end subroutine nsa_trabcs
