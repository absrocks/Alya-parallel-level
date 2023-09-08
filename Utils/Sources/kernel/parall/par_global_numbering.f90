!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_global_numbering.f90
!> @author  houzeaux
!> @date    2019-01-02
!> @brief   Global numbering
!> @details Global numbering
!> @} 
!-----------------------------------------------------------------------

subroutine par_global_numbering()

  use def_master
  use def_domain
  use def_parall
  use mod_memory
  use mod_parall,         only : par_memor
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_ALLGATHER
  use mod_messages,       only : messages_live
  implicit none

  integer(ip)          :: ipoin,ipart,npoin_offset
  integer(ip), pointer :: npoin_gat(:)

  nullify(npoin_gat)

  if( kfl_global_numbering_par == 0 ) then
     !
     ! Do nothing, take the global numbering given by the mesh generator
     !
  else if( kfl_global_numbering_par == 1 ) then
     !
     ! Lexical numbering
     !
     call messages_live('COMPUTE GLOBAL NUMBERING')
     call memory_alloca(par_memor,'NPOIN_GAT','par_global_numbering',npoin_gat,npart+1_ip,LBOUN=0_ip)
     call PAR_ALLGATHER(npoin_own,npoin_gat,1_4,'IN MY CODE')
     npoin_offset = 0
     do ipart = 1,kfl_paral-1
        npoin_offset = npoin_offset + npoin_gat(ipart)
     end do
     if(INOTMASTER ) lninv_loc = 0_ip
     do ipoin = 1,npoin_own
        lninv_loc(ipoin) = ipoin + npoin_offset
     end do
     call memory_deallo(par_memor,'NPOIN_GAT','par_global_numbering',npoin_gat)
     call PAR_INTERFACE_NODE_EXCHANGE(lninv_loc,'SUM','IN MY CODE')

  end if

end subroutine par_global_numbering
