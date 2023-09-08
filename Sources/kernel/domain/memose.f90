subroutine memose(order)
  !-----------------------------------------------------------------------
  !****f* Domain/memose
  ! NAME
  !    memose
  ! DESCRIPTION
  !    Allocate memory for set arrays
  ! USED BY
  !    reaset
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_inpout
  use mod_memory
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: order

  select case (order) 

  case (1_ip)
     if( neset > 0 ) then
        call memory_alloca(memor_dom,'LESET','memose',leset,nelem)
     end if

  case (-1_ip)
     if( neset > 0 ) then
        call memory_deallo(memor_dom,'LESET','memose',leset)
     end if

  case (2_ip)
     if( nbset > 0 ) then
        call memory_alloca(memor_dom,'LBSET','memose',lbset,nboun)
     end if

  case (-2_ip)
     if( nbset > 0 ) then
        call memory_deallo(memor_dom,'LBSET','memose',lbset)
     end if

  case (3_ip)
     if( nnset > 0 ) then
        call memory_alloca(memor_dom,'LNSEC','memose',lnsec,nnset)
     end if

  case (4_ip)
     if( nbset > 0 ) then
        call memory_alloca(memor_dom,'LBSEC','memose',lbsec,nbset)
     end if

  case (5_ip)
     if( neset > 0 ) then
        call memory_alloca(memor_dom,'LESEC','memose',lesec,neset)
     end if

  case (-5_ip)
     if( neset > 0 ) then
        call memory_deallo(memor_dom,'LESEC','memose',lesec)
     end if

  case (6_ip)
     if( mwitn > 0 ) then
        call memory_alloca(memor_dom,'SHWIT','memose',shwit,mnode,mwitn)
        call memory_alloca(memor_dom,'DEWIT','memose',dewit,ndime,mnode,mwitn)
        call memory_alloca(memor_dom,'LEWIT','memose',lewit,mwitn)
     end if

  case (7_ip)
     if( mwitn > 0 ) then
        call memory_alloca(memor_dom,'COWIT','memose',cowit,3_ip,mwitn)
     end if

  case (10_ip)
     !
     ! Deallocate memory
     !
     if( neset > 0 ) then
        call memory_deallo(memor_dom,'LESET','memose',leset)
     end if
     if( nbset > 0 ) then
        call memory_deallo(memor_dom,'LBSET','memose',lbset)
     end if
     if( nnset > 0 ) then
        call memory_deallo(memor_dom,'LNSET','memose',lnset)
     end if

  case (11_ip)
     !
     ! LNWIT
     !
     call memory_alloca(memor_dom,'LNWIT','memose',lnwit,nwitn)

  case (12_ip)
     !
     ! LNSET
     !
     if( nnset > 0 ) then
        call memory_alloca(memor_dom,'LNSET','memose',lnset,npoin)
     end if

  case (-12_ip)
     !
     ! LNSET
     !
     if( nnset > 0 ) then
        call memory_deallo(memor_dom,'LNSET','memose',lnset)
     end if

  end select

end subroutine memose
