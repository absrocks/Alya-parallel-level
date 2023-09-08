subroutine membcs(order)
  !-----------------------------------------------------------------------
  !****f* Domain/membcs
  ! NAME
  !    membcs
  ! DESCRIPTION
  !    Allocate/Deallocate memory for node and boundary codes
  ! USED BY
  !    reaset
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip
  use def_domain, only : kfl_codno,kfl_codbo
  use def_domain, only : kfl_geobo,kfl_geono
  use def_domain, only : kfl_icodn,kfl_icodb
  use def_domain, only : kfl_coded
  use def_domain, only : npoin,nboun,nboun_2
  use def_domain, only : nbopo,mcono,memor_dom
  use def_domain, only : meshe
  use def_kermod, only : ndivi
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none
  integer(ip), intent(in) :: order

  select case (order) 

  case (1_ip)

     if( kfl_icodn > 0 ) then 
        call memory_alloca(memor_dom,'KFL_CODNO','membcs',kfl_codno,mcono,npoin)
        !do ipoin = 1,npoin
        !   kfl_codno(:,ipoin) = mcodb+1
        !end do
     end if

  case (-1_ip)

     if( kfl_icodn > 0 ) then
        call memory_deallo(memor_dom,'KFL_CODNO','membcs',kfl_codno)
     end if

  case (2_ip)

     if( kfl_icodb > 0 ) then
        call memory_alloca(memor_dom,'KFL_CODBO','membcs',kfl_codbo,max(1_ip,nboun))
     end if

  case (-2_ip) 

     if( kfl_icodb > 0 ) then
        call memory_deallo(memor_dom,'KFL_CODBO','membcs',kfl_codbo)
     end if

  case (3_ip)

     if( kfl_icodn > 0 ) then
        call memory_deallo(memor_dom,'KFL_CODNO','membcs',kfl_codno)
     end if

     if( kfl_icodb > 0 ) then 
        call memory_deallo(memor_dom,'KFL_CODBO','membcs',kfl_codbo)
     end if

  case (4_ip)

     if( kfl_icodb > 0 ) then
        call memory_alloca(memor_dom,'KFL_GEOBO','membcs',kfl_geobo,nboun_2)
     end if

  case (5_ip)

     call memory_alloca(memor_dom,'KFL_GEONO','membcs',kfl_geono,nbopo)

  case (7_ip)

     if( kfl_icodb > 0 ) then
        call memory_alloca(memor_dom,'KFL_CODBO','membcs',kfl_codbo,max(1_ip,nboun))
     else
        call memory_alloca(memor_dom,'KFL_CODBO','membcs',kfl_codbo,1_ip)     
     end if

  case (-7_ip)

     call memory_deallo(memor_dom,'KFL_CODBO','membcs',kfl_codbo) 

  case ( 12_ip)
     !
     ! Allocate edge codes KFL_CODED
     !
     call memory_alloca(memor_dom,'KFL_CODED','membcs',kfl_coded,2_ip,meshe(ndivi) % nedge)

  case (-12_ip)
     !
     ! Deallocate edge codes KFL_CODED
     !
     call memory_deallo(memor_dom,'KFL_CODED','membcs',kfl_coded)

  end select

end subroutine membcs
