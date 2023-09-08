subroutine dod_inivar(order)
  !------------------------------------------------------------------------
  !****f* Dodeme/dod_inivar
  ! NAME
  !    dod_inivar
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Dodeme
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master 
  use def_domain
  use def_dodeme
  use def_inpout
  implicit none
  integer(ip),  intent(in)  :: order
  integer(ip)               :: isubd,ipoin,idime,jsubd

  select case ( order )

  case ( 1_ip )
     !
     ! Read and write files
     !
     lispa     = 0
     lisda     = lun_pdata_dod                    ! Reading file
     lisre     = lun_outpu_dod                    ! Writing file
  
  case ( 2_ip )
     !
     ! Define if there are chimera process
     !
     ihole_dod = 0
     ipatc_dod = 0
     ipres_dod = 0
     do isubd = 1,nsubd
        do jsubd = 1,nsubd
           if( intyp_dod(isubd,jsubd) /= 0 ) then
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_CHIMERA     ) ihole_dod(isubd) = 1
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_PATCH       ) ipatc_dod(isubd) = 1 
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_PRESCRIBED  ) ipres_dod(isubd) = 1 
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_HOLED_PATCH ) ipatc_dod(isubd) = 1 
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_HOLED_PATCH ) ihole_dod(isubd) = 1 
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_HOLED_PATCH ) ipatc_dod(isubd) = 1 
              if( ictop_dod(intyp_dod(isubd,jsubd)) == DOD_HOLED_PATCH ) ihole_dod(isubd) = 1 
           end if
        end do
     end do
     !
     ! Embedding boxes
     !
     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)
        do idime = 1,ndime
           current_subdomain % embox(idime,1) =  1e9_rp
           current_subdomain % embox(idime,2) = -1e9_rp
        end do
        do ipoin = 1,current_subdomain % npoin
           do idime = 1,ndime
              current_subdomain % embox(idime,1) = min( current_subdomain % embox(idime,1) , current_subdomain % coord(idime,ipoin) )
              current_subdomain % embox(idime,2) = max( current_subdomain % embox(idime,2) , current_subdomain % coord(idime,ipoin) )
           end do
        end do
     end do

  end select

end subroutine dod_inivar
 
