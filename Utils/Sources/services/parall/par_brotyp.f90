subroutine par_brotyp(itask)
  !-------------------------------------------------------------------------------
  !****f* Parall/par_brotyp
  ! NAME
  !    par_brotyp
  ! DESCRIPTION
  !    Broadcast PAI1P(PARD1) array of type(i1p)
  ! INPUT
  !    PAI1P from master
  ! OUTPUT
  ! USED BY
  !    par_partit
  !***
  !-------------------------------------------------------------------------------
  use def_parame 
  use def_domain 
  use def_parall
  use def_master
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipari,nsize_loc,ii,kk
  integer(4)              :: istat
  integer(ip), pointer    :: lsize_loc(:),lai1p_loc(:)

  nparr = 0
  npari = 0
  nparc = 0

  select case(itask)

  case(1_ip)

     !-------------------------------------------------------------------
     ! 
     ! Type(1ip)
     !
     !-------------------------------------------------------------------

     allocate( lsize_loc(pard1),stat=istat)
     call memchk( zero, istat, mem_servi(1:2,servi), 'LSIZE_LOC', 'par_brotyp', lsize_loc )
     ! 
     ! Master: LAI1P (int) <= PAI1P (type)
     !
     if( IMASTER ) then
        nsize_loc = 0
        do ipari = 1,pard1
           lsize_loc(ipari) = size(pai1p(ipari)%l)
           nsize_loc = nsize_loc + lsize_loc(ipari)
        end do

        allocate( lai1p_loc(nsize_loc),stat=istat)
        call memchk( zero, istat, mem_servi(1:2,servi), 'LAI1P_LOC', 'par_brotyp', lai1p_loc )
        kk = 0
        do ipari = 1,pard1
           do ii = 1,lsize_loc(ipari)
              kk = kk + 1
              lai1p_loc(kk) = pai1p(ipari)%l(ii) 
           end do
        end do
     end if

     npari =  pard1
     parin => lsize_loc
     strin =  'LSIZE_LOC'
     call par_broadc()  

     if( ISLAVE ) then
        nsize_loc = 0
        do ipari = 1,pard1
           nsize_loc = nsize_loc + lsize_loc(ipari)
        end do
        allocate( lai1p_loc(nsize_loc),stat=istat)
        call memchk( zero, istat, mem_servi(1:2,servi), 'LAI1P_LOC', 'par_brotyp', lai1p_loc )
     end if

     npari =  nsize_loc
     parin => lai1p_loc
     strin =  'LAI1P_LOC'
     call par_broadc()  

     if( ISLAVE ) then
        kk = 0
        do ipari = 1,pard1
           allocate( pai1p(ipari)%l(lsize_loc(ipari)),stat=istat)
           call memchk( zero, istat, mem_servi(1:2,servi), 'PAI1P', 'par_brotyp', pai1p(ipari)%l )
           do ii = 1,lsize_loc(ipari)
              kk = kk + 1
              pai1p(ipari)%l(ii) = lai1p_loc(kk) ! PAI1P <= LAI1P
           end do
        end do
     end if

     call memchk( two, istat, mem_servi(1:2,servi), 'LAI1P_LOC','par_brotyp', lai1p_loc)
     deallocate( lai1p_loc, stat=istat )
     if(istat/=0) call memerr( two, 'LAI1P_LOC', 'par_brotyp',0_ip)

     call memchk( two, istat, mem_servi(1:2,servi), 'LSIZE_LOC','par_brotyp', lsize_loc)
     deallocate( lsize_loc, stat=istat )
     if(istat/=0) call memerr( two, 'LSIZE_LOC', 'par_brotyp',0_ip)

  end select

end subroutine par_brotyp

