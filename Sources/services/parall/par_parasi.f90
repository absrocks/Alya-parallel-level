subroutine par_parasi(wtask,ndims,ivars,ndimr,ivarr)
  !------------------------------------------------------------------------
  !****f* Parall/par_parasi
  ! NAME
  !    par_parasi
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  npari,parin,npasi,paris
  use def_master, only     :  ISLAVE,kfl_paral
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ndims,ndimr
  integer(ip),  target     :: ivars(max(1_ip,ndims))
  integer(ip),  target     :: ivarr(max(1_ip,ndimr))

  if( ISLAVE ) then

     npari = 0
     npasi = 0 

     select case ( wtask )

     case ( 'SLX' ) 
        !
        ! par_slaves
        !
        npasi =  ndims
        npari =  ndimr
        paris => ivars
        parin => ivarr
        call par_slaves(1_ip)
        
     case default
        
        call runend('PARASI: WRONG CASE')
        
     end select

     npari = 0
     npasi = 0

  end if
  
end subroutine par_parasi
 
