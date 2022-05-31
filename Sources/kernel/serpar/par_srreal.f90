subroutine par_srreal(itask,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/par_srreal
  ! NAME
  !    par_srreal
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_master, only    :  npari,nparr,nparc,parin,parre
  implicit none
  integer(ip), intent(in) :: itask,ndimr
  real(rp),    target     :: rvarr(max(1_ip,ndimr))

  npari = 0
  nparc = 0

  select case ( itask )

  case ( 1_ip ) 

     nparr =  ndimr
     parre => rvarr
     call par_sendin()
     nullify(parre)

  case ( 2_ip )

     nparr =  ndimr
     parre => rvarr
     call par_receiv()
     nullify(parre)

  end select
  
end subroutine par_srreal
 
