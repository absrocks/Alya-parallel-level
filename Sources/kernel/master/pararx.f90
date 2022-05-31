subroutine pararx(wtask,ntype,ndimx,rvarx)
  !------------------------------------------------------------------------
  !****f* Parall/pararx
  ! NAME
  !    pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  npari,nparr,nparc,parin,parre,parx1,kfl_paral
  use def_master, only     :  party,parki,pardi,IPARALL,pard1
  use def_master, only     :  NPOIN_TYPE,NBOPO_TYPE,lntra
  use def_domain, only     :  npoin,nbopo
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimx
  complex(rp),  target     :: rvarx(ndimx)

  if( IPARALL ) then

     npari = 0
     nparc = 0
     nparr = 0

     select case ( wtask )

     case ( 'BCT' ) 

        call runend('PARARX: NOT CODE')

     case ( 'SND' ) 

        call runend('PARARX: NOT CODE')

     case ( 'RCV' )

        call runend('PARARX: NOT CODE')

     case ( 'MIN' ) 

        call runend('PARARX: NOT CODE')

     case ( 'SUM' ) 

        call runend('PARARX: NOT CODE')

     case ( 'MAX' )

        call runend('PARARX: NOT CODE')

     case ( 'S2M' ) 

        call runend('PARARX: NOT CODE')

     case ( 'GAT' ) 

        call runend('PARARX: NOT CODE')

     case ( 'IBI' ) 

        call runend('PARARX: NOT CODE')

     case ( 'SLX' )
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimx/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimx/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  4
              pardi =  1
           else
              parki =  7
              pardi =  1
           end if
           parx1 => rvarx
           call par_slexch()
        else
           call runend('PARARX: NOT CODED')
        end if

     case ( 'SLA' )
        !
        ! par_slexca
        !
        call runend('PARARX: NOT CODED')

     case default

        call runend('PARARX: WRONG CASE')

     end select

     !nparr = 0

  end if

end subroutine pararx
 
