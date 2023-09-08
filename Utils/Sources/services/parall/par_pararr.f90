subroutine par_pararr(wtask,ntype,ndimr,rvarr)
  !------------------------------------------------------------------------
  !****f* Parall/par_pararr
  ! NAME
  !    par_pararr
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : npari,nparr,nparc,parin,parre,parr1,kfl_paral
  use def_master,         only : party,parki,pardi,IPARALL,pard1,pari1,kfl_desti_par
  use def_master,         only : NPOIN_TYPE,NBOPO_TYPE,NBOUN_TYPE,NELEM_TYPE,lntra
  use def_domain,         only : npoin,nbopo,nboun,nelem
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimr
  real(rp),     target     :: rvarr(ndimr)

  if( IPARALL ) then

     npari = 0
     nparc = 0
     nparr = 0

     select case ( wtask )

     case ( 'BCT' ) 

        nparr =  ndimr
        parre => rvarr
        call par_broadc()

     case ( 'SND' ) 

        nparr =  ndimr
        parre => rvarr
        call par_sendin() 

     case ( 'RCV' )

        nparr =  ndimr
        parre => rvarr
        call par_receiv() 

     case ( 'MIN' ) 

        call PAR_MIN(ndimr,rvarr)
        !nparr =  ndimr
        !parre => rvarr
        !call par_operat(1_ip)

     case ( 'SUM' ) 

        call PAR_SUM(ndimr,rvarr)
        !nparr =  ndimr
        !parre => rvarr
        !call par_operat(3_ip)

     case ( 'MAX' )

        call PAR_MAX(ndimr,rvarr)
        !nparr =  ndimr
        !parre => rvarr
        !call par_operat(2_ip) 

     case ( 'S2M' ) 

        nparr =  ndimr
        parre => rvarr
        kfl_desti_par=0                                     ! Send to master
        call par_sendin()

     case ( 'GAT' ) 

        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE .or. ntype == NBOUN_TYPE .or. ntype == NELEM_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NELEM_TYPE ) pard1 = ndimr/nelem
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/max(1_ip,nbopo)
           if( ntype == NBOUN_TYPE ) pard1 = ndimr/max(1_ip,nboun)
           party = ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  6
              pardi =  1
           end if
        else
           party =  ntype
           parki =  2
           pardi =  1
           npari =  ndimr
        end if
        parr1 => rvarr
        call par_mygather()

     case ( 'IBI' ) 

        pari1 => lntra
        pard1 =  ndimr/npoin
        parr1 => rvarr(1:ndimr)
        call par_slexib(2_ip)

     case ( 'SLX' )
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call par_slexch()
        else
           call runend('PARARR: NOT CODED')
        end if

     case ( 'SLA' )
        !
        ! par_slexca
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimr/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimr/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  2
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
           parr1 => rvarr
           call par_slexca()
        else
           call runend('PARARR: NOT CODED')
        end if

     case default

        call runend('PARARR: WRONG CASE')

     end select

     nparr = 0
     nullify(parre)

  end if

end subroutine par_pararr
 
