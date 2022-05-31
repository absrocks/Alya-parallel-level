subroutine par_parari(wtask,ntype,ndimi,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/par_parari
  ! NAME
  !    par_parari
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : npari,npari,nparc,nparr,parin,parin,pari1
  use def_master,         only : party,parki,pardi,IPARALL,pard1,kfl_desti_par
  use def_master,         only : NPOIN_TYPE,NBOPO_TYPE,kfl_paral
  use def_domain,         only : npoin,nbopo
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimi
  integer(ip),  target     :: rvari(ndimi)
 
  if( IPARALL ) then

     npari = 0
     nparr = 0
     nparc = 0

     select case ( wtask )

     case ( 'BCT' ) 
        !
        ! par_broadc
        !
        npari =  ndimi
        parin => rvari
        call par_broadc()

     case ( 'SND' ) 
        !
        ! par_sendin
        !
        npari =  ndimi
        parin => rvari
        call par_sendin() 

     case ( 'S2M' ) 

        npari =  ndimi
        parin => rvari
        kfl_desti_par=0                                     ! Send to master
        call par_sendin()

     case ( 'RCV' )
        !
        ! par_receiv
        !
        npari =  ndimi
        parin => rvari
        call par_receiv() 

     case ( 'MIN' )
        !
        ! par_operat (minimum)
        !
        call PAR_MIN(ndimi,rvari)
        !npari =  ndimi
        !parin => rvari
        !call par_operat(1_ip)

     case ( 'SUM' )
        !
        ! par_operat (sum)
        !
        call PAR_SUM(ndimi,rvari)
        !npari =  ndimi
        !parin => rvari
        !call par_operat(3_ip)

     case ( 'MAX' )
        !
        ! par_operat (maximum)
        !
        call PAR_MAX(ndimi,rvari)
        !npari =  ndimi
        !parin => rvari
        !call par_operat(2_ip)

     case ( 'GAT' )
        !
        ! par_gather
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimi/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimi/nbopo
           party = ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  5
              pardi =  1
           end if
        else
           party =  ntype
           parki =  1
           pardi =  1
           npari =  ndimi
        end if
        pari1 => rvari
        call par_mygather() 

     case ( 'SLX' )
        !
        ! par_slexch
        !
        party =  ntype
        parki =  1
        pardi =  1
        npari =  ndimi
        pari1 => rvari
        call par_slexch() 

     case default
        
        call runend('PAR_PARARI: WRONG CASE')
        
     end select
     
     npari = 0
     nullify(parin)
     nullify(pari1)

  end if
  
end subroutine par_parari
 
