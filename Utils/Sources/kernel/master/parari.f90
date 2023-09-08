subroutine parari(wtask,ntype,ndimi,rvari)
  !------------------------------------------------------------------------
  !****f* Parall/parari
  ! NAME
  !    parari
  ! DESCRIPTION
  !    Works with arrays to deal with parallelization
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only : ip,rp
  use def_master,         only : npari,npari,nparc,nparr,parin,parin,pari1
  use def_master,         only : party,parki,pardi,IPARALL,pard1,npasr,paris
  use def_master,         only : NPOIN_TYPE,NFACE_TYPE,NBOPO_TYPE,npasi,nfacg
  use def_domain,         only : npoin,nbopo,ndime,npoin_2
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_SUM
  implicit none
  character(3), intent(in) :: wtask
  integer(ip),  intent(in) :: ntype,ndimi
  integer(ip),  target     :: rvari(ndimi)
 
  if( IPARALL ) then

     nparr = 0
     nparc = 0

     if( wtask == 'BCT' ) then 
        !
        ! par_broadc
        !
        npari =  ndimi
        parin => rvari
        call Parall(2_ip)

     else if( wtask =='SND' ) then  
        !
        ! par_sendin
        !
        npari =  ndimi
        parin => rvari
        call Parall(3_ip)

     else if( wtask =='S2M' ) then  

        npari =  ndimi
        parin => rvari
        call Parall(24_ip)

     else if( wtask =='RCV' ) then 
        !
        ! par_receiv
        !
        npari =  ndimi
        parin => rvari
        call Parall(4_ip)

     else if( wtask =='MIN' ) then 
        !
        ! par_operat (minimum)
        !
        if( ntype == NPOIN_TYPE ) then
           party = 3
           pardi = ndimi/npoin
           parki = 1
           npari =  ndimi
           pari1 => rvari           
           call Parall(424_ip)
        else
           call PAR_MIN(ndimi,rvari)
        end if

     else if( wtask =='SUM' ) then 
        !
        ! par_operat (sum)
        !
        call PAR_SUM(ndimi,rvari)

     else if( wtask =='MAX' ) then 
        !
        ! par_operat (maximum)
        !
        if( ntype == NPOIN_TYPE ) then
           party = 3
           pardi = ndimi/npoin
           parki = 1
           npari =  ndimi
           pari1 => rvari           
           call Parall(425_ip)
        else
           call PAR_MAX(ndimi,rvari)
        end if

     else if( wtask =='GAT' ) then 
        !
        ! par_gather
        !
        call runend('PARARI: GATHER NO LONGER HERE')

     else if( wtask == 'SLX' .or. wtask == 'SLG' ) then
        !
        ! par_slexch
        !
        if( wtask == 'SLG' ) then
           call Parall(611_ip)
        end if
        !
        ! par_slexch for vectors(ndimi)
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimi/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimi/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  6
              pardi =  2
           end if
           pari1 => rvari
           call Parall(400_ip)

        else if( ntype == NFACE_TYPE ) then

           pard1 = ndimi/nfacg
           party = ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  6
              pardi =  2
           end if
           pari1 => rvari
           call Parall(607_ip)
        else
           call runend('PARARR: NOT CODED')
        end if

        if( wtask == 'SLG' ) then
           call Parall(612_ip)
        end if

     else if( wtask =='SNR' ) then 
        !
        ! par_senrcv
        !
        npasr = 0
        nparr = 0
        call Parall(411_ip)
        npasi = 0
        npari = 0

     else if( wtask =='SRA' ) then 
        !
        ! par_senrcv: to any slave
        !
        npasr = 0
        nparr = 0
        call Parall(408_ip)
        npasi = 0
        npari = 0

     else if( wtask =='AGA' ) then 
        !
        ! Allgather
        !
        npari =  ndimi
        paris => rvari 
        call Parall(46_ip)

     else if( wtask =='FRI' ) then 
        !
        ! par_slexfr: take fringe geometry nodes
        !
        party = ntype
        if( ntype == NPOIN_TYPE ) then
           pard1 = ndimi/npoin_2
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              parki =  6
              pardi =  2
           end if
           pari1 => rvari
           call Parall(429_ip)
        end if

     else if( wtask =='ONL' ) then  
        !
        ! Only
        !
        party = ntype
        if( ntype == NPOIN_TYPE ) then           
           pard1 = ndimi/npoin
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           end if
           if( pardi /= 1 ) then
              call runend('PARARI: NOT CODED 3')
           end if
           npari =  ndimi
           pari1 => rvari 
           call Parall(430_ip)
        else
           call runend('PARARI: NOT CODED 2') 
        end if

     else if( wtask =='SMI' ) then 
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimi/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimi/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              call runend('PARARI: SMI NOT CODED')
              parki =  5
              pardi =  1
           end if
           pari1 => rvari
           call Parall(424_ip)
        else
           call runend('PARARI: SMI NOT CODED')
        end if

     else if( wtask =='SMA' ) then 
        !
        ! par_slexch
        !
        if( ntype == NPOIN_TYPE .or. ntype == NBOPO_TYPE ) then
           if( ntype == NPOIN_TYPE ) pard1 = ndimi/npoin
           if( ntype == NBOPO_TYPE ) pard1 = ndimi/nbopo
           party =  ntype
           if( pard1 == 1 ) then
              parki =  1
              pardi =  1
           else
              call runend('PARARI: SMA NOT CODED')
              parki =  5
              pardi =  1
           end if
           pari1 => rvari
           call Parall(425_ip)
        else
           call runend('PARARI: SMA NOT CODED')
        end if

     else
        
        call runend('PARARI: WRONG CASE') 
        
     end if
     
     npari = 0
     nullify(parin)
     nullify(pari1)

  end if
  
end subroutine parari
 
