subroutine par_outprt()
  !-------------------------------------------------------------------------------
  !****f* Parall/par_outprt
  ! NAME
  !    par_outprt
  ! DESCRIPTION
  !    Output Master's partition: Put here all the things needed by the master 
  !    when doing restart
  ! INPUT
  ! OUTPUT
  !
  ! USED BY
  !    par_partit
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_parall 
  use def_master
  use mod_memory
  implicit none
  integer(ip) :: iunit,ii,ji,npart_tmp

  iunit         = lun_aonlp_par + 0
  kfl_desti_par = 0
  npart_tmp     = npart_par

  if( PART_AND_WRITE() .or. READ_AND_RUN() ) then

     strin = 'readim_reastr_reageo_cderda'
     strre = 'readim_reastr_reageo_cderda'
     strch = 'readim_reastr_reageo_cderda'
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0

        call iexcha(npart_tmp)

        !-------------------------------------------------------------------
        !
        !  Read in readim
        !
        !-------------------------------------------------------------------
 
        call iexcha(kfl_autbo)
        call iexcha(npoin)
#ifdef NDIMEPAR

#else
        call iexcha(ndime)
#endif         
        call iexcha(nelem)
        call iexcha(nboun)
        call iexcha(nperi)
        call iexcha(nfiel)
        call iexcha(ngrou_dom)

        !-------------------------------------------------------------------
        !
        ! Read in reaset
        !
        !-------------------------------------------------------------------

        call iexcha(neset)
        call iexcha(nbset)
        call iexcha(nnset)

        !-------------------------------------------------------------------
        !
        ! Read in reabcs
        !
        !-------------------------------------------------------------------

        call iexcha(kfl_icodn) 
        call iexcha(kfl_icodb)
        !
        ! Send or receive data
        !
        if(parii==1) then
           call memory_alloca(mem_servi(1:2,servi),'PARIN','par_outprt',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_servi(1:2,servi),'PARRE','par_outprt',parre,nparr,'DO_NOT_INITIALIZE')
           if( READ_AND_RUN() ) call par_receiv()
        end if
 
     end do

     if( PART_AND_WRITE() ) call par_sendin()

     call memory_deallo(mem_servi(1:2,servi),'PARIN','par_outprt',parin)
     call memory_deallo(mem_servi(1:2,servi),'PARRE','par_outprt',parre)
     !
     ! reafie arrays 
     !
     if( nfiel > 0 ) then
        if( READ_AND_RUN() ) then
           call memgeo(28_ip)
           call par_parari('RCV',0_ip,5_ip*nfiel,kfl_field)
        else
           call par_parari('SND',0_ip,5_ip*nfiel,kfl_field)   
        end if
     end if
     
     !-------------------------------------------------------------------
     !
     ! READ AND RUN is carried out with wrong number of subdomains
     !
     !-------------------------------------------------------------------

     if( npart_tmp /= npart_par ) then
        print*,'npart_tmp,npart_par',npart_tmp,npart_par
        call runend('PAR_OUTPRT:WRONG NUMBER OF SUBDOMAINS')
     end if

     !-------------------------------------------------------------------
     !
     ! Partition structure
     !
     !-------------------------------------------------------------------
     !
     ! Allocate memory for partition structure
     !
     if( READ_AND_RUN() ) then
        call par_memory(8_ip)
        call par_memset(1_ip)
        call par_memset(2_ip)
     end if 
     !
     ! Subdomain and partition dimensions, and sets
     !
     strin = 'Subdomain_and_partition_dimensions'
     strre = 'Subdomain_and_partition_dimensions'
     strch = 'Subdomain_and_partition_dimensions'
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        do ii=1,npart_par
           call iexcha(npoin_par(ii))
        end do
        do ii=1,npart_par
           call iexcha(nelem_par(ii))
        end do
        do ii=1,npart_par
           call iexcha(nboun_par(ii))
        end do
        do ii=1,npart_par
           call iexcha(slfbo_par(ii))
        end do
        call iexcha(npoin_total)
        call iexcha(nelem_total)
        call iexcha(nboun_total)
        call iexcha(gnb)
        call iexcha(gni)
        !
        ! Send or receive data
        !
        if(parii==1) then
           call memory_alloca(mem_servi(1:2,servi),'PARIN','par_outprt',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_servi(1:2,servi),'PARRE','par_outprt',parre,nparr,'DO_NOT_INITIALIZE')
           if( READ_AND_RUN() ) call par_receiv()
        end if

     end do

     if( PART_AND_WRITE() ) call par_sendin()

     call memory_deallo(mem_servi(1:2,servi),'PARIN','par_outprt',parin)
     call memory_deallo(mem_servi(1:2,servi),'PARRE','par_outprt',parre)

     !-------------------------------------------------------------------
     !
     ! Partition structure arrays
     !
     !-------------------------------------------------------------------

     npari = 0
     nparr = 0
     nparc = 0
     !
     ! GINDE_PAR
     !
     npari =  4*(npart_par+1)
     call par_senprt_int(npari,ginde_par)
     !
     ! LNEIG_PAR
     !
     npari =  npart_par
     parin => lneig_par
     strin =  'LNEIG_PAR'
     call par_senprt()

     goto 100
     if( kfl_outfo /= 50 ) then
        if( READ_AND_RUN() ) call par_memory(2_ip)
        if( READ_AND_RUN() ) call par_memory(20_ip) ! GrowSmarter
        !
        ! LNINV_LOC: Only needed when master does postprocess
        !
        npari =  npoin_total
        parin => lninv_loc
        strin =  'LNINV_LOC'
        call par_senprt()
        !
        ! XLNIN_LOC
        !
        npari =  npart_par+1
        parin => xlnin_loc
        strin =  'XLNIN_PAR'
        call par_senprt()
        !
        ! LBPER_PAR (GrowSmarter)
        !
        ! This boundary permutation array is needed in partition files
        ! to allow reread/send initial conditions through reabcs and
        ! par_reabcs when kfl_rread is enabled (REREAD tag in .dat)
        !
        if (nboun_total > 0) then
           npari =  nboun_total
           parin => lbper_par
           strin =  'LBPER_PAR'
           call par_senprt()
           !
           ! GrowSmarter TODO:
           ! We can save kfl_codbo which is not stored in partition file OR
           ! Reread .fix file which causes to recompute all the structures...
           !
        endif
     end if
100 continue
  end if

  nullify(parin)
  nullify(parre)

end subroutine par_outprt

subroutine par_senprt
  use def_master
  implicit none
  
  if( PART_AND_WRITE() ) then
     call par_sendin()
  else
     call par_receiv()
  end if

  nullify(parin)
  nullify(parre)
  
end subroutine par_senprt


