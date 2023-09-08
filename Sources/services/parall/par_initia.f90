subroutine par_initia()
  !------------------------------------------------------------------------
  !****f* Parall/par_initia
  ! NAME
  !    par_initia
  ! DESCRIPTION
  !    This routine initialize MPI and open files
  ! OUTPUT
  !    nproc_par ... number of processes
  !    iproc_par ... my PID
  !    kfl_paral ... =iproc_par: tell Master if the process was initiated
  !                  by MPI
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use mod_parall,         only : PAR_COMM_MY_CODE4
  use mod_parall,         only : PAR_COMM_WORLD
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_parall,         only : PAR_CODE_SIZE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE

  implicit none
#ifdef MPI_OFF
#else
  !include  'mpif.h'
#endif

  integer(4)  :: istat,iproc_par4,nproc_par4,color
  integer(4)  :: iproc_par4_par,nproc_par4_par
  integer(ip) :: i

#ifdef MPI_OFF
  iproc_par = 0
  nproc_par = 1
#else

  istat = 0_4
 
  !call MPI_COMM_RANK(PAR_COMM_MY_CODE4,iproc_par4,istat)
  !call MPI_COMM_SIZE(PAR_COMM_MY_CODE4,nproc_par4,istat)
  call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE4,iproc_par4,nproc_par4)
  iproc_par = int(iproc_par4,ip)
  nproc_par = int(nproc_par4,ip)

#endif
  !
  ! if(there is more than one process) then
  !   KFL_PARAL= 0 .... I am the master
  !   KFL_PARAL> 0 .... I am slave KFL_PARAL
  ! else if(there is only one process) then
  !   KFL_PARAL=-1 .... No parallelization
  ! end if
  !
  if( nproc_par > 1 ) then
     kfl_paral        = iproc_par
     PAR_MY_CODE_RANK = iproc_par4
     PAR_CODE_SIZE    = nproc_par4
     call vocabu(-1_ip,0_ip,0_ip)
  else
     PAR_MY_CODE_RANK = -1
     PAR_CODE_SIZE    =  1
  end if
  !
  ! Initialize communication parameters
  !
  npari = 0
  nparr = 0
  nparc = 0
  nparx = 0
  strre = 'NULL'
  strin = 'NULL'
  strch = 'NULL'
  strcx = 'NULL'
  !
  ! Initialize variables
  !
  cpu_paral     = 0.0_rp
  !
  ! Partition postprocess variables
  !
  wopos_par(1,1) = 'NEIGB'
  wopos_par(1,2) = 'INTER'
  wopos_par(1,3) = 'BOUND'

  wopos_par(2,1) = 'SCALA'
  wopos_par(2,2) = 'SCALA'
  wopos_par(2,3) = 'SCALA'
  !
  ! My rank in current code
  !
  !
  ! Others
  !
  ipass_par      = 0        ! Number of passes for asynchronous
  !
  ! Nullify pointers
  !
  nullify(padja_par)    
  nullify(ladja_par)
  nullify(lepar_par)    
  nullify(lnpar_par)    
  nullify(lbpar_par)    
  nullify(leper_par)    
  nullify(lbper_par)    
  nullify(lneig_par)    
  nullify(ginde_par)  
  nullify(lsubz_par)  
  nullify(lcomm_par)  
  nullify(nskew_par)    
  nullify(ngive_par)    
  nullify(slfbo_par)    
  nullify(leind_par)    
  nullify(lbind_par)    
  nullify(lnods_par)  
  nullify(nhang_par)  
  nullify(xlnin_loc)    
  nullify(exnpe_loc)
  nullify(neighDom)    
  nullify(xadjDom)    
  nullify(adjDom)      
  nullify(translDual)   
  nullify(iaDual)      
  nullify(jaDual)       
  nullify(colours)
  nullify(badj)         
  nullify(bdom)         
  nullify(bpoin)
  nullify(leinv_par)    
  nullify(lbinv_par)    
  nullify(lnper_par)    
  nullify(lninv_par)    
  nullify(ireq4)       
  nullify(parws)       
  nullify(lowns_par) 
  nullify(lownr_par)   
  nullify(lnsec_par) 
  nullify(nnset_par)   
  nullify(nwitn_par)   

end subroutine par_initia
