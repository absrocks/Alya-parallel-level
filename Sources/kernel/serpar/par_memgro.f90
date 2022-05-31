subroutine par_memgro(itask)
  !-------------------------------------------------------------------------------
  !****f* parall/par_memgro
  ! NAME
  !    par_memgro
  ! DESCRIPTION
  !    Allocate memgro for partition dimensions and arrays
  ! INPUT
  ! OUTPUT
  ! USED BY
  !    par_create_graph_arrays
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_parall 
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Allocate in domgro
     !
     allocate(comle(icoml)%neighDom((npart_par*(npart_par+1))/2),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%neighDom','par_memgro',comle(icoml)%neighDom)

     allocate(comle(icoml)%lneig_par(npart_par),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%LNEIG_PAR','par_memgro',comle(icoml)%lneig_par)

     allocate(comle(icoml)%ngrou_par(npart_par),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%NGROU_PAR','par_memgro',comle(icoml)%ngrou_par)

     allocate(comle(icoml)%lnpar_par(comle(icoml)%ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%LNPAR_PAR','par_memgro',comle(icoml)%lnpar_par)

     allocate(comle(icoml)%ndomi(comle(icoml)%ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%ndomi','par_memgro',comle(icoml)%ndomi)

     allocate(comle(icoml)%domli(npart_par,comle(icoml)%ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%domli','par_memgro',comle(icoml)%domli)

     allocate(comle(icoml)%xadjDom(npart_par+1),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%xadjDom','par_memgro',comle(icoml)%xadjDom)

  case(2_ip)
     !
     ! Allocate in domgro
     !
     allocate(comle(icoml)%adjDom(comle(icoml)%xadjDom(npart_par+1)-1),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%adjDom','par_memgro',comle(icoml)%adjDom)

  case(3_ip)
     !
     ! Allocate in arrgro
     !
     allocate( comle(icoml)%badj(comle(icoml)%gnb+1),stat=istat)
     call memchk(zero,istat,par_memor, 'badj', 'par_memgro', comle(icoml)%badj )    
    
     allocate( comle(icoml)%bdom(comle(icoml)%ngrou_total-comle(icoml)%gni),stat=istat)
     call memchk(zero,istat,par_memor, 'bdom', 'par_memgro', comle(icoml)%bdom )  
      
     allocate( comle(icoml)%bpoin(comle(icoml)%ngrou_total-comle(icoml)%gni),stat=istat)
     call memchk(zero,istat,par_memor, 'bpoin', 'par_memgro', comle(icoml)%bpoin )

  case(4_ip)
     !
     ! Allocate in comgro
     !
     allocate(comle(icoml)%lcomm_par(comle(icoml)%nbcol,npart_par),stat=istat)
     call memchk(zero,istat,par_memor,'comle(icoml)%lcomm_par','par_memgro',comle(icoml)%lcomm_par)  

  case(10_ip)
     !
     ! Deallocate in grogro
     !
     call memchk(two,istat,par_memor, 'comle(icoml)%bpoin', 'par_memgro', comle(icoml)%bpoin )
     deallocate( comle(icoml)%bpoin,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%bpoin','par_memgro',0_ip)

     call memchk(two,istat,par_memor, 'comle(icoml)%bdom', 'par_memgro', comle(icoml)%bdom )  
     deallocate( comle(icoml)%bdom,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%bdom','par_memgro',0_ip)

     call memchk(two,istat,par_memor, 'comle(icoml)%badj', 'par_memgro', comle(icoml)%badj )    
     deallocate( comle(icoml)%badj,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%badj','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%domli','par_memgro',comle(icoml)%domli)
     deallocate(comle(icoml)%domli,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%domli','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%ndomi','par_memgro',comle(icoml)%ndomi)
     deallocate(comle(icoml)%ndomi,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%ndomi','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%LNPAR_PAR','par_memgro',comle(icoml)%lnpar_par)
     deallocate(comle(icoml)%lnpar_par,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%lnpar_par','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%NGROU_PAR','par_memgro',comle(icoml)%ngrou_par)
     deallocate(comle(icoml)%ngrou_par,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%ngrou_par','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%LNEIG_PAR','par_memgro',comle(icoml)%lneig_par)
     deallocate(comle(icoml)%lneig_par,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%lneig_par','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%neighDom','par_memgro',comle(icoml)%neighDom)
     deallocate(comle(icoml)%neighDom,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%neighDom','par_memgro',0_ip)

     call memchk(zero,istat,par_memor,'comle(icoml)%lcomm_par','par_memgro',comle(icoml)%lcomm_par)
     deallocate(comle(icoml)%lcomm_par,stat=istat)
     if(istat/=0) call memerr(two,'comle(icoml)%lcomm_par','par_memgro',0_ip)

  end select

end subroutine par_memgro
