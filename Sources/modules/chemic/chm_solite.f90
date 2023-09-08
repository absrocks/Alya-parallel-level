subroutine chm_solite()
  !-----------------------------------------------------------------------
  !****f* partis/chm_solite
  ! NAME 
  !    chm_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the species equations.
  ! USES
  !    chm_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper 
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: izrhs,izmat
  real(rp)    :: cpu_refe1,cpu_refe2,time1,time2

  reset: do
  call cputim(cpu_refe1)
  !
  ! Update inner iteration counter
  !
  itinn(modul)  =  itinn(modul) + 1_ip
  ittot_chm     =  ittot_chm    + 1_ip
  rtpts_chm     =  0.0_rp
  comin_chm     =  1.0e9_rp
  comax_chm     = -1.0e9_rp
  kfl_goit2_chm =  0_ip 
  !
  ! Solve the algebraic system
  !
  if( kfl_coupl_chm == 0_ip ) then
     !
     ! Gauss-Seidel: Assemble A^i, b^i and solve x^i
     !
     call livinf(160_ip,' ',1_ip)
     !
     ! In combustion, update initial values of specific heat and viscosity with temperature coming from outside
     !
     if (kfl_model_chm == 4_ip .and. kfl_gauss_chm < 2_ip) then
        call chm_upcpmu(1_ip) 
        call chm_upwmea(6_ip) 
        call ker_updpro()                ! We need to update the global viscosity,density, and Cp
        call chm_omegak(1_ip,nspec_chm)  ! Update all the mass source terms in combustion
     endif

     do iclas_chm = 1,nclas_chm
        iclai_chm = iclas_chm
        iclaf_chm = iclas_chm
        kfl_gocla_chm =  1
        if (INOTSLAVE .and. (kfl_model_chm==4 .or. kfl_model_chm==5)) then
           call livinf(165_ip,speci(iclas_chm)%name,0_ip)
           if (iclas_chm < nclas_chm) call livinf(165_ip,'-',0_ip)
        endif
        do while(kfl_gocla_chm<=kfl_spite_chm)
           call chm_updunk(7_ip)
           !
           ! Strong Gauss-Seidel we update everything for each species and iteration
           !
           if (kfl_model_chm == 4_ip .and. kfl_gauss_chm == 2_ip) then
              call chm_upcpmu(1_ip)
              call chm_upwmea(6_ip) 
              call ker_updpro()
              call chm_omegak(iclai_chm,iclaf_chm)
           endif
           call chm_matrix()    
           call cputim(time1)
           if (kfl_reset == 1_ip) exit reset
           !call chm_edgeba(rhsid,unkno,amatr,pmatr)
           call solver(rhsid,unkno,amatr,pmatr)
           call cputim(time2)
           cputi_chm(3) = cputi_chm(3) + (time2-time1)       
           if( kfl_normc_chm == 3_ip ) ripts_chm(iclas_chm) = solve(1)%resin        
           call chm_endite(3_ip)    ! Convergence and update solution
           kfl_gocla_chm=kfl_gocla_chm+1_ip
        enddo
     enddo

     call livinf(164_ip,' ',1_ip)
     call livinf(56_ip,' ',modul)        

  else if( kfl_coupl_chm == 1 ) then
     !
     ! Jacobi
     !
     call chm_matrix() ! Assemble all A^i, b^i
     if (kfl_reset == 1) exit reset
     izrhs = 1         ! Solve A^i x^i = b^i for all i successively
     izmat = 1
     do iclas_chm = 1,nclas_chm
        call cputim(time1)
        call solver(rhsid(izrhs),unkno(izrhs),amatr(izmat),pmatr)
        call cputim(time2)
        cputi_chm(3) = cputi_chm(3) + (time2-time1)       
        if( kfl_normc_chm == 3 ) ripts_chm(iclas_chm) = solve(1)%resin
        izrhs = izrhs + solve(1)%nzrhs
        izmat = izmat + solve(1)%nzmat
     end do

  else if( kfl_coupl_chm == 2 ) then
     !
     ! Monolithic
     !
     !
     ! In combustion, update initial values of specific heat and viscosity with temperature coming from outside
     !
     if (kfl_model_chm == 4) then
        call chm_upcpmu(1_ip) 
        call chm_upwmea(6_ip) 
        call ker_updpro() ! We need to update the global viscosity,density, and Cp
        call chm_omegak(1_ip,nspec_chm)   ! Update all the mass source terms in combustion
     endif
     iclai_chm = 1
     iclaf_chm = nspec_chm
     if (INOTSLAVE .and. kfl_model_chm==4) then
        call livinf(165_ip,' -- MONOLITHIC',0_ip)
     endif
     call chm_updunk(7_ip)
     call chm_matrix()    
     if (kfl_reset == 1) exit reset
     call cputim(time1)
     call solver(rhsid,unkno,amatr,pmatr)
     call cputim(time2)
     cputi_chm(3) = cputi_chm(3) + (time2-time1)       
     if( kfl_normc_chm == 3 ) ripts_chm = solve(1)%resin        
     call chm_endite(3_ip) ! Convergence and update solution

  end if

  call cputim(cpu_refe2)
  cpu_modul(3,modul) = cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 
  exit reset
enddo reset

end subroutine chm_solite
