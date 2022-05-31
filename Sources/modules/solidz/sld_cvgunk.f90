!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_cvgunk.f90
!> @author  Gerard Guillamet
!> @date    February, 2017
!>           - Subroutine written
!> @date    February, 2018
!>           - Add energies
!> @brief   Convergence checks for Solidz module
!> @details
!>          Convergence criteria implemented:
!>           - Residual error criterion (Belytschko p354)
!>           - Displ. increment error criterion (Belytschko p354)
!>           - Energy convergence criterion (Belytschko p354)
!>
!>          \verbatim
!>          ITASK = 0 ... Writting and timings output
!>                  1 ... Convergence inner iterations
!>                  2 ... Convergence outer iterations
!>                  3 ... Check residual time iterations
!>          \endverbatim
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures.
!>
!> @todo    To do list::\n
!>           - Extract into another subroutine outputs displacements
!>             reaction forces and exmedi stuff.
!>           - Revision rhsid with and without fixity for force and energy
!>             criterions.
!>           - Include natural forces (from fsi problems) in the energy
!>             criterion.
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_cvgunk(itask)

  use def_kintyp,            only : ip, rp
  use def_master,            only : zeror, isnain, INOTSLAVE, ITER_K, TIME_N
  use def_master,            only : ITASK_ENDINN, ITASK_ENDITE, ITASK_ENDSTE
  use def_master,            only : unkno, displ, rhsid, modul, kfl_rstar
  use def_master,            only : cpu_initi, momod
  use def_master,            only : itinn, itcou, ittim, cutim
  use def_domain,            only : ndime, npoin
  use mod_outfor,            only : outfor
  use mod_array_operations,  only : array_operations_residual_norm
  use mod_array_operations,  only : array_operations_norm2
  use def_solidz,            only : last_iters_sld
  use def_solidz,            only : kfl_resid_sld, kfl_goite_sld
  use def_solidz,            only : kfl_timei_sld, kfl_stead_sld
  use def_solidz,            only : fint2_sld, fext2_sld, fine2_sld, fnatu_sld
  use def_solidz,            only : allie_sld, allwk_sld, allke_sld, etota_sld
  use def_solidz,            only : resid_sld, eener_sld
  use def_solidz,            only : miinn_sld, dtcri_sld
  use def_solidz,            only : cotol_sld, sstol_sld
  use def_solidz,            only : resdispl_sld, resforce_sld
  use def_solidz,            only : SLD_STATIC_PROBLEM
  use def_solidz,            only : kfl_volca_sld
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_write_cvg

  implicit none

  integer(ip), intent(in) :: itask !< What to do
  integer(ip), save       :: kpass = 0
  real(rp),    save       :: cpuit_sld=0.0_rp
  real(rp),    save       :: risld,rfsld,resld,resti=1.0_rp
  real(rp)                :: time1
  real(rp)                :: fmax2, rhsn2, emax
  real(rp)                :: cpu_ass_ave
  real(rp)                :: cpu_ass_max
  real(rp)                :: cpu_sol_ave
  real(rp)                :: cpu_sol_max
  real(rp)                :: displ_max

  select case ( itask )

  case ( 0_ip )

     !-------------------------------------------------------------------
     !
     ! Timings and output
     !
     !-------------------------------------------------------------------

     !
     ! Compute some minmax values
     !
     call sld_minmax(cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,displ_max)
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        !
        ! Write heading
        !
        if( kpass == 0 .and. kfl_rstar /= 2 ) then
           write(momod(modul) % lun_conve,100)
           if( kfl_volca_sld > 0_ip ) call sld_cardiac_cycle_write_cvg(1_ip)
        end if

        if( kpass == 1 ) then
           time1 = time1 - cpuit_sld
        else
           time1 = time1 - cpu_initi
        end if
        !
        ! Write results
        !
        write(momod(modul) % lun_conve,101) &
             !   1     2            3              4     5         6
             ittim,itcou,itinn(modul),last_iters_sld,cutim,dtcri_sld,&
             !   7     8     9        10                11                12
             risld,rfsld,resld,displ_max,allie_sld(ITER_K),allwk_sld(ITER_K),&
             !              13                14
             allke_sld(ITER_K),etota_sld(ITER_K),&
             !  15          16          17          18          19
             time1,cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,&
             !              20                 21                 22
             resforce_sld(1,1), resforce_sld(2,1), resforce_sld(3,1),&
             !              23                 24                 25
             resdispl_sld(1,1), resdispl_sld(2,1), resdispl_sld(3,1)
             
        call cputim(cpuit_sld)
        flush(momod(modul) % lun_conve)
        !
        ! Write results cardiac cycle
        !
        if( kfl_volca_sld > 0_ip ) call sld_cardiac_cycle_write_cvg(2_ip)
        
     end if

     kpass = 1_ip

  case ( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the inner iterations (N-R solidz iterations)
     ! RISLD: Displacement residual
     ! RFSLD: Force residual
     !
     !-------------------------------------------------------------------

     !
     ! Displacement increment error (Belytschko p332)
     ! || d(n,i,j) - d(n,i,j-1)|| / ||d(n,i,j)|| (L2 norm)
     !
     risld = array_operations_residual_norm(2_ip,ndime,ndime,unkno,displ,0_ip,0_ip,ndime,1.0_rp)

     !
     ! Force residual error (Belytschko p332)
     ! Note that a meaningful value is added to avoid small relative
     ! differences between numerator and denominator typical in fully
     ! damaged elements. Then rhs < tol*1e-3
     ! ||rhs|| / max( ||fext||, ||fint||, ||M*a|| ) (L2 norm)
     !
     rhsn2 = array_operations_norm2(rhsid,DOFS=ndime)
     fmax2 = max(fint2_sld, fext2_sld, fine2_sld, fnatu_sld, 1.0e-3_rp)
     rfsld = rhsn2/fmax2

     if( INOTSLAVE ) then
       write (*,"(A,3E10.2)",advance="no") ' |RES|=' , rfsld
     end if

     !
     ! Energy convergence criterion (Belytschko p354)
     ! |du*rhs | / max( Wint, Wext, Wkin )
     !
     emax  = max(abs(allie_sld(ITER_K)), abs(allwk_sld(ITER_K)), allke_sld(ITER_K), 1.0e-3_rp)
     resld = eener_sld/emax
     !
     ! KFL_GOITE_SLD: Check convergence
     !
     if( isnain(risld) .or. isnain(rfsld) .or. isnain(resld) ) &
                         call runend("SOLIDZ: Nan or +/- Inf in residuals")
     if( risld > 1.0e10_rp ) call runend("SOLIDZ: Huge displacement residual")
     if( rfsld > 1.0e10_rp ) call runend("SOLIDZ: Huge force residual")
     !if (resld > 1.0e10_rp ) call runend("SOLIDZ: Huge energy residual")
     !
     ! Types of convergence criteria used to terminate iterations:
     !
     if (kfl_resid_sld == 1_ip) then
        !
        ! FORCE: Only force residual (Belytschko)
        !
        if ((rfsld <= cotol_sld) .or. &                       ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     else if (kfl_resid_sld == 2_ip) then
        !
        ! ENERGY: Only energy residual (Belytschko)
        !
        if (resld <= cotol_sld .or. &                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     else if (kfl_resid_sld == 3_ip) then
        !
        ! TOTAL: Displacement,force and energy residuals (Belytschko)
        !
        !if ((risld <= cotol_sld .and. rfsld <= cotol_sld &
        !     .and. resld <= cotol_sld) .or. &                 ! Inner tolerance achieved
        !     itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations
        if ((risld <= cotol_sld .and. rfsld <= cotol_sld &
             ) .or. &                                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations
     else
        !
        ! DISPL: Only displacement residual (Belytschko)
        !
        if (risld <= cotol_sld .or. &                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     end if

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the outer iterations
     ! External/global iterations in a coupling problem
     ! RESID_SLD: only displacement residual
     !
     ! ||u^{n+1,i+1}-u^{n+1,i-1}||/||u^{n+1,i+1}||
     !-------------------------------------------------------------------

     resid_sld = array_operations_residual_norm(2_ip,ndime,ndime,displ,displ,0_ip,npoin*ndime,ndime)

  case ( ITASK_ENDSTE )

     !-------------------------------------------------------------------
     !
     ! Check residual of the time iterations
     ! Only for dynamic (transient) problems reaching a steady state.
     ! RESTI: Time residual measured only using displacement
     !
     ! ||u^{n+1}-u^n|/||u^{n+1}||
     !-------------------------------------------------------------------

     if( kfl_timei_sld /= SLD_STATIC_PROBLEM ) then
        resti = array_operations_residual_norm(2_ip,ndime,ndime,displ,displ,0_ip,2_ip*npoin*ndime,ndime)
     end if
     !
     ! Some check to stop the code if things are going bad
     !
     if( isnain(resti)     ) kfl_stead_sld = 1_ip    ! NaN or +/- Inf
     if( resti > 1.0e10_rp ) kfl_stead_sld = 1_ip    ! Huge
     !
     ! Steady-state has been reached
     if( resti <= sstol_sld ) then
        kfl_stead_sld = 1_ip
        call outfor(28_ip,momod(modul) % lun_outpu,' ')
     end if

  end select

  !----------------------------------------------------------------------
  !
  ! Formats
  !
  !----------------------------------------------------------------------

100 format('# --| ALYA Convergence',/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --|                                                                      ',/,&
       &   '# --|  1. Time step           2. Global Iteration    3. Inner Iteration    ',/,&
       &   '# --|  4. Solver iterations   5. Current time        6. Stable time incr.  ',/,&
       &   '# --|  7. RSI Displacement    8. RSI Force           9. RSI Energy         ',/,&
       &   '# --| 10. Max. Displacement  11. Internal energy    12. External work      ',/,&
       &   '# --| 13. Kinetic energy     14. Total energy       15. Elapsed CPU time   ',/,&
       &   '# --| 16. Ass. ave cpu time  17. Ass. max cpu time  18. Sol. ave cpu time  ',/,&
       &   '# --| 19. Sol. max cpu time                                                ',/,&
       &   '# --|                                                                      ',/,&
       &   '# --| Reaction forces and displacements at specified Dirichelet BCs        ',/,&
       &   '# --|                                                                      ',/,&
       &   '# --| 20. RFX                21. RFY                22. RFZ                ',/,&
       &   '# --| 23. UX                 24. UY                 25. UZ                 ',/,&
       &   '# --|                                                                      ')

101 format(4x,i9,2x,i9,2x,i9,2x,i9,50(2x,e16.8e3))

end subroutine sld_cvgunk
