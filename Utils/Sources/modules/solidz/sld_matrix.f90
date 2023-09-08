!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_matrix.f90
!> @author  Solidz Team
!> @date    November, 2017-Adds doxygen
!> @brief   This routine computes the matrix and right hand side
!> @details
!>      OUTPUT
!>      USES
!>        sld_elmope
!>        sld_bouope
!>
!>      USED BY
!>        mod_solution_methods_sld
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_matrix(itask)

  use def_kintyp,         only : ip, rp
  use def_master,         only : INOTMASTER, coupling, IMASTER
  use def_master,         only : npoi1, npoi2, npoi3, mem_modul
  use def_master,         only : rhsid, amatr
  use def_master,         only : stateland, stretlam, gdedet, gdeinv, troponin
  use def_master,         only : ITASK_MATRIX, solve_sol, modul
  use def_master,         only : CPU_ASSEMBLY,CPU_MINI_ASSEMBLY
  use def_master,         only : CPU_MAXI_ASSEMBLY,CPU_COUNT_ASSEMBLY
  use def_master,         only : cutim, timei, timef
  use def_master,         only : gdepo
  use def_master,         only : kfl_cellmod,CELL_TT2006_EXMEDI,CELL_OHARA_EXMEDI,CELL_FITZHUGH_EXMEDI,CELL_NOMOD_EXMEDI
  use def_domain,         only : ndime, npoin, nmate, mgaus, mnode
  use def_domain,         only : vmass, kfl_elcoh
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo, memory_alloca_min
  use def_solver,         only : SOL_MATRIX_HAS_CHANGED
  use def_solidz,         only : kfl_penal_sld, kfl_limit_sld, kfl_vofor_sld, kfl_conta_sld
  use def_solidz,         only : fexte_sld, finte_sld, macce_sld, fcont_sld, SLD_PDN_UNILATERAL
  use def_solidz,         only : fext2_sld, fint2_sld, fine2_sld
  use def_solidz,         only : kfl_cohes_sld, kfl_xfeme_sld
  use def_solidz,         only : kfl_pseud_sld, kfl_gdepo
  use def_solidz,         only : vofor_sld, cpu_ass_sol_sld
  use def_solidz,         only : ndofn_sld, frxid_sld, vdiag_sld
  use def_solidz,         only : nzrhs_sld, vmasx_sld, vmass_sld
  use def_solidz,         only : fnatu_sld
  use mod_ker_timeline,   only : ker_timeline
  use mod_communications, only : PAR_BARRIER
  use mod_timings,        only : timings_assembly
  use mod_alya2dlb,       only : alya2dlb_DLB_Enable
  use mod_alya2dlb,       only : alya2dlb_DLB_Disable

  implicit none

  integer(ip), intent(in)     :: itask                !< What to do
  integer(ip)                 :: izrhs,ipoin,idime,itott
  integer(ip)                 :: mephy,pmate,ierr
  real(rp)                    :: time1,time2,time3,time4
  real(rp)                    :: dummr(mgaus*mnode)
  real(rp),    pointer        :: frnat(:,:)

  call ker_timeline('INI_ASSEMBLY')

  !-----------------------------------------------------------------
  !
  ! Initializations
  !
  !-----------------------------------------------------------------

  time1     =  0.0_rp
  time2     =  0.0_rp
  time3     =  0.0_rp

  nullify(frnat)

  if( INOTMASTER ) then
     !
     ! Initialize solver
     !
     if (itask == 1_ip) then        ! explicit

        do izrhs=1,nzrhs_sld
           rhsid(izrhs)=0.0_rp
        end do
        call inisol()

     else if(itask == 2_ip) then    ! implicit

        do izrhs=1,nzrhs_sld
           rhsid(izrhs)=0.0_rp
        end do
        call inisol()

     end if

  end if
  !
  ! Initialize arrays
  !
  if( itask <= 2_ip .and. INOTMASTER ) then
     !
     ! Mass matrices
     !
     vdiag_sld(:) = 0.0_rp                          ! Diagonal mass matrix
     vmass_sld(:) = 0.0_rp                          ! Lumped mass matrix
     if( kfl_xfeme_sld == 1 ) vmasx_sld(:) = 0.0_rp ! Lumped mass matrix (X-FEM)
     !
     ! Force vectors
     !
     frxid_sld(:) = 0.0_rp                          ! Reaction forces
     finte_sld(:) = 0.0_rp                          ! Internal force vector
     fexte_sld(:) = 0.0_rp                          ! External force vector
     macce_sld(:) = 0.0_rp                          ! Inertial force vector
     if( kfl_conta_sld /= 0 ) fcont_sld(:) = 0.0_rp ! Contact force (reaction) vector
     !
     ! Nodal push back of the deformation gradient
     !
     if( kfl_gdepo /= 0 ) then
        gdepo(:,:,:) = 0.0_rp
        gdeinv(:,:,:) = 0.0_rp
        gdedet(:) = 0.0_rp
     end if
     !
     ! Coupling with exmedi: 1. computing the activation vector
     !
     mephy = 0
     troponin = 0.0_rp
     stretlam = 0.0_rp
     if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
        do pmate = 1,nmate
          if (kfl_cellmod(pmate) == CELL_FITZHUGH_EXMEDI) mephy = 1_ip
        end do
        if( mephy == 1_ip ) then
           call sld_exmedi(0_ip,0_ip,0_ip,dummr,dummr,0_ip)
        end if
        call sld_parall(4_ip)
     end if

  end if

  !-----------------------------------------------------------------
  !
  ! Matrix, preconditioner, RHS and elment assemblies
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) ierr = alya2dlb_DLB_Enable()

  call cputim(time1)
  if( itask <= 2 ) call sld_elmope(itask)
  call cputim(time2)

#ifdef ALYA_DLB
  call PAR_BARRIER()
  if( INOTMASTER ) ierr = alya2dlb_DLB_Disable()
#endif

  !-----------------------------------------------------------------
  !
  ! Solve push forward (slaves)
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) then

     if( kfl_gdepo /= 0 ) then
        call rhsmod(ndime*ndime,gdepo)
        call rhsmod(ndime*ndime,gdeinv)
        call rhsmod(       1_ip,gdedet)
        do ipoin=1, npoin
           gdepo(1:ndime,1:ndime,ipoin)  = gdepo( 1:ndime,1:ndime,ipoin) / vmass(ipoin)
           gdeinv(1:ndime,1:ndime,ipoin) = gdeinv(1:ndime,1:ndime,ipoin) / vmass(ipoin)
           gdedet(ipoin)                 = gdedet(ipoin) / vmass(ipoin)
        end do
     end if

  end if

  !-----------------------------------------------------------------
  !
  ! Pseudo time step and EXMEDI coupling
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! Explicit case: Interchange vdiag_sld among the slaves and set it up properly (DONE ONLY IF INOTMASTER)
     ! and pseudotime
     !
     if (itask == 1 .and. kfl_pseud_sld == 1) then
        call runend('mariano hay que dividir por la matriz de masa aqui')
        call rhsmod(1_ip,vdiag_sld)
     end if
     !
     ! Interchange exmedi-related state variables
     !
     if (itask <= 2) then
        if(( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) then
           ! Parallel exchange of stretlam
           call rhsmod(1_ip,stretlam)
           call rhsmod(     1_ip,troponin)
           if (size(stateland) > 12) then

!!!              call rhsmod(1_ip,vdiag_sld)
!!!!              acaaaaaaaaaaaaaaaaaa ver esto para intercambiar la del estado

           end if
           if (size(troponin) > 1_ip) then
              call rhsmod(     1_ip,troponin)
              do ipoin=1,npoin
                 troponin(ipoin) = troponin(ipoin)/vmass(ipoin)
              end do
           end if
        end if
     end if
     !
     ! XFEM - Parall exchange of mass matrix
     !
     if (kfl_xfeme_sld == 1) call rhsmod(1_ip,vmasx_sld)
  end if

  !-----------------------------------------------------------------
  !
  ! Boundary assembly
  !
  !-----------------------------------------------------------------

  call cputim(time3)
  if( INOTMASTER ) then
     call sld_bouope(itask)
  end if
  call cputim(time4)

  !-----------------------------------------------------------------
  !
  ! Other stuff
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! Cohesive laws
     !
     if (kfl_cohes_sld > 0 .and. kfl_elcoh > 0) call sld_cohesi(itask)
     if (kfl_cohes_sld > 0 .and. kfl_elcoh > 0) call sld_cohele(itask)
     !
     ! Coupling
     !
     call sld_coupli(ITASK_MATRIX)
     !
     ! Limiter
     !
     if (kfl_limit_sld==1) call sld_limite()
     !
     ! Add external volume forces
     !
     if( kfl_vofor_sld > 0 ) then
        !
        ! Internal nodes
        !
        do ipoin = 1,npoin
           if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
              itott = (ipoin-1) * ndofn_sld
              do idime = 1,ndime
                 itott = itott + 1
                 rhsid(itott) = rhsid(itott) + vofor_sld(idime,ipoin)*(cutim/(timef-timei))
                 fexte_sld(itott) = fexte_sld(itott) + vofor_sld(idime,ipoin)*(cutim/(timef-timei))
              end do
           end if
        end do

     end if
     !
     ! Penalize system
     !
     if (kfl_penal_sld == 1) call sld_pensys(amatr)

  end if
  !
  ! Parallel exchage: rhsmod for force vectors and mass matrix
  !
  if( INOTMASTER ) then
     call rhsmod(ndime,frxid_sld)
     call rhsmod(ndime,finte_sld)
     call rhsmod(ndime,fexte_sld)
     call rhsmod(ndime,macce_sld)
     call rhsmod(1_ip, vmass_sld)
     if ( kfl_conta_sld == SLD_PDN_UNILATERAL ) call rhsmod(ndime, fcont_sld)
  end if
  !
  ! L2 norms of the force vectors required for the converge criteria
  !
  call norm2x(ndime,finte_sld,fint2_sld)
  call norm2x(ndime,fexte_sld,fext2_sld)
  call norm2x(ndime,macce_sld,fine2_sld)
  !
  ! Timings
  !
  cpu_ass_sol_sld(1) = time2 - time1
  cpu_ass_sol_sld(4) = time4 - time3
  call timings_assembly(cpu_ass_sol_sld(1),cpu_ass_sol_sld(4))
  !
  ! L2 norm of natural force imposed in solver
  !
  call memory_alloca(mem_modul(1:2,modul),'FRNAT','sld_matrix',frnat,ndime,npoin,'DO NOT INITIALIZE')
  if( IMASTER ) call memory_alloca_min(frnat)
  frnat = 0.0_rp
  fnatu_sld = 0.0_rp
  if( solve_sol(1) % kfl_bvnat == 1_ip .and. associated(solve_sol(1) % bvnat) ) then
     if( INOTMASTER ) then
        frnat = solve_sol(1) % bvnat
        call rhsmod(ndime,frnat)
     end if
  end if
  call norm2x(ndime,frnat,fnatu_sld)
  call memory_deallo(mem_modul(1:2,modul),'FRNAT','sld_matrix',frnat)
  !
  ! Timeline
  !
  call ker_timeline('END_ASSEMBLY')
  !
  ! Matrix has been assembled
  !
  solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED ! Matrix has been assembled

end subroutine sld_matrix
