!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @addtogroup CPU_Time
!> @{
!> @file    nsi_outcpu.f90
!> @author  houzeaux
!> @date    2018-12-30
!> @brief   Output Nastin CPU time
!> @details his routine writes a summary of spent computing time. The
!>          quantities computed by the code are:
!>          CPUTI_NSI(1) .... Element assembly 
!>          CPUTI_NSI(2) .... Boundary assembly
!>          CPUTI_NSI(3) .... SGS solver
!>          CPUTI_NSI(4) .... Correction
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_outcpu()

  use def_master
  use def_domain
  use def_nastin
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
  use mod_outfor,         only : outfor
  implicit none
  real(rp) :: xfact,cpuot
  real(rp) :: cpmax
  !
  ! Guillume explained that he does not use max because if you do max of each component it is not ok.
  ! One option would be to do max of the total and then print the values for that process.
  ! But this is not easy.
  ! One simple solution would be to leave everything averge but at least print additionally the max for the total.   
!  cpmax = sum(cputi_assembly_nsi)
!  call PAR_MAX    (1_ip,cpmax         ,'IN MY CODE')
  !  write(*,*) 'max(cputi_assembly_nsi)', cpmax         ! MISING SEND TO  NSI_OUTCPU  
  
  call PAR_MAX    (10_ip,cputi_nsi         ,'IN MY CODE')
  call PAR_AVERAGE(10_ip,cputi_assembly_nsi,'IN MY CODE')
  
  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     routp(1) = cpu_modul(30,modul)
     call outfor(29_ip,momod(modul)%lun_outpu,' ')

     if( cpu_modul(30,modul) > 0.0_rp ) then
        xfact = 100.0_rp / routp(1)
     else
        xfact = 1.0_rp
     end if
     cpuot = 0.0_rp

     if( NSI_MONOLITHIC ) then
        !
        ! MONOLITHIC
        !
        coutp(1)  = 'Solver mom.+con.'
        routp(1)  = solve(1) % cputi(1)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')

     else if( NSI_SCHUR_COMPLEMENT ) then
        !
        ! SCHUR COMPLEMENT
        !
        coutp(1)  = 'Solver momentum='
        routp(1)  = solve(1) % cputi(1)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')
        coutp(1)  = 'Solver continuity='
        routp(1)  = solve(2) % cputi(1)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')       

     else if( NSI_FRACTIONAL_STEP ) then
        !
        ! FRACTIONAL STEP
        !
        coutp(1)  = 'Solver momentum='
        routp(1)  = cputi_nsi(1)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')
        coutp(1)  = 'Solver continuity='
        routp(1)  = cputi_nsi(2)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')
        
     end if

     if( kfl_corre_nsi /= 0 ) then
        coutp(1)  = 'Velocity correction='
        routp(1)  = cputi_nsi(4)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')
     end if
     if( kfl_sgsti_nsi /= 0 .or. kfl_sgsco_nsi /= 0 .or. kfl_stabi_nsi > 0 ) then
        coutp(1)  = 'SGS + projection=' 
        routp(1)  = cputi_nsi(3)
        routp(2)  = xfact * routp(1)
        cpuot     = cpuot + routp(1)
        call outfor(30_ip,momod(modul)%lun_outpu,' ')
     end if

     cpuot = cpu_modul(30,modul) - cpuot
     !
     ! Others
     !
     coutp(1)  = 'OTHERS (ASSEMBLY...)'
     routp(1)  = cpuot
     routp(2)  = xfact*routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')
     !
     ! Timing for assembly
     !
     ioutp(1)    = cpu_modul(CPU_COUNT_ASSEMBLY,modul)
     routp(1:10) = cputi_assembly_nsi(1:10) 
     call outfor(83_ip,momod(modul)%lun_outpu,' ')
     
  end if

end subroutine nsi_outcpu
