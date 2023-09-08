!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_memall.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Allocates memory for the arrays for porous
!> @details Allocates memory for the arrays for porous
!> @} 
!------------------------------------------------------------------------
subroutine por_memall()
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_porous
  use mod_memchk
  use mod_memory
  implicit none
  !
  ! Solver
  !
  solve_sol => solve(1:)
  call soldef(4_ip)
  solve(1) % bvess     => bvess_por(:,:,1)
  solve(1) % kfl_fixno => kfl_fixno_por(:,:,1)
  solve(2) % bvess     => bvess_por(:,:,2)
  solve(2) % kfl_fixno => kfl_fixno_por(:,:,2)  
  !
  ! Problem unknowns PRESS, WASAT and solver initialization
  !
  if( INOTMASTER ) then
     !
     ! POROUS: Porous unknowns  - Pressure & Saturation &  Velocity  
     ! 
     call memory_alloca(mem_modul(1:2,modul),'VELOC','por_memall',veloc,ndime,npoin,ncomp_por)
     call memory_alloca(mem_modul(1:2,modul),'PRESS','por_memall',press,npoin,ncomp_por)
     call memory_alloca(mem_modul(1:2,modul),'WASAT','por_memall',wasat,npoin,ncomp_por)
     !
     ! Other arrays
     !
     call memory_alloca(mem_modul(1:2,modul),'SATIN_POR','por_memall',satin_por,nelem)
     call memory_alloca(mem_modul(1:2,modul),'PORO0_POR','por_memall',poro0_por,nelem)
     call memory_alloca(mem_modul(1:2,modul),'NODPO_POR','por_memall',nodpo_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'NODPE_POR','por_memall',nodpe_por,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'PERME_POR','por_memall',perme_por,ndime,nelem)
     call memory_alloca(mem_modul(1:2,modul),'IWELL_POR','por_memall',iwell_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'WINDE_POR','por_memall',winde_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'WMASS_POR','por_memall',wmass_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'IPWEL_POR','por_memall',ipwel_por,mheiw_por,nwell_por)
     call memory_alloca(mem_modul(1:2,modul),'DATAW_POR','por_memall',dataw_por,mheiw_por,nwell_por)
     call memory_alloca(mem_modul(1:2,modul),'IHEIP_POR','por_memall',iheip_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'XWELL_POR','por_memall',xwell_por,npoin) 

  else

     !
     ! Master: allocate minimum memory
     !
     call memory_alloca(mem_modul(1:2,modul),'VELOC','por_memall',veloc,1_ip,1_ip,3_ip)
     call memory_alloca(mem_modul(1:2,modul),'PRESS','por_memall',press,1_ip,3_ip)
     call memory_alloca(mem_modul(1:2,modul),'WASAT','por_memall',wasat,1_ip,3_ip)
     !
     ! Other arrays
     !
     call memory_alloca(mem_modul(1:2,modul),'SATIN_POR','por_memall',satin_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'PORO0_POR','por_memall',poro0_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'NODPO_POR','por_memall',nodpo_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'NODPE_POR','por_memall',nodpe_por,ndime,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'PERME_POR','por_memall',perme_por,1_ip,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'IWELL_POR','por_memall',iwell_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'WINDE_POR','por_memall',winde_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'WMASS_POR','por_memall',wmass_por,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'IPWEL_POR','por_memall',ipwel_por,mheiw_por,nwell_por)
     call memory_alloca(mem_modul(1:2,modul),'DATAW_POR','por_memall',dataw_por,mheiw_por,nwell_por)
     call memory_alloca(mem_modul(1:2,modul),'IHEIP_POR','por_memall',iheip_por,npoin)
     call memory_alloca(mem_modul(1:2,modul),'XWELL_POR','por_memall',xwell_por,1_ip) 

  end if

end subroutine por_memall
      
