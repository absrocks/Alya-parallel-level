!-----------------------------------------------------------------------
!> @addtogroup NeutroTurnon
!> @{
!> @file    neu_memall.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory 
!> @details Allocate memory 
!> @} 
!-----------------------------------------------------------------------
subroutine neu_memall()
  use def_kintyp, only : ip
  use def_master, only : neutr
  use def_master, only : solve
  use def_master, only : INOTMASTER
  use def_master, only : mem_modul,modul
  use def_domain, only : npoin,ndime
  use def_solver, only : solve_sol
  use def_neutro, only : num_energies_neu
  use def_neutro, only : num_directions_neu
  use def_neutro, only : ncomp_neu
  use def_neutro, only : direc_neu
  use def_neutro, only : weigd_neu
  use def_neutro, only : scattering_neu
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_alloca_min
  implicit none

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------

  solve_sol => solve(1:1)
  call soldef(4_ip)

  !----------------------------------------------------------------------
  !
  ! Arrays
  !
  !----------------------------------------------------------------------
  !
  ! Directions
  !
  call memory_alloca(mem_modul(1:2,modul),'DIREC_NEU'     ,'neu_memall',direc_neu     ,3_ip,num_directions_neu)
  call memory_alloca(mem_modul(1:2,modul),'WEIGD_NEU'     ,'neu_memall',weigd_neu     ,num_directions_neu)
  call memory_alloca(mem_modul(1:2,modul),'SCATTERING_NEU','neu_memall',scattering_neu,num_directions_neu,num_directions_neu)
  
  if( INOTMASTER ) then
     !
     ! RADIA
     !
     call memory_alloca(mem_modul(1:2,modul),'RADIA','neu_memall',neutr,num_energies_neu,num_directions_neu,npoin,ncomp_neu)
 
  else
     !
     ! RADIA
     ! 
     call memory_alloca_min(mem_modul(1:2,modul),'NEUTR','neu_memall',neutr)

  end if
  
end subroutine neu_memall
