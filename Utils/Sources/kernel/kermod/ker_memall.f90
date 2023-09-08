subroutine ker_memall()
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_memall
  ! NAME 
  !    ker_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    NS equations
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !    ker_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_memory
  use mod_maths,          only : maths_equalize_arrays
  use mod_chktyp,         only : check_type
  use mod_communications, only : PAR_SUM
  use mod_ker_subdomain,  only : ker_subdomain_motion_exists
  implicit none
  integer(ip) :: kfl_value,ipoin
  
  
  if( INOTMASTER ) then
     !
     ! WALLD: wall distance
     ! 
     if( kfl_walld /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD'    ,'ker_memall',walld,    npoin)
        if (kfl_walld == 2 .or. kfl_walld == 3) then 
          call memory_alloca(mem_modul(1:2,modul),'WALLO'    ,'ker_memall',wallo,    npoin)
          call memory_alloca(mem_modul(1:2,modul),'WALLCOOR'    ,'ker_memall',wallcoor,  ndime,  npoin)
        endif
        call memory_alloca(mem_modul(1:2,modul),'UWALL_KER','ker_memall',uwall_ker,npoin)
        if( kfl_delta == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'UWAL2_KER','ker_memall',uwal2_ker,npoin)
        end if
     end if
     !
     ! WALLN: wall normal
     ! 
     if( kfl_walln /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLN'    ,'ker_memall',walln,ndime,npoin)
     end if
     !
     ! ROUGH: Allocate roughness
     !
     if( kfl_rough == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'ROUGH','ker_memall',rough,npoin)
     end if
     !
     ! CANHE: Allocate canopy heigh 
     !
     if( kfl_canhe == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'CANHE','ker_memall',canhe,npoin)
     end if
     !
     ! Support boundary
     !
     if( kfl_suppo == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DISPL_KER','ker_memall',displ_ker,ndime,npoin)
     end if
     !
     ! adjoint optimization variables
     !
     if( kfl_adj_prob == 1 .and. kfl_dvar_type /= 5 .and. kfl_dvar_type /= 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS','ker_memall',sens,kfl_ndvars_opt)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 5) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,npoin)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     end if
     
!      call memory_alloca(mem_modul(1:2,modul),'UNTUR_FIX','ker_memall',untur_fix,npoin)
     !
     ! Advection vector
     ! KFL_VEFUN = 0 ... ADVEC points to VELOC 
     !           < 0 ... Copy values from a field
     !           > 0 ... User-defined function
     !
     if( kfl_vefun == 0 ) then
        advec => veloc
     else
        call memory_alloca(mem_modul(1:2,modul),'ADVEC','ker_memall',advec,ndime,npoin,3_ip)
        if( kfl_vefun < 0 ) then 
           kfl_value = -kfl_vefun
           call check_type(xfiel,kfl_value,ndime,npoin,1_ip) 
           call maths_equalize_arrays(ndime,npoin,xfiel(kfl_value) % a(:,:,1),advec)
           call maths_equalize_arrays(ndime,npoin,advec(1:ndime,1:npoin,1),advec(1:ndime,1:npoin,2))
           call maths_equalize_arrays(ndime,npoin,advec(1:ndime,1:npoin,1),advec(1:ndime,1:npoin,3))
        end if
     end if
     !
     ! Force based on residual for rigid body motion
     !
     if( kfl_forca_res == 1_ip  )then
        call memory_alloca(mem_modul(1:2,modul),'FORCA','ker_memall',forca,ndime,npoin,3_ip)
     end if
     !
     ! DISPM and VELOM
     !
     if( ker_subdomain_motion_exists() ) then
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , dispm , ndime , npoin , 1_ip )
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , velom , ndime , npoin )
     end if
     !
     ! no slip wall
     !
     if( kfl_noslw_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'EL_NSW_VISC'    ,'ker_memall',el_nsw_visc    ,nelem)
        call memory_alloca(mem_modul(1:2,modul),'FACT_NSW_KER'   ,'ker_memall',fact_nsw_ker   ,npoin)
     end if

     
  else
     !
     ! Wall distance
     !
     if( kfl_walld /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD'    ,'ker_memall',walld,    1_ip)
        if (kfl_walld == 2 .or. kfl_walld == 3) then
          call memory_alloca(mem_modul(1:2,modul),'WALLO'    ,'ker_memall',wallo,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'WALLCOOR'    ,'ker_memall',wallcoor,  ndime,  1_ip)
        endif  
        call memory_alloca(mem_modul(1:2,modul),'UWALL_KER','ker_memall',uwall_ker,1_ip)
        if( kfl_delta == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'UWAL2_KER','ker_memall',uwal2_ker,1_ip)
        end if
     end if
     !
     ! ROUGH: Allocate roughness
     !
     if( kfl_rough == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'ROUGH','ker_memall',rough,1_ip)
     end if
     !
     ! CANHE: Allocate canopy heigh 
     !
     if( kfl_canhe == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'CANHE','ker_memall',canhe,1_ip)
     end if
     !
     ! Support boundary
     !
     if( kfl_suppo == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DISPL_KER','ker_memall',displ_ker,1_ip,1_ip)
     end if
     !
     ! Fields prescribed by master
     !
     if( kfl_vefun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'ADVEC','ker_memall',advec,1_ip,1_ip,3_ip)
     end if
     !
     ! Force based on residual for rigid body motion
     !
     if( kfl_forca_res == 1_ip  )then
        call memory_alloca_min(forca)
     end if
     !
     ! adjoint optimization variables
     !
     if( kfl_adj_prob == 1 .and. kfl_dvar_type /= 5 .and. kfl_dvar_type /= 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS','ker_memall',sens,1_ip)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 5) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     end if
     !
     ! no slip wall
     !
     if( kfl_noslw_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'EL_NSW_VISC'    ,'ker_memall',el_nsw_visc    ,1_ip)   ! guess his is the correect for master        
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Memory for solver
  !
  !----------------------------------------------------------------------

  if( kfl_extro /= 0 ) then  
     solve_sol => solve(1:1)  ! Roughness extension
     call soldef(4_ip)
  end if
  if( kfl_walld /= 0 ) then
     solve_sol => solve(2:2)  ! Roughness extension
     call soldef(4_ip)
  end if
  if( kfl_suppo /= 0 ) then
     solve_sol => solve(3:3)  ! Support geometry for mesh multiplication 
     call soldef(4_ip)
  end if
  if( kfl_walln /= 0 ) then
     solve_sol => solve(5:5)  ! Wall normal
     call soldef(4_ip)
  end if
 ! if( kfl_defor /= 0 ) then
 !    solve_sol => solve(4:4)  ! Mesh deformation
 !    call soldef(4_ip)
 ! end if

end subroutine ker_memall
