subroutine tur_maxlen()

  !*****************************************************
  !tur_maxlen:
  ! Used by tur_endste
  ! This routine calculates de maximum mixing length using Mellor Yamada
  ! Max mixing length = int (sqrt(k)z dz) /int (sqrt(k)dz) 
  !*****************************************************
! use of global variables and subroutines
  use def_parame
  use def_master,         only : rhsid, amatr,INOTMASTER,IPARALL, nparr, parre, unkno
  use def_master,         only : solve, mem_modul, modul, untur
  use def_domain,         only : meshe, elmar, coord, npoin,ndime
  use mod_ADR,            only : ADR_assemble_convective, ADR_assemble_extension
  use mod_memory,         only : memory_alloca, memory_deallo
  use def_solver,         only : solve_sol
  use mod_solver,         only : solver_solve
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use def_turbul,         only : kfl_fixno_tur, TUR_FAMILY_K_EPS, tur_max_mixlen
  use def_turbul,         only : numer, denom, z_wall, gravi_tur, kfl_cotem_tur
  use def_kermod,         only : ndivi
  use mod_communications, only : PAR_MAX, PAR_MIN
  implicit none  
  ! local variables
  integer(ip)             :: ipoin, idime
  integer(ip), save       :: ipass=0
  real(rp),    save       :: z_top, z_min         ! maximum height
  real(rp),    target     :: rmaxi(1)  
  integer(ip), pointer    :: kfl_fixno_tmp(:,:)
  real(rp),    pointer    :: bvess_tmp(:,:)
  real(rp) ,   pointer    :: conve_tmp(:,:)
!  real(rp),    pointer    :: numer(:), denom(:), z_wall(:)
  real(rp),    pointer    :: force(:)

  integer(ip)            ::  vcoor 
  vcoor = 3 ! vertical coordinate

  if (kfl_cotem_tur .gt. 0) then !if temper coupling, set vertical
                                 !coordinate in terms of gravity
     do idime =1, ndime
        if (abs(gravi_tur(idime)).gt.0.5_rp) vcoor = idime    
     end do
  end if


  ! calculate z_top of domain (only once)
  if(ipass==0) then ! to be done only once ! calculates maximum height
     z_top = -1.0e25_rp
     z_min = 1e20_rp
       do ipoin=1, npoin
           z_top = max(coord(vcoor,ipoin),z_top)
           z_min = min(coord(vcoor,ipoin),z_min)
        end do
     call PAR_MAX(z_top) 
     call PAR_MIN(z_min)
 ! end if
  !print *, 'z_top', z_top
  !print *, 'z_min', z_min
  
  !
  ! Allocate structures 
  !  
     nullify(z_wall)
     nullify(numer)
     nullify(denom)
     call memory_alloca(mem_modul(1:2,modul),'Z_WALL ','tur_endste',z_wall,npoin)
     call memory_alloca(mem_modul(1:2,modul),'NUMER','tur_endste',numer ,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DENOM','tur_endste',denom ,npoin)
  end if ! ipass ==0 ! first step

  nullify(kfl_fixno_tmp)
  nullify(bvess_tmp)
  nullify(conve_tmp)
  nullify(force)
  call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_TMP','tur_endste',kfl_fixno_tmp,1_ip ,npoin)
  call memory_alloca(mem_modul(1:2,modul),'BVESS_TMP','tur_endste',bvess_tmp    ,1_ip ,npoin)
  call memory_alloca(mem_modul(1:2,modul),'CONVE_TMP','tur_endste',conve_tmp ,ndime,npoin)
  call memory_alloca(mem_modul(1:2,modul),'FORCE','tur_endste',force ,npoin)

 ! if (INOTMASTER) then
     ! solver 
  solve_sol                => solve(3:)
  solve_sol(1) % kfl_fixno => kfl_fixno_tmp 
  solve_sol(1) % bvess     => bvess_tmp
  solve_sol(1) % kfl_iffix =  1  
  if(ipass==0) then 
     
     do ipoin =1,npoin
        conve_tmp(1:3,ipoin) = 0.0_rp ! convection in vertical direction
        conve_tmp(vcoor,ipoin) = 1.0_rp ! from bottom to top
        force(ipoin) = 1.0_rp
        kfl_fixno_tmp(1,ipoin) =  0
        unkno(ipoin) = 0.0_rp ! coord(3, ipoin)
        if(kfl_fixno_tur(1,ipoin,2) == 3) then
           bvess_tmp(1,ipoin) = 0.0_rp  !coord(3, ipoin)
           kfl_fixno_tmp(1,ipoin) =  1
        end if
     end do
     call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)  
     if (INOTMASTER) &
          call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
     call solver_solve(solve_sol,amatr,rhsid,unkno)
     !  end if
     if (INOTMASTER) z_wall(1:npoin) = unkno(1:npoin)
  end if
  
  ! integrate from wall to top
  do ipoin = 1,npoin
     conve_tmp(1:3,ipoin) = 0.0_rp ! convection in vertical direction
     conve_tmp(vcoor,ipoin) = 1.0_rp ! from bottom to top
     kfl_fixno_tmp(1,ipoin) = 0      ! boundary condition free
     unkno(ipoin) = numer(ipoin) ! 0.0_rp
     !  numer(ipoin) = 0.0_rp  ! Unknown sqrt(k) z dz 
     !  denom(ipoin) = 0.0_rp  ! Unknown sqrt(k) dz
     if(kfl_fixno_tur(1,ipoin,2) == 3) then ! if wall law for eps
        kfl_fixno_tmp(1,ipoin) =  1
        bvess_tmp(1,ipoin)     =  0.0_rp  ! boundary value
     end if
  end do

  if (INOTMASTER) force(1:npoin) = sqrt(untur(1,1:npoin,1))*z_wall(1:npoin) !  coord(3,1:npoin)


  ! assemble matrix and rhs for numer
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid) 
  call solver_solve(solve_sol,amatr,rhsid,unkno)
  if (INOTMASTER) numer (1:npoin) =  unkno(1:npoin)

  ! assemble matrix and rhs  for denom

  if (INOTMASTER) then
     force(1:npoin) = sqrt(untur(1,1:npoin,1))
     unkno(1:npoin) = denom(1:npoin) !0.0_rp
  end if
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
  call solver_solve(solve_sol,amatr,rhsid,unkno)
  if (INOTMASTER) denom (1:npoin)= unkno(1:npoin)



  ! topography height ! integral from botton with rhs 0, only done once

  if (ipass==0) then
     do ipoin =1, npoin
        if (kfl_fixno_tur(1,ipoin,2) /= 3) then ! not wall
           tur_max_mixlen(ipoin) = 0.075_rp*(numer(ipoin)/denom(ipoin) )
        else
           tur_max_mixlen(ipoin) = 0.0_rp
        end if
     end do
  end if


  if (INOTMASTER) unkno(1:npoin) =  tur_max_mixlen(1:npoin)
  !     if (kfl_fixno_tur(1,ipoin,2) /= 3) & ! initial condition, denom .gt. 0
  !          tur_max_mixlen(ipoin) = 0.075_rp*(numer(ipoin)/denom(ipoin) ) !-z_wall(ipoin))     



  !
  ! project solution on the top vertically to all domain
  !
  if (INOTMASTER) then
     conve_tmp(1:3,1:npoin) = 0.0_rp
     conve_tmp(vcoor,1:npoin) = - 1.0_rp  ! vertical convection from top to bottom
  end if
  do ipoin =1, npoin
     kfl_fixno_tmp(1,ipoin) =  0
     force(ipoin) = 0.0_rp
     if (abs(coord(vcoor, ipoin)-z_top)<0.1_rp ) then ! prescribe top of domain
        kfl_fixno_tmp(1,ipoin) =  1 
        unkno(ipoin)           =   0.075_rp*numer(ipoin)/denom(ipoin)
        bvess_tmp(1,ipoin)     =  unkno(ipoin) !tur_max_mixlen(ipoin)
     end if
  end do
  call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)  
  if (INOTMASTER) &
       call ADR_assemble_convective(1_ip,meshe(ndivi),elmar,conve_tmp,amatr,force,rhsid)
  call solver_solve(solve_sol,amatr,rhsid,unkno)

  if (INOTMASTER) tur_max_mixlen(1:npoin) =  unkno(1:npoin)
  !  if (INOTMASTER) tur_max_mixlen(1:npoin) =  z_wall(1:npoin)
  !
  ! Deallocate structures 
  !  
  nullify(solve_sol(1) % kfl_fixno) 
  nullify(solve_sol(1) % bvess    )
  if (INOTMASTER) then
     deallocate (kfl_fixno_tmp)
     deallocate (bvess_tmp)
     !  deallocate (z_wall)
     deallocate (conve_tmp)
     !  deallocate (numer)
     !  deallocate (denom)
     deallocate (force)
  end if


  ipass = 1
end subroutine tur_maxlen
