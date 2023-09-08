!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis memory
!! @file    pts_memall.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine allocate memory for particles
!! @details Allocate memory and initialize particle type
!> @}
!------------------------------------------------------------------------

subroutine pts_memall()
  use def_parame
  use def_domain
  use def_master
  use def_kermod
  use def_partis
  use def_solver
  use def_inpout
  use mod_memory
  use mod_messages,       only : livinf
  implicit none
  integer(ip) :: itype,ivari,ivar2
  integer(8)  :: memorysize_for_particles  ! Type integer and of the system-dependent kind C_SIZE_T (from the ISO_C_BINDING module)
  type(latyp) :: particle_structure        ! Used to calculate the amount of memory allocated
  integer(4)  :: istat
  !
  ! Allocate memory: can only have MLAGR living at the same time
  !
  allocate( lagrtyp(mlagr), stat=istat )
  if(istat/=0) then
     !memorysize_for_particles = int(sizeof(particle_structure),8)*int(mlagr,8)
     memorysize_for_particles = 0_8 !size(particle_structure,KIND=8)*int(mlagr,8)
     call livinf(-7_ip,"Partis failed to allocate storage for particles (lagrtyp), Number of particles:",mlagr)
     call livinf(-7_ip,"Requested memory in bytes (lagrtyp):",int(memorysize_for_particles,ip))
     call runend("pts_memall: Particle storage allocation failure")
  end if
  
  lagrtyp(1:mlagr) = lagrtyp_init
  !
  ! Number of used types
  !
  ntyla_pts = 0
  do itype = 1,mtyla
     if( parttyp(itype) % kfl_exist == 1 ) ntyla_pts = max(ntyla_pts,itype)
  end do


  !
  ! Velocity field if nastin is off and the field is provided in the dat
  !
  if(veloc_field_id>0) then
     if(INOTMASTER) then
        if (memory_size(veloc)==0) then
           call memory_alloca(mem_modul(1:2,modul),'VELOC','pts_memall',veloc,ndime,npoin, pts_advec_narrays) !uses only 3 veloc fields
        end if
     else
        if (memory_size(veloc)==0) then
           call memory_alloca(mem_modul(1:2,modul),'VELOC','pts_memall',veloc,1_ip,1_ip, pts_advec_narrays)
        end if
     end if
  end if
  

  if( INOTMASTER ) then
     !
     ! Velocity deformation tensor
     !
     itype_loop: do itype = 1,ntyla_pts
        if( parttyp(itype) % kfl_exist /= 0 .and. parttyp(itype) % kfl_saffm /= 0 ) then
           call memory_alloca(mem_modul(1:2,modul),'DEFOR_PTS','pts_memall',defor_pts,ntens,npoin)
           exit itype_loop
        end if
     end do itype_loop
     !
     ! Wall element
     !
     call memory_alloca(mem_modul(1:2,modul),'LBOUE_PTS','pts_memall',lboue_pts,nelem)
     !
     ! Element natural length
     !
     call memory_alloca(mem_modul(1:2,modul),'HLENG_PTS','pts_memall',hleng_pts,nelem)
     !
     ! Distance to slip boundaries and friction coefficient
     !
     if( kfl_slip_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_SLIP_PTS','pts_memall',walld_slip_pts,npoin)
        call memory_alloca(mem_modul(1:2,modul),'FRICTION_PTS'  ,'pts_memall',friction_pts,npoin)
     end if
     !
     ! Distance to bouncing boundaries and friction coefficient
     !
     if( kfl_bouncing_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_BOUNCING_PTS','pts_memall',walld_bouncing_pts,npoin)
     end if
     !
     ! Two-way coupling: momentum is accumulated
     !
     if( kfl_coupl(ID_NASTIN,ID_PARTIS) /= 0  ) then
        call memory_alloca(mem_modul(1:2,modul),'MOMEN','pts_memall',momen,ndime,npoin_2)
     end if

  else

     if( kfl_slip_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_SLIP_PTS','pts_memall',walld_slip_pts,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'FRICTION_PTS'  ,'pts_memall',friction_pts,1_ip)
     end if
     if( kfl_bouncing_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_BOUNCING_PTS','pts_memall',walld_bouncing_pts,1_ip)
     end if

  end if
  !
  ! Residence time
  !
  ivari = 5
  call posdef(25_ip,ivari)
  if( ivari > 0 ) then
     kfl_resid_pts = 1
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'RESID_PTS','pts_memall',resid_pts,ntyla_pts,nelem_2)
     end if
  end if
  !
  ! Deposition
  !
  ivari = 2
  call posdef(25_ip,ivari)
  ivar2 = 7
  call posdef(25_ip,ivar2)
  if( ivari > 0 .or. ivar2 > 0 ) then
     kfl_depos_pts = 1
     if( ntyla_pts > 0 ) then
        if( INOTMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'DEPOE_PTS','pts_memall',depoe_pts,ntyla_pts,nelem_2)
           call memory_alloca(mem_modul(1:2,modul),'DEPOB_PTS','pts_memall',depob_pts,ntyla_pts,nboun_2)
        else
           call memory_alloca(mem_modul(1:2,modul),'DEPOE_PTS','pts_memall',depoe_pts,ntyla_pts,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'DEPOB_PTS','pts_memall',depob_pts,ntyla_pts,1_ip)
        end if
     end if
  end if
  !
  ! Deposition surface
  !
  if( kfl_depos_pts == 0 .and. kfl_depos_surface_pts == 1 .and. kfl_oudep_pts /= 0 ) then
     if( ntyla_pts > 0 ) then
        if( INOTMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'DEPOB_PTS','pts_memall',depob_pts,ntyla_pts,nboun_2)
        else
           call memory_alloca(mem_modul(1:2,modul),'DEPOB_PTS','pts_memall',depob_pts,ntyla_pts,1_ip)
        end if
     end if
  end if
  !
  ! Postprocess
  !  
  nvarp_pts = count(postprocess_var_pts)
  call memory_alloca(mem_modul(1:2,modul),'POSTPROCESS_LIST_PTS','pts_output',postprocess_list_pts,nvarp_pts)
  ivar2 = 0
  do ivari = 1,mvarp_pts
     if( postprocess_var_pts(ivari) ) then
        ivar2 = ivar2 + 1
        postprocess_list_pts(ivar2) = ivari
     end if
  end do
  nvard_pts = count(deposition_var_pts)
  call memory_alloca(mem_modul(1:2,modul),'DEPOSITION_LIST_PTS','pts_output',deposition_list_pts,nvard_pts)
  ivar2 = 0
  do ivari = 1,mvarp_pts
     if( deposition_var_pts(ivari) ) then
        ivar2 = ivar2 + 1
        deposition_list_pts(ivar2) = ivari
     end if
  end do
  !
  ! Solver
  !
  solve_sol                => solve(1:1)
  solve_sol(1) % kfl_fixno => kfl_fixno_walld_slip_pts
  call soldef(4_ip)

  solve_sol                => solve(2:2)
  solve_sol(1) % kfl_fixno => kfl_fixno_walld_bouncing_pts
  call soldef(4_ip)
  
end subroutine pts_memall
