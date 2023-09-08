subroutine rad_memall()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_memall
  ! NAME 
  !    rad_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    heat transfer equation
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_radiat
  use mod_memchk
  implicit none
  integer(ip) :: ielem,pelty,pgaus,dummi
  integer(4)  :: istat
  !
  ! Problem unknowns Heat Radiation Source, and whatever method is used to obtain it
  !   (In first version the method is the P1 equation so that the radiation average is
  !    allocated, not the intensity). 
  !    Also some solver initialization is done
  !
  if( INOTMASTER ) then
     !
     ! RADSO: Heat radiation source 
     ! 
     allocate(radso(npoin,ncomp_rad),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RADIAT','rad_memall',radso)
     !
     ! AVRAD_RAD: Average radiation intensity
     !
!     if(postp(1) % npp_stepi(5)>0) then
        allocate(radav_rad(npoin,ncomp_rad),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'radav_RAD','rad_memall',radav_rad)        
!     end if
     !
     ! RASGS: Subgrid scale radiation source term 
     !
     if(postp(1) % npp_stepi(4) /= 0 .or.  kfl_sgsno_rad == 1 ) then     
        allocate(rasgs_rad(nelem),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RASGS','rad_memall',rasgs_rad)
        do ielem=1,nelem
           pelty=ltype(ielem)
           pgaus=ngaus(pelty)
           allocate(rasgs_rad(ielem)%a(mgaus,2),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'RASGS','rad_memall',rasgs_rad(ielem)%a)
        end do
     end if
     !
     ! RAPRO: Projections for orthogonal SGS
     !
!!$     if( kfl_ortho_rad /= 0 ) then
!!$        allocate(rapro_rad(npoin),stat=istat)
!!$        call memchk(zero,istat,mem_modul(1:2,modul),'RAPRO_TEM','rad_memall',rapro_rad)    
!!$     end if     
     !
     ! GRTEM: Radiation gradients
     !
     if(kfl_ellen_rad==-1) then
        allocate(grrad_rad(ndime,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRRAD_RAD','rad_memall',grrad_rad)        
     end if
     !
     ! RAOLD_RAD: Old radiation
     !
     if(postp(1) % npp_stepi(8)>0) then
        allocate(raold_rad(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RAOLD_RAD','rad_memall',raold_rad)        
     end if
     !
     !  Couple to other modules to get density, temperature, and concentrations
     !
     call rad_couple()
     !
     !  If this is a test case we must allocate our own memory
     !
     if(kfl_atest_rad/=0_ip) then
        allocate(tempe_rad(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RADIAT','tempe_rad',tempe_rad)
     end if
     !
     ! Solver memory
     !
     solve_sol => solve(1:)
     call soldef(4_ip)

  else

     allocate(radso(1,3),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'RADSO','rad_memall',radso)
     if(postp(1) % npp_stepi(4)/=0.or.kfl_sgsno_rad==1) then     
        allocate(rasgs_rad(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RASGS','rad_memall',rasgs_rad)
     end if
!     if(postp(1) % npp_stepi(5)>0) then
        allocate(radav_rad(1,3),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RADAV_RAD','rad_memall',radav_rad)        
!     end if
     if(kfl_ellen_rad==-1) then
        allocate(grrad_rad(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRRAD_RAD','rad_memall',grrad_rad)        
     end if
     if(postp(1) % npp_stepi(8)>0) then
        allocate(raold_rad(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RAOLD_RAD','rad_memall',raold_rad)        
     end if
     !
     !  Couple to other modules to get density, temperature, and concentrations
     !
     call rad_couple()
     !
     !  Wether or not we are in a test case
     !
     if(kfl_atest_rad/=0_ip) then
        allocate(tempe_rad(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'RADIAT','tempe_rad',tempe_rad)
     endif
  end if

end subroutine rad_memall
      
