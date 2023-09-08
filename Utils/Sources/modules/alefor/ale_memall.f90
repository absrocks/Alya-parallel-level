!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_memall.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Arrays allocation subroutine
!> @details Arrays allocation subroutine
!> @} 
!-----------------------------------------------------------------------
subroutine ale_memall
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      def_alefor
  use      mod_memchk
  implicit none
  integer(ip)    :: ipoin,idime
  integer(4)     :: istat
  !
  ! ALE unknown
  !
  if( INOTMASTER ) then
     !
     ! Unknowns
     !  
     allocate(velom(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOM','ale_memall',velom) 
     ! -----------------------------------------------------------------------
     ! Change to keep the values of previous iterations and time step of dispm
     allocate(dispm(ndime,npoin,3),stat=istat)
     !      allocate(dispm(ndime,npoin,2),stat=istat)
     ! -----------------------------------------
     call memchk(zero,istat,mem_modul(1:2,modul),'DISPM','ale_memall',dispm)
     ! ---------------------------------------------------------------------------
     ! Change to keep the values of previous iterations and time step of coord_ale
     allocate(coord_ale(ndime,npoin,3),stat=istat)
     !      allocate(coord_ale(ndime,npoin,2),stat=istat)
     ! ---------------------------------------------
     call memchk(zero,istat,mem_modul(1:2,modul),'COORD_ALE','ale_memall',coord_ale)
     ! ---------------------------------------------------------------------------
     ! Addition to keep the values of previous iterations and time step of coord_ale
     allocate(coord_ori(ndime,npoin),stat=istat)
     ! ---------------------------------------------
     call memchk(zero,istat,mem_modul(1:2,modul),'COORD_ORI','ale_memall',coord_ori)
     allocate(bvess_ref(ndime,npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_REF','ale_memall',bvess_ref)

     !
     ! Save initial coordinates
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           coord_ale(idime,ipoin,1) = coord(idime,ipoin)
           coord_ale(idime,ipoin,2) = coord(idime,ipoin)
           ! -----------------------------------------------------------------------------
           ! Addition to keep the values of previous iterations and time step of coord_ale
           coord_ale(idime,ipoin,3) = coord(idime,ipoin)
           coord_ori(idime,ipoin)   = coord(idime,ipoin)
        end do
     end do

  else

     allocate(velom(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOM','ale_memall',velom) 
     ! -----------------------------------------------------------------------
     ! Change to keep the values of previous iterations and time step of dispm
     allocate(dispm(1,1,3),stat=istat)
     !      allocate(dispm(1,1,2),stat=istat)
     ! -------------------------------
     call memchk(zero,istat,mem_modul(1:2,modul),'DISPM','ale_memall',dispm)
     ! ---------------------------------------------------------------------------
     ! Change to keep the values of previous iterations and time step of coord_ale
     allocate(coord_ale(ndime,1,3),stat=istat)
     !      allocate(coord_ale(ndime,1,2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'COORD_ALE','ale_memall',coord_ale)
     ! ---------------------------------------------
     ! Addition to keep the values of previous iterations and time step of coord_ale
     allocate(coord_ori(1,1),stat=istat)
     ! ---------------------------------------------
     call memchk(zero,istat,mem_modul(1:2,modul),'COORD_ORI','ale_memall',coord_ori)
     allocate(bvess_ref(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_REF','ale_memall',bvess_ref)

  end if
  !
  ! Solver memory
  !
  solve_sol => solve(1:1)
  call soldef(4_ip)  
  solve(1) % bvess     => bvess_ale
  solve(1) % kfl_fixno => kfl_fixno_ale

  !
  ! Save initial displacement directions (for motion in multiple directions)
  !
  if( INOTMASTER ) then 
     do ipoin = 1,npoin
        do idime = 1,ndime
           bvess_ref(idime,ipoin) = bvess_ale(idime,ipoin)
        end do
     end do
  endif

end subroutine ale_memall
      
