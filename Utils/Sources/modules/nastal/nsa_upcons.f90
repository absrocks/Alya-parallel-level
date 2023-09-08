subroutine nsa_upcons
!-----------------------------------------------------------------------
!****f* Nastal/nsa_upcons
! NAME 
!    nsa_upcons
! DESCRIPTION
!    Update conservative set
! USES
!    nsa_...
! USED BY
!    nsa_gocons
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_solver
  use      mod_solver,  only : solver_solve
  use mod_commdom_alya, only: CPLNG
  implicit none
  integer(ip) :: iauxi
  real(rp)    :: xdumy,elsou(nevat_nsa)

  call nsa_updunk(8_ip) ! unkno <-- u(,ITER_K), dunkn <--- 0.0_rp

!!  if (kfl_timul_nsa /= 0) then     ! update and correct residuals according to the time 
!!   call nsa_uptimi(4_ip)         ! integration used
!!  end if


  if (kfl_resmo_nsa == 1) then
     !
     ! Residual smoothing
     !    
     call smoot4(ndofn_nsa,rhsid,vmass,resmo_nsa)
     call nsa_updunk(12_ip) ! rhsid <-- resmo_nsa

  end if
  !
  ! Call the solver
  !
  if (kfl_timet_nsa == 1) then


     if (kfl_mod_elmop_nsa == 0)  then
        !
        ! OLD EXPLICIT SCHEMES
        !
        solve(1)%xdiag =  1.0_rp/dtinv_nsa
        iauxi= kfl_dttyp_nsa(1)+kfl_dttyp_nsa(2)+kfl_dttyp_nsa(3)
        if (kfl_diagi_nsa == 1 .or. iauxi > 0) then
           call solver(rhsid,unkno,xdumy,vdiag_nsa)
        else if (kfl_diagi_nsa == 0) then        
           call solver(rhsid,unkno,xdumy,xdumy)        
        end if
        
        if (kfl_unkse_nsa == 1) then
           !
           ! Compute cons variables from primitive ones, storing cons in unkno
           !
           call nsa_setvar(4_ip , 1_ip)
        end if
        
        !
        ! Compute prim variables from cons ones, taking into account boundary conditions
        !
        call nsa_setvar(3_ip , 1_ip)
        
     else
        !
        ! NEW EXPLICIT SCHEMES
        !
        solve(1)%xdiag =  1.0_rp/dtinv_nsa  
        call solver_solve(momod(modul) % solve, amatr, rhsid, dunkn_nsa, vdiag_nsa)

        !
        ! Rotate (generalized) back boundary conditions when needed
        call nsa_roback
        !
        ! Correct unkno, because implicit is in the delta form
        call nsa_updunk(11_ip)  ! Correct unkno, because new explicit is in the delta form
        
        ! final conversions are done now in nsa_computephysical

        
     end if


  else if (kfl_timet_nsa == 2) then
     
     ! IMPLICIT (AND MATRIX-EXPLICIT) SCHEMES
     
     !
     ! Solve
     !

!     if (kfl_delun_nsa == 0) then
!        call runend('NSA_UPCONS: IMPLICIT SCHEMES ARE IN DELTA FORM')
!        call solver(rhsid,unkno,amatr,pmatr)
!        !
!        ! Rotate (generalized) back boundary conditions when needed
!        !
!        call nsa_roback
!
!     else if (kfl_delun_nsa == 1) then
        !
        ! dunkn is the unknown in the delta form
!        call solver(rhsid,dunkn_nsa,amatr,pmatr)

     call solver_solve(momod(modul) % solve, amatr, rhsid, dunkn_nsa, pmatr)
     
     !
     ! Rotate (generalized) back boundary conditions when needed
     call nsa_roback
     !
     ! Correct unkno, because implicit is in the delta form
     call nsa_updunk(11_ip)
     !
     

     !
     ! Compute prim variables from cons ones
     !

!!!! probando nsa_computephysical, llamada en endite
!!     call nsa_setvar(5_ip , 1_ip)


  end if

  ! acaaaa ojooooo para shofino
!  write(194,1000) ittim,ittim,itinn(modul),unkno(5),unkno(7),unkno(8),cutim
!  write(194,1000) ittim,ittim,itinn(modul),unkno(398*4+1),unkno(399*4+1),unkno(400*4+1),cutim
!1000 format(4x,i9,2x,i9,2x,i9,20(2x,e12.6))



! para minicucu
!  if (itinn(modul)==1) then
!     write(6677,100) dunkn_nsa(5),unkno(5),densi(2,TIME_N),densi(2,ITER_K)
!  end if
!100 format(10(2x,e))


end subroutine nsa_upcons
