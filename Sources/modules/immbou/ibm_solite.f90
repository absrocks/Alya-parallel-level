subroutine ibm_solite(itask)
  !-----------------------------------------------------------------------
  !****f* ibm_solite/ibm_solite
  ! NAME
  !    ibm_solite
  ! DESCRIPTION
  !    This routines solves the Euler and Newton equations for rigid bodies
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ittib
  integer(ip)             :: iimbo,iwaib,idime

  if( itask == 1 ) then 
     !
     ! Counters
     !
     itinn(modul) = itinn(modul) + 1
     !
     ! Begin iteration
     ! 
     if( itcou == 1 ) call ibm_tistep()     
     !
     ! Add forces: gravity, magnetic: F^n+1 and T^n+1
     !     
     call ibm_forces()
     !
     ! Solver Euler's equations at t^n+x
     !     
     call ibm_eulers(dtime)     
     !-------------------------------------------------------------------
     !
     ! Collisions: Internal time advance
     !
     !-------------------------------------------------------------------
     do iimbo = 1, nimbo
        imbou(iimbo) % cotim = -1.0e10_rp
        do idime = 1,3
           imbou(iimbo) % accel(idime,3) = imbou(iimbo) % accel(idime,2)
           imbou(iimbo) % accea(idime,3) = imbou(iimbo) % accea(idime,2)
        end do
     end do
     
     if( kfl_colli_ibm /= 0 ) then
        call livinf(-4_ip,'IMMBO: START COLLISIONS',0_ip)
        cutim_ibm = cutim - dtime ! Go back in time
        dtime_ibm = dtime

        ittib     = 0
        do iwaib = 1,nwaib
           twall_ibm(iwaib) % cotim = -1.0e10_rp
        end do
        
        do while( cutim_ibm < cutim )
           !
           ! Time step size for prevent collision: dtime_ibm
           !
           call ibm_coldet(ittib)
           !
           ! Newmark scheme
           !
           call ibm_newmak(dtime_ibm)
           !
           ! Resolve a collision if there is one
           !
           call ibm_colres(ittib)
           !
           ! Update the angular and linear displacements and velocities of the particles
           !        
           call ibm_updunk(5_ip)
           !
           ! Update time
           !
           cutim_ibm = cutim_ibm + dtime_ibm
           dtime_ibm = cutim - cutim_ibm     
           ittib     = ittib + 1
           routp(1)  = cutim_ibm
           if( kfl_colli_ibm /= 0 ) call livinf(98_ip,'TIME= ',0_ip)
        end do
        
        call livinf(-5_ip,'IMMBO: END COLLISIONS',0_ip)

     !-------------------------------------------------------------------
     !
     ! Collisions: End internal time advance
     !
     !-------------------------------------------------------------------
     else
        !
        ! Newmark scheme
        !
        call ibm_newmak(dtime)
        !
        ! Update the angular and linear displacements and velocities of the particles
        !        
        call ibm_updunk(5_ip)   
     end if
  end if
  !
  ! Coupling with the mesh
  !
  call ibm_insout()
  !
  ! Interpolation
  !
  call ibm_linter()
  !
  ! IB: Identify IB Gauss points host elements
  !   
  !call ibm_chgaib()
  !
  ! Mass conservation arrays
  !
  if( kfl_diric_ibm > 0 ) call ibm_massma()

  call ibm_coupli(ITASK_ENDITE)
   
  if( ittim /= 0 ) call livinf(164_ip,' ',1_ip) 
  !
  ! Compute mass matrices as new holes have been created
  ! Problem with EXNOR and LPOTY. Local axes could change because
  ! of the presence of holes. Therefore we would have to recompute
  ! all NBOPO depending arrays: we do not do this NOW
  !
  if( kfl_diric_ibm > 0 ) then
     call domarr(2_ip)
  end if

end subroutine ibm_solite
