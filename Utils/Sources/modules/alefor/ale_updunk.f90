!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_updunk.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Updates.
!> @details Updates.
!> @} 
!-----------------------------------------------------------------------
subroutine ale_updunk(itask)
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use def_kermod,         only : kfl_adj_prob
  implicit none
  integer(ip), intent(in) :: itask  !> where the subrutine is called
  integer(ip)             :: ipoin,idime,iimbo,kpoin,itotn

  real(rp),    pointer    :: force(:,:),accel(:,:),velol(:,:),posil(:,:)          ! Linear motion
  real(rp),    pointer    :: torqu(:,:),accea(:,:),veloa(:,:),posia(:,:),quate(:,:),q_dot(:,:) ! Angular motion
  real(rp),    pointer    :: vpfor(:,:),vptor(:,:)
  !
  ! For the adjoint case, not to enter here
  !
  if (kfl_adj_prob == 1_ip) return

  iimbo = 1_ip  ! for the moment only one RB

  select case (itask)

  case(1_ip)
     !
     ! Initial guess for dispm: d(n,0,*) <-- d(n-1,*,*)
     !
     if ( INOTMASTER ) then

        do ipoin = 1, npoin

           do idime = 1, ndime

              dispm(idime,ipoin,2) = dispm(idime,ipoin,3)

           end do

        end do
        ! print*, "DEBUG: case(1) Initialization outer loop"
     endif

  case(2_ip)
     !
     ! Initial guess for dispm to be used in the inner iterations 
     !
     if ( INOTMASTER ) then

        do ipoin = 1, npoin

           do idime = 1, ndime

              dispm(idime,ipoin,1) = dispm(idime,ipoin,2)

           end do

        end do
        ! print*, "DEBUG: case(2) Initialization inner loop"
     endif

  case(3_ip)

     if( INOTMASTER ) then

        !-------------------------------------------------------------------
        ! 
        ! Update displacement, mesh velocity and new mesh coordinate
        ! FMALE: 1. does not update coordinate
        !        2. Invert mesh velocity
        !
        !-------------------------------------------------------------------

        do ipoin = 1,npoin
           itotn = (ipoin-1_ip) * ndime ! Before in deform_deform

           do idime = 1,ndime

              if( kfl_fixno_ale(idime,ipoin) == -1 .or. kfl_fixno_ale(idime,ipoin) == 3 ) then

                 !
                 ! FMALE type
                 !
                 ! ---------------------------------------------------------------------------
                 ! Change to use the indexes of coord_ale and dispm in the same way as in the 
                 ! other modules
                 !
                 itotn                    = itotn + 1_ip
                 dispm(idime,ipoin,1)     = unkno(itotn)
                 coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
                 velom(idime,ipoin)       =-dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
                 !                  dispm(idime,ipoin,1)     =  coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,2)
                 !                  velom(idime,ipoin)       = -dtinv * ( coord_ale(idime,ipoin,1) - coord(idime,ipoin) )
                 bvess_ale(idime,ipoin)   =  0.0_rp                 
              else
                 !
                 ! Normal type
                 !
                 itotn                    = itotn + 1_ip
                 dispm(idime,ipoin,1)     = unkno(itotn)
                 coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,1) + dispm(idime,ipoin,1)
                 velom(idime,ipoin)       = dtinv * ( coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,3) )
                 coord(idime,ipoin)       = coord_ale(idime,ipoin,1)
                 !                  dispm(idime,ipoin,1)   =  coord_ale(idime,ipoin,1) - coord_ale(idime,ipoin,2)
                 !                  velom(idime,ipoin)     =  dtinv * ( coord_ale(idime,ipoin,1) - coord(idime,ipoin) )
              end if

           end do

        end do
        ! print*, "DEBUG: case(3) Doiter"
     end if

  case(4_ip)
     !
     ! Look case 24
     !

  case(5_ip)

     if( INOTMASTER ) then

        if ( kfl_solve_ale == 1_ip ) then

           do ipoin = 1, npoin

              do idime = 1, ndime

                 dispm(idime,ipoin,3) = dispm(idime,ipoin,1)
                 coord_ale(idime,ipoin,3) = coord_ale(idime,ipoin,1)

              end do

           end do
           ! print*, "DEBUG: case(5) End step"
        else

           !-------------------------------------------------------------------
           !
           ! When ALE is not solved ( kfl_solve_ale /= 1 ) put VELOM and DISPM to zero
           !
           !-------------------------------------------------------------------

           do ipoin = 1,npoin
              do idime = 1,ndime
                 dispm(idime,ipoin,1)  =  coord(idime,ipoin) - coord_ale(idime,ipoin,3)
                 velom(idime,ipoin)    =  0.0_rp                 
              end do
           end do
           ! print*, "DEBUG: case(5) End step not solved"           
        end if

     end if

  case(6_ip)   ! For global iterations in coupling schemes, FSI 
     if( INOTMASTER ) then

        do ipoin = 1, npoin

           do idime = 1, ndime

              ! dispm(idime,ipoin,1)     = dispm(idime,ipoin,3)
              coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,3)
              coord(idime,ipoin)       = coord_ale(idime,ipoin,1)

           end do

        end do

     end if
     !print*, "DEBUG: case(6) global iterations"

  case(21_ip)    !  In (20 + i) I will put the 'equivalent' to i from nsi
     !
     ! Assign a(n,0,*) <-- a(n-1,*,*),  RB initial guess for outer iterations
     ! For the Force and Moment this might be a good place to extrapolate, force( ,nprev_ale) and torqu( ,nprev_ale)
     ! must have been set in nsi after solving if in the previous step - This supposes ale is solved before
     ! nsi if this is not the case some other strategy must be thought
     !
     if ( kfl_rigid_ale == 0 ) return

     if ( kfl_genco_ale == 1_ip ) then
        !
        ! For the moment, only angular coordinates are being used as generalized coordinates
        !
        posia =>  rbbou(iimbo) % posia
        veloa =>  rbbou(iimbo) % veloa
        vpfor =>  rbbou(iimbo) % vpfor
        
        do idime = 1_ip, 3_ip
           
           posia(idime,2_ip)   = posia(idime,nprev_ale)
           veloa(idime,2_ip)   = veloa(idime,nprev_ale)
           vpfor(idime,2_ip)   = vpfor(idime,nprev_ale)
           
        end do

     else

        accel =>  rbbou(iimbo) % accel
        velol =>  rbbou(iimbo) % velol
        posil =>  rbbou(iimbo) % posil
        accea =>  rbbou(iimbo) % accea
        veloa =>  rbbou(iimbo) % veloa
        posia =>  rbbou(iimbo) % posia
        quate =>  rbbou(iimbo) % quate
        q_dot =>  rbbou(iimbo) % q_dot
        vpfor =>  rbbou(iimbo) % vpfor
        vptor =>  rbbou(iimbo) % vptor
        !
        quate(1,2) = quate(1,nprev_ale)
        q_dot(1,2) = q_dot(1,nprev_ale)
        do idime = 1,3
           accel(idime,2)   = accel(idime,nprev_ale) 
           velol(idime,2)   = velol(idime,nprev_ale) 
           posil(idime,2)   = posil(idime,nprev_ale)
           accea(idime,2)   = accea(idime,nprev_ale) 
           veloa(idime,2)   = veloa(idime,nprev_ale) 
           posia(idime,2)   = posia(idime,nprev_ale)                 
           quate(idime+1,2) = quate(idime+1,nprev_ale)
           q_dot(idime+1,2) = q_dot(idime+1,nprev_ale)
        end do

        if ( ( kfl_foexo_ale == 1_ip ) .or. ( kfl_crist_ale == 1_ip ) ) then ! Force ( & Torque) extrapolation order
           do idime = 1,3
              vpfor(idime,2)   = 1.0_rp * vpfor(idime,nprev_ale)
              vptor(idime,2)   = 1.0_rp * vptor(idime,nprev_ale)
           end do
        else
           do idime = 1,3
              vpfor(idime,2)   = 1.5_rp * vpfor(idime,nprev_ale) - 0.5_rp * vpfor(idime,nprev_ale+1_ip)
              vptor(idime,2)   = 1.5_rp * vptor(idime,nprev_ale) - 0.5_rp * vptor(idime,nprev_ale+1_ip)
           end do
        end if

     end if  ! kfl_genco_ale

     ! print*, "DEBUG: case(21) endite"
  case(22_ip)
     !
     ! Assign a(n,i,0) <-- a(n,i-1,*), RB initial guess for inner iterations
     !
     if ( kfl_rigid_ale == 0 ) return

     if ( kfl_genco_ale == 1_ip ) then
        !
        ! For the moment, only angular coordinates are being used as generalized coordinates
        !
        posia =>  rbbou(iimbo) % posia
        veloa =>  rbbou(iimbo) % veloa
        vpfor =>  rbbou(iimbo) % vpfor
        
        do idime = 1_ip, 3_ip

           posia(idime,1_ip)  = posia(idime,2_ip)
           veloa(idime,1_ip)  = veloa(idime,2_ip)
           vpfor(idime,1_ip)  = vpfor(idime,2_ip)
           
        end do

     else

        accel =>  rbbou(iimbo) % accel
        velol =>  rbbou(iimbo) % velol
        posil =>  rbbou(iimbo) % posil
        accea =>  rbbou(iimbo) % accea
        veloa =>  rbbou(iimbo) % veloa
        posia =>  rbbou(iimbo) % posia
        quate =>  rbbou(iimbo) % quate
        q_dot =>  rbbou(iimbo) % q_dot
        vpfor =>  rbbou(iimbo) % vpfor
        vptor =>  rbbou(iimbo) % vptor
        !
        quate(1,1) = quate(1,2)
        q_dot(1,1) = q_dot(1,2)
        do idime = 1,3
           accel(idime,1)   = accel(idime,2) 
           velol(idime,1)   = velol(idime,2) 
           posil(idime,1)   = posil(idime,2)
           accea(idime,1)   = accea(idime,2) 
           veloa(idime,1)   = veloa(idime,2) 
           posia(idime,1)   = posia(idime,2)                 
           quate(idime+1,1) = quate(idime+1,2)
           q_dot(idime+1,1) = q_dot(idime+1,2)
           vpfor(idime,1)   = vpfor(idime,2) 
           vptor(idime,1)   = vptor(idime,2) 
        end do

     end if ! kfl_genco_ale

  case(23_ip)
     !
     ! Assign a(n,i,j-1) <-- a(n,i,j), update of RB unknown
     ! this is not used - the values of posil(,1)  are set directly in ale_solrbo
     !
     if ( kfl_rigid_ale == 0 ) return
     ! print*, "DEBUG: case(23)"

  case(24_ip)
     !
     ! keep the last dispm for the next coupling iteration
     !
     if( INOTMASTER ) then  ! added this because in the master dispm(1,1,2)       
        do ipoin = 1,npoin

           do idime = 1,ndime

              dispm(idime,ipoin,2) = dispm(idime,ipoin,1)
              !
              ! Added to use the indexes in coor_ale as the equivalent variables in other modules.
              ! This will be used in the strong coupling schemes
              !
              coord_ale(idime,ipoin,2) = coord_ale(idime,ipoin,1)

           end do

        end do
        ! print*, "DEBUG: case(24) endite"
     end if

     !
     ! Assign a(n,i-1,*) <-- a(n,i,*)
     !
     if ( kfl_rigid_ale == 0 ) return


     if ( kfl_genco_ale == 1_ip ) then
        !
        ! For the moment, only angular coordinates are being used as generalized coordinates
        !
        posia =>  rbbou(iimbo) % posia
        veloa =>  rbbou(iimbo) % veloa
        force =>  rbbou(iimbo) % force
        
        do idime = 1_ip, 3_ip

           posia(idime,2_ip)  = posia(idime,1_ip)
           veloa(idime,2_ip)  = veloa(idime,1_ip)
           force(idime,2_ip)  = force(idime,1_ip)
           
        end do

     else

        accel =>  rbbou(iimbo) % accel
        velol =>  rbbou(iimbo) % velol
        posil =>  rbbou(iimbo) % posil
        accea =>  rbbou(iimbo) % accea
        veloa =>  rbbou(iimbo) % veloa
        posia =>  rbbou(iimbo) % posia
        quate =>  rbbou(iimbo) % quate
        q_dot =>  rbbou(iimbo) % q_dot
        force =>  rbbou(iimbo) % force
        torqu =>  rbbou(iimbo) % torqu
        !
        quate(1,2) = quate(1,1)
        q_dot(1,2) = q_dot(1,1)
        do idime = 1,3
           accel(idime,2)   = accel(idime,1) 
           velol(idime,2)   = velol(idime,1) 
           posil(idime,2)   = posil(idime,1)
           accea(idime,2)   = accea(idime,1) 
           veloa(idime,2)   = veloa(idime,1) 
           posia(idime,2)   = posia(idime,1)                 
           quate(idime+1,2) = quate(idime+1,1)
           q_dot(idime+1,2) = q_dot(idime+1,1)
           force(idime,2)   = force(idime,1) 
           torqu(idime,2)   = torqu(idime,1) 
        end do

     end if ! kfl_genco_ale

  case(25_ip)
     !
     ! a(n-1,*,*) <-- a(n,*,*)
     !        

     !veloc(idime,ipoin,4) = veloc(idime,ipoin,3)
     !veloc(idime,ipoin,3) = veloc(idime,ipoin,1)
     !
     ! Assign a(n,i-1,*) <-- a(n,i,*)
     !
     if ( kfl_rigid_ale == 0 ) return

     if ( kfl_genco_ale == 1_ip ) then
        !
        ! For the moment, only angular coordinates are being used as generalized coordinates
        !
        posia =>  rbbou(iimbo) % posia
        veloa =>  rbbou(iimbo) % veloa
        force =>  rbbou(iimbo) % force
        
        do idime = 1_ip, 3_ip

           posia(idime,4_ip)  = posia(idime,3_ip)
           posia(idime,3_ip)  = posia(idime,1_ip)

           veloa(idime,4_ip)  = veloa(idime,3_ip)
           veloa(idime,3_ip)  = veloa(idime,1_ip)
           
           force(idime,4_ip)  = force(idime,3_ip)
           force(idime,3_ip)  = force(idime,1_ip)
           
        end do

     else

        accel =>  rbbou(iimbo) % accel
        velol =>  rbbou(iimbo) % velol
        posil =>  rbbou(iimbo) % posil
        accea =>  rbbou(iimbo) % accea
        veloa =>  rbbou(iimbo) % veloa
        posia =>  rbbou(iimbo) % posia
        quate =>  rbbou(iimbo) % quate
        q_dot =>  rbbou(iimbo) % q_dot
        force =>  rbbou(iimbo) % force
        torqu =>  rbbou(iimbo) % torqu
        !
        quate(1,4) = quate(1,3)
        q_dot(1,4) = q_dot(1,3)
        do idime = 1,3
           accel(idime,4)   = accel(idime,3) 
           velol(idime,4)   = velol(idime,3) 
           posil(idime,4)   = posil(idime,3)
           accea(idime,4)   = accea(idime,3) 
           veloa(idime,4)   = veloa(idime,3) 
           posia(idime,4)   = posia(idime,3)                 
           quate(idime+1,4) = quate(idime+1,3)
           q_dot(idime+1,4) = q_dot(idime+1,3)
           force(idime,4)   = force(idime,3) 
           torqu(idime,4)   = torqu(idime,3) 
        end do
        !
        quate(1,3) = quate(1,1)
        q_dot(1,3) = q_dot(1,1)
        do idime = 1,3
           accel(idime,3)   = accel(idime,1) 
           velol(idime,3)   = velol(idime,1) 
           posil(idime,3)   = posil(idime,1)
           accea(idime,3)   = accea(idime,1) 
           veloa(idime,3)   = veloa(idime,1) 
           posia(idime,3)   = posia(idime,1)                 
           quate(idime+1,3) = quate(idime+1,1)
           q_dot(idime+1,3) = q_dot(idime+1,1)
           if ( kfl_crist_ale /=1 ) force(idime,3)   = force(idime,1)   ! for cristobal's case this is not needed
           if ( kfl_crist_ale /=1 ) torqu(idime,3)   = torqu(idime,1) 
        end do
        ! print*, "DEBUG: case(25) ",posil(:,1) 

     end if  ! kfl_genco_ale

  case(31_ip)
     !
     ! Assign a(n,i,*)  <-- u(n-1,*,*), initial guess after reading restart
     !
     if ( kfl_rigid_ale == 0 ) return

     if ( kfl_genco_ale == 1_ip ) then
        !
        ! For the moment, only angular coordinates are being used as generalized coordinates
        !
        posia =>  rbbou(iimbo) % posia
        veloa =>  rbbou(iimbo) % veloa
        vpfor =>  rbbou(iimbo) % vpfor
       
        do idime = 1_ip, 3_ip

           posia(idime,1_ip)  = posia(idime,nprev_ale)
           posia(idime,2_ip)  = posia(idime,nprev_ale)
           
           veloa(idime,1_ip)  = veloa(idime,nprev_ale)
           veloa(idime,2_ip)  = veloa(idime,nprev_ale)
           
           vpfor(idime,1_ip)  = vpfor(idime,nprev_ale)
           vpfor(idime,2_ip)  = vpfor(idime,nprev_ale)
           
           !vpfor(idime,2_ip)   = 1.0_rp * vpfor(idime,nprev_ale)
           !vpfor(idime,1_ip)   = vpfor(idime,2)
                      
        end do

     else

        accel =>  rbbou(iimbo) % accel
        velol =>  rbbou(iimbo) % velol
        posil =>  rbbou(iimbo) % posil
        accea =>  rbbou(iimbo) % accea
        veloa =>  rbbou(iimbo) % veloa
        posia =>  rbbou(iimbo) % posia
        quate =>  rbbou(iimbo) % quate
        q_dot =>  rbbou(iimbo) % q_dot
        vpfor =>  rbbou(iimbo) % vpfor
        vptor =>  rbbou(iimbo) % vptor
        !
        quate(1,1) = quate(1,nprev_ale)
        quate(1,2) = quate(1,nprev_ale)
        q_dot(1,1) = q_dot(1,nprev_ale)
        q_dot(1,2) = q_dot(1,nprev_ale)

        ! print*, "DEBUG: case(31) posil", posil(:,1_ip)

        do idime = 1,3
           accel(idime,1)   = accel(idime,nprev_ale) 
           velol(idime,1)   = velol(idime,nprev_ale) 
           posil(idime,1)   = posil(idime,nprev_ale)
           accea(idime,1)   = accea(idime,nprev_ale) 
           veloa(idime,1)   = veloa(idime,nprev_ale) 
           posia(idime,1)   = posia(idime,nprev_ale)                 
           quate(idime+1,1) = quate(idime+1,nprev_ale)
           q_dot(idime+1,1) = q_dot(idime+1,nprev_ale)

           accel(idime,2)   = accel(idime,nprev_ale) 
           velol(idime,2)   = velol(idime,nprev_ale) 
           posil(idime,2)   = posil(idime,nprev_ale)
           accea(idime,2)   = accea(idime,nprev_ale) 
           veloa(idime,2)   = veloa(idime,nprev_ale) 
           posia(idime,2)   = posia(idime,nprev_ale)                 
           quate(idime+1,2) = quate(idime+1,nprev_ale)
           q_dot(idime+1,2) = q_dot(idime+1,nprev_ale)
        end do

        if ( ( kfl_foexo_ale == 1_ip ) .or. ( kfl_crist_ale == 1_ip ) ) then ! Force ( & Torque) extrapolation order
           do idime = 1,3
              vpfor(idime,2)   = 1.0_rp * vpfor(idime,nprev_ale)
              vptor(idime,2)   = 1.0_rp * vptor(idime,nprev_ale)
              vpfor(idime,1)   = vpfor(idime,2)
              vptor(idime,1)   = vptor(idime,2)
           end do
        else
           do idime = 1,3
              vpfor(idime,2)   = 1.5_rp * vpfor(idime,nprev_ale) - 0.5_rp * vpfor(idime,nprev_ale+1_ip)
              vptor(idime,2)   = 1.5_rp * vptor(idime,nprev_ale) - 0.5_rp * vptor(idime,nprev_ale+1_ip)
              vpfor(idime,1)   = vpfor(idime,2)
              vptor(idime,1)   = vptor(idime,2)
           end do
        end if

     end if  ! kfl_genco_ale
     ! print*, "DEBUG: case(31)", posil(:,1_ip)
  end select

end subroutine ale_updunk

