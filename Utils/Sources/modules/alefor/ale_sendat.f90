subroutine ale_sendat(itask)
  !-----------------------------------------------------------------------
  !****f* master/ale_sendat
  ! NAME
  !    ale_sendat
  ! DESCRIPTION
  !    This routines creates the IB structure
  ! for the moment only itask=0 & 1 are ready
  ! y borre algo de 7 
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ii,jj,inoib,ipoib,iboib,idime,jdime

  if( itask == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     !-------------------------------------------------------------------

     call ale_mealrb(2_ip)
     do iimbo = 1,nrbod         

        rbbou(iimbo) % massa      = -1.0_rp
        rbbou(iimbo) % densi      = -1.0_rp 
        rbbou(iimbo) % volum      =  0.0_rp
        rbbou(iimbo) % momin      = -1.0_rp
        rbbou(iimbo) % posgr      = -1.0e12_rp

        rbbou(iimbo) % npoib      =  0
        rbbou(iimbo) % nboib      =  0

        rbbou(iimbo) % posil      =  -1.0e12_rp
        rbbou(iimbo) % velol      =  0.0_rp
        rbbou(iimbo) % accel      =  0.0_rp
        rbbou(iimbo) % force      =  0.0_rp
        rbbou(iimbo) % vpfor      =  0.0_rp
        rbbou(iimbo) % vforce     =  0.0_rp
        rbbou(iimbo) % pforce     =  0.0_rp

        rbbou(iimbo) % torqu      =  0.0_rp
        rbbou(iimbo) % vptor      =  0.0_rp
        rbbou(iimbo) % vtorqu     =  0.0_rp
        rbbou(iimbo) % ptorqu     =  0.0_rp

        rbbou(iimbo) % posia      =  0.0_rp
        rbbou(iimbo) % veloa      =  0.0_rp
        rbbou(iimbo) % accea      =  0.0_rp
        rbbou(iimbo) % rotac      =  0.0_rp
        rbbou(iimbo) % quate(1,1) =  1.0_rp
        rbbou(iimbo) % quate(1,2) =  1.0_rp
        rbbou(iimbo) % quate(1,3) =  1.0_rp
        rbbou(iimbo) % quate(1,4) =  1.0_rp
        rbbou(iimbo) % q_dot(1,1) =  0.0_rp
        rbbou(iimbo) % q_dot(1,2) =  0.0_rp
        rbbou(iimbo) % q_dot(1,3) =  0.0_rp
        rbbou(iimbo) % q_dot(1,4) =  0.0_rp

        do idime = 1,3
           rbbou(iimbo) % quate(idime+1,1) = 0.0_rp
           rbbou(iimbo) % quate(idime+1,2) = 0.0_rp                      
           rbbou(iimbo) % quate(idime+1,3) = 0.0_rp                      
           rbbou(iimbo) % quate(idime+1,4) = 0.0_rp                      
           rbbou(iimbo) % q_dot(idime+1,1) = 0.0_rp
           rbbou(iimbo) % q_dot(idime+1,2) = 0.0_rp                      
           rbbou(iimbo) % q_dot(idime+1,3) = 0.0_rp                      
           rbbou(iimbo) % q_dot(idime+1,4) = 0.0_rp                      
           do jdime = 1,3
              rbbou(iimbo) % rotac(idime,jdime) = 0.0_rp
           end do
        end do

        do idime = 1,3
           rbbou(iimbo) % rotac(idime,idime) = 1.0_rp      
        end do

     end do

  else if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------
     !
     ! Read in ale_reaphy
     !
     call iexcha( kfl_rigid_ale )
     do ii = 1,3
        call iexcha( nstli_ale(ii) )
        call iexcha( nstro_ale(ii) )
     end do
     call iexcha( kfl_mvext_ale )
     call iexcha( kfl_ralei_ale )
     call iexcha( nstra_ale )
     call iexcha( nenra_ale )
     call iexcha( kfl_catfo_ale )             ! Catamaran force
     call iexcha( kfl_grafo_ale )             ! Gravity force
     call iexcha( kfl_forca_res )             ! Residual based forces for rigid body motion
     call iexcha( kfl_sprin_ale )             ! Spring activated for rigid body motion
     call iexcha( kfl_topmo_ale )             ! Top model activated for rigid body motion
     call iexcha( kfl_pertu_ale )             ! Use perturbation of generalized forces
     call iexcha( kfl_genco_ale )             ! Generalized coordinates solver for rigid body

     call rexcha( ralei_ale )
     !
     ! Spring constants for rigid body motion
     !
     do ii = 1_ip, 3_ip
        call rexcha( sprin_ale(ii) )
     end do

  else if( itask == 2 ) then   ! This was previously in itask 1 but now I treat it separatelly so that everybody has kfl_rigid_ale    
      
     if ( kfl_rigid_ale == 1 ) then

        do iimbo = 1,nrbod
           call iexcha( rbbou(iimbo) % npoib )
           call iexcha( rbbou(iimbo) % nboib )
           call iexcha( rbbou(iimbo) % nrbse )
           do ii = 1,10
              call iexcha( rbbou(iimbo) % lrbse(ii) )
           end do

           call rexcha( rbbou(iimbo) % massa )
           call rexcha( rbbou(iimbo) % densi )
           call rexcha( rbbou(iimbo) % volum )
           do ii = 1,6
              call rexcha( rbbou(iimbo) % momin(ii) )
           end do
           do ii = 1,3
              call rexcha( rbbou(iimbo) % posgr(ii) ) 
           end do
        end do

     end if

  else if ( itask == 3_ip ) then   ! I split in 2 part what before was in itask=2  here I put what is also used in restart


     if ( kfl_rigid_ale == 1 ) then

        do iimbo = 1,nrbod
           !
           ! Linear motion
           !
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % posil(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % velol(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % accel(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % force(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % vpfor(ii,jj) )
              end do
           end do
           do ii = 1,3
              call rexcha( rbbou(iimbo) % vforce(ii) )
           end do
           do ii = 1,3
              call rexcha( rbbou(iimbo) % pforce(ii) )
           end do
           !
           ! Angular motion
           !
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % posia(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % veloa(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % accea(ii,jj) )
              end do
           end do
           do jj = 1,3
              do ii = 1,3
                 call rexcha( rbbou(iimbo) % rotac(ii,jj) )
              end do
           end do
           do ii = 1,4
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % quate(jj,ii) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % torqu(ii,jj) )
              end do
           end do
           do ii = 1,3
              do jj = 1,4
                 call rexcha( rbbou(iimbo) % vptor(ii,jj) )
              end do
           end do
           do ii = 1,3
              call rexcha( rbbou(iimbo) % vtorqu(ii) )
           end do
           do ii = 1,3
              call rexcha( rbbou(iimbo) % ptorqu(ii) )
           end do
        end do

     end if


  else if ( itask == 7_ip ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

!     call Parall(20_ip)    ! so that everybody has kfl_rigid_ale idea mia pero iba mal
     if ( kfl_rigid_ale == 1 ) then

        do iimbo = 1,nrbod
           !
           ! Boundary mesh
           !
           do ipoib = 1,rbbou(iimbo) % npoib
              do idime = 1,ndime
                 call rexcha( rbbou(iimbo)%cooin(idime,ipoib) ) ! COOIN
              end do
           end do
           do ipoib = 1,rbbou(iimbo) % npoib
              do idime = 1,ndime
                 call rexcha( rbbou(iimbo)%cooib(idime,ipoib) ) ! COOIB
              end do
           end do
           !        do ipoib = 1,rbbou(iimbo) % npoib
           !           do idime = 1,ndime
           !              call rexcha( rbbou(iimbo)%cooi2(idime,ipoib) ) ! COOI2  ! not used in ale
           !           end do
           !        end do
           do iboib = 1,rbbou(iimbo) % nboib
              do inoib = 1,mnoib
                 call iexcha( rbbou(iimbo)%lnoib(inoib,iboib) ) ! LNOIB
              end do
           end do
           do iboib = 1,rbbou(iimbo) % nboib
              call iexcha( rbbou(iimbo)%ltyib(iboib) )          ! LTYIB
           end do


        end do

     end if

  end if

end subroutine ale_sendat
