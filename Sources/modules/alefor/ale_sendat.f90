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

  if( itask == 1 ) then

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

end subroutine ale_sendat
