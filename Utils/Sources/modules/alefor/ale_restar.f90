!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_restar.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Restarting.
!> @details Restarting.
!>          ITASK = 1 ... Reads the initial values from the restart file
!>                  2 ... Writes restart file
!> @} 
!-----------------------------------------------------------------------
subroutine ale_restar(itask)
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_postpr
  use mod_memchk
  
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,kfl_gores
  integer(ip)             :: ipoin,idime
  
  !
  ! Check if restart file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0_ip ) return

  !----------------------------------------------------------------------
  !
  ! Displacement & Velocity
  !
  !----------------------------------------------------------------------
  gevec => dispm(:,:,1)
  call postpr(gevec,postp(1)%wopos(1:3,1),ittim,cutim)

  gevec => velom(:,:)
  call postpr(gevec,postp(1)%wopos(1:3,2),ittim,cutim)

  call postpr(coord,postp(1)%wopos(1:3,3),ittim,cutim)
  if( INOTMASTER )then
     do ipoin = 1,npoin
        do idime = 1,ndime
           coord_ale(idime,ipoin,1) = coord(idime,ipoin)
           coord_ale(idime,ipoin,2) = coord(idime,ipoin)
           ! -----------------------------------------------------------------------------
           ! Addition to keep the values of previous iterations and time step of coord_ale
           coord_ale(idime,ipoin,3) = coord(idime,ipoin)
           ! coord_ori(idime,ipoin)   = coord(idime,ipoin)
        end do
     end do
  end if

!!! Ojo sacar !!!
! coord(:,:) = coord(:,:) + dispm(:,:,1)

  kfl_domar = 1

  !----------------------------------------------------------------------
  !
  ! Rigid Body 
  !
  !----------------------------------------------------------------------

  if ( kfl_rigid_ale == 1 ) then
     if( INOTSLAVE ) then
        !
        ! Open files
        !
        call ale_openfi(0_ip)
        if( itask == READ_RESTART_FILE) then
           read(lun_resta_ale) rbbou(1) % posil(1:3,1:4) 
           read(lun_resta_ale) rbbou(1) % velol(1:3,1:4) 
           read(lun_resta_ale) rbbou(1) % accel(1:3,1:4) 
           read(lun_resta_ale) rbbou(1) % force(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % vpfor(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % pforce(1:3)  
           read(lun_resta_ale) rbbou(1) % vforce(1:3)  
           read(lun_resta_ale) rbbou(1) % posia(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % veloa(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % accea(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % rotac(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % torqu(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % vptor(1:3,1:4)
           read(lun_resta_ale) rbbou(1) % ptorqu(1:3)
           read(lun_resta_ale) rbbou(1) % vtorqu(1:3)
           read(lun_resta_ale) rbbou(1) % quate(1:4,1:4)
           read(lun_resta_ale) rbbou(1) % q_dot(1:4,1:4)

!           print*,'A:ale_restar:rbbou(1) % posil(:,:)', rbbou(1) % posil(:,:)
!           print*,'A:ale_restar:', rbbou(1) % velol(:,:) 
!           print*,'A:ale_restar:', rbbou(1) % accel(:,:) 
!           print*,'A:ale_restar:', rbbou(1) % force(:,:)
!           print*,'A:ale_restar:', rbbou(1) % vpfor(:,:)
!           print*,'A:ale_restar:', rbbou(1) % pforce(:)  
!           print*,'A:ale_restar:', rbbou(1) % vforce(:)  
!           print*,'A:ale_restar:', rbbou(1) % posia(:,:)
!           print*,'A:ale_restar:', rbbou(1) % veloa(:,:)
!           print*,'A:ale_restar:', rbbou(1) % accea(:,:)
!           print*,'A:ale_restar:', rbbou(1) % rotac(:,:)
!           print*,'A:ale_restar:', rbbou(1) % torqu(:,:)
!           print*,'A:ale_restar:', rbbou(1) % vptor(:,:)
!           print*,'A:ale_restar:', rbbou(1) % ptorqu(:)
!           print*,'A:ale_restar:', rbbou(1) % vtorqu(:)
!           print*,'A:ale_restar:', rbbou(1) % quate(:,:)
!           print*,'A:ale_restar:', rbbou(1) % q_dot(:,:)

        else

!           print*,'B:ale_restar:', rbbou(1) % posil(:,:)
!           print*,'B:ale_restar:', rbbou(1) % velol(:,:) 
!           print*,'B:ale_restar:', rbbou(1) % accel(:,:)

           write(lun_resta_ale) rbbou(1) % posil(1:3,1:4)  
           write(lun_resta_ale) rbbou(1) % velol(1:3,1:4) 
           write(lun_resta_ale) rbbou(1) % accel(1:3,1:4) 
           write(lun_resta_ale) rbbou(1) % force(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % vpfor(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % pforce(1:3)  
           write(lun_resta_ale) rbbou(1) % vforce(1:3)  
           write(lun_resta_ale) rbbou(1) % posia(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % veloa(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % accea(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % rotac(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % torqu(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % vptor(1:3,1:4)
           write(lun_resta_ale) rbbou(1) % ptorqu(1:3)
           write(lun_resta_ale) rbbou(1) % vtorqu(1:3)
           write(lun_resta_ale) rbbou(1) % quate(1:4,1:4) 
           write(lun_resta_ale) rbbou(1) % q_dot(1:4,1:4) 
  
!           print*,'B:ale_restar:', rbbou(1) % force(:,:)
!           print*,'B:ale_restar:', rbbou(1) % vpfor(:,:)
!           print*,'B:ale_restar:', rbbou(1) % pforce(:)  
!           print*,'B:ale_restar:', rbbou(1) % vforce(:)  
!           print*,'B:ale_restar:', rbbou(1) % posia(:,:)
!           print*,'B:ale_restar:', rbbou(1) % veloa(:,:)
!           print*,'B:ale_restar:', rbbou(1) % accea(:,:)
!           print*,'B:ale_restar:', rbbou(1) % rotac(:,:)
!           print*,'B:ale_restar:', rbbou(1) % torqu(:,:)
!           print*,'B:ale_restar:', rbbou(1) % vptor(:,:)
!           print*,'B:ale_restar:', rbbou(1) % ptorqu(:)
!           print*,'B:ale_restar:', rbbou(1) % vtorqu(:)
!           print*,'B:ale_restar:', rbbou(1) % quate(:,:)
!           print*,'B:ale_restar:', rbbou(1) % q_dot(:,:)

        end if
        !
        !close files
        !
        call ale_openfi(3_ip)
     end if
  end if
  !
  ! Finish
  !
  call respre(3_ip,kfl_gores)

end subroutine ale_restar
 
