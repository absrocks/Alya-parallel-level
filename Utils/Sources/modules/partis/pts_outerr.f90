subroutine pts_outerr()
  !------------------------------------------------------------------------
  !****f* master/pts_outerr
  ! NAME
  !    pts_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    Turnon
  !***
  !------------------------------------------------------------------------

  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0,itype,injec
  character(20) :: messa

  ierro = 0
  iwarn = 0
  !
  ! Particle injection
  !
  do itype = 1,mtyla
     if( parttyp(itype) % kfl_exist /= 0 ) then
        if(  parttyp(itype) % kfl_grafo == 1 .or. &
             parttyp(itype) % kfl_buofo == 1 ) then
           if(  grnor == 0.0_rp .and. gravi(1) == 0.0_rp .and. &
                gravi(2) == 0.0_rp .and. gravi(3) == 0.0_rp ) then
              ierro = ierro + 1
              call outfor(1_ip,lun_outpu,&
                   'GRAVITY IS MISSING FOR PARTICLE INJECTION')
           end if
        end if
        if(  parttyp(itype) % kfl_modla /= 1 .and. parttyp(itype) % diame <= 0.0_rp ) then
           ierro = ierro + 1
           call outfor(1_ip,lun_outpu,&
                'ZERO OR NEGATIVE DIAMETER PARTICLE')
        end if
     end if
     if( parttyp(itype) % kfl_brown /= 0 ) then
        if( kfl_prope == 0 .and. parttyp(itype) % diffu == 0.0_rp ) then
           ierro = ierro + 1
           call outfor(1_ip,lun_outpu,&
                'FOR BROWNIAN MOTION, VISCOSITY SHOULD BE DEFINED IN KERMOD')
        end if
     end if
  end do
  if( kfl_prope == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,lun_outpu,&
          'PROPERTIES NOT DEFINED')
  end if
  !
  ! Element bin
  !
  if( kfl_usbin_pts /= 0 .and. kfl_element_bin == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,lun_outpu,&
          'ACTIVATE BIN_ELEMENT OPTION IN KER.DAT FILE')
  end if
  !
  ! If there are sliup walls and normals are not allocated
  !
  if(IMASTER) then
     if( kfl_slip_wall_pts > 0 .or. kfl_bouncing_wall_pts > 0 ) then
        if( kfl_walln/=1 ) then
           ierro = ierro + 1
           call outfor(1_ip,lun_outpu,&
                'For slip walls in partis, add WALL_NORMALS to ker.dat to NUMERICAL_TREATMENT section')
        end if
     end if
  end if  
  !
  !
  !
  do injec = 1,size(kfl_injla_pts)
     if( kfl_injla_pts(injec) /= 0 ) then
        if( kfl_injty_pts(injec) /= 0 ) then
           itype =  kfl_injty_pts(injec)
           if( parttyp(itype) % kfl_exist == 0 ) then
              ierro = ierro + 1
              call outfor(1_ip,lun_outpu,'NON-DEFINED PARTICLE TYPE IS INJECTED BY INJECTOR '//trim(intost(injec)))
           end if
        end if
     end if
  end do
  !
  ! Check if there is a wall distance criterion but distance to the wall does not exist
  !
  !if( dimin_pts > 0.0_rp .and. kfl_walld == 0 ) then
  !   ierro = ierro + 1
  !   call outfor(1_ip,lun_outpu,&
  !        'DISTANCE TO THE WALL SHOULD BE COMPUTED BY KERMOD TO USE WALL-DISTANCE BASED DEPOSITION CRITERION')
  !end if
  !
  ! Stop
  !
  call errors(3_ip,ierro,iwarn,'')

  !messa = intost(ierro)
  !if( ierro == 1 ) then
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERROR HAS BEEN FOUND')
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERROR HAS BEEN FOUND')
  !else if( ierro >= 2 ) then
  !   call runend('PTS_OUTERR: '//trim(messa)//' ERRORS HAVE BEEN FOUND')
  !end if

end subroutine pts_outerr
