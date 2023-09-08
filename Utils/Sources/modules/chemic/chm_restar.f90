subroutine chm_restar(itask)
  !------------------------------------------------------------------------
  !****f* chemic/chm_restar
  ! NAME 
  !    chm_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_postpr
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores,ispec
  integer(ip)             :: dummi(10), nreal, ipoin
  character(5)            :: wopos(3), outna
  !
  ! Check if restrt file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = 3
  else
     icomp = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Concentration
  !
  !----------------------------------------------------------------------

  do ispec = 1,nspec_chm
     iwopo = 1
     wopos(1) = speci(ispec) % name
     wopos(2) = postp(1) % wopos(2,iwopo)
     wopos(3) = postp(1) % wopos(3,iwopo)
     if( INOTMASTER ) gesca => conce(:,ispec,icomp)
     call postpr(gesca,wopos,ittim,cutim)

     if (kfl_tisch_chm==2) then  ! BDF
        if (ncomp_chm > 3) then  ! BDF2 and bigger
           if( ispec < 10 ) then
             write(outna,'(a,i1,a)') 'C0', ispec, 'T2'
           else
             write(outna,'(a,i2,a)') 'C', ispec, 'T2'
           end if
           wopos(1) = outna
           wopos(2) = postp(1) % wopos(2,iwopo)
           wopos(3) = postp(1) % wopos(3,iwopo)
           if( INOTMASTER ) gesca => conce(:,ispec,4)
           call postpr(gesca,wopos,ittim,cutim)
 
           if (itask == READ_RESTART_FILE) then
              !
              ! Broadcast of file_opened
              !
              if( IPARALL ) then
                 if (IMASTER) then
                    if (file_opened) then
                       dummi(1) = 1
                    else
                       dummi(1) = 0
                    end if
                 end if
 
                 nreal = 1
                 call parari('BCT',0_ip,nreal,dummi)
                 if (dummi(1)==1) then
                    file_opened = .true.
                 else if (dummi(1) == 0) then
                    file_opened = .false.
                 end if
 
              end if
              kfl_rsta2_chm = file_opened ! last run was BDF2
           end if
   
           if (itask == READ_RESTART_FILE .and. .not. kfl_rsta2_chm ) then
              neule_chm = ittim + 1
              if (INOTMASTER) then
                 do ipoin =1, npoin
                     conce(ipoin,ispec,4) = conce(ipoin,ispec,icomp)
                 end do
              end if
           end if
        end if
     end if

  end do

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine chm_restar
 
