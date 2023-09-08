subroutine rad_outset()

  !-----------------------------------------------------------------------
  !****f* Radiat/rad_outset
  ! NAME 
  !    rad_outset
  ! DESCRIPTION
  !    Compute and write results on sets
  ! USES
  ! USED BY
  !    rad_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip) :: ieset,ibset,inset,dummi

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if(maxval(postp(1) % npp_setse)>0) then

     if( INOTMASTER ) then
        do ieset=1,neset        
           call rad_elmset(lesec(ieset),postp(1) % veset(postp(1) % nvaes+1,ieset),&
                postp(1) % veset(1,ieset))
        end do
     end if
     call posdef(21_ip,dummi)
     if(postp(1) % npp_setse(1)/=0) then
        do ieset=1,neset
           if(postp(1) % veset(postp(1) % nvaes+1,ieset)/=0.0_rp)&
                postp(1) % veset(1,ieset)=postp(1) % veset(1,ieset)&
                /postp(1) % veset(postp(1) % nvaes+1,ieset)
        end do
     end if
 
  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if(maxval(postp(1) % npp_setsb)>0) then

     if( INOTMASTER ) then
        do ibset=1,nbset
           call rad_bouset(lbsec(ibset),postp(1) % vbset(postp(1) % nvabs+1,ibset),&
                postp(1) % vbset(1,ibset),postp(1) % vbset(2,ibset))
        end do
     end if
     call posdef(22_ip,dummi)
     if(postp(1) % npp_setsb(1)/=0) then
        do ibset=1,nbset
           if(postp(1) % vbset(postp(1) % nvabs+1,ibset)/=0.0_rp)&
                postp(1) % vbset(1,ibset)=postp(1) % vbset(1,ibset)&
                /postp(1) % vbset(postp(1) % nvabs+1,ibset)
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if(maxval(postp(1) % npp_setsn)>0) then

     if( INOTMASTER ) then
        do inset=1,nnset
           if(lnsec(inset)/=0) then
              if(postp(1) % npp_setsn(1)/=0) postp(1) % vnset(1,inset)=radav_rad(lnsec(inset),1)
           end if
        end do
     end if
     call posdef(23_ip,dummi)

  end if

end subroutine rad_outset
