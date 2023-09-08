subroutine chm_outset()
  !-----------------------------------------------------------------------
  !****f* partis/chm_outset
  ! NAME 
  !    chm_outset
  ! DESCRIPTION
  !    Compute and write results on sets
  ! USES
  ! USED BY
  !    chm_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ieset,ibset,inset,iclas,dummi

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if(maxval(postp(1) % npp_setse)>0) then

     if( INOTMASTER ) then
        do ieset = 1,neset
           call chm_elmset(lesec(ieset),ieset)
        end do
     end if
     call posdef(21_ip,dummi)
     !
     ! Normalize results of necessary
     !
     if( INOTSLAVE ) then
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then

     do ibset = 1,nbset 
        call chm_bouset(lbsec(ibset),ibset)
     end do
     call posdef(22_ip,dummi)
     !
     ! Normalize results
     !
     if( INOTSLAVE ) then
        do ibset = 1,nbset
           do iclas=2,nclas_chm 
              !
              ! <Yk> = int_S rho*u* Y_k dS / int_S rho*u dS, for k = 1,...,8
              !
              if(postp(1) % npp_setsb(2)/=0.and.vbset(1,ibset) /= 0.0_rp) then 
                 vbset( iclas,ibset) = vbset( iclas,ibset) / vbset( 1,ibset)
              end if

           end do
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
              iclas=1
              if(postp(1) % npp_setsn(1)/=0) postp(1) % vnset(1,inset)=conce(lnsec(inset),iclas,1)
           end if
        end do
     end if
     call posdef(23_ip,dummi)

  end if

end subroutine chm_outset
