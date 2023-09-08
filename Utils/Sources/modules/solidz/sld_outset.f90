!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outset.f90
!> @author  Solidz Team
!> @date    April, 2010
!>          - Subroutine creation
!> @brief   Compute and write results on sets
!> @details Sets can be defined at element, boundary and node levels.
!> @note    The keyflag kfoten is only activated for specific variables
!>          (node set) which need it or always (element set).
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_outset()

  use def_kintyp, only : ip, rp
  use def_master, only : INOTMASTER, INOTSLAVE
  use def_master, only : postp, vbset
  use def_domain, only : nnset, neset, nbset
  use def_domain, only : lnsec, lbsec, lesec
  use def_solidz, only : kfote_sld

  implicit none

  integer(ip)         :: ieset, ibset, inset, dummi

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse) > 0 ) then

     if( INOTMASTER ) then
        !
        ! Call fotens if necessary
        if( any(postp(1) % npp_setse(2:14) /= 0_ip) ) then
           kfote_sld = kfote_sld + 1_ip
           if (kfote_sld == 1_ip) call sld_fotens
        end if

        do ieset = 1,neset
           call sld_elmset(lesec(ieset),ieset)
        end do
     end if
     !
     ! Parall
     !
     call posdef(21_ip,dummi)

  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if ( maxval(postp(1) % npp_setsb) > 0 ) then

     if ( INOTMASTER ) then
        do ibset = 1,nbset
           call sld_bouset(lbsec(ibset),ibset)
        end do
     end if
     !
     ! Parall
     !
     call posdef(22_ip,dummi)
     !
     ! Averaged/Norm of some variables
     !
     if ( INOTSLAVE ) then
        do ibset = 1,nbset
           ! DIBOX
           if ( postp(1) % npp_setsb(5) /= 0 .and. vbset(postp(1) % nvabs+1,ibset)>0 ) then
              vbset(5,ibset) = vbset(5,ibset)/vbset(postp(1) % nvabs+1,ibset)
           end if
           ! DIBOY
           if ( postp(1) % npp_setsb(6) /= 0 .and. vbset(postp(1) % nvabs+1,ibset)>0 ) then
              vbset(6,ibset) = vbset(6,ibset)/vbset(postp(1) % nvabs+1,ibset)
           end if
           ! DIBOZ
           if ( postp(1) % npp_setsb(7) /= 0 .and. vbset(postp(1) % nvabs+1,ibset)>0 ) then
              vbset(7,ibset) = vbset(7,ibset)/vbset(postp(1) % nvabs+1,ibset)
           end if
           ! DIBOU
           if ( postp(1) % npp_setsb(8) /= 0 .and. vbset(postp(1) % nvabs+1,ibset)>0 ) then
              vbset(8,ibset) = vbset(8,ibset)/vbset(postp(1) % nvabs+1,ibset)
           end if
           ! FCONT
           if ( postp(1) % npp_setsb(16) /= 0 ) then
              vbset(16,ibset) = sqrt(sum(vbset(13:15,ibset)**2))
           end if
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then

     if( INOTMASTER ) then
        !
        ! Call fotens if any of these variables are called:
        ! SIGMA, EPSIL, LNEPS, SEQVM
        if( any(postp(1) % npp_setsn(13:30) /= 0_ip) .or. postp(1) % npp_setsn(34) /= 0_ip ) then

           kfote_sld = kfote_sld + 1_ip
           if ( kfote_sld == 1_ip ) call sld_fotens

        end if

        do inset = 1,nnset
           call sld_nodset(lnsec(inset),inset)
        end do

     end if
     !
     ! Parall
     !
     call posdef(23_ip,dummi)

  end if

end subroutine sld_outset
