subroutine exm_outset()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_outset
  ! NAME 
  !    exm_outset
  ! DESCRIPTION
  !    Compute and write results on sets:
  !    - Element, boundary and node sets:
  !      1. INTRA
  !      2. EXTRA
  !      3. RECOV
  ! USES
  ! USED BY
  !    exm_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_iofile
  implicit none
  integer(ip) :: ipoin,dummi
  integer(ip) :: inset
  integer(ip) :: ieset

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then

     if( INOTMASTER ) then

        do inset = 1,nnset
          call exm_nodset(lnsec(inset),inset)
        end do
     end if

     call posdef(23_ip,dummi)

  end if


   !----------------------------------------------------------------------
   !
   ! Element sets
   !
   !----------------------------------------------------------------------

   if( maxval(postp(1) % npp_setse) > 0 ) then

      if( INOTMASTER ) then
         do ieset = 1,neset
            call exm_elmset(lesec(ieset),ieset)
         end do
      end if
      !
      ! Parall
      !
      call posdef(21_ip,dummi)

   end if



end subroutine exm_outset
