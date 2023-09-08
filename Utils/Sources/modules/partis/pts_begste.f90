!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_begste.f90
!> @author  houzeaux
!> @date    2018-09-22
!> @brief   Start a time step
!> @details Inject particles if required
!> @} 
!-----------------------------------------------------------------------

subroutine pts_begste()
  
  use def_parame
  use def_master
  use def_domain
  use def_partis
  use mod_pts_injection, only : pts_injection_injectors
  implicit none
  integer(ip) :: idime,kk
  real(rp) :: xx



  !Update fields, read if necessary
  !There is another call in master/Begste.f90, but it happens after begste(), i.e. after this function.
  !however I need to have it called before to have correct values for kk, xx and preload the files in time
  call calc_kx_tran_fiel()
  !
  ! Read velocity form file
  !
  if (INOTMASTER) then
     if ( veloc_field_id>0 ) then
        if  (kfl_field(4,veloc_field_id) > 1) then !if there is more tan 1 timestep
            !call read_veloc()

            !print *,'Loading field:',veloc_field_id,', timestep:', k_tran_fiel(veloc_field_id)
            advec(:,:,3) = advec(:,:,1) ! <-- in theory kermod does this, apparently not

            kk = k_tran_fiel(veloc_field_id)
            xx = x_tran_fiel(veloc_field_id)

            !print *, cutim,': Partis: kk=', kk
            !print *, cutim,': Partis: xx=', xx
            do idime = 1,ndime
               advec(idime,:,1)  = xfiel(veloc_field_id) % a(idime,:,kk) * xx + &
                                xfiel(veloc_field_id) % a(idime,:,kk+1) * (1.0_rp-xx)
            end do
         
            !print *,cutim,': veloc='
            !do idime = 1,size(advec,2,ip)
            !    print *, veloc(:,idime,1)
            !end do
            !print *,'===='
            !print *,'Size of veloc', size(veloc,1),'x', size(veloc,2)
            !print *,'Size of xfiel', size(xfiel(veloc_field_id) % a,1),'x', size(xfiel(veloc_field_id) % a,2)
        end if
        !else -- see pts_inivar(2_ip)
     end if
  end if


  if( kfl_timco == 1 ) dtinv = max(dtinv,1.0_rp/dtime)
  !
  ! Inject particles
  !
  nlagr = nlacc_pts
  call pts_injection_injectors()

  ! injetion happens before timestep is updated, so save injection into the old file still
  !
  ! If pts.res splitting is requested, open the new file
  !
  if(kfl_ptsres_split>0) then
     if( mod(ittim, kfl_oufre_pts) == 0 ) then
        call pts_openfi(2_ip)
     end if
  end if

  
end subroutine pts_begste


!subroutine read_veloc()
  !unfinished, read the velocity on each iteration from the files
!  use def_master, only : ittim, advec, IMASTER
!  use def_domain, only : npoin, ndime
!  use def_partis, only : veloc_field_prefix, veloc_field_Ndigits, veloc_field_Nsteps
!  use mod_mpio_seq_io, only : SEQ_FILE_READ_HEADER
!  use def_mpio, only : mpio_header
!  use def_kintyp, only : ip
!  implicit none

  


!  type(mpio_header) :: header
!  CHARACTER(LEN=255) :: filename
!  character(LEN=10)  :: str_ndigits
!  integer(ip)        :: iter2load

!  if(IMASTER) then

      !do modulo here to prediodically load the steps
!      iter2load = mod(ittim-1, veloc_field_Nsteps)+1

!      write(str_ndigits , '(i10)') veloc_field_Ndigits


!      WRITE (filename, '(a, i'//trim(adjustl(str_ndigits))//'.'//trim(adjustl(str_ndigits))//',a)')  trim(adjustl(veloc_field_prefix)), iter2load,'.mpio.bin'
!      print *, "<"//trim(adjustl(filename))//">"

!      print *,'Reading ',trim(adjustl(filename))
!      call SEQ_FILE_READ_HEADER(trim(adjustl(filename)), header)

!      if (header % dimension /= ndime) then
!         call runend('Number of dimensions in the file does not coincide with the dimensionality of the problem')
!      end if

!      if (header % lines /= npoin) then
!         call runend('Number of points in the file does not coincide with the number of points in the domain')
!      end if


!     call SEQ_FILE_READ_ALL(advec, filename, ndime, npoin)
!  end if
