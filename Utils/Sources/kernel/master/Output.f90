subroutine Output(itask)
  !-----------------------------------------------------------------------
  !****f* master/Output
  ! NAME
  !    Output
  ! DESCRIPTION
  !    This routine output and postprocess the solution
  ! USES
  !    moduls
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master

#ifdef CATA
  use def_domain, only : coord, npoin, lnods, nelem, ltype
  use tcp
#endif

  implicit none

  integer(ip), intent(in) :: itask
  integer(ip)             :: dummi,ivari,ivarp
  !
  ! Where we are
  ! 1 = initial solution
  ! 2 = end of a time step
  ! 3 = end of the run
  !
!!$  if( itask == 1 ) then
!!$     ittyp = ITASK_INITIA
!!$  else if( itask == 2 ) then
!!$     ittyp = ITASK_ENDTIM
!!$  else if( itask == 3 ) then
!!$     ittyp = ITASK_ENDRUN
!!$  else if( itask == 4 ) then
!!$     ittyp = ITASK_ENDINN
!!$  end if
  ittyp = itask
  !
  ! Set header
  !
  if( ittyp ==  ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     do modul = 1,mmodu
        if( kfl_modul(modul) /= 0 ) then
           postp => momod(modul) % postp
           call posdef(3_ip,dummi)
        end if
     end do
     modul = 0
  end if
  !
  ! Output and postprocess of modules
  !
  call Kermod(-ITASK_OUTPUT)
  do iblok = 1,nblok
     call moduls(ITASK_OUTPUT)
  end do
  !
  ! Set writing
  !
  if( ittyp ==  ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     do modul = 1,mmodu
        if( kfl_modul(modul) /= 0 ) then
           postp => momod(modul) % postp
           call posdef(4_ip,dummi)
        end if
     end do
     modul = 0
  end if
  !
  ! Close H5 file for current timestep
  !
  call Hdfpos( 7_ip )
  !
  !VTK
  !
  !FLUSH AND CLOSE for current timestep
  !
  !#ifdef VTK
  !  if (ittim > 0.and.ISLAVE.and.vtk_id/=0) then
  !     call vtkXMLWriterF_Write(vtk_id, ierr)
  !     call vtkXMLWriterF_Delete(vtk_id)
  !  endif
  !#endif
  !
  !CATALYST
  !
#ifdef CATA
  !write(kfl_paral+700,*)'npoin=',npoin
  !write(kfl_paral+700,*)'coord=',coord(1:3,npoin)
  !if (ittim==10)write(kfl_paral+700,*)'veloc=',veloc(1:3,1:1,1)
  !if (ittim==10)write(kfl_paral+700,*)'press=',press(1:1,1)
  !write(kfl_paral+700,*)'------------------------------'
  !write(kfl_paral+700,*)'veloc2=',veloc(1:3,1,1)
  call testcoprocessor(ittim,cutim,npoin,coord,nelem,lnods,ltype,kfl_paral,veloc,press)
  !call testcoprocessor(ittim,cutim,npoin,coord,nelem,lnods,ltype,kfl_paral)
#endif

end subroutine Output
