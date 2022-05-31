!------------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    vtkmov.f90
!> @date    15/01/2015
!> @author  Hadrien Calmet
!> @brief   moving vtk data from scratch to gpfs
!> @details CASE 1
!!               naked coping data
!!          CASE 2
!!               tar and copy data
!!          CASE 3
!!               targz and copy data
!> @}
!------------------------------------------------------------------------
subroutine vtkmov(itask)
  !
  use def_parame
  use def_master
  use mod_std
  !--------------------------------------
  implicit none
  integer(ip), intent(in) :: itask       !> option for moving data
  integer(ip)             :: ipart
  character(150)          :: TMPDIR,my_tmpdir,nunam_pos1,nunam_pos2,nunam_pos3
  integer(ip)             :: values(8),auxi1,auxi2,auxi3,auxi_tar,auxi_cptar,auxi
  integer(ip)             :: auxi4,auxi5,auxi_targz,auxi_cptargz
  !--------------------------------------
  if (kfl_outfo == 40 .OR. kfl_outfo == 41 ) then
     !
     if (kfl_vtk == 0 ) then
        !
        if (ISLAVE) then
           !
           select case ( itask )
              !
              ! 3 options : 1 - naked coping data
              !             2 - tar and copy data
              !             3 - targz and copy data
           case ( 1_ip )
              !
              ! naked coping data (example :cp /scratch/ *_1.vtu /gpfs/file_name_vtk )
              !
              call date_and_time(VALUES=values)
              !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
              auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              !
              call get_environment_variable("TMPDIR", my_tmpdir)
              ipart = kfl_paral - 1
              if      (ipart<10)then
                 write(nunam_pos2,'(a,i1,a)') '*_',ipart,'.vtu'! composition of part of data name (example:*_1.vtu)
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
              else if (ipart<100)then
                 write(nunam_pos2,'(a,i2,a)') '*_',ipart,'.vtu'! composition of part of data name (example:*_1.vtu)
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
              else if (ipart<1000)then
                 write(nunam_pos2,'(a,i3,a)') '*_',ipart,'.vtu'! composition of part of data name (example:*_1.vtu)
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
              else if (ipart<10000)then
                 write(nunam_pos2,'(a,i4,a)') '*_',ipart,'.vtu'! composition of part of data name (example:*_1.vtu)
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
              else if (ipart<100000)then
                 write(nunam_pos2,'(a,i5,a)') '*_',ipart,'.vtu'! composition of part of data name (example:*_1.vtu)
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
              endif
              !
              ! naked coping data duration
              !
              call date_and_time(VALUES=values)
              !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
              auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              auxi=auxi2-auxi1
              !write (*,*),'kfl_paral=',kfl_paral,'copy naked data',auxi
              !
           case ( 2_ip )
              !
              ! tar and copy data (example :tar -cf /scratch/fensap_1.tar  /scratch/ *_1.vtu
              !                             cp /scratch/fensap_1.tar /gpfs/file_name_vtk/ )
              !
              call get_environment_variable("TMPDIR", my_tmpdir)
              !
              ipart = kfl_paral - 1
              if      (ipart<10)then
                 write(nunam_pos1,'(a,a,i1,a)') trim(title),'_',ipart,'.tar'    ! composition of tar name (example:fensap_1.tar)
                 write(nunam_pos3,'(a,i1,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! tar step
                 call execute_command_line ('tar -cf '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  tar step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi4=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi4=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<100)then
                 write(nunam_pos1,'(a,a,i2,a)') trim(title),'_',ipart,'.tar'    ! composition of tar name (example:fensap_1.tar)
                 write(nunam_pos3,'(a,i2,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! tar step
                 call execute_command_line ('tar -cf '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  tar step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi4=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi4=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<1000)then
                 write(nunam_pos1,'(a,a,i3,a)') trim(title),'_',ipart,'.tar'    ! composition of tar name (example:fensap_1.tar)
                 write(nunam_pos3,'(a,i3,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! tar step
                 call execute_command_line ('tar -cf '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  tar step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi4=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi4=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<10000)then
                 write(nunam_pos1,'(a,a,i4,a)') trim(title),'_',ipart,'.tar'    ! composition of tar name (example:fensap_1.tar)
                 write(nunam_pos3,'(a,i4,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! tar step
                 call execute_command_line ('tar -cf '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  tar step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi4=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi4=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<100000)then
                 write(nunam_pos1,'(a,a,i5,a)') trim(title),'_',ipart,'.tar'    ! composition of tar name (example:fensap_1.tar)
                 write(nunam_pos3,'(a,i5,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! tar step
                 call execute_command_line ('tar -cf '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  tar step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos1)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi4=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi4=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              endif
              !
              ! tar and copy data duration
              !
              auxi_tar=auxi2-auxi1     ! tar time
              auxi_cptar=auxi4-auxi2   ! copy tar time
              !write (*,*),'kfl_paral=',kfl_paral,'tar copy data',auxi_cptar+auxi_tar
              !
           case ( 3_ip )
              !
              ! targz and copy data (example :tar cfz /scratch/fensap_1.targz  /scratch/ *_1.vtu
              !                               cp /scratch/fensap_1.targz /gpfs/file_name_vtk/ )
              !
              call get_environment_variable("TMPDIR", my_tmpdir)
              !
              ipart = kfl_paral - 1
              if      (ipart<10)then
                 write(nunam_pos2,'(a,a,i1,a)') trim(title),'_',ipart,'.tar.gz' ! composition of tar name (example:fensap_1.tar.gz)
                 write(nunam_pos3,'(a,i1,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! targz step
                 call execute_command_line ('tar cfz '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi3=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi3=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  targz step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi5=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi5=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<100)then
                 write(nunam_pos2,'(a,a,i2,a)') trim(title),'_',ipart,'.tar.gz' ! composition of tar name (example:fensap_1.tar.gz)
                 write(nunam_pos3,'(a,i2,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! targz step
                 call execute_command_line ('tar cfz '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi3=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi3=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  targz step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi5=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi5=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<1000)then
                 write(nunam_pos2,'(a,a,i3,a)') trim(title),'_',ipart,'.tar.gz' ! composition of tar name (example:fensap_1.tar.gz)
                 write(nunam_pos3,'(a,i3,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! targz step
                 call execute_command_line ('tar cfz '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi3=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi3=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  targz step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi5=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi5=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<10000)then
                 write(nunam_pos2,'(a,a,i4,a)') trim(title),'_',ipart,'.tar.gz' ! composition of tar name (example:fensap_1.tar.gz)
                 write(nunam_pos3,'(a,i4,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! targz step
                 call execute_command_line ('tar cfz '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi3=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi3=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  targz step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi5=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi5=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              else if (ipart<100000)then
                 write(nunam_pos2,'(a,a,i5,a)') trim(title),'_',ipart,'.tar.gz' ! composition of tar name (example:fensap_1.tar.gz)
                 write(nunam_pos3,'(a,i5,a)') '*_',ipart,'.vtu'                 ! composition of part of data name (example:*_1.vtu)
                 !initiale time
                 call date_and_time(VALUES=values)
                 !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! targz step
                 call execute_command_line ('tar cfz '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(my_tmpdir)//'/'//trim(title)//'-'//trim(nunam_pos3))
                 call date_and_time(VALUES=values)
                 !auxi3=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi3=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
                 ! cp  targz step
                 call execute_command_line ('cp '//trim(my_tmpdir)//'/'//trim(nunam_pos2)//' '//trim(title)//'_vtk/')
                 call date_and_time(VALUES=values)
                 !auxi5=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
                 auxi5=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
              endif
              !
              ! targz and copy data duration
              !
              auxi_targz=auxi3-auxi1   ! targz time
              auxi_cptargz=auxi5-auxi3 ! copy targz time
              !write (*,*),'kfl_paral=',kfl_paral,'targz copy data',auxi_cptargz+auxi_targz
              !
           end select

        endif
     endif
  endif

end subroutine vtkmov !vtkmov
