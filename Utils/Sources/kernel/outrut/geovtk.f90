!------------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    geovtk.f90
!> @date    01/07/2013
!> @author  Hadrien Calmet
!> @brief   Alya output in vtk format
!> @details CASE sequencial:
!!            Alya in serial with VTK lib
!!          CASE master:
!!            Alya in parallel with master ascii creata metadata
!!          CASE slave:
!!            Alya in parallel with slave with VTK lib
!!
!> @}
!------------------------------------------------------------------------
subroutine geovtk(bridge,wopos,itste,ttime,pdime)

  use def_kintyp
  use def_domain
  use def_postpr
  use mod_iofile
  use def_master
  use mod_memory
  use mod_memchk
  use mod_std
  use mod_postpr, only : postpr_at_current_time_step


  implicit none
  character(*), intent(in)    :: wopos(*) !> array of variable name
  real(rp),     intent(inout) :: bridge(*)!> array of variable data
  integer(ip),  intent(in)    :: itste    !> integer time step
  real(rp),     intent(in)    :: ttime    !> real time step
  integer(ip),  intent(in)    :: pdime    !> dimension
  !-----------------------------------------------
  !
  !small remark the connectivity had to begin by 0
  !
  !-----------------------------------------------
  integer(ip)             :: mastervtk = 666
  integer(ip)             :: masterpvd = 777
  logical                 :: dir_e
  integer(ip)             :: ipart,i
  integer(ip)             :: inode,ielem
  integer(ip)             :: offset,istat,icount,pnode
  !-----------------------------------------------
  !ojo
  !lnods_tmp and offset_tmp had to be integer double precision 8 bytes
  !ltype_tmp had to be integer simple precison 2 bytes
  !
  !-----------------------------------------------
  integer(8), pointer       :: lnods_tmp(:),offset_tmp(:)
  integer(4), pointer       :: ltype_tmp(:)

  integer(ip)               :: ierr,nvar
  character(150)            :: filsa
  character(150)            :: nunam_pos1,nunam_pos2
  character(150)            :: TMPDIR,my_tmpdir

  integer(ip)               :: values(8),auxi1,auxi2,auxi
  integer(ip),save          :: total_auxi=0
  !---------------------------------------------------


  nullify(lnods_tmp)
  nullify(ltype_tmp)
  nullify(offset_tmp)

  offset=0
  icount=1
  !
  ! no data writing at T=0
  !
  if (itste == 0) return
  !
  ! setting of the output name
  ! Save the name
  !
  !
  write(nunam_pos1,'(i7)') ittim
  if( ittim < 10 ) then
     write(nunam_pos1,'(a,i1)') '000000',ittim
  else if( ittim < 100 ) then
     write(nunam_pos1,'(a,i2)') '00000',ittim
  else if( ittim < 1000 ) then
     write(nunam_pos1,'(a,i3)') '0000',ittim
  else if( ittim < 10000 ) then
     write(nunam_pos1,'(a,i4)') '000',ittim
  else if( ittim < 100000 ) then
     write(nunam_pos1,'(a,i5)') '00',ittim
  else if( ittim < 1000000 ) then
     write(nunam_pos1,'(a,i6)') '0',ittim
  end if
  filsa = trim(title)//'-'//trim(nunam_pos1)
  !
  ! nvar = number of variable to post process
  !
  nvar=postpr_at_current_time_step()

  if (ISEQUEN) then
     !
     ! Serial
     !
     ! test lib phillipe  prob with integer x4 x8 -----> alya est av integer x4(simple) et real x8 (double)
     !                                            -----> vtk est par defaut int x8 et real x8
     !check the length of lnods_tmp (check +1 is pnode+1)
     !
     ! Check if result folder exist
     !
     inquire(file=trim(title)//'_vtk', exist=dir_e)
     if ( dir_e ) then
        !write(*,*)"exist file"
     else
        call execute_command_line ('mkdir '//trim(title)//'_vtk 2> /dev/null')
     end if
     !
     !
     icount=0
     do ielem=1,nelem
        if       ( ltype(ielem) == 30 ) then ! TETRA
           icount=icount+4+1
        else if  ( ltype(ielem) == 37 ) then ! HEXA
           icount=icount+8+1
        else if  ( ltype(ielem) == 32 ) then ! PYRA
           icount=icount+5+1
        else if  ( ltype(ielem) == 34 ) then ! PENTA
           icount=icount+6+1
        endif
     enddo
     !
     ! allocate memory temporary
     !
     allocate (lnods_tmp(icount),stat=istat)
     !call memchk(zero,istat,memor_dom,'lnods_tmp','geovtk',lnods_tmp) !pb with alya test
     allocate (offset_tmp(nelem),stat=istat)
     !call memchk(zero,istat,memor_dom,'offset_tmp','geovtk',offset_tmp)!pb with alya test
     allocate (ltype_tmp(nelem),stat=istat)
     !call memchk(zero,istat,memor_dom,'ltype_tmp','geovtk',ltype_tmp)!pb with alya test
     !
     ! CONNEC
     !
     icount=0
     do ielem=1,nelem
        if       ( ltype(ielem) == 30 ) then ! TETRA
           pnode = 4
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 10
        else if  ( ltype(ielem) == 37 ) then ! HEXA
           pnode = 8
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 12
        else if  ( ltype(ielem) == 32 ) then ! PYRA
           pnode = 5
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 14
        else if  ( ltype(ielem) == 34 ) then ! PENTA
           pnode = 6
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 13
        end if
     enddo
     !
#ifdef VTK_DEBUG
     write(*,*)'nelem=',nelem
     write(*,*)'npoin=',npoin
     write(*,*)'icount=',icount
     write(*,*)'size(lype_tmp)=',size(ltype_tmp,KIND=ip)
     write(*,*)'lype_tmp=',ltype_tmp
     write(*,*)'size(lnods_tmp)=',size(lnods_tmp,KIND=ip)
     write(*,*)'lnods_tmp=',lnods_tmp
     write(*,*)'size(offset_tmp)=',size(offset_tmp,KIND=ip)
     write(*,*)'offset_tmp=',offset_tmp
#endif
     !
#ifdef VTK
     if (vtk_time /= ittim) then ! open vtk file
        !
        call vtkXMLWriterF_New(vtk_id)
        !
        ! 4 ---> unstruct type (see vtkdatatype.h)
        !
        call vtkXMLWriterF_SetDataObjectType(vtk_id, 4)
        !
        ! case of ASCII (0=ascii,1,2=binary)
        !
        if (kfl_outfo == 41) then
           call vtkXMLWriterF_SetDataModeType(vtk_id, 0)
        endif
        !
        ! Initialize
        !
        call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//trim(filsa)//'.vtu')
        !
        !
     endif
     !
     !POINTDATA
     !
     if (nvar /= 0) then
        if (wopos(3) /= 'VORTX') then
           if      (wopos(2)=='SCALA') then
              call vtkXMLWriterF_SetPointData(vtk_id,wopos(1),11_4,bridge,INT(npoin,8),1,wopos(2))
           else if (wopos(2)=='VECTO') then
              call vtkXMLWriterF_SetPointData(vtk_id,wopos(1),11_4,bridge,INT(npoin,8),3,wopos(2))
           else
              call runend('TYPE OF VARIABLE NOT CODED')
           end if
        else
           call vortwr()
        end if
        ncoun_pos=ncoun_pos+1
     end if
     !
     !
     if (ncoun_pos == nvar ) then ! close the vtk file
        !
        !COORD (11_4-->FLOAT64 see vtktype.h)
        !
        call vtkXMLWriterF_SetPoints(vtk_id, 11_4, coord, INT(npoin,8))
        !
        !CONNEC
        !
        call vtkXMLWriterF_SetCellsWithTypes(vtk_id, ltype_tmp , INT(nelem,8), lnods_tmp, INT(icount,8))
        !
        !FLUSH
        !
        call vtkXMLWriterF_Write(vtk_id, ierr)
        !
        !CLOSE
        !
        Call vtkXMLWriterF_Delete(vtk_id)
        !
     endif
     !
#endif
     !
     ! Deallocate memory temporary
     !
     deallocate(lnods_tmp)
     deallocate(offset_tmp)
     deallocate(ltype_tmp)
     !
     ! writing of PVD (in case of time serie)
     !
     if (kfl_rstar == 0 ) then !(no restart)
        if (vtk_time == -1.0_rp) then ! write the first time step (no restart)
           vtk_time=ittim
           open (unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',status='unknown')
           write(masterpvd,'(a)')'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
           write(masterpvd,'(a)')' <Collection>'
           !
           write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.vtu"/>'
           !
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           close(unit=masterpvd)
        else if (ittim /= vtk_time) then ! write the new time step (no restart)
           vtk_time=ittim
           open(unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',position='append',status='old')
           backspace(masterpvd)
           backspace(masterpvd)
           write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.vtu"/>'
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           close(unit=masterpvd)
        end if
     else if (kfl_rstar == 2 .and. ittim /= vtk_time)  then !(restart)
        vtk_time=ittim
        !
        !Check if pvd file exist
        !
        inquire(file=trim(title)//'_vtk/'//trim(title)//'.pvd', exist=dir_e)
        if ( dir_e ) then
           !write(*,*)"exist file"
        else
           open (unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',status='unknown')
           write(masterpvd,'(a)')'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
           write(masterpvd,'(a)')' <Collection>'
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           !
        end if
        !
        open(unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',position='append',status='old')
        backspace(masterpvd)
        backspace(masterpvd)
        write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.vtu"/>'
        write(masterpvd,'(a)')' </Collection>'
        write(masterpvd,'(a)')'</VTKFile>'
        close(unit=masterpvd)
     end if

  end if

  if(IMASTER)then
     !
     ! Master PVTU (this function does not exist for fortran in VTK)
     !
     !
     !Check if result folder exist
     !
     inquire(file=trim(title)//'_vtk', exist=dir_e)
     if ( dir_e ) then
        !write(*,*)"exist file"
     else
        call execute_command_line ('mkdir '//trim(title)//'_vtk 2> /dev/null')
     end if
     !
     !
     if (vtk_time /= ittim) then ! open vtk file (writing of pvtu file)
        !
        !
        open (unit=mastervtk,file= trim(title)//'_vtk/'//trim(filsa)//'.pvtu',status='unknown')
        write(mastervtk,'(a)')'<?xml version="1.0"?>'
        write(mastervtk,'(a)')'<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
        write(mastervtk,'(a)')'  <PUnstructuredGrid GhostLevel="0">'
        write(mastervtk,'(a)')'    <PPointData>'
     endif
     !
     !
     if (nvar /= 0 ) then
        if (wopos(3) /= 'VORTX') then
           if      (wopos(2)=='SCALA') then
              write(mastervtk,'(a,a,a)')'    <PDataArray type="Float64" Name="',wopos(1),'" NumberOfComponents="1"/>'
           else if (wopos(2)=='VECTO') then
              write(mastervtk,'(a,a,a)')'    <PDataArray type="Float64" Name="',wopos(1),'" NumberOfComponents="3"/>'
           end if
        else
           call vortwr()
        endif
        ncoun_pos=ncoun_pos+1
     end if
     !
     !
     if (ncoun_pos == nvar ) then ! close the vtk file
        write(mastervtk,'(a)')'    </PPointData>'
        write(mastervtk,'(a)')'    <PPoints>'
        write(mastervtk,'(a)')'      <PDataArray type="Float64" NumberOfComponents="3"/>'
        write(mastervtk,'(a)')'    </PPoints>'
        do ipart = 0,npart-1
           if      (ipart<10)then
              write(mastervtk,'(a,a,a,i1,a)')'    <Piece Source="',trim(filsa),'_',ipart,'.vtu"/>'
           else if (ipart<100)then
              write(mastervtk,'(a,a,a,i2,a)')'    <Piece Source="',trim(filsa),'_',ipart,'.vtu"/>'
           else if (ipart<1000)then
              write(mastervtk,'(a,a,a,i3,a)')'    <Piece Source="',trim(filsa),'_',ipart,'.vtu"/>'
           else if (ipart<10000)then
              write(mastervtk,'(a,a,a,i4,a)')'    <Piece Source="',trim(filsa),'_',ipart,'.vtu"/>'
           else if (ipart<100000)then
              write(mastervtk,'(a,a,a,i5,a)')'    <Piece Source="',trim(filsa),'_',ipart,'.vtu"/>'
           end if
        end do
        write(mastervtk,'(a)')'  </PUnstructuredGrid>'
        write(mastervtk,'(a)')'</VTKFile>'
        close(unit=mastervtk)
     endif
     !
     ! writing of PVD file (case of time serie)
     !
     if (kfl_rstar == 0 ) then !(no restart)
        if (vtk_time == -1.0_rp) then ! write the first time step
           vtk_time=ittim
           open (unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',status='unknown')
           write(masterpvd,'(a)')'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
           write(masterpvd,'(a)')' <Collection>'
           !
           write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.pvtu"/>'
           !
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           close(unit=masterpvd)
           !
           !
        else if (ittim /= vtk_time) then ! write the new time step
           vtk_time=ittim
           open(unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',position='append',status='old')
           backspace(masterpvd)
           backspace(masterpvd)
           write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.pvtu"/>'
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           close(unit=masterpvd)
        end if
     else if (kfl_rstar == 2 .and. ittim /= vtk_time)  then !(restart)
        vtk_time=ittim
        !
        !Check if pvd file exist
        !
        inquire(file=trim(title)//'_vtk/'//trim(title)//'.pvd', exist=dir_e)
        if ( dir_e ) then
           !write(*,*)"exist file"
        else
           open (unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',status='unknown')
           write(masterpvd,'(a)')'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
           write(masterpvd,'(a)')' <Collection>'
           write(masterpvd,'(a)')' </Collection>'
           write(masterpvd,'(a)')'</VTKFile>'
           !
        end if
        !
        open(unit=masterpvd,file= trim(title)//'_vtk/'//trim(title)//'.pvd',position='append',status='old')
        backspace(masterpvd)
        backspace(masterpvd)
        write(masterpvd,'(a,i4,a,a,a)')'  <DataSet timestep="',ittim,'" group="" part="0" file="./',trim(filsa),'.pvtu"/>'
        write(masterpvd,'(a)')' </Collection>'
        write(masterpvd,'(a)')'</VTKFile>'
        close(unit=masterpvd)
     end if

  end if


  if (ISLAVE) then
     !
     ! Slave Binary (vtk lib)
     !
     ! test lib phillipe  prob with integer x4 x8 -----> alya est av integer x4(simple) et real x8 (double)
     !                                            -----> vtk etait pas defaut int x8 et real x8
     !check the length of lnods_tmp (check +1 is pnode+1)
     !
     !
     ! Check if result folder exist
     !
     inquire(file=trim(title)//'_vtk', exist=dir_e)
     if ( dir_e ) then
        !write(*,*)"exist file"
     else
        call execute_command_line ('mkdir '//trim(title)//'_vtk 2> /dev/null')
     end if
     !
     icount=0
     do ielem=1,nelem
        if       ( ltype(ielem) == 30 ) then ! TETRA
           icount=icount+4+1
        else if  ( ltype(ielem) == 37 ) then ! HEXA
           icount=icount+8+1
        else if  ( ltype(ielem) == 32 ) then ! PYRA
           icount=icount+5+1
        else if  ( ltype(ielem) == 34 ) then ! PENTA
           icount=icount+6+1
        endif
     enddo
     !
     ! allocate memory temporary
     !
     allocate (lnods_tmp(icount),stat=istat)
     !call memchk(zero,istat,memor_dom,'lnods_tmp','geovtk',lnods_tmp)
     allocate (offset_tmp(nelem),stat=istat)
     !call memchk(zero,istat,memor_dom,'offset_tmp','geovtk',offset_tmp)
     allocate (ltype_tmp(nelem),stat=istat)
     !call memchk(zero,istat,memor_dom,'ltype_tmp','geovtk',ltype_tmp)
     !
     ! CONNEC
     !
     icount=0
     do ielem=1,nelem
        if       ( ltype(ielem) == 30 ) then ! TETRA
           pnode = 4
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 10
        else if  ( ltype(ielem) == 37 ) then ! HEXA
           pnode = 8
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 12
        else if  ( ltype(ielem) == 32 ) then ! PYRA
           pnode = 5
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 14
        else if  ( ltype(ielem) == 34 ) then ! PENTA
           pnode = 6
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1
           end do
           offset = offset + pnode
           offset_tmp(ielem) = offset
           ltype_tmp(ielem) = 13
        end if
     enddo
     !
     !
#ifdef VTK_DEBUG
     if (kfl_paral==2) then
        write(*,*)'nelem=',nelem
        write(*,*)'npoin=',npoin
        write(*,*)'icount=',icount
        write(*,*)'size(lype_tmp)=',size(ltype_tmp,KIND=ip)
        write(*,*)'lype_tmp=',ltype_tmp
        write(*,*)'size(lnods_tmp)=',size(lnods_tmp,KIND=ip)
        write(*,*)'lnods_tmp=',lnods_tmp
        write(*,*)'size(offset_tmp)=',size(offset_tmp,KIND=ip)
        write(*,*)'offset_tmp=',offset_tmp
     end if
#endif
     !
     !
#ifdef VTK
     if (ncoun_pos == 0) then ! open vtk file
        !
        !
        call vtkXMLWriterF_New(vtk_id)
        !
        ! 4 ---> unstruct type
        !
        call vtkXMLWriterF_SetDataObjectType(vtk_id, 4)
        !
        ! case of ASCII (0=ascii,1,2=binary)
        !
        if (kfl_outfo == 41) then
           call vtkXMLWriterF_SetDataModeType(vtk_id, 0)
        endif
        !
        ! Initialize
        !
        !
        ! test writing data
        !
        call date_and_time(VALUES=values)
        !auxi1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
        auxi1=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
        !
        !
        !
        call get_environment_variable("TMPDIR", my_tmpdir)
        !
        !  the case of kfl_vtk=0 doesn't work so kfl_vtk=1 by default
        !  = 0 when making test performance for writing the data on scratch
        !
        kfl_vtk=1
        ipart = kfl_paral - 1
        if      (ipart<10)then
           write(nunam_pos2,'(a,a,i1,a)') trim(filsa),'_',ipart,'.vtu'
           if ( kfl_vtk == 0 ) then
              call vtkXMLWriterF_SetFileName(vtk_id, trim(my_tmpdir)//'/'//nunam_pos2 ) ! write on the scratch
           else
              call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//nunam_pos2) ! write on the gpfs
           endif
        else if (ipart<100)then
           write(nunam_pos2,'(a,a,i2,a)') trim(filsa),'_',ipart,'.vtu'
           if ( kfl_vtk == 0 ) then
              call vtkXMLWriterF_SetFileName(vtk_id, trim(my_tmpdir)//'/'//nunam_pos2 ) ! write on the scratch
           else
              call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//nunam_pos2) ! write on the gpfs
           endif
        else if (ipart<1000)then
           write(nunam_pos2,'(a,a,i3,a)') trim(filsa),'_',ipart,'.vtu'
           if ( kfl_vtk == 0 ) then
              call vtkXMLWriterF_SetFileName(vtk_id, trim(my_tmpdir)//'/'//nunam_pos2 ) ! write on the scratch
           else
              call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//nunam_pos2) ! write on the gpfs
           endif
        else if (ipart<10000)then
           write(nunam_pos2,'(a,a,i4,a)') trim(filsa),'_',ipart,'.vtu'
           if ( kfl_vtk == 0 ) then
              call vtkXMLWriterF_SetFileName(vtk_id, trim(my_tmpdir)//'/'//nunam_pos2 ) ! write on the scratch
           else
              call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//nunam_pos2) ! write on the gpfs
           endif
        else if (ipart<100000)then
           write(nunam_pos2,'(a,a,i5,a)') trim(filsa),'_',ipart,'.vtu'
           if ( kfl_vtk == 0 ) then
              call vtkXMLWriterF_SetFileName(vtk_id, trim(my_tmpdir)//'/'//nunam_pos2 ) ! write on the scratch
           else
              call vtkXMLWriterF_SetFileName(vtk_id, trim(title)//'_vtk/'//nunam_pos2) ! write on the gpfs
           endif
        end if
     end if
     !
     !POINTDATA
     !
     if (nvar /= 0) then
        if (wopos(3) /= 'VORTX') then
           if      (wopos(2)=='SCALA') then
              call vtkXMLWriterF_SetPointData(vtk_id,wopos(1),11_4,bridge,INT(npoin,8),1,wopos(2))
           else if (wopos(2)=='VECTO') then
              call vtkXMLWriterF_SetPointData(vtk_id,wopos(1),11_4,bridge,INT(npoin,8),3,wopos(2))
           else
              call runend('TYPE OF VARIABLE NOT CODED')
           end if
        else
           call vortwr()
        endif
        ncoun_pos=ncoun_pos+1
     end if
     !
     !
     if (ncoun_pos == nvar ) then ! close the vtk file
        !
        !COORD
        !
        call vtkXMLWriterF_SetPoints(vtk_id, 11_4, coord, INT(npoin,8))
        !
        !CONNEC
        !
        call vtkXMLWriterF_SetCellsWithTypes(vtk_id, ltype_tmp, int(nelem,8), lnods_tmp, int(icount,8))
        !
        !FLUSH
        !
        call vtkXMLWriterF_Write(vtk_id, ierr)
        !
        !CLOSE
        !
        call vtkXMLWriterF_Delete(vtk_id)
        !
     endif
     !
     ! test writing data
     !
     call date_and_time(VALUES=values)
     !auxi2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     auxi2=(values(5)*3600+values(6)*60+values(7))*1000 + values(8)   ! miniseconds
     auxi=auxi2-auxi1
     total_auxi=total_auxi + auxi
     !write (*,*),'kfl_paral=',kfl_paral,'writing data',total_auxi
     !
     !
     !
#endif
     !
     ! Deallocate memory temporary
     !
     deallocate(lnods_tmp)
     deallocate(offset_tmp)
     deallocate(ltype_tmp)

  end if

end subroutine geovtk
