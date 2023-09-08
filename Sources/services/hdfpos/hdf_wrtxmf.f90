subroutine hdf_wrtxmf(itask)
  !-----------------------------------------------------------------------
  !****f* hdfpos/hdf_wrtxmf
  ! NAME
  !    hdf_
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... 
  !    ITASK = 2 ... 
  ! USED BY
  !    Hdfpos
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_hdfpos
  use def_domain
  use def_postpr
  use def_parall
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif
  integer(ip),    intent(in) :: itask
  integer(ip)                :: ielem
  integer(ip)                :: ierr
  character(150), save       :: filsa
  logical                    :: dir_e

  !
  ! Save the name
  !
  filsa = adjustl(trim(namda))


  select case ( itask )

  case ( 1_ip )

     if (IMASTER) then
        hdf5_connedim = 0
     endif

#ifdef MPI_OFF
#else
     call MPI_REDUCE( hdf5_connedim, hdf5_condifin, 1, MPI_INTEGER8, &
          MPI_SUM, 0, PAR_COMM_MY_CODE, ierr )
#endif

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_wrtxmf: hdf5_connedim:', hdf5_connedim
     write(*,*) kfl_paral, 'hdf_wrtxmf: hdf5_condifin:', hdf5_condifin
#endif



     !-------------------------------------------------------------------
     !
     ! geo mesh ()
     !
     !-------------------------------------------------------------------
     !
     ! 
     write(666,*)''
     write(666,*)'   <Grid Name="',trim(filsa),'" Type="Uniform" >'
     write(666,'(a,e16.8E3,a)')'    <Time Type="Single" Value="',cutim,'" />'
     write(666,'(a,i12,a)')'    <Topology Type="Mixed" NumberOfElements="',nelem_total,'">'
     write(666,'(a,i12,a,a,a)')'     <DataStructure Dimensions="',hdf5_condifin,&
          '" NumberType="Int" Format="HDF" >',trim(fil_postp),':/CONNE </DataStructure>'
     write(666,*)'    </Topology>'
     write(666,*)'    <Geometry GeometryType="XYZ">'
     write(666,'(a,i12,a,a,a)')'     <DataStructure Dimensions="3',npoin_total, &
          '" NumberType="Float" Presicion="8" Format="HDF" >',trim(fil_postp), &
          ':/COORD </DataStructure>'
     write(666,*)'    </Geometry>'
     write(666,*)'   </Grid>'





  case ( 2_ip )
     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for scalar output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !
     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     hdf5_i=ittim
     !
     ! Compose file name: example-1234567.h5 ! cutim !
     !
     write(666,*)'    <Attribute Name="',wopos_hdf(1),&
          '" Center="Node" AttributeType="Scalar"> '
     write(666,'(a,i12,a,a,a,a,a)') &
          '     <DataStructure Format="HDF" DataType="Float" Precision="8" Dimensions="' &
          ,npoin_total,'">',trim(filsa)//'-'//adjustl(trim(nunam_pos)),'.h5:/',wopos_hdf(1), &
          ' </DataStructure>'
     write(666,*)'    </Attribute>'



  case ( 3_ip )
     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for vector output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !
     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     hdf5_i=ittim
     !
     ! Compose file name: example-1234567.h5
     !
     write(666,*)'    <Attribute Name="',wopos_hdf(1),&
          '" Center="Node" AttributeType="Vector"> '
     write(666,'(a,i12,a,a,a,a,a)') &
          '     <DataStructure Format="HDF" DataType="Float" Precision="8" Dimensions="3',   &
          npoin_total,'">',trim(filsa)//'-'//adjustl(trim(nunam_pos)),'.h5:/',wopos_hdf(1), & 
          ' </DataStructure>'
     write(666,*)'    </Attribute>'


  case ( 4_ip )
     !-------------------------------------------------------------------
     !   close the file
     !
     !------------------------------------------------------------------
     backspace(666)
     backspace(666)
     backspace(666)
     backspace(666)

     write(666,*)''
     write(666,*)'  </Grid>'
     write(666,*)' </Domain>'
     write(666,*)'</Xdmf>'

     close(unit=666)


  case ( 5_ip )
     !-------------------------------------------------------------------
     !   open the file and write the entete 
     !   initialize hdf5_connedim to 0
     !------------------------------------------------------------------ 

     if (IMASTER) then
        hdf5_connedim = 0
     endif
#ifdef MPI_OFF
#else
     call MPI_REDUCE( hdf5_connedim, hdf5_condifin, 1, MPI_INTEGER8, &
          MPI_SUM, 0, PAR_COMM_MY_CODE, ierr )
#endif

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_wrtxmf: hdf5_connedim:', hdf5_connedim
     write(*,*) kfl_paral, 'hdf_wrtxmf: hdf5_condifin:', hdf5_condifin
#endif
     !
     ! case of restart 
     ! Check if result folder exist
     !
     inquire(file=trim(filsa)//'.xmf', exist=dir_e)
     if ( dir_e ) then !  "exist file"
        if (IMASTER)then
           open(unit=666,file=trim(filsa)//'.xmf',position='append',status='old') 
        else   
           return
        endif
     else
        if (IMASTER) then
           open(unit=666,file=trim(filsa)//'.xmf',status='unknown') 

           write(666,'(a)')'<?xml version="1.0" ?>'
           write(666,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" ['
           write(666,'(a)')'<!ENTITY InitialHeavyData "data.mesh.h5">'
           write(666,'(a)')']>'
           write(666,*)''
           write(666,*)'<Xdmf>'
           write(666,*)' <Domain>'
           write(666,*)'   <Grid Name="',trim(filsa),'" GridType="Collection" CollectionType="Temporal" >'
           write(666,*)''
           write(666,*)''
           write(666,*)''
           write(666,*)''
        endif
     endif

  case ( 6_ip ) 
     !-------------------------------------------------------------------
     !   write the mesh for each new time step 
     !
     !------------------------------------------------------------------ 

     if (hdf5_i /= ittim ) then
        backspace(666)
        backspace(666)
        backspace(666)
        backspace(666)
        write(666,*)''
        write(666,*)'   <Grid Name="',trim(filsa),'" Type="Uniform" >'
        write(666,'(a,e16.8E3,a)')'    <Time Type="Single" Value="',cutim,'" />'
        write(666,'(a,i12,a)')'    <Topology Type="Mixed" NumberOfElements="',nelem_total,'">'
        write(666,'(a,i12,a,a,a)')'     <DataStructure Dimensions="',hdf5_condifin,&
             '" NumberType="Int" Format="HDF" >',trim(fil_postp),':/CONNE </DataStructure>'
        write(666,*)'    </Topology>'
        write(666,*)'    <Geometry GeometryType="XYZ">'
        write(666,'(a,i12,a,a,a)')'     <DataStructure Dimensions="3',npoin_total,&
             '" NumberType="Float" Presicion="8" Format="HDF" >',trim(fil_postp),':/COORD </DataStructure>'
        write(666,*)'    </Geometry>'
     end if


  case ( 7_ip )
     !-------------------------------------------------------------------
     !  close the <> Grid 
     !
     !------------------------------------------------------------------ 

     write(666,*)'    </Grid>'
     write(666,*)''
     write(666,*)'  </Grid>'
     write(666,*)' </Domain>'
     write(666,*)'</Xdmf>'


  case ( 8_ip )  
     !-------------------------------------------------------------------
     !  Open the <> Grid vortex 
     !
     !------------------------------------------------------------------ 
     if (IMASTER) then
        open(unit=667,file=trim(filsa)//'-vort.xmf',status='unknown') 

        write(667,'(a)')'<?xml version="1.0" ?>'
        write(667,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" ['
        write(667,'(a)')'<!ENTITY InitialHeavyData "data.mesh.h5">'
        write(667,'(a)')']>'
        write(667,*)''
        write(667,*)'<Xdmf>'
        write(667,*)' <Domain>'
        write(667,*)'   <Grid Name="',trim(filsa)//'-vort','" GridType="Collection" CollectionType="Temporal" >'
     end if

  case ( 9_ip )

     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for vortex output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !
     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     !
     ! Compose file name: example-1234567.h5
     !
     write(667,*)''
     write(667,*)'   <Grid Name="',wopos_hdf(1),'" Type="Uniform" >'
     write(667,'(a,e16.8E3,a)')'    <Time Type="Single" Value="',cutim,'" />'
     write(667,*)'    <Topology Type="Polyvertex" NumberOfElements="',nvort_total,'"/>'
     write(667,*)'     <Geometry GeometryType="XYZ">' 
     write(667,'(a,a,i12,a,a,a,a,a)')'      <DataStructure Format="HDF" ', &
          'DataType="Float" Precision="8" Dimensions="3', &
          nvort_total,'">',trim(filsa)//'-'//adjustl(trim(nunam_pos)) &
          ,'.h5:/',wopos_hdf(1),' </DataStructure>'
     write(667,*)'     </Geometry>'
     write(667,*)'    </Grid>'


  case ( 10_ip )  
     !-------------------------------------------------------------------
     !  Close the <> Grid vortex 
     !
     !------------------------------------------------------------------
     write(667,*)''
     write(667,*)'    </Grid>'
     write(667,*)' </Domain>'
     write(667,*)'</Xdmf>'

     close(unit=667)

  case ( 11_ip )  
     !-------------------------------------------------------------------
     !  Open the <> Grid filter 
     !
     !------------------------------------------------------------------ 
     if (IMASTER) then
        open(unit=668,file=trim(filsa)//'-filter.xmf',status='unknown') 
        write(668,'(a)')'<?xml version="1.0" ?>'
        write(668,'(a)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" ['
        write(668,'(a)')'<!ENTITY InitialHeavyData "data.mesh.h5">'
        write(668,'(a)')']>'
        write(668,*)''
        write(668,*)'<Xdmf>'
        write(668,*)' <Domain>'
        write(668,*)'   <Grid Name="',trim(filsa)//'-filter','" GridType="Collection" CollectionType="Temporal" >'
        !write(668,*)''
        !write(668,*)''
        !write(668,*)''
        !write(668,*)''
     end if

  case ( 12_ip )
     !-------------------------------------------------------------------
     !   write the mesh filtered for each new time step 
     !
     !------------------------------------------------------------------ 
     if (IMASTER) then
        !if (hdf5_i /= ittim ) then 
        !backspace(668)
        !backspace(668)
        !backspace(668)
        !backspace(668)
        write(668,*)''
        write(668,*)'   <Grid Name="',trim(filsa),'" Type="Uniform" >'
        write(668,'(a,e16.8E3,a)')'    <Time Type="Single" Value="',cutim,'" />'
        write(668,'(a,i12,a)')'    <Topology Type="Mixed" NumberOfElements="',nelem_total_filt,'">'
        write(668,'(a,i12,a,a,a)')'     <DataStructure Dimensions="',hdf5_condifin,&
             '" NumberType="Int" Format="HDF" >',trim(fil_postp),':/CONNE </DataStructure>'
        write(668,*)'    </Topology>'
        write(668,*)'    <Geometry GeometryType="XYZ">'
        write(668,'(a,i12,a,a,a)')'     <DataStructure Dimensions="3',npoin_total_filt,&
             '" NumberType="Float" Presicion="8" Format="HDF" >',trim(fil_postp),':/COORD </DataStructure>'
        write(668,*)'    </Geometry>'

        !end if
     end if


  case ( 13_ip )
     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for scalar filter output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !

     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     hdf5_i=ittim
     !
     ! Compose file name: example-1234567.h5 ! cutim !
     !
     write(668,*)'    <Attribute Name="',wopos_hdf(1),&
          '" Center="Node" AttributeType="Scalar"> '
     write(668,'(a,i12,a,a,a,a,i1,a,a,a)') &
          '     <DataStructure Format="HDF" DataType="Float" Precision="8" Dimensions="' &
          ,npoin_total_filt,'">',trim(filsa)//'-'//adjustl(trim(nunam_pos)),'.h5:', &
          '/FILTER',kfl_filte,'/',wopos_hdf(1),' </DataStructure>'
     write(668,*)'    </Attribute>'

  case ( 14_ip )
     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for scalar filter output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !

     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     hdf5_i=ittim
     !
     ! Compose file name: example-1234567.h5 ! cutim !
     !
     write(668,*)'    <Attribute Name="',wopos_hdf(1),&
          '" Center="Node" AttributeType="Vector"> '
     write(668,'(a,i12,a,a,a,a,i1,a,a,a)') &
          '     <DataStructure Format="HDF" DataType="Float" Precision="8" Dimensions="' &
          ,npoin_total_filt,'">',trim(filsa)//'-'//adjustl(trim(nunam_pos)),'.h5:', &
          '/FILTER',kfl_filte,'/',wopos_hdf(1),' </DataStructure>'
     write(668,*)'    </Attribute>'



  case ( 15_ip )  
     !-------------------------------------------------------------------
     !  Close the <> Grid filter output
     !
     !------------------------------------------------------------------
     write(668,*)''
     write(668,*)'    </Grid>'
     write(668,*)' </Domain>'
     write(668,*)'</Xdmf>'

     close(unit=668)



  end select

  !
  ! Format
  !
  !200 format('FIELD<double> ',a,'("',a,'",',i12,',',i12,');')
  !201 format('FIELD<float>  ',a,'("',a,'",',i12,',',i12,');')
  !205 format('TEXTE Time(" ',e13.6,'");')
  !210 format('SOLUTION Solution( ) =',/,'{')
  !220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
  !221 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,',',a,');')
  !230 format('};')


end subroutine hdf_wrtxmf

