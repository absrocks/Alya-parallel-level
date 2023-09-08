subroutine Cgns24(order)
!-----------------------------------------------------------------------
!****f* Services/cgns24
! NAME
!    cgns24
! DESCRIPTION
!    This routine
! USED BY
!    outdom
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      mod_memchk
  implicit none
  include  'cgnslib_f.h'
  integer(ip),  intent(in) :: order
  character(50)            :: basename,zonename,elname
  integer(4)               :: icelldim,iphysdim
  integer(4)               :: isize(3),ier
  integer(4)               :: nelem_start,nelem_end,nbdyelem
  integer(4)               :: index_section,index_base
  integer(4)               :: index_file,index_zone,index_coord,index_flow
  integer(4)               :: kfl_eltyp,nele2,ielty,inode,ielem,istat
  integer(4)               :: one4=1,two4=2,three4=3,five4=5
  integer(4)               :: idime,ipoin,ibopo
  integer(4), allocatable  :: lnod2(:,:)
  real(8)                  :: exno2(ndime,npoin)

  select case(order)

  case(1)

     call cg_open_f(trim(fil_outpu_dom),MODE_WRITE,index_file,ier)
     if(ier/=0) call cgnser

  case (2)
     

     basename='Alya problem: '//trim(title)
     icelldim=ndime
     iphysdim=ndime
     call cg_base_write_f(index_file,basename,icelldim,iphysdim,&
          index_base,ier) 
     if(ier/=0) call cgnser

     call cg_simulation_type_write_f(index_file,index_base,UserDefined,ier) 
     if(ier/=0) call cgnser

     
     !  create zone
     zonename = 'Zone 1'                     !  define zone name    
     isize(1)=npoin                          !  vertex size  
     isize(2)=nelem                          !  cell size
     isize(3)=0                              !  boundary vertex size (zero if elements not sorted)
       
     call cg_zone_write_f(index_file,index_base,zonename,isize,& 
          Unstructured,index_zone,ier)
     if(ier/=0) call cgnser

     ! put DataClass and DimensionalUnits under Base
      !call cg_goto_f(index_file,index_base,ier,'end')
      !call cg_dataclass_write_f(Dimensional,ier)
      !call cg_units_write_f(Kilogram,Meter,Second,Kelvin,Degree,ier)
          
     ! write grid coordinates (user must use SIDS-standard names here)
     if(ier/=0) call cgnser
     call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordinateX',coord(1,1:npoin),one4,ier)
     if(ier/=0) call cgnser
     call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordinateY',coord(2,1:npoin),two4,ier)
     if(ier/=0) call cgnser
     if(ndime==3)&
          call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordinateZ',coord(3,1:npoin),three4,ier)
     if(ier/=0) call cgnser
 
     nbdyelem    = 0         !  unsorted boundary elements
     nelem_start = 1
     nelem_end   = 0
     
     ! Loop over element types
     do ielty=1,nelty

        allocate(lnod2(nnode(ielty),nelem),stat=istat)
        call memchk(zero,istat,memor_dom,'CGNS24','cgns24',lnod2)
        nele2=0
        do ielem=1,nelem 
           if(ltype(ielem)==ielty) then
              nele2=nele2+1
              do inode=1,nnode(ielty)
                 lnod2(inode,nele2)=lnods(inode,ielem)
              end do
           end if
        end do 

        if(ndime==2) then
           if(nnode(ielty)==3) then 
              kfl_eltyp=TRI_3
              elname='TRI_3'
           else if(nnode(ielty)==4) then 
              kfl_eltyp=QUAD_4
              elname='QUAD_4'
           else if(nnode(ielty)==6) then
              kfl_eltyp=TRI_6
              elname='TRI_6'
           else if(nnode(ielty)==9) then
              kfl_eltyp=QUAD_9
              elname='QUAD_9'
           end if
        else
           if(nnode(ielty)==4) then 
              kfl_eltyp=TETRA_4
              elname='TETRA_4'
           else if(nnode(ielty)==6) then
              kfl_eltyp=PENTA_6
              elname='PENTA_6'
           else if(nnode(ielty)==8) then
              kfl_eltyp=HEXA_8
              elname='HEXA_8'
           else if(nnode(ielty)==27) then 
              kfl_eltyp=HEXA_27
              elname='HEXA_27'
           end if
        end if

        !nelem_start=1
        !nelem_end=nele2
        nelem_end=nelem_end+nele2

        !  write element connectivity
        call cg_section_write_f(index_file,index_base,index_zone,&
             trim(elname),kfl_eltyp,nelem_start,nelem_end,nbdyelem,&
             lnod2,index_section,ier)
        if(ier/=0) call cgnser

        nelem_start=nelem_start+nele2

        call memchk(two,istat,memor_dom,'LNOD2','cgns24',lnod2)
        deallocate(lnod2,stat=istat)
        if(istat/=0) call memerr(two,'LNOD2','ccgns24',0_ip)

     end do

     ! create flow solution node

     index_flow=1
     call cg_sol_write_f(index_file,index_base,index_zone,'Flow solution',&
          Vertex,index_flow,ier)

      call cg_dataclass_write_f(Dimensional,ier)
      call cg_units_write_f(Kilogram,Meter,Second,Kelvin,Degree,ier)

    do ipoin=1,npoin
        ibopo=lpoty(ipoin)
        if(ibopo>0) then
           do idime=1,ndime
              exno2(idime,ipoin)=exnor(idime,1,ibopo)
           end do
        else
           do idime=1,ndime
              exno2(idime,ipoin)=0.0_rp
           end do
        end if
     end do
     isize(1)=ndime
     isize(2)=npoin
     elname='EXNOR'
     call cg_field_write_f(index_file,index_base,index_zone,index_flow,&
          RealDouble,elname,exno2,isize,ier)     
     if(ier/=0) call cgnser
     !call cg_goto_f(index_file,index_base,zonename,1,'end',ier)
     !call cg_dataclass_write_f(Dimensional,ier)
     !if(ier/=0) call cgnser
     
     !  close CGNS file
     call cg_close_f(index_file,ier)
     if(ier/=0) call cgnser     
     
  end select

end subroutine Cgns24
