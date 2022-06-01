program periodic

  implicit none

  real*8,    parameter :: pi=3.141592653589793238462643383279502884197_8
  real*8,    parameter :: twopi=6.283185307179586476925286766559005768394_8
  character*5          :: c
  character*40         :: fil_meanv, fil_coord, fil_fix, fil_per, fil_inp, types, lineb
  character( Len=5000) :: line
  logical              :: duplicate, periodicity, node_is_periodic, back, iex, hex08, pen06
  logical, allocatable :: nboub(:,:), new_wall(:,:,:), node_belongs(:)
  real*8, allocatable  :: coord_node(:), coord_bound(:,:,:), coord_unique(:,:), coord(:,:,:), coord_all(:,:), veloc(:,:)
  real*8               :: start_time, stop_time, length, timestep, bcoor(4,3), cog(3), bnorm(3), bcog(3), produ, &
                          bnorm_aux(3), zvec(3), znorm, norm_aux, cosine, anglex, angley, temp1, temp2
  real*8               :: dot_prod
  integer              :: i,j,k,itask,inp_set,iset,io,ii,m,n,l,ielem,icoun
  integer              :: nlines, nboun, nelem, nperi, npoin, number_points, iper1, ilength, &
                          iper2, nper, nset, isteps, npoin_inlet, iboun, ibound(4), iboundary, &
                          nextra, ndime, inode
  integer, allocatable :: id_bound(:), node(:), nodb(:), node_bound(:,:), node_new(:,:), &
                          node_unique(:), node_comp(:), node_per(:,:), &
                          new_per(:,:), isetb(:), nbset(:), ibous(:), &
                          node_wall(:,:), nbcod(:), nnode(:,:), nnodb(:,:,:), &
                          node_per_bound(:), node_count(:), node_bset(:), bound_all(:,:)

  hex08 = .false. 
  pen06 = .false.

  open (unit=101,file='mesh.inp')
  read(101,*) itask

  if (itask == 0) then ! Inlet mesh generation

     call cpu_time(start_time)
     periodicity = .false.

     !
     !  Read input variables
     !
   
     read(101,*) i
     if (i==1) periodicity = .true.
   
     read(101,*) inp_set   ! inflow set
     read(101,*) length
     read(101,*) isteps
     read(101,*) nset      ! number of sets of non­periodic boun that share nodes with the inlet bound - guessed from p3 inlet.pdf
     allocate(isetb(nset))
     read(101,*) isetb(1:nset)
   
     timestep = length/isteps
   
     !
     !  Get the problem name (for obtaining the geometry files)
     !
   
     call getarg(1, fil_meanv)
   
     !
     !  Obtain geo, fix and per files
     !
   
     fil_coord = trim(fil_meanv)//'.geo.dat'
   
     inquire(file=fil_coord, EXIST=iex)
     if (.not.iex) then
        fil_coord = trim(fil_meanv)//'.dom.geo'
        inquire(file=fil_coord, EXIST=iex)
        if (.not.iex) then
           fil_coord = trim(fil_meanv)//'.geo'
           inquire(file=fil_coord, EXIST=iex)
           if (.not.iex) then 
              print*, 'Geometry file (.geo.dat, .dom.geo or .geo) not found'
              stop
           end if
        end if
     end if
   
     fil_fix = trim(fil_meanv)//'.fix'
   
     inquire(file=fil_fix, EXIST=iex)
     if (.not.iex) then
        fil_fix = trim(fil_meanv)//'.fix.dat'
        inquire(file=fil_fix, EXIST=iex)
        if (.not.iex) then
           fil_fix = trim(fil_meanv)//'.dom.fix'
           inquire(file=fil_fix, EXIST=iex)
           if (.not.iex) then
              print*, 'Boundary file (.fix.dat, .fix or .dom.fix) not found'
              stop
           end if
        end if
     end if
   
     if (periodicity) then 
   
        fil_per = trim(fil_meanv)//'.per'
     
        inquire(file=fil_per, EXIST=iex)
        if (.not.iex) then
           fil_per = trim(fil_meanv)//'.per.dat'
           inquire(file=fil_per, EXIST=iex)
           if (.not.iex) then
              print*, 'Periodicity file (.per or .per.dat) not found'
              stop
           end if
        end if
   
     end if
   
     !
     !  Open input files
     !
   
     open (unit=11,file=fil_coord)
     open (unit=12,file=fil_fix)
     if (periodicity) open (unit=13,file=fil_per)
   
     !
     !  Open output files
     !
   
     open (unit=14,file='inlet.geo.dat')
     open (unit=15,file='inlet.fix')
     open (unit=16,file='inlet.per')
     open (unit=17,file='BouCoords.dat')
     open (unit=18,file='boundaries.dat')
     open (unit=19,file='inlet.dom.dat')
     open (unit=20,file='angle.dat')
     open (unit=500,file='gen_inlet.log')
   
     !
     !  Obtain the total number of nodes, the total number of periodic nodes,  
     !  the total number of boundaries and the number of boundaries at the inlet                          
     !                                                 
     print*, 'Number of elements in the direction of the extrusion', isteps
     write(500,*) 'Number of elements in the direction of the extrusion', isteps
   
     print*, 'Size of elements in the direction of the extrusion', timestep
     write(500,*) 'Size of elements in the direction of the extrusion', timestep
     
     read(11,*) c
     do while(c/='COORD')
        read(11,*) c
     end do
   
     number_points = 0
     do while(c/='END_C')
        number_points = number_points + 1
        read(11,*) c
     end do
     number_points = number_points - 1
     rewind(11)
   
     print*, 'Number of nodes of original mesh:', number_points
     write(500,*) 'Number of nodes of original mesh:', number_points
   
     if (periodicity) then
        nper = 0
        do
           read(13,*,iostat=io)
           if (io/=0) exit
           nper = nper + 1
        end do
        print*, 'Number of periodic nodes:', nper
        write(500,*) 'Number of periodic nodes:', nper
        rewind(13)
     end if
   
     nlines = 0
     do
        read(12,*,iostat=io)
        if (io/=0) exit
        nlines = nlines + 1
     end do
     print*, 'Number of boundaries:', nlines
     write(500,*) 'Number of boundaries:', nlines
     rewind(12)
   
     nboun = 0
     do i=1, nlines
        read(12,*) iboun, iset
        if (iset == inp_set) then
           nboun = nboun + 1
        end if
     end do
     print*, 'Number of boundaries at the inlet:', nboun
     write(500,*) 'Number of boundaries at the inlet:', nboun
     rewind(12)
   
     !
     !  Allocate structures
     !
   
     allocate(id_bound(nboun))
     allocate(node_count(nlines))
     allocate(node_per_bound(nboun))
     allocate(node_bound(nboun,4))
     allocate(bound_all(nlines,4))
     allocate(node(4))
     allocate(coord_all(number_points,3))
     allocate(coord_bound(nboun,4,3))
     allocate(coord_node(3))
     allocate(node_new(nboun,4))
     allocate(node_unique(nboun*4))
     allocate(node_comp(nboun*4))
     if (periodicity) allocate(node_per(nper,2))
     allocate(node_belongs(number_points))
     allocate(coord_unique(nboun*4,3))
     allocate(nbset(nset))
     allocate(nodb(4))
     
     !
     !  Storing node coordinates
     !
   
     read(11,*) c
     do while(c/='COORD')
        read(11,*) c
     end do
     do i=1, number_points
        read(11,*) ii, coord_all(i,1:3)
     end do
     rewind(11)
   
     !
     !  Read boundaries at the inlet and boundary nodes
     !
   
     print*, 'Reading inlet boundaries'
     write(500,*) 'Reading inlet boundaries'
   
     k = 0
     do i=1, nlines
        read(12,*) iboun, iset
        if (iset == inp_set) then
           k = k + 1
           id_bound(k) = iboun
        end if
     end do
     rewind(12)
   
     k=1
     read(11,*) c
     do while(c/='BOUND')
        read(11,*) c
        k=k+1
     end do
     rewind(11)
   
     do i=1,k
        read(11,'(A)') lineb
     end do
   
     nextra = 1
     if (SCAN(lineb, "L")>0) then
        if (SCAN(lineb, "M")>0) then
           nextra = 2
        end if
     end if
     print*, 'nextra is', nextra
     rewind(11)
   
     read(11,*) c
     do while(c/='BOUND')
        read(11,*) c
     end do
     do i=1, nlines
        read(11, '(a)') line
        back = .true.
        ilength = len_trim(line)
        k = index(line(1:ilength), ' ', back)
        if (k == 0) then
           node_count(i) = 0
!           return  !error #6353: A RETURN statement is invalid in the main program.
        end if
   
        node_count(i) = 1
        do
           do
              if (k .le. 0) exit
              if (line(k:k) == ' ') then
                 k = k - 1
                 cycle
              else
                 node_count(i) = node_count(i) + 1
                 exit
              end if
           end do
           do
              if (k .le. 0) exit
              if (line(k:k) /= ' ') then
                 k = k - 1
                 cycle
              end if
              exit
           end do
           if (k .le. 0) exit
        end do
        node_count(i) = node_count(i) - nextra
     end do
     rewind(11)
   
     print*, 'Reading boundary nodes'
     write(500,*) 'Reading boundary nodes'
   
     read(11,*) c
     do while(c/='BOUND')
        read(11,*) c
     end do
     do i=1, nlines
        read(11,*) ii, bound_all(i,1:node_count(i))
     end do
     rewind(11)
   
     do j=1, nboun
        node_per_bound(j) = node_count(id_bound(j))
        node_bound(j,1:node_per_bound(j)) = bound_all(id_bound(j),1:node_count(id_bound(j)))
        if (node_per_bound(j) == 4) then
           hex08 = .true.
        else if (node_per_bound(j) == 3) then
           pen06 = .true.
        end if
        do k=1, node_per_bound(j)
           coord_bound(j,k,1:3) = coord_all(node_bound(j,k),1:3)
        end do
     end do
   
     !
     !  Renumber nodes, remove duplicates
     !
   
     print*, 'Renumbering nodes'
     write(500,*) 'Renumbering nodes'
   
     duplicate = .false.
     npoin_inlet = 0
     do j=1, node_per_bound(1)
        npoin_inlet = npoin_inlet + 1
        node_new(1,j) = npoin_inlet
        node_comp(npoin_inlet) = node_bound(1,j)
        node_unique(npoin_inlet) = node_new(1,j)
        coord_unique(npoin_inlet,1:3) = coord_bound(1,npoin_inlet,1:3)
     end do
   
     do i=2, nboun
        do j=1, node_per_bound(i)
           k = 1
           duplicate = .false.
           do while(k.le.npoin_inlet)
              if (node_bound(i,j)==node_comp(k)) then
                 duplicate = .true.
                 node_new(i,j) = k
              end if
              k = k + 1
              if (duplicate) k = npoin_inlet + 1
           end do
           if (.not.duplicate) then
              npoin_inlet = npoin_inlet + 1
              node_new(i,j) = npoin_inlet
              node_comp(npoin_inlet) = node_bound(i,j) 
              node_unique(npoin_inlet) = node_new(i,j)
              coord_unique(npoin_inlet,1:3) = coord_bound(i,j,1:3)
           end if 
        end do
     end do
   
     do i=1, npoin_inlet
        write(17,20) node_unique(i), coord_unique(i,1:3)
        write(18,*) node_unique(i), node_comp(i)
     end do
   
     close(17)
     close(18)
   
     print*, 'Number of inlet nodes:', npoin_inlet
     write(500,*) 'Number of inlet nodes:', npoin_inlet
   
     npoin = npoin_inlet * (isteps+1)          ! total number of nodes
     nelem = nboun * isteps
   
     allocate(nnode(npoin_inlet,isteps+1))
     allocate(nnodb(nboun,4,isteps+1))
     allocate(coord(npoin_inlet,isteps+1,3))
   
     !
     !  Coordinate transformation
     !
   
     print*, 'Rotating the inlet plane for extrusion along the Z axis'
     write(500,*) 'Rotating the inlet plane for extrusion along the Z axis'
      
     zvec(1:2) = 0.0d0
     zvec(3) = -1.0d0
     znorm = 1.0d0
   
     !  
     !  Rotating around the Y axis
     !
     angley = 0.0d0
   
     bnorm(1) = (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,3)-coord_bound(1,2,3)) - &
                (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,2)-coord_bound(1,2,2))
     bnorm(2) = (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,1)-coord_bound(1,2,1)) - &
                (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,3)-coord_bound(1,2,3))
     bnorm(3) = (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,2)-coord_bound(1,2,2)) - &
                (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,1)-coord_bound(1,2,1))
   
     print*, 'Inlet normal vector:', bnorm(1:3)
     write(500,*) 'Inlet normal vector:', bnorm(1:3)
   
     bnorm_aux(1) = bnorm(1)
     bnorm_aux(2) = 0.0d0
     bnorm_aux(3) = bnorm(3)
   
     if (abs(bnorm_aux(1)) .gt. 1e-12) then 
        dot_prod = bnorm_aux(1) * zvec(1) + bnorm_aux(2) * zvec(2) + bnorm_aux(3) * zvec(3)
        norm_aux = sqrt(bnorm_aux(1)**2.0d0+bnorm_aux(2)**2.0d0+bnorm_aux(3)**2.0d0)
   
        cosine = dot_prod / (znorm * norm_aux)
   
        if (bnorm_aux(1) .le. 0.0d0) then
           angley = acos(cosine)
        else
           angley = twopi - acos(cosine)
        end if
   
        print*, '1st angle is:', angley*180.0d0/pi
        write(500,*) '1st angle is:', angley*180.0d0/pi
   
        do i=1, nboun
           do j=1, node_per_bound(i)
              temp1 = coord_bound(i,j,1)
              temp2 = coord_bound(i,j,3)
              coord_bound(i,j,1) = temp1 * cos(angley) - temp2 * sin(angley)
              coord_bound(i,j,3) = temp1 * sin(angley) + temp2 * cos(angley)
           end do
        end do
   
        do i=1, npoin_inlet
           temp1 = coord_unique(i,1)
           temp2 = coord_unique(i,3)
           coord_unique(i,1) = temp1 * cos(angley) - temp2 * sin(angley)
           coord_unique(i,3) = temp1 * sin(angley) + temp2 * cos(angley)
        end do
     end if
    
     !  
     !  Rotating around the X axis
     !
     anglex = 0.0d0
   
     bnorm(1) = (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,3)-coord_bound(1,2,3)) - &
                (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,2)-coord_bound(1,2,2))
     bnorm(2) = (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,1)-coord_bound(1,2,1)) - &
                (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,3)-coord_bound(1,2,3))
     bnorm(3) = (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,2)-coord_bound(1,2,2)) - &
                (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,1)-coord_bound(1,2,1))
   
     print*, 'Inlet normal vector:', bnorm(1:3)
     write(500,*) 'Inlet normal vector:', bnorm(1:3)
   
     bnorm_aux(1) = 0.0d0
     bnorm_aux(2) = bnorm(2)
     bnorm_aux(3) = bnorm(3)
   
     if (abs(bnorm_aux(2)) .gt. 1e-12) then 
        dot_prod = bnorm_aux(1) * zvec(1) + bnorm_aux(2) * zvec(2) + bnorm_aux(3) * zvec(3)
        norm_aux = sqrt(bnorm_aux(1)**2.0d0 + bnorm_aux(2)**2.0d0 + bnorm_aux(3)**2.0d0)
   
        cosine = dot_prod / (znorm * norm_aux)
     
        if (bnorm_aux(2) .le. 0.0d0) then
           anglex = acos(cosine)
        else
           anglex = twopi - acos(cosine)
        end if
   
        print*, '2nd angle is:', anglex*180.0d0/pi
        write(500,*) '2nd angle is:', anglex*180.0d0/pi
   
        do i=1, nboun
           do j=1, node_per_bound(i)
              temp1 = coord_bound(i,j,2)
              temp2 = coord_bound(i,j,3)
              coord_bound(i,j,2) =   temp1 * cos(anglex) - temp2 * sin(anglex)
              coord_bound(i,j,3) =   temp1 * sin(anglex) + temp2 * cos(anglex)
           end do
        end do
   
        do i=1, npoin_inlet
           temp1 = coord_unique(i,2)
           temp2 = coord_unique(i,3)
           coord_unique(i,2) = temp1 * cos(anglex) - temp2 * sin(anglex)
           coord_unique(i,3) = temp1 * sin(anglex) + temp2 * cos(anglex)
        end do
     end if
    
     bnorm(1) = (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,3)-coord_bound(1,2,3)) - &
                (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,2)-coord_bound(1,2,2))
     bnorm(2) = (coord_bound(1,2,3)-coord_bound(1,1,3))*(coord_bound(1,3,1)-coord_bound(1,2,1)) - &
                (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,3)-coord_bound(1,2,3))
     bnorm(3) = (coord_bound(1,2,1)-coord_bound(1,1,1))*(coord_bound(1,3,2)-coord_bound(1,2,2)) - &
                (coord_bound(1,2,2)-coord_bound(1,1,2))*(coord_bound(1,3,1)-coord_bound(1,2,1))
   
     print*, 'Inlet normal vector:', bnorm(1:3)
     write(500,*) 'Inlet normal vector:', bnorm(1:3)
   
     !
     !  Calculate center of gravity of an element
     !  for correcting the boundary orientation
     !
   
     cog(1:3) = 0.0d0
     do i= 1, node_per_bound(1)
        cog(1:2) = cog(1:2) + coord_bound(1,i,1:2)
     end do
     cog(1:2) = cog(1:2)/dfloat(node_per_bound(1))
     cog(3) = coord_bound(1,1,3)+timestep/2d0
      
     print*, 'Center of gravity coordinates:', cog(1:3)
     write(500,*) 'Center of gravity coordinates:', cog(1:3)
    
     ! 
     !  Write element connectivity
     !
     
     print*, 'Writing geo file'
     write(500,*) 'Writing geo file'
   
     print*, 'Writing elements'
     write(500,*) 'Writing elements'
   
     do i = 1, isteps+1
        do j = 1, nboun
           nnodb(j,1:node_per_bound(j),i) = node_new(j,1:node_per_bound(j)) + npoin_inlet*(i-1)
        end do
     end do
    
     write(14,*) 'TYPES'
     ielem = 0
     do i = 1, isteps
        do j = 1, nboun
           if (node_per_bound(j) == 4) then
              ielem = ielem + 1
              write(14,*) ielem, '37'
           else
              ielem = ielem + 1
              write(14,*) ielem, '34'
           end if
        end do
     end do
     write(14,*) 'END_TYPES'

     write(14,*) 'ELEMENTS'
     ielem = 0
     do i = 1, isteps
        do j = 1, nboun
           znorm = (coord_bound(j,2,1)-coord_bound(j,1,1))*(coord_bound(j,3,2)-coord_bound(j,2,2)) - &
                   (coord_bound(j,2,2)-coord_bound(j,1,2))*(coord_bound(j,3,1)-coord_bound(j,2,1))
           if (znorm .gt. 0.0d0) then
              if (node_per_bound(j) == 4) then
                 ielem = ielem + 1
                 write(14,'(30(i10,1x))') ielem, nnodb(j,1,i), nnodb(j,2,i), nnodb(j,3,i),  nnodb(j,4,i), &
                                    nnodb(j,1,i+1), nnodb(j,2,i+1), nnodb(j,3,i+1), nnodb(j,4,i+1)
              else
                 ielem = ielem + 1
                 write(14,'(30(i10,1x))') ielem, nnodb(j,1,i), nnodb(j,2,i), nnodb(j,3,i), &
                                    nnodb(j,1,i+1), nnodb(j,2,i+1), nnodb(j,3,i+1)
              end if
           else
              if (node_per_bound(j) == 4) then
                 ielem = ielem + 1
                 write(14,'(30(i10,1x))') ielem, nnodb(j,1,i+1), nnodb(j,2,i+1), nnodb(j,3,i+1),  nnodb(j,4,i+1), &
                                    nnodb(j,1,i), nnodb(j,2,i), nnodb(j,3,i), nnodb(j,4,i)
              else
                 ielem = ielem + 1
                 write(14,'(30(i10,1x))') ielem, nnodb(j,1,i+1), nnodb(j,2,i+1), nnodb(j,3,i+1), &
                                    nnodb(j,1,i), nnodb(j,2,i), nnodb(j,3,i)
              end if
           end if
        end do
     end do
     write(14,*) 'END_ELEMENTS'
   
     !
     !  Define boundary intersections and nodes
     !
   
   
     print*, 'Writing boundaries'
     write(500,*) 'Writing boundaries'
     write(14,*) 'BOUNDARIES'
   
     allocate(nbcod(npoin_inlet))            ! node boundary code
     allocate(new_wall(nboun,4,nset))
     allocate(nboub(nboun,nset))
   
     nbcod(1:npoin_inlet) = 0
     iboundary = 0
   
     do i=1, nset
        nbset(i) = 0
        do j=1, nlines
           read(12,*) iboun, iset
           if (iset == isetb(i)) then
              nbset(i) = nbset(i) + 1
           end if
        end do
        print*, 'Number of boundaries at set', isetb(i), ':', nbset(i)
        write(500,*) 'Number of boundaries at set', isetb(i), ':', nbset(i)
        rewind(12)
     end do
   
     do i=1,nset
        node_belongs = .false.
        new_wall = .false.
        k=0
        allocate(ibous(nbset(i)))
        allocate(node_bset(nbset(i)))
        allocate(node_wall(nbset(i),4))
        do j=1, nlines
           read(12,*) iboun, iset
           if (iset == isetb(i)) then
              k=k+1
              ibous(k) = iboun
           end if
        end do
        rewind(12)
   
        do j=1, nbset(i)
           node_bset(j) = node_count(ibous(j))

!          node_wall(j,1:node_per_bound(j)) = bound_all(ibous(j),1:node_count(ibous(j)))  ! esto esta mal:
                                                                                          ! node_per_bound(j) no tiene sentido es de tamano nboun no nbset(i)
                                                                                          ! node_count(ibous(j)) esta bien  pero ya que tenemos definido 
                                                                                          ! node_bset(j)   mejor usarlo directamente
           node_wall(j,1:node_bset(j)) = bound_all(ibous(j),1:node_bset(j))
           do k=1, node_bset(j)
              node_belongs(node_wall(j,k)) = .true.
           end do
        end do
   
        do m=1, nboun
           do n=1, node_per_bound(m)
              if (node_belongs(node_bound(m,n))) then
                  new_wall(m,n,i) = .true.
              end if
           end do
        end do
   
        deallocate(ibous)
        deallocate(node_bset)
        deallocate(node_wall)
   
        !
        ! Define boundaries and correct orientation
        !
   
        do j=1, nboun
           l = 0
           do k=1, node_per_bound(j)
              if (new_wall(j,k,i)) then
                 l = l + 1
              end if
           end do
           if (l .ge. 2) then
              icoun = 0
              do k=1, node_per_bound(j)
                 if (new_wall(j,k,i)) then
                    if (icoun == 0) then
                       ibound(1) = nnodb(j,k,1)                      ! the 1 means plane 1 -inlet face
                       bcoor(1,1:3) = coord_bound(j,k,1:3)
                       ibound(2) = nnodb(j,k,2)
                       bcoor(2,1:2) = coord_bound(j,k,1:2)           ! 2 second plane
                       bcoor(2,3) = coord_bound(j,k,3)+timestep
                       icoun = 1
                    else
                       ibound(3) = nnodb(j,k,2)
                       bcoor(3,1:2) = coord_bound(j,k,1:2)
                       bcoor(3,3) = coord_bound(j,k,3)+timestep
                       ibound(4) = nnodb(j,k,1)
                       bcoor(4,1:3) = coord_bound(j,k,1:3)
                    end if
                 end if
              end do
              bnorm(1) = (bcoor(2,2)-bcoor(1,2))*(bcoor(3,3)-bcoor(2,3)) - (bcoor(2,3)-bcoor(1,3))*(bcoor(3,2)-bcoor(2,2))
              bnorm(2) = (bcoor(2,3)-bcoor(1,3))*(bcoor(3,1)-bcoor(2,1)) - (bcoor(2,1)-bcoor(1,1))*(bcoor(3,3)-bcoor(2,3))
              bnorm(3) = (bcoor(2,1)-bcoor(1,1))*(bcoor(3,2)-bcoor(2,2)) - (bcoor(2,2)-bcoor(1,2))*(bcoor(3,1)-bcoor(2,1))
              bcog(1:3) = bcoor(1,1:3) - cog(1:3)
   
              produ = bcog(1)*bnorm(1)+bcog(2)*bnorm(2)+bcog(3)*bnorm(3)
    
              do n=1, isteps
                 iboundary = iboundary + 1
                 if (produ .ge. 0.0d0) then
                    write(14,'(30(i10,1x))') iboundary, ibound(1)+npoin_inlet*(n-1), ibound(2)+npoin_inlet*(n-1), &
                                       ibound(3)+npoin_inlet*(n-1), ibound(4)+npoin_inlet*(n-1)
                    write(15,*) iboundary, isetb(i)
                 else         
                    write(14,'(30(i10,1x))') iboundary, ibound(4)+npoin_inlet*(n-1), ibound(3)+npoin_inlet*(n-1), &
                                       ibound(2)+npoin_inlet*(n-1), ibound(1)+npoin_inlet*(n-1)
                    write(15,*) iboundary, isetb(i)
                 end if
              end do
           end if
        end do
     end do
   
     write(14,*) 'END_BOUNDARIES'
   
     !
     !  Write node coordinates
     !
   
     print*, 'Writing coordinates'
     write(500,*) 'Writing coordinates'
     write(14,*) 'COORDINATES'
     do i = 1, isteps+1
        do j = 1, npoin_inlet
             nnode(j,i) = node_unique(j) + npoin_inlet*(i-1)
             coord(j,i,1:2) = coord_unique(j,1:2)
             coord(j,i,3) = coord_unique(j,3) + timestep*(i-1)
             write(14,20) nnode(j,i), coord(j,i,1:3)
        end do
     end do
     write(14,*) 'END_COORDINATES'
   
     !
     !  Periodicity
     !
   
     print*, 'Writing periodicity file'
     write(500,*) 'Writing periodicity file'
   
     nperi = 0
     if (periodicity) then
        l = 0
        do k=1,nper
           read(13,*) iper1, iper2
           do i=1, npoin_inlet
              if (node_comp(i) == iper1) then
                 l = l + 1
                 node_per(l,1) = iper1
                 node_per(l,2) = iper2
              end if
           end do
        end do
   
        allocate(new_per(l,2))
   
        do i=1, l
           do j=1, nboun
              do k=1 ,node_per_bound(j)
                 if (node_per(i,1) == node_bound(j,k)) then
                    new_per(i,1) = node_new(j,k)
                 end if
                 if (node_per(i,2) == node_bound(j,k)) then
                    new_per(i,2) = node_new(j,k)
                 end if
              end do
           end do
        end do
   
        do i=1, l
           write(16,*) new_per(i,1), new_per(i,2)
           write(16,*) new_per(i,1), new_per(i,1)+npoin_inlet*isteps
           write(16,*) new_per(i,1), new_per(i,2)+npoin_inlet*isteps
           nperi = nperi + 3
        end do
   
        do i=1, npoin_inlet
           node_is_periodic = .false.
           do j=1, l
              do k=1, 2
                 if (nnode(i,1) == new_per(j,k)) then
                    node_is_periodic = .true.
                 end if
              end do
           end do
           if (.not.node_is_periodic) then
              write(16,*) nnode(i,1), nnode(i,isteps+1)
              nperi = nperi + 1
           end if
        end do
     
        do i=2, isteps
           do j=1, l
              write(16,*) new_per(j,1)+npoin_inlet*(i-1), new_per(j,2)+npoin_inlet*(i-1)
              nperi = nperi + 1
           end do
        end do
   
     else
        do i=1, npoin_inlet
           write(16,*) nnode(i,1), nnode(i,isteps+1)
           nperi = nperi + 1
        end do
     end if
   
     ! 
     !  Writing dom.dat file
     !
   
     print*, 'Writing dom.dat file'
     write(500,*) 'Writing dom.dat file'
   
     if (hex08 .and. pen06) then
        types = 'PEN06, HEX08'
     else if (hex08) then
        types = 'HEX08'
     else if (pen06) then
        types = 'PEN06'
     end if
   
     npoin = npoin_inlet*(isteps+1)
     nelem = ielem
     nboun = iboundary
     write(19,100) npoin, &
           nelem, &
           types, &
           nboun, &
           nperi, &
           './inlet.geo.dat', &
           './inlet.per', &
           './inlet.fix'
         !
   100   format( &
          '$------------------------------------------------------------',/, &
          'DIMENSIONS',/,&
          '  NODAL_POINTS  ',i8,/,&
          '  ELEMENTS      ',i8,/,&
          '  SPACE_DIMENSIONS   3',/,&
          '  TYPES_OF_ELEMENTS ',a,/,&
          '  BOUNDARIES    ' ,i8,/,&
          '  PERIODIC_NODES=',i8,/,&
          'END_DIMENSIONS',/,&
          '$------------------------------------------------------------',/,&
          'STRATEGY',/,&
          '  DOMAIN_INTEGRATION_POINTS  0',/,&
          '  PERIODICITY:    MATRIX',/,&
          'END_STRATEGY',/,&
          '$-------------------------------------------------------------',/,&
          'GEOMETRY',/,&
          '  GROUPS = 500',/,&
          '  INCLUDE ',a,/,&
          '  PERIODIC ',/,&
          '    INCLUDE ',a,/,&
          '  END_PERIO ',/,&
          'END_GEOMETRY',/,&
          '$-------------------------------------------------------------',/,&
          'SETS',/,&
          'END_SETS',/,&
          '$-------------------------------------------------------------',/,&
          'BOUNDARY_CONDITIONS, EXTRAPOLATE',/,&
          ' ON_BOUNDARIES',/,&
          '  INCLUDE ',a,/,&
          ' END_ON_BOUNDARIES',/,&
          'END_BOUNDARY_CONDITIONS',/,&
          '$-------------------------------------------------------------')
   
   
     write(20,*) '#npoin, Anglex, Angley'
     write(20,*) npoin, anglex, angley
   
     call cpu_time(stop_time)
     write(*,"(A,f9.2,x,A)") "Elapsed time:", stop_time - start_time, 'seconds'
     write(500,"(A,f9.2,x,A)") "Elapsed time:", stop_time - start_time, 'seconds'
   
     close(14)
     close(15)
     close(16)
     close(19)
     close(20)
     close(500)

  else ! velocity field rotation to match the original mesh

     call getarg(1, fil_inp)
     ndime = 3

     open (unit=201,file='angle.dat')
     open (unit=202,file=fil_inp)
     open (unit=203,file='veloc.dat')

     read(201,*)
     read(201,*) npoin, anglex, angley

     allocate(veloc(npoin,ndime))

     do i = 1, npoin
        read(202,*) inode, veloc(i,1:3)
     end do

     anglex = twopi - anglex
     angley = twopi - angley

     !  
     !  Rotating around the X axis
     !

     do i=1, npoin
        temp1 = veloc(i,2)
        temp2 = veloc(i,3)
        veloc(i,2) = temp1 * cos(anglex) - temp2 * sin(anglex)
        veloc(i,3) = temp1 * sin(anglex) + temp2 * cos(anglex)
     end do

    !
    !  Rotating around the Y axis
    !

    do i=1, npoin
       temp1 = veloc(i,1)
       temp2 = veloc(i,3)
       veloc(i,1) = temp1 * cos(angley) - temp2 * sin(angley)
       veloc(i,3) = temp1 * sin(angley) + temp2 * cos(angley)
    end do

    do i=1, npoin
       write(203,20) i, veloc(i,1:3)
    end do

    close(201)
    close(202)
    close(203)

  end if

  close(101)

20 format(i10,3f22.10)

end program periodic
