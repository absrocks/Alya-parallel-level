      program Block3d
!-----------------------------------------------------------------------
!
!  This programs constructs the .fix and .geo file for Fantom for the
!  domain  V = [0,a] x [0,b] x [0,c] using P1, P2, Q1 and Q2 elements
!      
!-----------------------------------------------------------------------

      USE IFPORT   
      implicit none
      real*8,    parameter :: pi=3.141592653589793238462643383279502884197_8
      real*8,    parameter :: twopi=6.283185307179586476925286766559005768394_8
      integer &
           lugeo,lufix,luper, lumsh,nex,ney,nez,nnod,npx,npy,npz, &
           np,i,j,k, auxin, ipoin, tflag, funve, funte, luint, &
           code, codin, codup, codlo, ludom, npoin, nperi, iboun, & 
           ie, je, ke, no1, no2, no3, no4, nelem, nboun, codle, codri, &
           lugrp
      real*8 &
           a,b,c,x,y,z,dx,dy,dz,vx,vy,vz, trati,delte, gamma, force, yy, alpha,ynn, &
           invis, delta, ranam, tempe, veloc, ranve, ranmo, mu, theta, ranvx, ranvy, ranvz, &
           tecol, tehot
      character*20 &
           file_geo, file_fix, file_msh, label, file, file_per, file_grp
      logical &
           wallx,wally,wallz

      external ipoin
!
!*** Initializations
!      
      lugeo = 8                 ! geo.dat ! geometry
      lufix = 9                 ! fix.dat ! boundary conditions
      luper = 92                ! .per ! periodicity
      lugrp = 21

      lumsh = 13                ! .msh  gid format
      funve = 15
      funte = 17
      luint = 19
     
!!$c
!!$c*** This program is tuned to generate the mesh, initial and boundary conditions for a turbulent channel [2,2*pi,2*pi/3]]
!!$c*** Initial conditions are generated for a laminar channel with 30% random fluctuations in velocity and homogeneous
!!$c*** temperature
!!$c    
      force = 1.0d0;
      invis = 5.0d0

      CALL SRAND(123)  ! rand numbers seed
!!$c
!!$c*** gamma, hyperbolic gamma factor for mesh in y direction
!!$c
      gamma = 2.50d0
      alpha = 0.9d0
      ynn = dtanh(gamma*(2.0d0*alpha-1.0))/dtanh(gamma)
      
!!$c
!!$c*** Temperature ratios
!!$c      
      
      nex = 1001
      tflag = -1      
      do while (nex.gt.1000)
         write(*,'(a,$)') '>>> Elements in the x,y,z-direction : '
         read (*,*) nex
      end do      
      do while (tflag.lt.1.or.tflag.gt.3)
         write(*,'(a,$)') '>>> Temperature ratio (Tr) :Press 1: Tr=1.01',' Press 2: Tr = 2.00 ','Press 3: Tr = 4.00' ,'Press 4: Tr= 8.00 :'
         read (*,*) tflag 
      end do
         write(*,'(a,$)') '>>> Random amplitud of velocity: (%) '
         read (*,*) ranam
	ranam = ranam*0.01d0
      ! temperature ratio
      if (tflag == 1) trati = 1.01d0
      if (tflag == 2) trati = 2.0d0
      if (tflag == 3) trati = 4.0d0
      if (tflag == 4) trati = 8.0d0
      ! upper and lower temperature 

      delte = (trati-1.0d0)/(trati+1.0d0)      
      tehot = trati
      tecol = 1.0d0
      
      write(label,'(i3)') nex

      file = 'chan'
      file = trim(file)//trim(adjustl(label))//'T'
      write(label,'(i1)') tflag

      file = trim(file)//trim(adjustl(label))
      file_geo = trim(file)//'.geo.dat' 
      open(lugeo,file=file_geo,status='unknown')

      file_fix = trim(file)//'.fix' 
      open(lufix,file=file_fix,status='unknown')

      file_grp = trim(file)//'.grp' 
      open(lugrp,file=file_grp,status='unknown')

      file_msh = trim(file)//'.post.msh'          
      open(lumsh,file=file_msh,status='unknown')

      file_per = trim(file)//'.per'          
      open(luper,file=file_per,status='unknown')

      file_per=  trim(file)// '.vel.ini' 
      open(funve,file=file_per,status='unknown')

      file_per=  trim(file)// '.tem.ini' 
      open(funte,file=file_per,status='unknown')

      file_per=  trim(file)// '.cod_inter' 
      open(luint,file=file_per,status='unknown')

      file_per=  trim(file)// '.dom.dat' 
      open(ludom,file=file_per,status='unknown')

      write(*,'(a,$)') '>>> Length in the x-direction : '
!      read (*,*) a
      a = 4.0d0*pi
      write(*,*) 'a=', a
      write(*,'(a,$)') '>>> Length in the y-direction : '
!      read (*,*) b
      b = 2.0d0
      write(*,*) 'b=', b
      write(*,'(a,$)') '>>> Length in the z-direction : '
!      read (*,*) c
      c = (4.0d0/3.0d0)*pi
      write(*,*) 'c=', c

      ! elements in x and z directions
      ney=nex
      nez=nex      

      nnod=8 ! nodes per element


      ! nodes in x and z directions
      npx = nex + 1
      npy = ney + 1
      npz = nez + 1
      

      dx = a/float(npx-1)
      dy = b/float(npy-1)
      dz = c/float(npz-1)
      
      write(lumsh,'(a)') 'MESH "chan32" dimension 3 ElemType Hexahedra Nnode 8'
      write(lumsh,'(a)') ' '

!!$c
!!$c***  Coordinates and fixity conditions
!!$c  
    
      write(lugeo,'(a)') 'COORDINATES'
      write(lumsh,'(a)') 'COORDINATES'

            
      codin = 3 ! code of interior node
      codlo = 1 ! code of lower wall
      codup = 2 ! code of upper wall
      codle = 4 ! code of left wall
      codri = 5 ! code of right wall
      vx = 0.0
      vy = 0.0
      vz = 0.0
      np =  0
      x  = -dx
      do i = 1,npx
        wallx = .false.
        if((i.eq.1).or.(i.eq.npx)) wallx = .true.
        x = x + dx
        y = -dy
        do j = 1,npy
          tempe = 1.0d0
          wally = .false.
          code = codin
!          if(i.eq.1) code = codle
!          if(i.eq.npx) code = codri
          if(j.eq.npy) then 
             tempe =tehot  !upper wall (hot)
             code = codup
          else  if(j.eq.1) then
             tempe =tecol  !lower wall (cold)
             code = codlo
          end if
          if((j.eq.1).or.(j.eq.npy)) wally = .true.
!------     hyperbolic z distribution
          if (trati.lt.1.1) then   
            y = dtanh(gamma*(2.0d0*dfloat(j-1)/dfloat(ney)-1.0))/dtanh(gamma)

! nonsymmetric hyperbolic distribution as in Nicoud for greaters temp ratio
          else 
            yy = dtanh(gamma*(2.0d0*alpha*dfloat(j-1)/dfloat(ney)-1.0))/dtanh(gamma)
            y = 2.0d0*(yy+1.0)/(ynn+1.0)-1.0          
          end if
          tempe = 0.5d0*(tehot + tecol) + 0.5d0*(tehot - tecol)* y
          z = -dz
          do k = 1,npz
            wallz = .false.
            if((k.eq.1).or.(k.eq.npz)) wallz = .true.
            z = z + dz
            np = np + 1
            write(lugeo,20) np, x, y, z
            write(lumsh,20) np, x, y, z

!---------- To be modified according to the test case            
            if(wally) then !
              write(funve,'(i8,x,3(f17.10,x))') np, 0.0d0, 0.0d0, 0.0d0
              write(funte,'(i8,x,1(f17.10,x))') np, tempe                  
            else                             
!              veloc = 0.5d0*invis*(y+1.0d0)*(1.0d0-y) ! *force          
!              veloc = 0.5d0*invis*(1.0d0-y**(10.0) ) ! *force            
!              veloc = veloc*0.25d0 ! maximu 25% of laminar solution
              veloc = 0.5d0*invis*(1.0d0-y*y ) ! *force           	    
              ranve = ranam*veloc ! random velocity amplitud
              ranmo = ranam*RAND(0)*veloc
              mu = 2.0d0*(RAND(0)-0.5) ! -1 < mu < +1
              theta = twopi*RAND(0)
              ranvx = ranmo*sqrt(1.0d0-mu*mu)*cos(theta)
              ranvy = ranvx*tan(theta)
              ranvz = ranmo*mu
              ! linear temperature distribution, no! cause we want to start w incompressible
!              tempe = 0.5d0*((trati-1.0d0)*y+trati+1.0d0) 

              write(funve,'(i8,x,3(f17.10,x))') np, veloc+ ranvx,ranvy,ranvz
              write(funte,'(i8,x,1(f17.10,x))') np, tempe
!              if (.not.wallx) &
                   write(luint,'(2i8)') np, code ! writes all internal codes 
            end if
!----------------------------------------------------            
          end do
        end do
      end do
      write(lugeo,'(a)') 'END_COORDINATES'
      write(lumsh,'(a)') 'END COORDINATES'


!!$      
!!$c
!!$c***  Nodal connections
!!$c
      write(lugeo,'(a8)') 'ELEMENTS'
      write(lumsh,'(a8)') 'ELEMENTS'
!      if(nnod.eq. 4) then
!        call elemp1(nex,ney,nez,npy,npz,lugeo)
!      else if(nnod.eq. 8) then
        call elemq1(nex,ney,nez,npy,npz,lugeo,lumsh)
!      else if(nnod.eq.10) then
!        call elemp2(nex,ney,nez,npy,npz,lugeo)
!      else if(nnod.eq.27) then
!        call elemq2(nex,ney,nez,npy,npz,lugeo)
!      else
!        write(*,'(a)') 'Wrong number of nodes'
!        stop
!      end if
      write(lugeo,'(a12)') 'END_ELEMENTS'
      write(lumsh,'(a12)') 'END ELEMENTS'

      write(lugeo,'(a12)') 'BOUNDAR'
      iboun = 0
      do ie = 1, nex
         i = ie
         do ke =1, nez             
            k = ke 
            iboun = iboun +1
            no1 = ipoin(i  ,1  ,k   ,npy,npz)
            no2 = ipoin(i+1,1  ,k   ,npy,npz)
            no3 = ipoin(i+1,1  ,k +1,npy,npz)
            no4 = ipoin(i  ,1  ,k +1,npy,npz)
            write(lugeo, '(5(i8,x))') iboun, no1, no2, no3, no4
            write(lufix, '(2i8)') iboun, codlo
            write(lugrp, '(2i8)') iboun, codlo
         end do
      end do
      do ie = 1, nex
         i = ie
         do ke =1, nez             
            k = ke 
            iboun = iboun +1
            no1 = ipoin(i  ,npy  ,k   ,npy,npz)
            no2 = ipoin(i+1,npy  ,k   ,npy,npz)
            no3 = ipoin(i+1,npy  ,k+1,npy,npz)
            no4 = ipoin(i  ,npy  ,k +1,npy,npz)
            write(lugeo, '(5(i8,x))') iboun, no1, no2, no3, no4
            write(lufix, '(2i8)') iboun, codup
            write(lugrp, '(2i8)') iboun, codup
         end do
      end do
      ! left and right wall (inflow and outflow)
!!$      do je = 1, ney
!!$         j = je
!!$         do ke =1, nez             
!!$            k = ke 
!!$            iboun = iboun +1
!!$            no1 = ipoin(1 ,j  ,k   ,npy,npz)
!!$            no2 = ipoin(1 ,j+1,k   ,npy,npz)
!!$            no3 = ipoin(1 ,j+1,k +1,npy,npz)
!!$            no4 = ipoin(1 ,j  ,k +1,npy,npz)
!!$            write(lugeo, '(5(i8,x))') iboun, no1, no2, no3, no4
!!$            write(lugrp, '(2i8)') iboun, codle
!!$         end do
!!$      end do
!!$      do je = 1, ney
!!$         j = je
!!$         do ke =1, nez             
!!$            k = ke 
!!$            iboun = iboun +1
!!$            no1 = ipoin(npx ,j  ,k   ,npy,npz)
!!$            no2 = ipoin(npx ,j+1,k   ,npy,npz)
!!$            no3 = ipoin(npx ,j+1,k +1,npy,npz)
!!$            no4 = ipoin(npx ,j  ,k +1,npy,npz)
!!$            write(lugeo, '(5(i8,x))') iboun, no1, no2, no3, no4
!!$            write(lugrp, '(2i8)') iboun, codri
!!$         end do
!!$      end do
!!$
     
      write(lugeo,'(a12)') 'END_BOUNDARIES'

!c      
!c***  PERIODIC BOUNDARY CONDITIONS
!c
!      write(lugeo,'(a)') 'SLAVE_NODES'
      
!!$c
!!$c*** Periodic in z direction
!!$c      

      nperi = 0       ! # periodic nodes
      do i = 2, npx-1
        do j =1, npy
        auxin = (i-1)*npy*npz+(j-1)*npz
        write(luper,'(2(i10))') auxin + npz,  auxin +1
        nperi = nperi +1
        end do
      end do
!!$c
!!$c*** Periodic in x direction
!!$c      
      do j = 1, npy
        do k = 2, npz-1
        auxin = (j-1)*npz +k
        write(luper,'(2(i10))') auxin + npy*npz*(npx-1),  auxin
        nperi = nperi +1
        end do
      end do
!!$c
!!$c*** Periodic points with more than one master, only one can a be master, all the others, slaves
!!$c*** This happen to all points in the 4 edges in y direction in the inlet and output(x=z=0&x=a,z=b&x=0,z=b&x=a,z=0) 
!!$c*** For j = 1 (lower wall)  

!      write(lugeo,'(2(i10))') ipoin(npx,  1,  1, npy,npz), ipoin(1,1,1,npy,npz)
!      write(lugeo,'(2(i10))') ipoin(npx,  1,npz, npy,npz), ipoin(1,1,1,npy,npz)
!      write(lugeo,'(2(i10))') ipoin(1  ,  1,npz, npy,npz), ipoin(1,1,1,npy,npz)
      do j =1, npy
        write(luper,'(2(i10))') ipoin(1,j,1,npy,npz),  ipoin(npx,  j,  1, npy,npz)
        write(luper,'(2(i10))') ipoin(1,j,1,npy,npz),  ipoin(npx,  j,npz, npy,npz)
        write(luper,'(2(i10))') ipoin(1,j,1,npy,npz),  ipoin(1  ,  j,npz, npy,npz)
        nperi = nperi +3
      end do
      
!c*** For j = npy (upper wall)      
!      write(lugeo,'(2(i10))') ipoin(npx,npy,  1, npy,npz), ipoin(1,npy,1,npy,npz)
!      write(lugeo,'(2(i10))') ipoin(npx,npy,npz, npy,npz), ipoin(1,npy,1,npy,npz)
!      write(lugeo,'(2(i10))') ipoin(1  ,npy,npz, npy,npz), ipoin(1,npy,1,npy,npz)
         
      
!      write(lugeo,'(a)') 'END_SLAVE_NODES'

      ! WRITE DOM.DAT
      npoin = npx* npy* npz
      nelem = nex* ney* nez
      nboun = iboun
      write(ludom,100) npoin, &
           nelem, &       
           nboun, & 
           nperi, & 
           './'//TRIM(file_geo), &
           './'//TRIM(file)//'.per', &
           './'//TRIM(file_grp), &
           './'//TRIM(file_fix), &
           './'//TRIM(file)//'.cod_inter', &
           './'//trim(file)//'.vel.ini' ,&
           './'//trim(file)//'.tem.ini'          
      !
100   format( &
       '$------------------------------------------------------------',/, &
       'DIMENSIONS',/,&
       '  NODAL_POINTS  ',i8,/,&
       '  ELEMENTS      ',i8,/,&
       '  SPACE_DIMENSIONS   3',/,&
       '  TYPES_OF_ELEMENTS HEX08',/,&
       '  BOUNDARIES    ' ,i8,/,&
       '  PERIODIC_NODES=',i8,/,&
       'END_DIMENSIONS',/,&
       '$------------------------------------------------------------',/,&
       'STRATEGY',/,&
       '  INTEGRATION_RULE           Open',/,&  
       '  DOMAIN_INTEGRATION_POINTS  0',/,&
       '  PERIODICITY:    MATRIX',/,&
       'END_STRATEGY',/,&
       '$-------------------------------------------------------------',/,&
       'GEOMETRY',/,&
       '  GROUPS = 500',/,&
       '  TYPES, ALL =HEX08',/,&
       '  END_TYPES',/,&
       '  INCLUDE ',a,/,&
       '  PERIODIC ',/,&
       '    INCLUDE ',a,/,&
       '  END_PERIO ',/,&
       'END_GEOMETRY',/,&  
       '$-------------------------------------------------------------',/,&
       'SETS',/,&
       ' BOUNDARIES',/,&
       '  INCLUDE ',a,/,&
       ' END_BOUNDARIES',/,&
       'END_SETS',/,&
       '$-------------------------------------------------------------',/,&
       'BOUNDARY_CONDITIONS, EXTRAPOLATE',/,&
       ' ON_BOUNDARIES',/,&
       '  INCLUDE ',a,/,&
       ' END_ON_BOUNDARIES',/,&         
       ' ON_NODES',/,&
       '   INCLUDE ',a,/,&
       ' END_ON_NODES',/,&
       '  VALUES, FUNCTION=1, DIMENSION=3',/, &
       '     INCLUDE ',a,/,&
       '  END_VALUES',/,&
       '  VALUES, FUNCTION=2, DIMENSION=1',/,&
       '     INCLUDE ',a,/,&
       '  END_VALUES',/,&      
       'END_BOUNDARY_CONDITIONS',/,&
       '$-------------------------------------------------------------')
      
      close(luper)
      close(lufix)
      close(lumsh)
      close(funve)
      close(funte)
      close(luint)
      close(lugeo)
      close(ludom)

   20 format(i10,3f17.10)
   30 format(i8,2x,a3,3(1x,f17.10),2(i4))
   40 format(i8,2x,a1,1x,f17.10,i4)   

      end program
!c
!c***  Function ipoin
!c      
      integer function ipoin(i,j,k,ny,nz)
      integer i, j, k, ny, nz
      ipoin = (i-1)*ny*nz + (j-1)*nz + k
      end
!c
!c***  P1 elements
!c
      subroutine elemp1(nex,ney,nez,npy,npz,lugeo)
      implicit none
      integer   nex,ney,nez,npy,npz,lugeo
      integer   ne,ie,je,ke,ip,jp,kp, &
           no1,no2,no3,no4,no5,no6,no7,no8, &
           ipoin
      external ipoin

      ne = 0
      do ie = 1,nex
        ip = ie
        do je = 1,ney
          jp = je
          do ke = 1,nez
            kp = ke
            no1 = ipoin(ip  ,jp  ,kp  ,npy,npz)
            no2 = ipoin(ip+1,jp  ,kp  ,npy,npz)
            no3 = ipoin(ip+1,jp+1,kp  ,npy,npz)
            no4 = ipoin(ip  ,jp+1,kp  ,npy,npz)
            no5 = ipoin(ip  ,jp  ,kp+1,npy,npz)
            no6 = ipoin(ip+1,jp  ,kp+1,npy,npz)
            no7 = ipoin(ip+1,jp+1,kp+1,npy,npz)
            no8 = ipoin(ip  ,jp+1,kp+1,npy,npz)
            ne = ne + 1
            write(lugeo,10) ne, no1, no2, no4, no5
            ne = ne + 1
            write(lugeo,10) ne, no3, no6, no7, no8
            ne = ne + 1
            write(lugeo,10) ne, no2, no3, no4, no6
            ne = ne + 1
            write(lugeo,10) ne, no3, no4, no6, no8
            ne = ne + 1
            write(lugeo,10) ne, no2, no4, no5, no6
            ne = ne + 1
            write(lugeo,10) ne, no4, no5, no6, no8
          end do
        end do
      end do
   10 format(5i6)
	
      end
!c
!c***  Q1 elements
!c
      subroutine elemq1(nex,ney,nez,npy,npz,lugeo,lumsh)
      implicit none
      integer    nex,ney,nez,npy,npz,lugeo,lumsh
      integer  ne,ie,je,ke,ip,jp,kp, &
           no1,no2,no3,no4,no5,no6,no7,no8, &
           ipoin
      external ipoin

      ne = 0
      do ie = 1,nex
        ip = ie
        do je = 1,ney
          jp = je
          do ke = 1,nez
            kp = ke
            no1 = ipoin(ip  ,jp  ,kp  ,npy,npz)
            no2 = ipoin(ip+1,jp  ,kp  ,npy,npz)
            no3 = ipoin(ip+1,jp+1,kp  ,npy,npz)
            no4 = ipoin(ip  ,jp+1,kp  ,npy,npz)
            no5 = ipoin(ip  ,jp  ,kp+1,npy,npz)
            no6 = ipoin(ip+1,jp  ,kp+1,npy,npz)
            no7 = ipoin(ip+1,jp+1,kp+1,npy,npz)
            no8 = ipoin(ip  ,jp+1,kp+1,npy,npz)
            ne = ne + 1
            write(lugeo,10) &
                 ne, no1, no2, no3, no4, no5, no6, no7, no8
            write(lumsh,10) &
                 ne, no1, no2, no3, no4, no5, no6, no7, no8     
          end do
        end do
      end do
   10 format(9i8)

      end
!c
!c***  P2 elements
!c
      subroutine elemp2(nex,ney,nez,npy,npz,lugeo)
      implicit none
      integer &
       nex,ney,nez,npy,npz,lugeo
      integer &
       ne,ie,je,ke,ip,jp,kp,no(27),&
       ipoin
      external ipoin

      ne = 0
      do ie = 1,nex
        ip = 2*ie - 1
        do je = 1,ney
          jp = 2*je - 1
          do ke = 1,nez
            kp = 2*ke - 1
            no( 1) = ipoin(ip  ,jp  ,kp  ,npy,npz)
            no( 2) = ipoin(ip+2,jp  ,kp  ,npy,npz)
            no( 3) = ipoin(ip+2,jp+2,kp  ,npy,npz)
            no( 4) = ipoin(ip  ,jp+2,kp  ,npy,npz)
            no( 5) = ipoin(ip  ,jp  ,kp+2,npy,npz)
            no( 6) = ipoin(ip+2,jp  ,kp+2,npy,npz)
            no( 7) = ipoin(ip+2,jp+2,kp+2,npy,npz)
            no( 8) = ipoin(ip  ,jp+2,kp+2,npy,npz)
            no( 9) = ipoin(ip+1,jp  ,kp  ,npy,npz)
            no(10) = ipoin(ip+2,jp+1,kp  ,npy,npz)
            no(11) = ipoin(ip+1,jp+2,kp  ,npy,npz)
            no(12) = ipoin(ip  ,jp+1,kp  ,npy,npz)
            no(13) = ipoin(ip  ,jp  ,kp+1,npy,npz)
            no(14) = ipoin(ip+2,jp  ,kp+1,npy,npz)
            no(15) = ipoin(ip+2,jp+2,kp+1,npy,npz)
            no(16) = ipoin(ip  ,jp+2,kp+1,npy,npz)
            no(17) = ipoin(ip+1,jp  ,kp+2,npy,npz)
            no(18) = ipoin(ip+2,jp+1,kp+2,npy,npz)
            no(19) = ipoin(ip+1,jp+2,kp+2,npy,npz)
            no(20) = ipoin(ip  ,jp+1,kp+2,npy,npz)
            no(21) = ipoin(ip+1,jp+1,kp  ,npy,npz)
            no(22) = ipoin(ip+1,jp  ,kp+1,npy,npz)
            no(23) = ipoin(ip+2,jp+1,kp+1,npy,npz)
            no(24) = ipoin(ip+1,jp+2,kp+1,npy,npz)
            no(25) = ipoin(ip  ,jp+1,kp+1,npy,npz)
            no(26) = ipoin(ip+1,jp+1,kp+2,npy,npz)
            no(27) = ipoin(ip+1,jp+1,kp+1,npy,npz)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(1),no(2),no(4),no(5),no(9), &
                 no(21),no(12),no(13),no(22),no(25)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(3),no(6),no(7),no(8),no(23), &
                 no(18),no(15),no(24),no(26),no(19)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(2),no(3),no(4),no(6),no(10), &
                 no(11),no(21),no(14),no(23),no(27)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(2),no(4),no(5),no(6),no(21), &
                 no(25),no(22),no(14),no(27),no(17)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(4),no(5),no(6),no(8),no(25), &
                 no(17),no(27),no(16),no(20),no(26)
            ne = ne + 1
            write(lugeo,10) ne, &
                 no(3),no(4),no(6),no(8),no(11), &
                 no(27),no(23),no(24),no(16),no(26)
         end do
        end do
      end do
   10 format(11i6)

      end
!c
!c***  Q2 elements
!c
      subroutine elemq2(nex,ney,nez,npy,npz,lugeo)
      implicit none
      integer &
       nex,ney,nez,npy,npz,lugeo
      integer &
       ne,ie,je,ke,ip,jp,kp,no(27),i, &
       ipoin
      external ipoin

      ne = 0
      do ie = 1,nex
        ip = 2*ie - 1
        do je = 1,ney
          jp = 2*je - 1
          do ke = 1,nez
            kp = 2*ke - 1
            no( 1) = ipoin(ip  ,jp  ,kp  ,npy,npz)
            no( 2) = ipoin(ip+2,jp  ,kp  ,npy,npz)
            no( 3) = ipoin(ip+2,jp+2,kp  ,npy,npz)
            no( 4) = ipoin(ip  ,jp+2,kp  ,npy,npz)
            no( 5) = ipoin(ip  ,jp  ,kp+2,npy,npz)
            no( 6) = ipoin(ip+2,jp  ,kp+2,npy,npz)
            no( 7) = ipoin(ip+2,jp+2,kp+2,npy,npz)
            no( 8) = ipoin(ip  ,jp+2,kp+2,npy,npz)
            no( 9) = ipoin(ip+1,jp  ,kp  ,npy,npz)
            no(10) = ipoin(ip+2,jp+1,kp  ,npy,npz)
            no(11) = ipoin(ip+1,jp+2,kp  ,npy,npz)
            no(12) = ipoin(ip  ,jp+1,kp  ,npy,npz)
            no(13) = ipoin(ip  ,jp  ,kp+1,npy,npz)
            no(14) = ipoin(ip+2,jp  ,kp+1,npy,npz)
            no(15) = ipoin(ip+2,jp+2,kp+1,npy,npz)
            no(16) = ipoin(ip  ,jp+2,kp+1,npy,npz)
            no(17) = ipoin(ip+1,jp  ,kp+2,npy,npz)
            no(18) = ipoin(ip+2,jp+1,kp+2,npy,npz)
            no(19) = ipoin(ip+1,jp+2,kp+2,npy,npz)
            no(20) = ipoin(ip  ,jp+1,kp+2,npy,npz)
            no(21) = ipoin(ip+1,jp+1,kp  ,npy,npz)
            no(22) = ipoin(ip+1,jp  ,kp+1,npy,npz)
            no(23) = ipoin(ip+2,jp+1,kp+1,npy,npz)
            no(24) = ipoin(ip+1,jp+2,kp+1,npy,npz)
            no(25) = ipoin(ip  ,jp+1,kp+1,npy,npz)
            no(26) = ipoin(ip+1,jp+1,kp+2,npy,npz)
            no(27) = ipoin(ip+1,jp+1,kp+1,npy,npz)
            ne = ne + 1
            write(lugeo,10) ne, (no(i),i=1,27)
          end do
        end do
      end do
   10 format(10i8,' /',/,5x,9i8,' /',/,5x,9i8)

      end
      


