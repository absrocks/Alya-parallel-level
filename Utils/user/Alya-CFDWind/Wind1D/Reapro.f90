   subroutine reapro
     !------------
     ! **** Reapro
     ! this routines reads the main variables of the program 
     !-----------------------
     use def_master
     use Inpout
     implicit none
   

!
!*** Gets problemname from the call argument
!
     call GETARG(1,name) 
     !
     if(TRIM(name).eq.'') then
       call runend('Specify problemname') 
     endif  
     !
     !*** Gets filenames and the working path 
     !
     path ='./'
     !    Input files
     finp     = TRIM(path)//TRIM(name)
     
!     open(luinp,file=TRIM(finp),status='unknown')

     print *, 'input file', finp
!
!*** Reads the input file
!
     call readinp

     kfl_close =0   ! close integration rule (1 = activated)
     kfl_order =1   ! linear(1) or quadratic (2) elements
     !densi = 1.2290d0
     densi = 1.0d0     
     ngaus = 2      ! number of gauss integration points
     nnode = 2      ! nodes per element
     npoin = nelem + 1 ! number of grid points

     if (kfl_bouco_vel.eq.2.and.dwall.lt.1.0e-8) then
        write(*,*) 'WALL DISTANCE SHOULD BE GREATER THAN ZERO FOR WALL LAW'
        write(*,*) 'WALL DISTANCE turned EQUAL TO FIRST ELEMENT LENGTH'
        dwall= dz1
     else if(kfl_bouco_vel.eq.0) then ! non-slip velocity
        write(*,*) 'NON SLIP VELOCITY BOUNDARY CONDITION, WALL DISTANCE TURNED 0.0D0'
        dwall =0.0d0
     end if
     if (kfl_thmod.gt.0.and..not.kfl_abl2)  call wriwar('WARNING: temper without ABL2 wal model')
     

   end subroutine reapro
