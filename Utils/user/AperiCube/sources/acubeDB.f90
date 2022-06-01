subroutine acubeDB(&
     ndime,npoin,nelem,nboun,coord,ddref,xmaxi,ymaxi,zmaxi,npord,neord,mperm)
  implicit none

  integer :: ndime,npoin,nelem,nboun,npord(3),neord(3)
  integer :: mperm(npord(1),npord(2),npord(3)),lnode(8)
  real    :: coord(ndime,npoin),ddref(3),xmaxi,ymaxi,zmaxi,xfield(3)

  integer :: ipoin,ipoix,ipoiy,ipoiz,ielem,ielex,ieley,ielez,iboun,kcode(6)
  real    :: copoi(3),czero(3)
  
  ipoin= 0

!  czero(1)= -ddref(1)
!  czero(2)= -ddref(2)
!  czero(3)= -ddref(3)

!  czero(1)= - 0.5 * xmaxi 
!  czero(2)= - 0.5 * ymaxi
!  czero(3)= - 0.5 * zmaxi

  czero(1)= 0.0
  czero(2)= 0.0
  czero(3)= 0.0

  copoi(1)= czero(1)
  copoi(2)= czero(2)
  copoi(3)= czero(3)

  xfield(1)= 1.0
  xfield(2)= 0.0
  xfield(3)= 0.0

  write(100,'(a)') 'COORDINATES'
  do ipoix=1,npord(1)
     do ipoiy=1,npord(2)
        do ipoiz=1,npord(3)
           ipoin= ipoin+1
           coord(1,ipoin)= copoi(1)
           coord(2,ipoin)= copoi(2)
           coord(3,ipoin)= copoi(3)
           copoi(3)= copoi(3)+ddref(3)
           mperm(ipoix,ipoiy,ipoiz) = ipoin

           write(100,200) ipoin,coord(1:3,ipoin)
           write(500,200) ipoin,xfield(1:3)

        end do
        copoi(2)= copoi(2)+ddref(2)
        copoi(3)= czero(3)
     end do
     copoi(1)= copoi(1)+ddref(1)
     copoi(2)= czero(2)
     copoi(3)= czero(3)
  end do
  write(100,'(a)') 'ENDCOORDINATES'

  ielem=0
  iboun=0

  kcode(1)= 1 
  kcode(2)= 2 
  kcode(3)= 3 
  kcode(4)= 4 
  kcode(5)= 5 
  kcode(6)= 6 
  
  write(200,'(a)') 'ELEMENTS'
  write(300,'(a)') 'BOUNDARIES, ELEMENT'
  write(400,'(a)') 'ONBOUNDARIES'
  do ielex=1,neord(1)
     do ieley=1,neord(2)
        do ielez=1,neord(3)
           ielem=ielem+1
           lnode(1)=  mperm(ielex  ,ieley  ,ielez  ) 
           lnode(2)=  mperm(ielex+1,ieley  ,ielez  ) 
           lnode(3)=  mperm(ielex+1,ieley+1,ielez  ) 
           lnode(4)=  mperm(ielex  ,ieley+1,ielez  ) 
           lnode(5)=  mperm(ielex  ,ieley  ,ielez+1) 
           lnode(6)=  mperm(ielex+1,ieley  ,ielez+1) 
           lnode(7)=  mperm(ielex+1,ieley+1,ielez+1) 
           lnode(8)=  mperm(ielex  ,ieley+1,ielez+1) 

           write(200,100) ielem,lnode(1:8)           

           if (ielex == 1) then
              iboun=iboun+1
              write(300,100) iboun,lnode(1),lnode(5),lnode(8),lnode(4),ielem
              write(400,100) iboun,kcode(1)
           end if
           if (ieley == 1) then
              iboun=iboun+1
              write(300,100) iboun,lnode(1),lnode(2),lnode(6),lnode(5),ielem        
              write(400,100) iboun,kcode(2)
           end if
           if (ielez == 1) then
              iboun=iboun+1
              write(300,100) iboun,lnode(1),lnode(4),lnode(3),lnode(2),ielem                      
              write(400,100) iboun,kcode(3)
           end if
           if (ielex == neord(1)) then
              iboun=iboun+1
              write(300,100) iboun,lnode(2),lnode(3),lnode(7),lnode(6),ielem
              write(400,100) iboun,kcode(4)
           end if
           if (ieley == neord(2)) then
              iboun=iboun+1
              write(300,100) iboun,lnode(4),lnode(8),lnode(7),lnode(3),ielem        
              write(400,100) iboun,kcode(5)

              write(6,100) iboun,lnode(4),lnode(8),lnode(7),lnode(3),ielem        
              write(6,100) iboun,kcode(5)

           end if
           if (ielez == neord(3)) then
              iboun=iboun+1
              write(300,100) iboun,lnode(5),lnode(6),lnode(7),lnode(8),ielem                      
              write(400,100) iboun,kcode(6)
           end if
           
        end do
     end do
  end do
  write(200,'(a)') 'ENDELEMENTS'
  write(300,'(a)') 'ENDBOUNDARIES'
  write(400,'(a)') 'ENDONBOUNDARIES'

  write(50,'(a)') '$----Generated with Apericube--------------------------------------------------------'
  write(50,'(a)') 'DIMENSIONS'
  write(50, 150 ) '   NODAL_POINTS' , npoin
  write(50, 150 ) '   ELEMENTS    ' , nelem
  write(50, 150 ) '   BOUNDARIES  ' , nboun
  write(50,'(a)') '   SPACE_DIMENSIONS     3 '
  write(50,'(a)') '   TYPES_OF_ELEMS     HEX08'
  write(50,'(a)') 'END_DIMENSIONS'
  write(50,'(a)') '$------------------------------------------------------------'
  write(50,'(a)') 'STRATEGY'
  write(50,'(a)') '  INTEGRATION_RULE:          Open'
  write(50,'(a)') '  DOMAIN_INTEGRATION_POINTS: 0'
  write(50,'(a)') '$$$  SCALE XSCAL= 1.0E-3, YSCAL= 1.0E-3, ZSCAL= 1.0E-3'
  write(50,'(a)') 'END_STRATEGY'
  write(50,'(a)') '$-------------------------------------------------------------'
  write(50,'(a)') 'GEOMETRY'
  write(50,'(a)') '  INCLUDE  fort.100'
  write(50,'(a)') '  INCLUDE  fort.200'
  write(50,'(a)') '  INCLUDE  fort.300'
  write(50,'(a)') '  FIELDS, NUMBER = 1'
  write(50,'(a)') '    FIELD = 1, DIMENSION = 3, NODES '
  write(50,'(a)') '      INCLUDE  fort.500'
  write(50,'(a)') '    END_FIELD '
  write(50,'(a)') '  END_FIELDS '
  write(50,'(a)') '  SKEW_SYSTEMS'
  write(50,'(a)') '  END_SKEW_SYSTEMS'
  write(50,'(a)') 'END_GEOMETRY'
  write(50,'(a)') '$-------------------------------------------------------------'
  write(50,'(a)') 'SETS'
  write(50,'(a)') '$$$$  INCLUDE fort.400'
  write(50,'(a)') 'END_SETS'
  write(50,'(a)') '$-------------------------------------------------------------'
  write(50,'(a)') 'BOUNDARY_CONDITIONS  , EXTRAPOLATE'
  write(50,'(a)') '  INCLUDE  fort.400'
  write(50,'(a)') 'END_BOUNDARY_CONDITIONS'
  write(50,'(a)') '$-------------------------------------------------------------'


  
  100 format(20(2x,i8))
  150 format(a15,2x,i8)
  200 format(2x,i8,3(2x,f12.7))

end subroutine acubeDB
