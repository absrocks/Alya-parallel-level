subroutine domfac(ielty,pnodb,nface,lface,ltypf,nnodf)
  !-----------------------------------------------------------------------
  !****f* Domain/domfac
  ! NAME
  !    domfac
  ! DESCRIPTION
  !    Compute the faces of an element
  ! OUTPUT
  !    LFACE(:,IFACE) ... Face IFACE connectivity
  !    NNODF(IFACE) ..... Number of node of face IFACE
  !    LTYPF(IFACE) ..... Element type of face IFACE
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  nnode,cenam
  use def_master, only     :  lun_outpu
  use mod_elmgeo, only     :  list_faces_BAR02,type_faces_BAR02
  use mod_elmgeo, only     :  list_faces_TRI06,type_faces_TRI06
  use mod_elmgeo, only     :  list_faces_TRI03,type_faces_TRI03
  use mod_elmgeo, only     :  list_faces_QUA04,type_faces_QUA04
  use mod_elmgeo, only     :  list_faces_QUA08,type_faces_QUA08
  use mod_elmgeo, only     :  list_faces_QUA09,type_faces_QUA09
  use mod_elmgeo, only     :  list_faces_TET04,type_faces_TET04
  use mod_elmgeo, only     :  list_faces_TET10,type_faces_TET10
  use mod_elmgeo, only     :  list_faces_PYR05,type_faces_PYR05
  use mod_elmgeo, only     :  list_faces_PEN06,type_faces_PEN06
  use mod_elmgeo, only     :  list_faces_PEN15,type_faces_PEN15
  use mod_elmgeo, only     :  list_faces_PEN18,type_faces_PEN18
  use mod_elmgeo, only     :  list_faces_HEX08,type_faces_HEX08
  use mod_elmgeo, only     :  list_faces_HEX27,type_faces_HEX27
  use mod_outfor, only     :  outfor
  use def_elmtyp
  implicit none 
  integer(ip), intent(in)  :: ielty,pnodb,nface
  integer(ip), intent(out) :: lface(pnodb,nface),ltypf(nface),nnodf(nface)
  integer(ip)              :: iface,inodf,kface

  iface = 0

  if( ielty < 10 ) then
     !
     ! BAR elements: BAR02 to BAR04 
     !
     !   1       2
     !
     lface(1:1,1:2) = list_faces_BAR02(1:1,1:2)
     ltypf(1:2)     = type_faces_BAR02(1:2)

  else if( ielty == TRI03 ) then
     !
     ! TRI03
     !
     !   3              
     !                
     !                
     !   1       2      
     !
     lface(1:2,1:3) = list_faces_TRI03(1:2,1:3)
     ltypf(1:3)     = type_faces_TRI03(1:3)

  else if( ielty == QUA04 ) then
     !
     ! QUA04
     !
     !   4         3  
     !               
     !               
     !   1         2  
     !
     lface(1:2,1:4) = list_faces_QUA04(1:2,1:4)
     ltypf(1:4)     = type_faces_QUA04(1:4)

  else if( ielty == TRI06 ) then
     !
     ! TRI06
     !
     !   3               
     !                  
     !                  
     !                  
     !   6      5        
     !                  
     !                  
     !                  
     !                  
     !   1      4       2
     !
     lface(1:3,1:3) = list_faces_TRI06(1:3,1:3)
     ltypf(1:3)     = type_faces_TRI06(1:3)

  else if( ielty == QUA08 ) then
     !
     ! QUA08 
     !
     call runend('DOMFAC: NOT CODED')
     iface          = iface+1         !  4      7      3
     ltypf(iface)   = BAR03           !                 
     lface(1,iface) = 1               !                 
     lface(2,iface) = 2               !                 
     lface(3,iface) = 5               !  8      9      6
                                      !                 
     iface          = iface+1         !                 
     ltypf(iface)   = BAR03           !                 
     lface(1,iface) = 2               !  1      5      2
     lface(2,iface) = 3
     lface(3,iface) = 6

     iface          = iface+1
     ltypf(iface)   = BAR03
     lface(1,iface) = 3
     lface(2,iface) = 4
     lface(3,iface) = 7

     iface          = iface+1
     ltypf(iface)   = BAR03
     lface(1,iface) = 4
     lface(2,iface) = 1
     lface(3,iface) = 8


  else if( ielty == QUA09 ) then
     !
     ! QUA09
     !
     !  4      7      3
     !                 
     !                 
     !                 
     !  8      9      6
     !                 
     !                 
     !                 
     !  1      5      2
     lface(1:3,1:4) = list_faces_QUA09(1:3,1:4)
     ltypf(1:4)     = type_faces_QUA09(1:4)

  else if ( ielty == QUA16 ) then
     !
     ! QUA16
     !
     call runend('DOMFAC: NOT CODED')
     iface          = iface+1         ! 4    10    9    3
     ltypf(iface)   = BAR04           !                  
     lface(1,iface) = 1               !                  
     lface(2,iface) = 2               ! 11   16   15    8
     lface(3,iface) = 5               !                  
     lface(4,iface) = 6               !                  
                                      ! 12   13   14    7
     iface          = iface+1         !                  
     ltypf(iface)   = BAR04           !                  
     lface(1,iface) = 2               ! 1     5    6    2
     lface(2,iface) = 3
     lface(3,iface) = 7
     lface(4,iface) = 8

     iface          = iface+1
     ltypf(iface)   = BAR04     
     lface(1,iface) = 3
     lface(2,iface) = 4
     lface(3,iface) = 9
     lface(4,iface) = 10

     iface          = iface+1
     ltypf(iface)   = BAR04     
     lface(1,iface) = 4
     lface(2,iface) = 1
     lface(3,iface) = 11
     lface(4,iface) = 12

  else if( ielty == TET04 ) then
     !
     ! TET04
     !
     lface(1:3,1:4) = list_faces_TET04(1:3,1:4)
     ltypf(1:4)     = type_faces_TET04(1:4)

  else if( ielty == PYR05 ) then
     !
     ! PYR05
     !
     lface(1:4,1:5) = list_faces_PYR05(1:4,1:5)
     ltypf(1:5)     = type_faces_PYR05(1:5)

  else if( ielty == PEN06 ) then
     !
     ! PEN06
     !
     lface(1:4,1:5) = list_faces_PEN06(1:4,1:5)
     ltypf(1:5)     = type_faces_PEN06(1:5)

  else if( ielty == PEN15 ) then
     !
     ! PEN15
     !
     lface(1:8,1:5) = list_faces_PEN15(1:8,1:5)
     ltypf(1:5)     = type_faces_PEN15(1:5)

  else if( ielty == PEN18 ) then
     !
     ! PEN18
     !
     ! See GiD reference for node orientation!
     !
     lface(1:9,1:5) = list_faces_PEN18(1:9,1:5)
     ltypf(1:5)     = type_faces_PEN18(1:5)

  else if( ielty == HEX08 ) then
     !
     ! HEX08
     !
     lface(1:4,1:6) = list_faces_HEX08(1:4,1:6)
     ltypf(1:6)     = type_faces_HEX08(1:6)

  else if( ielty == TET10 ) then
     !
     ! TET10
     !
     lface(1:6,1:4) = list_faces_TET10(1:6,1:4)
     ltypf(1:4)     = type_faces_TET10(1:4)

!!$     iface          = iface+1
!!$     ltypf(iface)   = TRI06
!!$     lface(1,iface) = 1
!!$     lface(2,iface) = 2
!!$     lface(3,iface) = 4
!!$     lface(4,iface) = 5
!!$     lface(5,iface) = 9
!!$     lface(6,iface) = 8
!!$
!!$     iface          = iface+1
!!$     ltypf(iface)   = TRI06
!!$     lface(1,iface) = 2
!!$     lface(2,iface) = 3
!!$     lface(3,iface) = 4
!!$     lface(4,iface) = 6
!!$     lface(5,iface) = 10
!!$     lface(6,iface) = 9
!!$
!!$     iface          = iface+1
!!$     ltypf(iface)   = TRI06
!!$     lface(1,iface) = 3
!!$     lface(2,iface) = 1
!!$     lface(3,iface) = 4
!!$     lface(4,iface) = 7
!!$     lface(5,iface) = 8
!!$     lface(6,iface) = 10
!!$
!!$     iface          = iface+1
!!$     ltypf(iface)   = TRI06
!!$     lface(1,iface) = 1
!!$     lface(2,iface) = 2
!!$     lface(3,iface) = 3
!!$     lface(4,iface) = 5
!!$     lface(5,iface) = 6
!!$     lface(6,iface) = 7

  else if( ielty == HEX27 ) then
     !
     ! HEX27
     !
     lface(1:9,1:6) = list_faces_HEX27(1:9,1:6)
     ltypf(1:6)     = type_faces_HEX27(1:6)

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 1
!!$     lface(2,iface) = 4
!!$     lface(3,iface) = 3
!!$     lface(4,iface) = 2
!!$     lface(5,iface) = 12
!!$     lface(6,iface) = 11
!!$     lface(7,iface) = 10
!!$     lface(8,iface) = 9
!!$     lface(9,iface) = 21

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 1
!!$     lface(2,iface) = 2
!!$     lface(3,iface) = 6
!!$     lface(4,iface) = 5
!!$     lface(5,iface) = 9
!!$     lface(6,iface) = 14
!!$     lface(7,iface) = 17
!!$     lface(8,iface) = 13
!!$     lface(9,iface) = 22

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 2
!!$     lface(2,iface) = 3
!!$     lface(3,iface) = 7
!!$     lface(4,iface) = 6
!!$     lface(5,iface) = 10
!!$     lface(6,iface) = 15
!!$     lface(7,iface) = 18
!!$     lface(8,iface) = 14
!!$     lface(9,iface) = 23

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 3
!!$     lface(2,iface) = 4
!!$     lface(3,iface) = 8
!!$     lface(4,iface) = 7
!!$     lface(5,iface) = 11
!!$     lface(6,iface) = 16
!!$     lface(7,iface) = 19
!!$     lface(8,iface) = 15
!!$     lface(9,iface) = 24

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 1
!!$     lface(2,iface) = 5
!!$     lface(3,iface) = 8
!!$     lface(4,iface) = 4
!!$     lface(5,iface) = 13
!!$     lface(6,iface) = 20
!!$     lface(7,iface) = 16
!!$     lface(8,iface) = 12
!!$     lface(9,iface) = 25

!!$     iface          = iface+1
!!$     ltypf(iface)   = QUA09 
!!$     lface(1,iface) = 5
!!$     lface(2,iface) = 6
!!$     lface(3,iface) = 7
!!$     lface(4,iface) = 8
!!$     lface(5,iface) = 17
!!$     lface(6,iface) = 18
!!$     lface(7,iface) = 19
!!$     lface(8,iface) = 20
!!$     lface(9,iface) = 26

  else if ( ielty == SHELL ) then
     !
     ! SHELL
     !
     !iface          = iface+1         !  3        
     !ltypf(iface)   = BAR02           !           
     !lface(1,iface) = 1               !           
     !lface(2,iface) = 2               !  1       2
     !iface          = iface+1
     !ltypf(iface)   = BAR02
     !lface(1,iface) = 2
     !lface(2,iface) = 3
     !iface          = iface+1
     !ltypf(iface)   = BAR02
     !lface(1,iface) = 3
     !lface(2,iface) = 1
     lface(1:2,1:3) = list_faces_TRI03(1:2,1:3)
     ltypf(1:3)     = type_faces_TRI03(1:3)

  else if( ielty == BAR3D ) then
     !
     ! BAR3D
     !
     !iface          = iface+1         !  1------------2
     !ltypf(iface)   = POINT           
     !lface(1,iface) = 1              
     !iface          = iface+1
     !ltypf(iface)   = POINT           
     !lface(1,iface) = 2              
     lface(1:1,1:2) = list_faces_BAR02(1:1,1:2)
     ltypf(1:2)     = type_faces_BAR02(1:2)

  else

     call runend('FACES NOT PROGRAMMED')
  end if
  !
  ! NNODF: Number of node for each face
  !
  do kface = 1,nface
     nnodf(kface) = nnode(ltypf(kface))
  end do
  !
  ! Check errors
  !
  !if( iface /= nface )&
  !     call outfor(1_ip,lun_outpu,&
  !     'CHECK FACES OF ELEMENT '//trim(cenam(ielty)))
  do iface = 1,nface
     do inodf = 1,nnodf(iface)
        if( lface(inodf,iface) == 0 ) then
           call outfor(1_ip,lun_outpu,&
                'WRONG FACE LIST FOR ELEMENT '//trim(cenam(ielty)))
        end if
     end do
  end do

end subroutine domfac

