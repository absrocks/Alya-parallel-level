subroutine elsest_chkelm(&
     ndime,ptopo,pnode,elcod,shapf,deriv,&
     coglo,coloc,ifoun,lmini,lmaxi)
  !-----------------------------------------------------------------------
  !****f* Domain/elsest_chkelm
  ! NAME
  !    elsest_chkelm
  ! DESCRIPTION
  !    Check if a point belongs to an element
  ! INPUT
  !    NDIME ... Dimension
  !    PTOPO ... Element topology
  !              0 Quadrilateral  Hexahedra
  !              1 Triangle       Tetrahedra
  !              2    -           Pentahedra (wedge-shaped)
  !              3    -           Pyramid
  !    PNODE ... Number of element nodes
  !    ELCOD ... Element node coordinates
  !    COGLO ... Global coordinates of test point
  !    LMINI ... Minimum local coordinate (0.0)
  !    LMAXI ... Maximum local coordinate (1.0)
  ! OUTPUT
  !    IFOUN ... 1 if point is in element
  !              0 otherwise
  !    COLOC ... Local coordinates of test point
  !    SHAPF ... Shape function of test point in element
  !    DERIV ... Shape function derivatives of test point in element
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only     :  ip,rp
  use def_elmtyp, only     :  TET04
  use def_elmtyp, only     :  PYR05
  use def_elmtyp, only     :  PEN06
  use def_elmtyp, only     :  HEX08
  use mod_elmgeo, only     :  elmgeo_newrap
  use mod_elmgeo, only     :  elmgeo_inside_TRI03_QUA04
  use mod_elmgeo, only     :  elmgeo_inside_TET04
  use mod_elmgeo, only     :  elmgeo_inside_element_bounding_box
  use mod_elmgeo, only     :  elmgeo_inside_element_using_faces
  implicit none
  integer(ip), intent(in)  :: ndime,ptopo,pnode
  real(rp),    intent(in)  :: coglo(ndime),elcod(ndime,pnode)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(3),deriv(ndime,pnode),shapf(pnode)
  real(rp)                 :: ezzzt,lmaxi,lmini
  integer(ip)              :: i,j,k,pelty

  if( elmgeo_inside_element_bounding_box(ndime,pnode,elcod,coglo,abs(lmini)) ) then

     if(      pnode == 4 .and. ndime == 3 ) then
        pelty = TET04
     else if( pnode == 5 .and. ndime == 3 ) then
        pelty = PYR05
     else if( pnode == 6 .and. ndime == 3 ) then
        pelty = PEN06
     else if( pnode == 8 .and. ndime == 3 ) then
        pelty = HEX08
     end if

     if(  (ptopo==0.and.ndime==2.and.pnode==4).or.&
          (ptopo==1.and.ndime==2.and.pnode==3) ) then
        !
        ! Specific treatment: P1 and Q1 elements
        !
        call elmgeo_inside_TRI03_QUA04(&
             pnode,lmini,lmaxi,elcod,coglo,coloc,shapf,deriv,ifoun)

     else if( pelty == TET04 ) then
        !
        ! TETO04 
        !
        call elmgeo_inside_TET04(&
             lmini,lmaxi,elcod,coglo,coloc,&
             shapf,deriv,ifoun)

     else
        !
        ! PYR05 
        !
        if( elmgeo_inside_element_using_faces(pelty,elcod,coglo) ) then
           call elmgeo_newrap(coglo,coloc,ndime,pnode,elcod,shapf,deriv)
           ifoun = 1
        end if

     !   call elmgeo_inside_PYR05(&
     !        lmini,lmaxi,elcod,coglo,coloc,&
     !        shapf,deriv,ifoun)
     !else if( pnode == 6 .and. ndime == 3 ) then
     !   !
     !   ! PEN06
     !   !
     !   call elmgeo_inside_PEN06(&
     !        lmini,lmaxi,elcod,coglo,coloc,&
     !        shapf,deriv,ifoun)
!!$     else
!!$        !
!!$        ! Compute local from global coordinates
!!$        !
!!$        call elmgeo_newrap(coglo,coloc,ndime,pnode,elcod,shapf,deriv)
!!$
!!$        ifoun=0
!!$
!!$        if(ptopo==0.or.ptopo==-1) then
!!$           !
!!$           ! Quadrilateral and Hexahedra
!!$           !
!!$           if((coloc(1)>=-lmaxi).and.(coloc(1)<=lmaxi)) then
!!$              if((coloc(2)>=-lmaxi).and.(coloc(2)<=lmaxi)) then
!!$                 if((coloc(3)>=-lmaxi).and.(coloc(3)<=lmaxi)) then
!!$                    ifoun=1
!!$                 end if
!!$              end if
!!$           end if
!!$
!!$        else if(ptopo==1) then
!!$           !
!!$           ! Triangle and Tetrahedra
!!$           !
!!$           if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
!!$              if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
!!$                 if(ndime==2) then
!!$                    coloc(3) = 1.0_rp-coloc(1)-coloc(2)
!!$                    ezzzt    = 0.0_rp
!!$                    !else
!!$                    !   ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
!!$                 end if
!!$                 if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
!!$                    if(ndime==3) then
!!$                       ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
!!$                    end if
!!$                    if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
!!$                       ifoun=1
!!$                    end if
!!$                 end if
!!$              end if
!!$           end if
!!$
!!$        else if(ptopo==2) then
!!$           !
!!$           ! Pentahedra
!!$           !
!!$           if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
!!$              if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
!!$                 ezzzt = 1.0_rp - coloc(1) - coloc(2)
!!$                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
!!$                    if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
!!$                       ifoun=1
!!$                    end if
!!$                 end if
!!$              end if
!!$           end if
!!$
!!$        else if( ptopo == 3 ) then
!!$           !
!!$           ! Pyramids
!!$           !
!!$           if((coloc(3)>=-lmaxi).and.(coloc(3)<=lmaxi)) then
!!$              ezzzt = (1.0_rp-coloc(3))*0.5_rp + abs(lmini)
!!$              if((coloc(1)>=-ezzzt).and.(coloc(1)<=ezzzt)) then
!!$                 if((coloc(2)>=-ezzzt).and.(coloc(2)<=ezzzt)) then
!!$                    ifoun=1
!!$                 end if
!!$              end if
!!$           end if
!!$
!!$        end if
     end if

  end if

end subroutine elsest_chkelm
