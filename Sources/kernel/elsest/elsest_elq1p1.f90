subroutine elsest_elq1p1(&
     pnode,lmini,lmaxi,elcod,coglo,coloc,&
     shapf,deriv,ifoun,xjacm,xjaci)
  !-----------------------------------------------------------------------
  !****f* domain/elmder
  ! NAME
  !    elmder
  ! DESCRIPTION
  !    Check if point with global coordinates (x,y)=COGLO is inside
  !    a triangle P1 or a quadrilateral Q1. The Q1 element is
  !    divided into two P1 elements. Returns the local coordinates
  !    (s,t)=COLOC
  !
  !    For P1 triangles we have:
  !    x = (1-s-t)*x1 + s*x2 + t*x3
  !    y = (1-s-t)*y1 + s*y2 + t*y3
  !
  !    This linear problem is solved for (s,t):
  !         (x3-x1)(y-y1) -(y3-y1)(x-x1)
  !    s =  ----------------------------
  !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  !
  !         (x-x1)(y2-y1) -(y-y1)(x2-x1)
  !    t =  ----------------------------
  !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: lmini,lmaxi,elcod(2,pnode),coglo(2)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(*),shapf(pnode),deriv(2,pnode)
  real(rp),    intent(out) :: xjacm(*),xjaci(*)
  real(rp)                 :: deter,colo3,x2x1,y2y1,x3x1,y3y1,xx1,yy1

  ifoun=0
  !
  ! P1 and Q1: Check if point is in first triangle: nodes 1-2-3
  !
  x2x1     = elcod(1,2)-elcod(1,1)
  y2y1     = elcod(2,2)-elcod(2,1)
  x3x1     = elcod(1,3)-elcod(1,1)
  y3y1     = elcod(2,3)-elcod(2,1)
  xx1      = coglo(1)  -elcod(1,1)
  yy1      = coglo(2)  -elcod(2,1)
  deter    = 1.0_8/(x3x1*y2y1-y3y1*x2x1)
  coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
  coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
  if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
     if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
        colo3 = 1.0_rp-coloc(1)-coloc(2)
        if(colo3>=lmini.and.colo3<=lmaxi) ifoun=1
     end if
  end if

  if(pnode==4) then
     !
     ! Q1: Check if point is in second triangle: nodes 1-3-4
     !
     if(ifoun==0) then
        x2x1     = elcod(1,3)-elcod(1,1)
        y2y1     = elcod(2,3)-elcod(2,1)
        x3x1     = elcod(1,4)-elcod(1,1)
        y3y1     = elcod(2,4)-elcod(2,1)
        xx1      = coglo(1)  -elcod(1,1)
        yy1      = coglo(2)  -elcod(2,1)
        deter    = 1.0_8/(x3x1*y2y1-y3y1*x2x1)
        coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
        coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              colo3 = 1.0_rp-coloc(1)-coloc(2)
              if(colo3>=lmini.and.colo3<=lmaxi) ifoun=1
           end if
        end if
     end if
     if(ifoun==1) then
        call elsest_newrap(&
             coglo,coloc,2_ip,4_ip,elcod,xjacm,xjaci,shapf,deriv)
     end if

  else if(pnode==3.and.ifoun==1) then
     !
     ! P1: Compute shape function and derivatives
     !
     shapf(1)  = colo3
     shapf(2)  = coloc(1)
     shapf(3)  = coloc(2)
     deriv(1,1)=-1.0_rp
     deriv(1,2)= 1.0_rp
     deriv(1,3)= 0.0_rp
     deriv(2,1)=-1.0_rp
     deriv(2,2)= 0.0_rp
     deriv(2,3)= 1.0_rp

  end if

end subroutine elsest_elq1p1
