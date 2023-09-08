subroutine elsest_elq2p2(&
     pnode,lmini,lmaxi,elcod,coglo,coloc,&
     shapf,deriv,ifoun,xjacm,xjaci)
  !-----------------------------------------------------------------------
  !****f* domain/elmder
  ! NAME
  !    elmder
  ! DESCRIPTION
  !    Check if point with global coordinates (x,y)=COGLO is inside
  !    a triangle P1 or a quadrilateral Q2. The Q2 element is
  !    divided into five tetraeder elements. Returns the local coordinates
  !    (s,t)=COLOC
  !
  !    For tetraeder we have:
  !    x = (1-s-t-u)*x1 + s*x2 + t*x3 +u*x4
  !    y = (1-s-t-u)*y1 + s*y2 + t*y3 +u*x4
  !    z = (1-s-t-u)*z1 + s*z2 + t*z3 +u*z4
  !
  !    This linear problem is solved for (s,t,u):
  !         (x....
  !    s =  ----------------------------
  !         (x...
  !
  !         (x...
  !    t =  ----------------------------
  !         (x...
  !
  !         (x...
  !    u =  ----------------------------
  !         (x...
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
  real(rp),    intent(in)  :: lmini,lmaxi,elcod(3,pnode),coglo(3)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(*),shapf(pnode),deriv(3,pnode)
  real(rp),    intent(out) :: xjacm(*),xjaci(*)
  real(rp)                 :: deter,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,x,y,z,ezzzt

  ifoun=0
  if(pnode==4) then
     !
     ! P1 and Q1: Check if point is in first tetraeder: nodes 1-2-3-4
     !
     x1  = elcod(1,1)
     x2  = elcod(1,2)
     x3  = elcod(1,3)
     x4  = elcod(1,4)
     y1  = elcod(2,1)
     y2  = elcod(2,2)
     y3  = elcod(2,3)
     y4  = elcod(2,4)
     z1  = elcod(3,1)
     z2  = elcod(3,2)
     z3  = elcod(3,3)
     z4  = elcod(3,4)
     x   = coglo(1)
     y   = coglo(2)
     z   = coglo(3)
     deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
          & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
          & x1*y3*z4+x2*y4*z3+y2*x3*z4)
     coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
          &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
          &-z3*y*x4-z*x3*y4)
     coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
          & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
          & -y2*z*x4+x1*y2*z)
     coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
          &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
          &+z1*x*y2-x1*z3*y)
     if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
        if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
           if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
              ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
              if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                 ifoun=1
              end if
           end if
        end if
     end if
  else if(pnode==8) then
     !
     ! Q1: Check if point is first tetraeder: nodes 1-2-4-5
     !
     x1  = elcod(1,1)
     x2  = elcod(1,2)
     x3  = elcod(1,4)
     x4  = elcod(1,5)
     y1  = elcod(2,1)
     y2  = elcod(2,2)
     y3  = elcod(2,4)
     y4  = elcod(2,5)
     z1  = elcod(3,1)
     z2  = elcod(3,2)
     z3  = elcod(3,4)
     z4  = elcod(3,5)
     x   = coglo(1)
     y   = coglo(2)
     z   = coglo(3)
     deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
          & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
          & x1*y3*z4+x2*y4*z3+y2*x3*z4)
     coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
          &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
          &-z3*y*x4-z*x3*y4)
     coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
          & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
          & -y2*z*x4+x1*y2*z)
     coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
          &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
          &+z1*x*y2-x1*z3*y)
     if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
        if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
           if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
              ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
              if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                 ifoun=1
              end if
           end if
        end if
     end if
     if(ifoun==0) then
        !
        ! Check tetraeder: nodes 2-3-4-7
        !
        x1  = elcod(1,2)
        x2  = elcod(1,3)
        x3  = elcod(1,4)
        x4  = elcod(1,7)
        y1  = elcod(2,2)
        y2  = elcod(2,3)
        y3  = elcod(2,4)
        y4  = elcod(2,7)
        z1  = elcod(3,2)
        z2  = elcod(3,3)
        z3  = elcod(3,4)
        z4  = elcod(3,7)
        x   = coglo(1)
        y   = coglo(2)
        z   = coglo(3)
        deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
             & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
             & x1*y3*z4+x2*y4*z3+y2*x3*z4)
        coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
             &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
             &-z3*y*x4-z*x3*y4)
        coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
             & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
             & -y2*z*x4+x1*y2*z)
        coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
             &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
             &+z1*x*y2-x1*z3*y)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
                 end if
              end if
           end if
        end if
     end if
     if(ifoun==0) then
        !
        ! Check tetraeder: nodes 2-6-7-5
        !
        x1  = elcod(1,2)
        x2  = elcod(1,6)
        x3  = elcod(1,7)
        x4  = elcod(1,5)
        y1  = elcod(2,2)
        y2  = elcod(2,6)
        y3  = elcod(2,7)
        y4  = elcod(2,5)
        z1  = elcod(3,2)
        z2  = elcod(3,6)
        z3  = elcod(3,7)
        z4  = elcod(3,5)
        x   = coglo(1)
        y   = coglo(2)
        z   = coglo(3)
        deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
             & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
             & x1*y3*z4+x2*y4*z3+y2*x3*z4)
        coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
             &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
             &-z3*y*x4-z*x3*y4)
        coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
             & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
             & -y2*z*x4+x1*y2*z)
        coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
             &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
             &+z1*x*y2-x1*z3*y)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
                 end if
              end if
           end if
         end if
      end if
     if(ifoun==0) then
        !
        ! Check tetraeder: nodes 4-5-7-8
        !
        x1  = elcod(1,4)
        x2  = elcod(1,5)
        x3  = elcod(1,7)
        x4  = elcod(1,8)
        y1  = elcod(2,4)
        y2  = elcod(2,5)
        y3  = elcod(2,7)
        y4  = elcod(2,8)
        z1  = elcod(3,4)
        z2  = elcod(3,5)
        z3  = elcod(3,7)
        z4  = elcod(3,8)
        x   = coglo(1)
        y   = coglo(2)
        z   = coglo(3)
        deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
             & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
             & x1*y3*z4+x2*y4*z3+y2*x3*z4)
        coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
             &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
             &-z3*y*x4-z*x3*y4)
        coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
             & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
             & -y2*z*x4+x1*y2*z)
        coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
             &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
             &+z1*x*y2-x1*z3*y)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
                end if
              end if
           end if
        end if
     end if
     if(ifoun==0) then
        !
        ! Check tetraeder: nodes 2-4-5-7
        !
        x1  = elcod(1,2)
        x2  = elcod(1,4)
        x3  = elcod(1,5)
        x4  = elcod(1,7)
        y1  = elcod(2,2)
        y2  = elcod(2,4)
        y3  = elcod(2,5)
        y4  = elcod(2,7)
        z1  = elcod(3,2)
        z2  = elcod(3,4)
        z3  = elcod(3,5)
        z4  = elcod(3,7)
        x   = coglo(1)
        y   = coglo(2)
        z   = coglo(3)
        deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
             & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
             & x1*y3*z4+x2*y4*z3+y2*x3*z4)
        coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
             &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
             &-z3*y*x4-z*x3*y4)
        coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
             & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
             & -y2*z*x4+x1*y2*z)
        coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
             &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
             &+z1*x*y2-x1*z3*y)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
               end if
              end if
           end if
        end if
     end if
     if(ifoun==1) then
        call elsest_newrap(&
             coglo,coloc,3_ip,8_ip,elcod,xjacm,xjaci,shapf,deriv)
     end if
  else if(pnode==5) then
     !
     ! Q1: Check if point is first tetraeder: nodes 1-2-4-5
     !
     x1  = elcod(1,1)
     x2  = elcod(1,2)
     x3  = elcod(1,4)
     x4  = elcod(1,5)
     y1  = elcod(2,1)
     y2  = elcod(2,2)
     y3  = elcod(2,4)
     y4  = elcod(2,5)
     z1  = elcod(3,1)
     z2  = elcod(3,2)
     z3  = elcod(3,4)
     z4  = elcod(3,5)
     x   = coglo(1)
     y   = coglo(2)
     z   = coglo(3)
     deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 &
          & +z1*x2*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+  &
          & x1*y3*z4+x2*y4*z3+y2*x3*z4)
     coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 &
          &-z1*x4*y3+y1*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x &
          &-z3*y*x4-z*x3*y4)
     coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 &
          & -y1*z4*x+y2*z4*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 &
          & -y2*z*x4+x1*y2*z)
     coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 &
          &-x2*z*y3+x1*z*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 &
          &+z1*x*y2-x1*z3*y)
     if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
        if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
           if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
              ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
              if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                 ifoun=1
            end if
           end if
        end if
     end if
     if(ifoun==0) then
        !
        ! Check tetraeder: nodes 2-3-4-5
        !
        x1  = elcod(1,2)
        x2  = elcod(1,3)
        x3  = elcod(1,4)
        x4  = elcod(1,5)
        y1  = elcod(2,2)
        y2  = elcod(2,3)
        y3  = elcod(2,4)
        y4  = elcod(2,5)
        z1  = elcod(3,2)
        z2  = elcod(3,3)
        z3  = elcod(3,4)
        z4  = elcod(3,5)
        x   = coglo(1)
        y   = coglo(2)
        z   = coglo(3)
        deter    = 1.0_8/( z2*x4*y3-x1*z2*y3+x1*z2*y4-y1*z2*x4+y1*z2&
             &*x3-z2*x3*y4-x2*y3*z4-y2*x4*z3-z1*y2*x3+z1*y2*x4 +z1*x2&
             &*y3-z1*x2*y4-z1*x4*y3+z1*x3*y4-y1*x2*z3+y1*x2*z4+y1*x4&
             &*z3-y1*x3*z4+x1*y2*z3-x1*y2*z4-x1*y4*z3+ x1*y3*z4+x2*y4&
             &*z3+y2*x3*z4)
        coloc(1) = deter*(x1*z3*y+x3*z*y1-x3*z1*y-x1*z*y3-x*y1*z3+z1&
             &*x*y3-x1*y4*z3+z1*x3*y4-y1*x3*z4+x1*y3*z4 -z1*x4*y3+y1&
             &*x4*z3+y*z1*x4+y1*z4*x-z*y1*x4+y4*x1*z-z1*y4*x-y*x1*z4&
             &+z*x4*y3-z4*x*y3+z4*x3*y+z3*y4*x -z3*y*x4-z*x3*y4)
        coloc(2) = deter*(x*z2*y1-z2*y4*x+z2*y*x4+x1*z2*y4-x1*y*z2-y1&
             &*z2*x4+z1*y2*x4-z1*x2*y4+x2*y*z1-z1*x*y2 -y1*z4*x+y2*z4&
             &*x+z1*y4*x-x1*y2*z4+z*y1*x4+y4*x2*z-y4*x1*z-x2*z*y1+y&
             &*x1*z4+y1*x2*z4-y*x2*z4-y*z1*x4 -y2*z*x4+x1*y2*z)
        coloc(3) = deter*(y1*z2*x3-x*z2*y1-x1*z2*y3+z2*x*y3-z2*x3*y&
             &+x1*y*z2-z1*x*y3+x*y1*z3-x*y2*z3+x2*z*y1 -x2*z*y3+x1*z&
             &*y3-x1*y2*z+x3*z1*y-x3*z*y1+x3*y2*z+x2*z3*y-z1*y2*x3+z1&
             &*x2*y3-y1*x2*z3+x1*y2*z3-x2*y*z1 +z1*x*y2-x1*z3*y)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
                 end if
              end if
           end if
        end if
     end if
     if(ifoun==1) then
        call elsest_newrap(&
             coglo,coloc,3_ip,5_ip,elcod,xjacm,xjaci,shapf,deriv)
     end if

  end if

  if(pnode==4.and.ifoun==1) then
     !
     ! P1: Compute shape function and derivatives for tetraeder
     !
     shapf(1)  = ezzzt
     shapf(2)  = coloc(1)
     shapf(3)  = coloc(2)
     shapf(4)  = coloc(3)
     deriv(1, 1) =-1.0_rp
     deriv(2, 1) =-1.0_rp
     deriv(3, 1) =-1.0_rp
     deriv(1, 2) = 1.0_rp
     deriv(2, 3) = 1.0_rp
     deriv(3, 4) = 1.0_rp

  end if

end subroutine elsest_elq2p2
