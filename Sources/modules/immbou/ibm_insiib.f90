subroutine ibm_insiib(facoo,point,ndime,ifoun)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_insiib
  ! DESCRIPTION
  !    Determine if a point is inside a triangle using the same side technique
  !    point: point
  !    facoo(dimension,vertices): triangle coordinates
  !    ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  ! USED BY
  !    pofadi,faliib
  !----------------------------------------------------------------------- 

  use def_kintyp
  use def_immbou


  implicit   none
  integer(ip),   intent(in)     :: ndime
  integer(ip),   intent(out)    :: ifoun
  integer(ip)                   :: test1,test2,test3
  real(rp),      intent(in)     :: point(ndime)
  real(rp),      intent(in)     :: facoo(ndime,ndime)


  ! The same side technique is applied at each side of the triangle
  call ibm_samesi(facoo(:,1),facoo(:,2),point,facoo(:,3),ndime,test1)
  call ibm_samesi(facoo(:,1),facoo(:,3),point,facoo(:,2),ndime,test2)
  call ibm_samesi(facoo(:,2),facoo(:,3),point,facoo(:,1),ndime,test3)

  if ( test1 == 1 .and. test2 == 1 .and. test3 == 1) then
     ifoun = 1
  else
     ifoun = 0
  end if 

  
end subroutine ibm_insiib


subroutine ibm_samesi(orige,desti,poin1,poin2,ndime,istru)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_samesi
  ! DESCRIPTION
  !    Determine if a point is inside a triangle using the same side technique
  ! USED BY
  !    ibm_insiib
  !----------------------------------------------------------------------- 

  use def_kintyp
  
  implicit   none
  integer(ip), intent(in)   :: ndime
  integer(ip), intent(out)  :: istru
  integer(ip)               :: idime
  real(rp),    intent(in)   :: poin1(ndime),poin2(ndime),orige(ndime),desti(ndime)
  real(rp)                  :: side1(3),side2(3),side3(3),cp1(3),cp2(3),resul
  
  side1(3) = 0.0_rp
  side2(3) = 0.0_rp 
  side3(3) = 0.0_rp 
  do idime=1,ndime
     side1(idime) = desti(idime) - orige(idime)
     side2(idime) = poin1(idime) - orige(idime)
     side3(idime) = poin2(idime) - orige(idime)
  end do
  
  call vecpro(side1,side2,cp1,3)
  call vecpro(side1,side3,cp2,3)
  
  resul=0.0_rp
  do idime = 1,ndime
     resul = resul + (cp1(idime) * cp2(idime))
  end do
  if (resul >= 0) then
     istru = 1
  else
     istru = 0
  end if
  
end subroutine ibm_samesi
