subroutine bouder(nnodb,ndime,ndimb,deriv,bocod,baloc,eucta)
  !-----------------------------------------------------------------------
  !****f* Domain/bouder
  ! NAME
  !    bouder
  ! DESCRIPTION
  !    This routine calculates the baloc and eucta.
  ! USES
  !    vecnor
  !    vecpro
  ! USED BY
  !    nsm_bouope
  !    tem_bouope
  ! SOURCE
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: nnodb,ndime,ndimb
  real(rp),    intent(in)  :: deriv(max(1_ip,ndimb),nnodb)
  real(rp),    intent(in)  :: bocod(ndime,nnodb)
  real(rp),    intent(out) :: eucta,baloc(ndime,ndime)
 
  if( ndime == 1 ) then
     !
     ! 1D
     !
     baloc(1,1) =   1.0_rp
     eucta      =   baloc(1,1)*baloc(1,1)
     eucta      =   sqrt(eucta)
     
  else if( ndime == 2 ) then
     !
     ! 2D
     !
     call mbmabt(baloc,bocod,deriv,ndime,ndimb,nnodb)  ! Evaluates the tangent vectors
     baloc(1,2) =   baloc(2,1)
     baloc(2,2) = - baloc(1,1)
     eucta      =   baloc(1,2) * baloc(1,2) &
          &       + baloc(2,2) * baloc(2,2)
     eucta      =   sqrt(eucta)

  else if( ndime == 3 ) then
     !
     ! 3D
     !
     call mbmabt(baloc,bocod,deriv,ndime,ndimb,nnodb)  ! Evaluates the tangent vectors
     baloc(1,3) =   baloc(2,1) * baloc(3,2) - baloc(3,1) * baloc(2,2)
     baloc(2,3) =   baloc(3,1) * baloc(1,2) - baloc(1,1) * baloc(3,2)
     baloc(3,3) =   baloc(1,1) * baloc(2,2) - baloc(2,1) * baloc(1,2)
     eucta      =   baloc(1,3) * baloc(1,3) &
          &       + baloc(2,3) * baloc(2,3) &
          &       + baloc(3,3) * baloc(3,3)
     eucta      =   sqrt(eucta)

     ! recalculate t1 so that it is orthogonal to t2
     baloc(1,1) =   baloc(2,3) * baloc(3,2) - baloc(3,3) * baloc(2,2)
     baloc(2,1) =   baloc(3,3) * baloc(1,2) - baloc(1,3) * baloc(3,2)
     baloc(3,1) =   baloc(1,3) * baloc(2,2) - baloc(2,3) * baloc(1,2)

  end if

end subroutine bouder

