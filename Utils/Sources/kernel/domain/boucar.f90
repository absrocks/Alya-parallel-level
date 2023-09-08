subroutine boucar(&
     pnodb,pgaub,gbsha,deriv,weigb,bocod,gbsur,baloc)
  !-----------------------------------------------------------------------
  !****f* domain/boucar
  ! NAME
  !    boucar
  ! DESCRIPTION
  !    This routine calculates:
  !    GBSUR: Unit surface
  !    BALOC: Tangent system
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_parame, only     :  twopi
  use def_domain, only     :  ndime,ndimb,kfl_naxis
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnodb,pgaub
  real(rp),    intent(in)  :: gbsha(pnodb,pgaub)
  real(rp),    intent(in)  :: deriv(ndimb,pnodb,pgaub)
  real(rp),    intent(in)  :: weigb(pgaub)
  real(rp),    intent(in)  :: bocod(ndime,pnodb)
  real(rp),    intent(out) :: gbsur(pgaub)
  real(rp),    intent(out) :: baloc(ndime,ndime)
  integer(ip)              :: igaub,inodb
  real(rp)                 :: eucta,gpcod

  if(ndime==1) then
     gbsur(1)   = 1.0_rp
     baloc(1,1) = 1.0_rp
  else
     do igaub=1,pgaub
        call bouder(&
             pnodb,ndime,ndimb,deriv(1,1,igaub),&
             bocod,baloc,eucta)
        gbsur(igaub)=weigb(igaub)*eucta 
     end do
     !
     ! Cylindrical coordinates
     !  
     if(kfl_naxis==1) then
        do igaub=1,pgaub
           gpcod=0.0_rp
           do inodb=1,pnodb
              gpcod=gpcod+bocod(1,inodb)*gbsha(inodb,igaub)
           end do
           gbsur(igaub)=gbsur(igaub)*gpcod*twopi
        end do
     end if

  end if
  
end subroutine boucar
