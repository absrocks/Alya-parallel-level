subroutine nsa_elmchl(tragl,hleng,chave,chale,ndime,&
     pnode,hnatu,kfl_advec,zeroc)

!------------------------------------------------------------------------
!
! This routine computes the characteristic element length in the
! flow direction (if there are convective terms)
!
!------------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ndime,pnode,kfl_advec
  real(rp),    intent(in)  :: hnatu,zeroc
  real(rp),    intent(out) :: chale(2)
  real(rp),    intent(in)  :: tragl(ndime,ndime),hleng(ndime)
  real(rp),    intent(out) :: chave(ndime,2)
  integer(ip)              :: idime,ievab,inode,ivepo
  real(rp)                 :: elno1,elno2
!
! Initialization
!
  chale(1)=hleng(1) ! Maximum element length
  chale(2)=hleng(1)
!
! Length in the flow direction
!      
  if(kfl_advec==1) then 
     ! Characteristic element length
     call mbvab1(chave(1,2),tragl,chave(1,1),ndime,ndime,elno2,elno1)
     if(elno2>zeroc.and.elno1>zeroc) chale(1)=hnatu*elno1/elno2
  end if
!
! Divide h by 2 for quadratic elements
!
  if(ndime==2) then
     if(pnode>5) chale = 0.5_rp*chale
  else if(ndime==3) then
     if(pnode>9) chale = 0.5_rp*chale
  end if
  
end subroutine nsa_elmchl


