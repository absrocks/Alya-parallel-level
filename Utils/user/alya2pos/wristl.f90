subroutine wristl(kfl_stlbo,namda,ittim,lun,coord,kfl_displ,displ,npoin_total,lnodb,mnodb,ltypb,nboun_total)

  use def_kintyp
  use def_elmtyp

  implicit none

  character(150), intent(in)    :: namda
  integer(ip),    intent(in)    :: ittim,kfl_stlbo,kfl_displ,nboun_total,lun,npoin_total,mnodb
  real(rp),       intent(in)    :: coord(3,*),displ(3,*)
  integer(ip),    intent(in)    :: lnodb(mnodb,nboun_total),ltypb(nboun_total) !mnodb nboun_t
  character(8)                  :: chtim
  character(80)                 :: stl_header
  character(150)                :: filename
  integer(ip)                   :: iboun, p1,p2,p3
  real(rp)                      :: xx_coords(3,3)

  !
  ! Take time stamp and add to file name, nothing if time is negative (for single file output)
  !
  if(ittim<0) then
     chtim = 'ini'
  else if(ittim<10) then
     write(chtim,'(a,i1)') '0000000',ittim
  else if(ittim<100) then
     write(chtim,'(a,i2)') '000000',ittim
  else if(ittim<1000) then
     write(chtim,'(a,i3)') '00000',ittim
  else if(ittim<10000) then
     write(chtim,'(a,i4)') '0000',ittim
  else if(ittim<100000) then
     write(chtim,'(a,i5)') '000',ittim
  else if(ittim<1000000) then
     write(chtim,'(a,i6)') '00',ittim
  else if(ittim<10000000) then
     write(chtim,'(a,i7)') '0',ittim
  end if
  filename = trim(namda)//'.'//trim(chtim)//'.stl'

  if( kfl_stlbo == 2 ) then !Binary format
     open (unit=lun, FILE = trim(filename),ACTION='WRITE', STATUS = 'UNKNOWN', form='unformatted', access='stream' )
     !Header: 80 bytes (title) + 4 bytes (size)
     stl_header = '01234567012345670123456701234567012345670123456701234567012345670123456701234567'
     write(lun) stl_header
     write(lun) int(nboun_total,4)
  else
     open(unit=lun,file=trim(filename),form='formatted')
     write(lun,'(a)') 'solid '//trim(namda)
  endif

  do iboun = 1, nboun_total
     if( abs(ltypb(iboun)) == TRI03 ) then
        p1 = lnodb(1,iboun)
        p2 = lnodb(2,iboun)
        p3 = lnodb(3,iboun)
        call displace(p1,p2,p3,coord,kfl_displ,displ,xx_coords)
        call wristl_facet(lun,kfl_stlbo,xx_coords)
     else if(abs(ltypb(iboun)) == QUA04 ) then
        p1             = lnodb(1,iboun)
        p2             = lnodb(2,iboun)
        p3             = lnodb(3,iboun)
        call displace(p1,p2,p3,coord,kfl_displ,displ,xx_coords)
        call wristl_facet(lun,kfl_stlbo,xx_coords)
        p1             = lnodb(1,iboun)
        p2             = lnodb(3,iboun)
        p3             = lnodb(4,iboun)
        call displace(p1,p2,p3,coord,kfl_displ,displ,xx_coords)
        call wristl_facet(lun,kfl_stlbo,xx_coords)
     else
        write(*,*) 'STL output can only handle TRI03 and QUA04 boundaries'
     end if
  end do

  if( kfl_stlbo /= 2 ) then ! Ascii format needs tailer
     write(lun,'(a)') 'endsolid '//trim(namda)
  end if

end subroutine wristl

!
! Displace node coordinates if DISPL variable is available
!
subroutine displace(p1,p2,p3,coord,kfl_displ,displ,xx_coords)
  use def_kintyp
  implicit none
  integer(ip),    intent(in)    :: kfl_displ
  real(rp),       intent(in)    :: coord(3,*),displ(3,*)
  integer(ip),    intent(in)    :: p1,p2,p3
  real(rp),       intent(out)   :: xx_coords(3,3)
  integer  :: dim

  if (kfl_displ==0) then
     do dim=1,3
        xx_coords(dim,1) = coord(dim,p1)
        xx_coords(dim,2) = coord(dim,p2)
        xx_coords(dim,3) = coord(dim,p3)
     enddo
  else
     do dim=1,3
        xx_coords(dim,1) = coord(dim,p1)+displ(dim,p1)
        xx_coords(dim,2) = coord(dim,p2)+displ(dim,p2)
        xx_coords(dim,3) = coord(dim,p3)+displ(dim,p3)
     enddo
  endif
end subroutine displace

!
!   Write one STL face
!
subroutine wristl_facet(lun,kfl_stlbo,xx_coords)
  use def_kintyp
  implicit none
  integer(ip),intent(IN) :: lun, kfl_stlbo
  real(rp),intent(IN) :: xx_coords(3,3)
  real(rp) :: normal(3)

  call nortri(xx_coords,normal)
  if( kfl_stlbo == 2 ) then !Binary format
     call wristl_facet_bin(lun,xx_coords,normal)
  else
     call wristl_facet_asc(lun,xx_coords,normal)
  endif

end subroutine wristl_facet

subroutine wristl_facet_bin(lun,coords,normal)
  use def_kintyp
  implicit none

  real(rp),intent(IN) :: coords(3,3),normal(3)
  integer(ip),intent(IN) :: lun
  integer(2),parameter :: ii2 = 0

  write(lun) real(normal(1),4)
  write(lun) real(normal(2),4)
  write(lun) real(normal(3),4)
  write(lun) real(coords(1,1),4)
  write(lun) real(coords(2,1),4)
  write(lun) real(coords(3,1),4)
  write(lun) real(coords(1,2),4)
  write(lun) real(coords(2,2),4)
  write(lun) real(coords(3,2),4)
  write(lun) real(coords(1,3),4)
  write(lun) real(coords(2,3),4)
  write(lun) real(coords(3,3),4)
  write(lun) ii2 ! Garbage at end of STL record

end subroutine wristl_facet_bin

subroutine wristl_facet_asc(lun,coords,normal)
  use def_kintyp
  implicit none

  real(rp),intent(IN) :: coords(3,3),normal(3)
  integer(ip),intent(IN) :: lun

  write(lun,'(a,3(1x,e13.6))') 'facet normal ',normal(1),normal(2),normal(3)
  write(lun,'(a)') 'outer loop'
  write(lun,'(a,3(1x,e13.6))') 'vertex ',coords(1,1),coords(2,1),coords(3,1)
  write(lun,'(a,3(1x,e13.6))') 'vertex ',coords(1,2),coords(2,2),coords(3,2)
  write(lun,'(a,3(1x,e13.6))') 'vertex ',coords(1,3),coords(2,3),coords(3,3)
  write(lun,'(a)') 'endloop'
  write(lun,'(a)') 'endfacet'

end subroutine wristl_facet_asc

subroutine nortri(coord,normal)
!-----------------------------------------------------------------------
!****f* Domain/nortri
! NAME
!    nortri
! DESCRIPTION
!    This routine computes the boundary normals
! USES
! USED BY
!    bounor
!***
!-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  implicit none
  real(rp),    intent(in)  :: coord(3,3)
  real(rp),    intent(out) :: normal(3)
  real(rp)                 :: vec(3,3),absnormal

  vec(1,1) = coord(1,2) - coord(1,1)
  vec(2,1) = coord(2,2) - coord(2,1)
  vec(3,1) = coord(3,2) - coord(3,1)
  vec(1,2) = coord(1,3) - coord(1,1)
  vec(2,2) = coord(2,3) - coord(2,1)
  vec(3,2) = coord(3,3) - coord(3,1)
  call vecpro(vec(:,1),vec(:,2),vec(:,3))
  absnormal = sqrt(vec(1,3)*vec(1,3)+vec(2,3)*vec(2,3)+vec(3,3)*vec(3,3))
  if (absnormal.ne.0.0_rp) then
     normal(1) = -vec(1,3)/absnormal
     normal(2) = -vec(2,3)/absnormal
     normal(3) = -vec(3,3)/absnormal
  else
     print *, "NORMAL WITH ABS EQUAL TO ZERO..."
     normal(1) = 1.0
     normal(2) = 0.0
     normal(3) = 0.0
  endif

end subroutine nortri

subroutine vecpro(v1,v2,v3)

!-----------------------------------------------------------------------
!
! Tthree-dimensional vectorial product of two vectors  v3 = v1 x v2.
! The same pointer as for v1 or v2 may be used for v3.
!
!-----------------------------------------------------------------------
  use      def_kintyp, only : ip,rp
  implicit none
  real(rp),    intent(in)  :: v2(3),v1(3)
  real(rp),    intent(out) :: v3(3)
  real(rp)                 :: c1,c2,c3

  c1=v1(2)*v2(3)-v1(3)*v2(2)
  c2=v1(3)*v2(1)-v1(1)*v2(3)
  c3=v1(1)*v2(2)-v1(2)*v2(1)
  v3(1)=c1
  v3(2)=c2
  v3(3)=c3

end subroutine vecpro
