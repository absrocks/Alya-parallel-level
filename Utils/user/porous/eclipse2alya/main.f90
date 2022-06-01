program main_reagr
use def_kintyp
implicit none


  integer(ip)   :: ndi(3),ndime
  real(rp),allocatable     :: coord(:,:,:,:)
  real(rp),allocatable     :: zcorn(:,:,:,:)

  integer(ip)    :: nzerom
  integer(ip)    :: nz1,nz2
  integer(ip)    :: kfl_brugge

  kfl_brugge=1   ! 1 brugge   - 0 cube
  !
  ! For Brugge
  !
  if (kfl_brugge==1) then
     nzerom=204
     nz1=107
     nz2=97
     ndi(1) =139
     ndi(2) =48
     ndi(3) =9
  else
     nzerom=0
     nz1=0
     nz2=0
     ndi(1) =20
     ndi(2) =20
     ndi(3) =5
  end if

  ndime=3

  allocate(coord(ndime,ndi(1)+1_ip,ndi(2)+1_ip,2))
  allocate(zcorn(ndi(1),ndi(2),ndi(3),8))

  call reagrid2(coord,zcorn,ndi,ndime,nzerom,nz1,nz2)


end program main_reagr
