subroutine reastl()
  !-----------------------------------------------------------------------
  !****f* Domain/reastl
  ! NAME
  !    reastl
  ! DESCRIPTION
  !    Allocate the geometry arrays and read or define them
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_memchk
  use mod_iofile
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)          :: ivert,iface,iloca
  integer(ip), pointer :: lnodv(:,:) => null()
  real(rp),    pointer :: exwor(:,:) => null()
  real(rp),    pointer :: coorv(:,:) => null()

!!$  open(unit=77,file='romain.stl',status='old')
!!$  if(kfl_paral==-1.or.(kfl_paral==0.and.kfl_ptask/=2)) then
!!$
!!$     lispa = 0              ! 0 passes through ecoute
!!$     lisda = 77             ! Temporary data file
!!$     lisre = lun_outpu      ! Results file           
!!$     !
!!$     ! Read options and arrays
!!$     ! 
!!$     call ecoute('REASTL')
!!$     do while(words(1)/='SOLID')
!!$        call ecoute('REASTL')
!!$     end do
!!$
!!$     do while(words(1)/='ENDSO')
!!$        call ecoute('reastl')
!!$
!!$        if(words(1)=='FACET') then
!!$           !
!!$           ! Facet normal
!!$           !
!!$           iface=iface+1
!!$           iloca=0
!!$           exwor(1,iface)=param(2)
!!$           exwor(2,iface)=param(3)
!!$           exwor(3,iface)=param(4)
!!$
!!$           do while(words(1)/='OUTER') 
!!$              call ecoute('reastl')
!!$           end do
!!$           do while(words(1)/='ENDLO') 
!!$              if(words(1)=='VERTE') then
!!$                 ivert=ivert+1
!!$                 iloca=iloca+1
!!$                 lnodv(ivert,iface)=ivert
!!$                 coorv(1,ivert)=param(2)
!!$                 coorv(2,ivert)=param(3)
!!$                 coorv(3,ivert)=param(4)
!!$              end if
!!$              call ecoute('reastl')
!!$           end do
!!$        end if
!!$
!!$     end do
!!$
!!$  end if
!!$
!!$  return
!!$
!!$1 call runend('REASTL: ERROR')

end subroutine reastl
