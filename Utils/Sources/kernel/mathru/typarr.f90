subroutine typarr(itask,msize,nleng,iarra,itype,memor)
  !-----------------------------------------------------------------------
  !****f* mathru/typarr
  ! NAME
  !    typarr
  ! DESCRIPTION
  !    Converts array to type and type to array:
  !    IARRA(MSIZE,NLENG) <=> ITYPE(NLENG)%L(:)
  !
  !    ITASK=1 ... Converts a type to an array
  !         =2 ... Converts an array to a type (and allocate memory)
  !
  !    For example:
  !    LNODS(MNODE,NELEM) <=> LNODS(NELEM)%L(:)
  ! OUTPUT
  ! USED BY
  !    par_sengeo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use mod_memchk
  implicit none
  integer(ip), intent(in)    :: itask,msize,nleng
  integer(ip), intent(inout) :: iarra(msize,nleng)
  type(i1p),   intent(inout) :: itype(nleng)
  integer(8),  intent(out)   :: memor(2)
  integer(ip)                :: ileng,isize,nsize
  integer(4)                 :: istat

  select case(itask)

  case(1)
     !
     ! Convert an i1P type to an (overdimensioned) array
     !
     do ileng=1,nleng
        nsize=size(itype(ileng)%l,KIND=ip)
        do isize=1,nsize
           iarra(isize,ileng)=itype(ileng)%l(isize)
        end do
        do isize=nsize+1,msize
           iarra(isize,ileng)=0
        end do
     end do

  case(2)
     !
     ! Convert an array to a i1P type 
     !
     do ileng=1,nleng
        nsize=0
        nsize_loop: do while(iarra(nsize+1,ileng)/=0)
           nsize=nsize+1
           if(nsize==msize) exit nsize_loop 
        end do nsize_loop
        if(nsize>0.and.nsize<=msize) then
           allocate(itype(ileng)%l(nsize),stat=istat)
           call memchk(0_ip,istat,memor,'ITYPE(ILENG)%L','typarr',itype(ileng)%l)
           do isize=1,nsize
              itype(ileng)%l(isize)=iarra(isize,ileng)
           end do
        end if
     end do

  end select

end subroutine typarr
