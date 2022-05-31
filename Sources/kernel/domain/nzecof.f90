subroutine nzecof(&
     nlist,ncoef,nzdom,liste,touch)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the number of non-zero coefficients of a
  ! mesh graph stored in compressed sparse row (CSR) format 
  !                 
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  nelem,lnods,mnode,lnnod
  implicit none
  integer(ip), intent(in)    :: nlist,ncoef
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(inout) :: nzdom
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jelem,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kelem,knode,kpoin,nnodk,kposi,klist

  do jlist=1,nlist                                      ! Loop over those elements 
     jelem=liste(jlist)                                 ! where the point is
     nnodj=lnnod(jelem)
     !nnodj=nnode(ltype(jelem))
     do jnode=1,nnodj
        jpoin=lnods(jnode,jelem)
        jposi=(jlist-1)*mnode+jnode
        if(.not.touch(jposi)) then                      ! Position not touched           
           do klist=1,nlist                             ! Search other elements 
              kelem=liste(klist)                        ! where JPOIN is and 
              nnodk=lnnod(kelem)
              !nnodk=nnode(ltype(kelem))
              do knode=1,nnodk                          ! touch their position
                 kpoin=lnods(knode,kelem)
                 if(kpoin==jpoin) then
                    kposi=(klist-1)*mnode+knode
                    touch(kposi)=.true.
                 end if
              end do
           end do
           nzdom = nzdom+1
        end if
     end do
  end do

end subroutine nzecof
