subroutine ensdir(amatr,cmass,elmat,elmas,lnods,pnode,pevat,nunkn)
  !-----------------------------------------------------------------------
  !****f* eigdir/ensdir
  ! NAME
  !    csrase
  ! DESCRIPTION
  !    This routine assemble the matrix with a force brute method
  ! INPUT
  !    ELMAT
  !    NDOFN
  !    PNODE
  !    MEVAT
  !    LNODE
  !    IPROB
  ! OUTPUT
  !    AMATR
  ! USED BY
  !    ***_eigmat
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none 

  integer(ip), intent(in)    :: pnode,pevat,nunkn
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elmat(pevat,pevat),elmas(pevat,pevat)
  real(rp),    intent(inout) :: amatr(nunkn,nunkn),cmass(nunkn,nunkn)
  integer(ip)                :: ipoin,jpoin,inode,jnode

  
  do inode=1,pnode
    ipoin=lnods(inode)  
    do jnode=1,pnode
        jpoin=lnods(jnode)
        amatr(ipoin,jpoin) = amatr(ipoin,jpoin) + elmat(inode,jnode)   
	    cmass(ipoin,jpoin) = cmass(ipoin,jpoin) + elmas(inode,jnode)   
    enddo
  enddo


end subroutine ensdir
