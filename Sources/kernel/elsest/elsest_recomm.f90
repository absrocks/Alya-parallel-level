subroutine elsest_recomm(&
     mnode,ndime,npoin,nelem,nnode,lnods,ltype,&
     coord,recmeth)
  use def_elsest, only      : ip,rp,memor
  use mod_elsest
  implicit none
  integer(ip), intent(in) :: mnode,ndime,npoin,nelem
  integer(ip), intent(in) :: nnode(*),lnods(mnode,nelem),ltype(*)
  integer(ip), intent(out):: recmeth
  real(rp),    intent(in) :: coord(ndime,npoin)
  integer(ip)             :: istat,ielem,inode,ipoin,jnode,jpoin,idime,pnode
  real(rp)                :: dist(4),sedist(3),sesum(3),dista
  real(rp), allocatable   :: semax(:)
  !
  ! Compute characterstic number for anisotrophy of the mesh using max element length
  !
  allocate(semax(nelem),stat=istat)
  dist(:)=0.0_rp
  dist(1)=100000000.0_rp
  sedist(:) = 0
  sesum(:)=0
  !*OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ielem,pnode,inode,ipoin,jnode,jpoin,dista,idime)
  do ielem=1,nelem
     semax(ielem)=-1.0_rp
     pnode=nnode(ltype(ielem))
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        do jnode=inode+1,pnode
           jpoin=lnods(jnode,ielem)
           dista=0.0_rp
           do idime=1,ndime
              dista=dista&
                   +(coord(idime,jpoin)-coord(idime,ipoin))&
                   *(coord(idime,jpoin)-coord(idime,ipoin))
           end do
           !*OMP CRITICAL(sem)
           if(dista>semax(ielem)) semax(ielem)=dista
           !*OMP END CRITICAL(sem)
           !*OMP CRITICAL(d4)
           if(dista>dist(4)) dist(4)=dista
           !*OMP END CRITICAL(d4)
           !*OMP CRITICAL(d1)
           if(dista<dist(1)) dist(1)=dista
           !*OMP END CRITICAL(d1)
        end do
     end do
  end do
  !*OMP END PARALLEL DO
  dist(2)=(dist(4)-dist(1))/3.0_rp
  dist(3)=2.0_rp*(dist(4)-dist(1))/3.0_rp

  recmeth = int((real(nelem,rp)**(0.943_rp/real(ndime,rp)))*((dist(4)/dist(1))**(0.06_rp)))
  deallocate(semax,stat=istat)

end subroutine elsest_recomm
