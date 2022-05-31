subroutine strfun(vecto)
  !-----------------------------------------------------------------------
  !****f* Mathru/strfun
  ! NAME 
  !    strfun
  ! DESCRIPTION
  !    This routine computes the stream function
  ! USES
  !    strloc
  ! USED BY
  !    nsi_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use mod_memchk
  implicit none
  real(rp), intent(in)     :: vecto(ndime,npoin)
  real(rp)                 :: elvel(ndime,mnode),elcod(ndime,mnode)
  integer(ip)              :: ielem,kpoin,itouc,jnode,iesta,iesto
  integer(ip)              :: jpoin,ipoin,inode,pelty,pnode
  integer(8)               :: metmp(2)=0_8
  integer(4)               :: istat
  integer(ip), allocatable :: mpoin(:)
  integer(ip), allocatable :: nocha(:,:),noinv(:,:)
  real(rp),    allocatable :: westr(:,:)

  if(kfl_paral/=0) then
     !
     ! Allocate memory
     !
     allocate(nocha(16,nelty),stat=istat)
     call memchk(zero,istat,memor_dom,'NOCHA','strfun',nocha)
     allocate(noinv(16,nelty),stat=istat)
     call memchk(zero,istat,memor_dom,'NOINV','strfun',noinv)
     allocate(westr(16,nelty),stat=istat)
     call memchk(zero,istat,memor_dom,'WESTR','strfun',westr)     
     allocate(mpoin(npoin),stat=istat)
     call memchk(zero,istat,memor_dom,'MPOIN','strfun',mpoin)
     if(ndime==2) then
        iesta=10
        iesto=29
     else if(ndime==3) then
        iesta=30
        iesto=50
     end if
     do pelty=iesta,iesto
        if(lexis(pelty)/=0)&
             call chanum(&
             ltopo(pelty),nnode(pelty),nocha(1,pelty),&
             noinv(1,pelty),westr(1,pelty))
     end do
     !
     ! Compute stream function
     !    
     ielem=0
     kpoin=1
     mpoin(1)=1
     gesca(1)=0.0_rp
     do while(kpoin<npoin)
        ielem=mod(ielem+one,nelem)
        if(ielem==0) ielem=nelem
        pelty=ltype(ielem)
        pnode=nnode(pelty)     
        itouc=0
        do jnode=1,pnode
           jpoin=lnods(jnode,ielem)
           if(mpoin(jpoin)>=1) itouc=itouc+1 
        end do
        if(itouc>0.and.itouc<pnode) then 
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              elvel(1:ndime,inode) = vecto(1:ndime,ipoin)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
           call strloc(pnode,ndime,npoin,&
                lnods(1,ielem),mpoin,kpoin,elcod,elvel,gesca,&
                nocha(1,pelty),noinv(1,pelty),westr(1,pelty))
        end if
     end do
     do ipoin=1,npoin
        gesca(ipoin)=gesca(ipoin)/real(mpoin(ipoin),rp)
     end do
     !
     ! Deallocate memory
     !
     call memchk(two,istat,metmp, 'NOCHA','strfun',nocha)
     deallocate(nocha,stat=istat)
     if(istat/=0) call memerr(two,'NOCHA','strfun',0_ip)
     call memchk(two,istat,metmp, 'NOINV','strfun',noinv)
     deallocate(noinv,stat=istat)
     if(istat/=0) call memerr(two,'NOINV','strfun',0_ip)
     call memchk(two,istat,metmp, 'WESTR','strfun',westr)
     deallocate(westr,stat=istat)
     if(istat/=0) call memerr(two,'STREA','strfun',0_ip)
     call memchk(two,istat,metmp, 'MPOIN','strfun',mpoin)
     deallocate(mpoin,stat=istat)
     if(istat/=0) call memerr(two,'MPOIN','strfun',0_ip)

  end if

end subroutine strfun
