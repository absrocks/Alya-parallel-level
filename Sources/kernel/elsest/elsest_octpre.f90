subroutine elsest_octpre(&
     ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltype,coord)
  !------------------------------------------------------------------------
  !****f* elsest/elsest_octpre
  ! NAME 
  !    elsest_octpre
  ! DESCRIPTION
  !    Create the quad/oct tree
  ! USES
  ! USED BY
  !    elsest_octpre
  !***
  !------------------------------------------------------------------------
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)         :: ipara(*),imesh,ithre
  integer(ip), intent(in)         :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)         :: nnode(*)
  integer(ip), intent(in)         :: lnods(mnode,*),ltype(*)
  real(rp),    intent(in), target :: coord(ndime,npoin)
  integer(ip)                     :: idime,ipoin,istat,kthre
  real(rp)                        :: time1,time2,time4,time5

  call elsest_cputim(time1)
  oct_struc(imesh) % iallo = 1 
  !
  ! Allocate memory
  !
  !*OMP PARALLEL
  !*OMP SECTIONS
  !*OMP SECTION
  allocate(oct_struc(imesh)%cputi(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'CPUTI','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(oct_struc(imesh)%memor(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'MEMOR','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(oct_struc(imesh)%kstat(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSTAT','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(oct_struc(imesh)%ksear(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSEAR','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(oct_struc(imesh)%kfirs(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KFIRS','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(oct_struc(imesh)%kseco(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSECO','elsest_alloc',0_ip)
  !*OMP END SECTIONS
  !*OMP END PARALLEL
  !
  ! Point to current mesh (IMESH) structure
  !
  call elsest_octpoi(imesh)
  !
  ! Initialize parameters
  !
  !*OMP PARALLEL DO PRIVATE(kthre)
  do kthre=1,nthre
     cputi(:,kthre) = 0.0_rp
     memor(:,kthre) = 0_ip
     kstat(:,kthre) = 0
     kstat(1,kthre) = huge(1_ip)
     kstat(3,kthre) = huge(1_ip)
     ksear(kthre)   = 0
     kfirs(kthre)   = 0
     kseco(kthre)   = 0
  end do
  !*OMP END PARALLEL DO
  memax = 0_ip
  limit = ipara(9)
  !
  ! Compute bounding box
  !
  call elsest_boubox(ndime,npoin,coord,comin,comax)
  !
  ! Allocate memory for tree root
  !
  allocate(oct_struc(imesh) % tree_root,stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'TREE_ROOT','elsest_octpre',0_ip)
  oct_struc(imesh) % tree_root = octbox_init
 
  allocate(oct_struc(imesh) % tree_root % nodes(npoin),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'TREE_ROOT%NODES','elsest_octpre',oct_struc(imesh) % tree_root%nodes)

  allocate(oct_struc(imesh) % elcod(ndime,mnode,nthre),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_octpre',oct_struc(imesh) % elcod)
  !
  ! Init tree_root values
  !
  tree_root                     => oct_struc(imesh) % tree_root
  elcod                         => oct_struc(imesh) % elcod
  current(ithre) % o            => tree_root
  current(ithre) % o % npoinbox =  0
  current(ithre) % o % id       =  0
  current(ithre) % o % level    =  0
  current(ithre) % o % whoiam   =  0
  current(ithre) % o % childid  =  0
  current(ithre) % o % nelembox =  0
  current(ithre) % o % npoinbox =  npoin

  !*OMP PARALLEL DO PRIVATE(ipoin)
  do ipoin = 1,npoin
     current(ithre) % o % nodes(ipoin) = ipoin
  end do
  !*OMP END PARALLEL DO
  do idime = 1,ndime
     current(ithre) % o % minc(idime) = comin(idime)
     current(ithre) % o % maxc(idime) = comax(idime)
  end do
  nullify(current(ithre) % o % parent)
  call elsest_cputim(time2)
  cputi(1,ithre) = time2-time1
  !
  ! Generation of quad/oct tree
  !
  call elsest_octsub(ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode)
  !
  ! Postprocess of quad/oct tree
  !
  call elsest_octpos(ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)

  call elsest_cputim(time4)
  call elsest_cputim(time5)
  cputi(4,ithre)=time5-time4

  if(ipara(7)/=0) call elsest_statis(1_ip,imesh,ipara,ithre)

end subroutine elsest_octpre
