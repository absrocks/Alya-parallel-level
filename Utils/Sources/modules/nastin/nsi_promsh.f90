subroutine nsi_promsh
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_promsh
  ! NAME 
  !    nsi_promsh
  ! DESCRIPTION
  !    This routine projects data from the old to the new mesh
  !    NOT TESTED
  ! USES
  ! USED BY
  !    nsi_newmsh
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastin
  use      mod_memchk
  implicit none
  integer(ip) :: idime,ipoin,ielem,inode,igaus,icomp,iitsm
  integer(ip) :: izdom,jpoin,pnode,pelty,pgaus
  integer(4)  :: istat
  real(rp)    :: detjm,dvolu,mapre,mavel(ndime)
  real(rp)    :: cartd(ndime,mnode),elcod(ndime,mnode),gpcod(ndime)

  integer(ip) :: jnode,jelem,jfoun
  real(rp)    :: oshap(mnode),oderi(ndime,mnode)
  real(rp)    :: xjaci(ndime,ndime),xjacm(ndime,ndime),deltx(ndime)
  real(rp)    :: delts(ndime),inter(ndime),xbbox(2*ndime),tolel
  !
  tolel=1.0e-4_rp
  !
  ! Allocate lumped mass matrix (vmass), consistent mass matrix (amatr),
  ! generic vector (gevec) and genetric tensor (geten).
  !
  allocate(amatr(nzdom),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'AMATR','nsi_promsh',amatr)

  allocate(vmass(npoin),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'VMASS','nsi_promsh',vmass)

  allocate(gevec(npoin,ncomp_nsi),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'GEVEC','nsi_promsh',gevec)

  allocate(geten(ndime,npoin,ncomp_nsi),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'GETEN','nsi_promsh',geten)

  amatr=0.0_rp
  vmass=0.0_rp
  gevec=0.0_rp
  geten=0.0_rp
  !
  ! Loop over elements
  !
  elements: do ielem=1,nelem

     ! Element properties and dimensions
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)

     ! Gather
     elcod=0.0_rp
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        do idime=1,ndime
           elcod(idime,inode)=coord(idime,ipoin)
        end do
     end do

     gauss_points: do igaus=1,pgaus
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&           ! Cartesian derivative
             elcod,cartd,detjm,xjacm,xjaci)                        ! and Jacobian
        dvolu=elmar(pelty)%weigp(igaus)*detjm

        gpcod=0.0_rp
        do inode=1,nnode(ielem)
           do idime=1,ndime
              gpcod(idime)=gpcod(idime)+ elmar(pelty)%shape(inode,igaus)*elcod(idime,inode)
           end do
        end do

        ! Call to elsest
        !call  Elsest(1,2,0,1,jelem,jfoun,old_npoin,pnode,   &
        !             old_nelem,old_nelty,pnode,ndime,ip,    &
        !             rp,0,0,0,0,0,0,namod(modul),old_coord,old_lnods,&
        !             gpcod,oshap,oderi,xjacm,xjaci,deltx,   &
        !             delts,inter,tolel,xbbox)

        if(jfoun==0) then
           do inode=1,nnode(ielem)
              ipoin=lnods(inode,ielem)
              do icomp=1,ncomp_nsi
                 call runend('NOT CODED')
                 !do jnode=1,old_nnode(jelem)
                    !gevec(ipoin,icomp) = gevec(ipoin,icomp) &
                    !     +elmar(pelty)%shape(inode,igaus)&
                    !     *oshap(jnode)*old_press(jnode,icomp)*dvolu
                    !do idime=1,ndime
                    !   geten(idime,ipoin,icomp) = geten(idime,ipoin,icomp) &
                    !        +elmar(pelty)%shape(inode,igaus)&
                    !        *oshap(jnode)*old_veloc(idime,jnode,icomp)*dvolu
                    !end do
                 !end do
              end do
           end do
        else
           call runend('nsi_promsh: ERROR IN MESH PROJECTION')
        end if

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           vmass(ipoin) = vmass(ipoin) &
                + elmar(pelty)%shape(inode,igaus)*dvolu
        end do

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do jnode=1,pnode
              jpoin=lnods(jnode,ielem)
              do izdom=r_dom(ipoin),r_dom(ipoin+1)-1
                 if(c_dom(izdom)==jpoin)      &
                      amatr(izdom) = amatr(izdom)  &
                      +elmar(pelty)%shape(inode,igaus)&
                      *elmar(pelty)%shape(jnode,igaus)*dvolu
              end do
           end do
        end do


     end do gauss_points

  end do elements
  !
  ! Inverse of lumped mass matrix
  !
  do ipoin=1,npoin
     vmass(ipoin)=1.0_rp/vmass(ipoin)
  end do
  !
  ! New pressure and velocities
  !
  do iitsm = 1,mitsm
     do icomp = 1,ncomp_nsi
        do ipoin = 1,npoin
           mapre=0.0_rp
           mavel=0.0_rp
           do izdom=r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin=c_dom(izdom)
              mapre = mapre + amatr(izdom)*press(jpoin,icomp)
              do idime=1,ndime
                 mavel(idime) = mavel(idime) + amatr(izdom)*veloc(idime,jpoin,icomp)
              end do
           end do
           press(ipoin,icomp)=  vmass(ipoin)*gevec(ipoin,icomp) &
                + press(ipoin,icomp)*(1.0_rp - vmass(ipoin)*mapre)
           do idime = 1,ndime
              veloc(idime,ipoin,icomp) = vmass(ipoin)*geten(idime,ipoin,icomp) &
                   + veloc(idime,ipoin,icomp)*(1.0_rp - vmass(ipoin)*mavel(idime))
           end do
        end do
     end do
  end do
  !
  ! Call to elsest (memory deallocation)
  !call  Elsest(2,2,0,1,jelem,jfoun,old_npoin,pnode,   &
  !             old_nelem,old_nelty,pnode,ndime,ip,    &
  !             rp,0,0,0,0,0,0,namod(modul),old_coord,old_lnods,&
  !             gpcod,oshap,oderi,xjacm,xjaci,deltx,   &
  !             delts,inter,tolel,xbbox)
  !
  ! Deallocate memory
  !
  call memchk(two,istat,mem_modul(1:2,modul),'AMATR','nsi_promsh',amatr)
  deallocate(amatr,stat=istat)
  if(istat.ne.0)  call memerr(two,'AMATR','nsi_promsh',0_ip)

  call memchk(two,istat,mem_modul(1:2,modul),'VMASS','nsi_promsh',vmass)
  deallocate(vmass,stat=istat)
  if(istat.ne.0)  call memerr(two,'VMASS','nsi_promsh',0_ip)

  call memchk(two,istat,mem_modul(1:2,modul),'GEVEC','nsi_promsh',gevec)
  deallocate(gevec,stat=istat)
  if(istat.ne.0)  call memerr(two,'GEVEC','nsi_promsh',0_ip)

  call memchk(two,istat,mem_modul(1:2,modul),'GETEN','nsi_promsh',geten)
  deallocate(geten,stat=istat)
  if(istat.ne.0)  call memerr(two,'GETEN','nsi_promsh',0_ip)

end subroutine nsi_promsh
