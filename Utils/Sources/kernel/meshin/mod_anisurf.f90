module mod_anisurf

  use mod_cart
  use mod_srftol

contains

  subroutine anisurf(lface,nnofa,nface,ndim,nblay,coor,npoin,rblay,&
       lside,nside,lsmark,tolscal,nnosi,coorold,eltoelold,lfold,nfold,&
       npold,rnofaold,lelemold,lsurfold,nsurf,rsize,rsuni,rnopo,&
       lcell,ncell,lcart,lpsur,lptri,lmark,lfnew,nfnew,lptype,&
       tolcol,tolcolr,tolref,tolrefr,lline,nboup)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)   :: nnofa,ndim,nblay,nnosi,nfold,npold,nsurf,ncell,nboup
    integer(ip), intent(inout):: nface,npoin,nside,nfnew
    real(rp), intent(in)      :: rblay(nblay),tolscal,rsuni 
    integer(ip),pointer       :: lface(:,:),lside(:,:),lfnew(:,:),lline(:)
    integer(ip),pointer       :: lpoin(:),lsmark(:),lmark(:),lptype(:,:)
    real(rp),pointer          :: coor(:,:),rsize(:)
    real(rp),pointer          :: rnopo(:,:)
    integer(4)                :: istat
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold),lsurfold(nfold)
    type(cell),intent(in)     :: lcell(ncell)
    integer(ip),intent(inout) :: lelemold(nfold)
    integer(ip),pointer       :: lpsur(:),lptri(:),lfmark(:),lcart(:)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold)
    real(rp),intent(in)       :: tolcol,tolcolr,tolref,tolrefr
    integer(ip)               :: isurf
    !
    !     This subroutine is responsible of the generation of the whole anisotropic
    !     surface mesh 
    !
    do isurf=1,nsurf

       call anisurfl(lface,nnofa,nface,ndim,nblay,coor,npoin,rblay,&
            lside,nside,lsmark,tolscal,nnosi,coorold,eltoelold,lfold,nfold,&
            npold,rnofaold,lelemold,lsurfold,isurf,rsize,rsuni,rnopo,&
            lcell,ncell,lcart,lpsur,lptri,lfnew,nfnew,lptype,&
            tolcol,tolcolr,tolref,tolrefr,lline,nboup)

    enddo

  end subroutine anisurf

  subroutine anisurfl(lface,nnofa,nface,ndim,nblay,coor,npoin,rblay,&
       lside,nside,lsmark,tolscal,nnosi,coorold,eltoelold,lfold,nfold,&
       npold,rnofaold,lelemold,lsurfold,isurf,rsize,rsuni,rnopo,lcell,ncell,&
       lcart,lpsur,lptri,lfnew,nfnew,lptype,&
       tolcol,tolcolr,tolref,tolrefr,lline,nboup)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)   :: nnofa,ndim,nblay,nnosi,nfold,npold,isurf,ncell,nboup
    integer(ip), intent(inout):: nface,npoin,nside,nfnew
    real(rp), intent(in)      :: rblay(nblay),tolscal,rsuni 
    integer(ip),pointer       :: lface(:,:),lside(:,:),lfnew(:,:),lptype(:,:)
    integer(ip),pointer       :: lpoin(:),lsmark(:),lcart(:),lpsur(:),lline(:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),eltoel(:,:),lplay(:)
    integer(ip),pointer       :: ptoel1new(:),ptoel2new(:),eltoelnew(:,:)
    real(rp),pointer          :: coor(:,:),rsize(:)
    real(rp),pointer          :: rnopo(:,:),rtapo(:,:)
    integer(4)                :: istat
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold)
    integer(ip),intent(in)    :: lfold(nnofa,nfold),lsurfold(nfold)
    type(cell),intent(in)     :: lcell(ncell)
    integer(ip),intent(inout) :: lelemold(nfold)
    integer(ip),pointer       :: lptri(:),lfmark(:),lfsid(:)
    integer(ip),pointer       :: lbsid(:,:),lsconn(:),lfath(:)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold)
    real(rp),intent(in)       :: tolcol,tolcolr,tolref,tolrefr
    integer(ip)               :: nbsid,iface,nsconn,ipoin
    !
    !     This subroutine is responsible of the generation of the anisotropic
    !     surface mesh for surface isurf 
    !
    !     On input:
    !                 - lface: the current mesh
    !                 - lptri: a face number of the current mesh for each point
    !                 - lpsur: a face number of the old mesh for each point
    !                 - lfold: the old mesh
    !                 - lside: the sides of this patch
    !
    !     On output:
    !                 - lface: the current mesh modified
    !                 - lfnew: the anisotropic mesh
    !
    nullify(ptoel1,ptoel2,eltoel,lsconn)

    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','anisurfl',lpoin) 
    allocate(rtapo(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RTAPO','anisurfl',rtapo) 
    allocate(lfath(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LFATH','anisurfl',lfath) 
    !
    !     Initialize lfath
    ! 
    do ipoin=1,npoin
       lfath(ipoin)=ipoin
    enddo 
    !
    !     Output lface && lside
    !
    !call outaniface(nnofa,nface,npoin,ndim,lface,coor,lside,nnosi,nside,rnopo)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    !
    !     Reorder the sides
    !
    call reorside(nside,lside,nnosi,eltoel,lface,nnofa,nface,npoin,&
         lline,ptoel1,ptoel2,lsmark,nsconn,lsconn,coor,rnopo,rtapo,&
         lpoin,lpsur,lptri,lcart,lptype,lfath)
    !
    !     Get the wetted sides and point normals
    !
    call mptside(lside,lsmark,nside,npoin,nnosi,&
         tolscal,coor,ndim,rnopo,rtapo,lpoin,lpsur,&
         lcart,lptri,lptype,lline,nsconn,lsconn,lfath)
    !
    !     Allocate lplay
    !   
    allocate(lplay(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPLAY','anisurfl',lplay) 
    !
    !     Generate the surface mesh
    !
    call gtanisurf(nfnew,nnofa,lfnew,lside,lsmark,nside,npoin,lpoin,&
         coor,nnosi,rblay,nblay,rtapo,lplay,lfsid)
    !
    !     Resize lpoin
    !
    call memrea(npoin,memor_msh,'LPOIN','anisurfl',lpoin)
    call memrea(npoin,memor_msh,'LCART','anisurfl',lcart)
    call memrea(npoin,memor_msh,'LPSUR','anisurfl',lpsur)
    call memrea(npoin,memor_msh,'LPTRI','anisurfl',lptri)
    call memrea(npoin,memor_msh,'RSIZE','anisurfl',rsize)
    call memrea(npoin,memor_msh,'RNOPO','anisurfl',rnopo)
    call memrea(npoin,memor_msh,'LPTYPE','anisurfl',lptype)
    !
    !     Output lface && lside
    !
    call outaniface(nnofa,nfnew,npoin,ndim,lfnew,coor,lside,nnosi,nside,rtapo)
    !
    !     Advance one layer at a time and check intersection
    !
    call aniadv(nfnew,nnofa,lfnew,lside,lsmark,nside,npoin,ndim,nblay,&
         coor,rnopo,rblay,coorold,eltoelold,lfold,nfold,&
         npold,rnofaold,lelemold,lsurfold,isurf,&
         rsize,rsuni,lcell,ncell,lcart,lpsur,lptri,lfmark,&
         rtapo,lpoin,lfsid,nnosi,lfath)
    !
    !     Output the mesh
    !   
    !call outaniface3(nnofa,nfnew,npoin,ndim,lfnew,coor,lfmark)
    !
    !     Take out one more level
    !
    call onemore(nnofa,nfnew,lfnew,lfmark,lplay,npoin)
    !
    !     Output the mesh
    !   
    !call outaniface3(nnofa,nfnew,npoin,ndim,lfnew,coor,lfmark)
    !
    !     Allow for one level jump only
    !
    call rmvjmp(lfnew,nfnew,npoin,nnofa,lplay,lfmark,nside,nblay)
    !
    !     Output the mesh
    !   
    !call outaniface3(nnofa,nfnew,npoin,ndim,lfnew,coor,lfmark)
    !
    !     Add elements to make a smooth transition
    !
    call addcover(lfnew,nfnew,nnofa,lfmark,lplay,npoin)
    !
    !     Output the mesh
    !   
    call outaniface3(nnofa,nfnew,npoin,ndim,lfnew,coor,lfmark)
    !
    !     Compress lfnew
    !
    call compresf(lfnew,nnofa,nfnew,lfmark)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lfnew,nfnew,npoin,nnofa,ptoel1new,ptoel2new)
    !
    !     Get the faces surrounding faces
    !
    call trtotr(lfnew,nnofa,nfnew,ptoel1new,ptoel2new,npoin,eltoelnew)
    !
    !     Get the bl boundary
    !
    call blbound(lfnew,nfnew,nnofa,npoin,lbsid,nbsid,eltoelnew)
    !
    !     Force the bl boundary
    !
    call forcbl(lface,nnofa,nface,lbsid,nbsid,eltoel,lside,nside,&
         npoin,coor,ndim,lptri,lsmark,lpoin,isurf,rnopo,lcell,ncell,&
         lcart,rsize,rsuni)
    !
    !     Take out the old mesh delimited by the bl mesh boundaries
    ! 
    call cleansurf(lface,nnofa,nface,npoin,lside,nside,nnosi,eltoel,&
         coor,ndim,lbsid,nbsid,lptri)
    !
    !     Collapse the short edges and optimize the mesh
    !   
    call anicol(coor,ndim,npoin,lface,nface,nnofa,lptype,ptoel1,ptoel2,&
         rsize,rnopo,eltoel,nbsid,&
         nnosi,tolcol,tolcolr,tolref,tolrefr,lbsid,lsmark,lptri,nboup)
    !
    !     Resize lface
    !
    !call memrea(nface+nfnew,memor_msh,'LFACE','anisurfl',lface)
    !
    !     Add the new mesh to the old mesh
    !
    !do iface=1,nfnew
    !   nface=nface+1
    !   lface(1,nface)=lfnew(1,iface)
    !   lface(2,nface)=lfnew(2,iface)
    !   lface(3,nface)=lfnew(3,iface)
    !enddo
    !
    !     Output the mesh
    !
    !call outaniface2(nnofa,nface,npoin,ndim,lface,coor)
    !
    !     Compress the points
    !
    !call compresp(lface,nnofa,nface,lfnew,nfnew,npoin,coor,lptri,ndim,lpoin)

    call memchk(2_ip,istat,memor_msh,'LPLAY','anisurf',lplay)
    deallocate(lplay,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPLAY','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'RTAPO','anisurf',rtapo)
    deallocate(rtapo,stat=istat)
    if(istat/=0) call memerr(2_ip,'RTAPO','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPOIN','anisurf',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOELNEW','anisurfl',eltoelnew)
    deallocate(eltoelnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOELNEW','anisurfl',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','anisurfl',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','anisurfl',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','anisurf',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','anisurf',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1NEW','anisurf',ptoel1new)
    deallocate(ptoel1new,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1NEW','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2NEW','anisurf',ptoel2new)
    deallocate(ptoel2new,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2NEW','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFSID','anisurf',lfsid)
    deallocate(lfsid,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFSID','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSCONN','anisurf',lsconn)
    deallocate(lsconn,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSCONN','anisurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFATH','anisurf',lfath)
    deallocate(lfath,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFATH','anisurf',0_ip)

  end subroutine anisurfl

  subroutine mptside(lside,lsmark,nside,npoin,nnosi,&
       tolscal,coor,ndim,rnopo,rtapo,lpoin,lpsur,&
       lcart,lptri,lptype,lline,nsconn,lsconn,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use mod_mshtol, only        :  ptoelm,sitosa
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)    :: nnosi,ndim,nsconn
    integer(ip), intent(inout) :: npoin,nside
    integer(ip), intent(inout) :: lsconn(nsconn+1)
    real(rp), intent(in)       :: tolscal
    real(rp), pointer          :: coor(:,:),rnopo(:,:),rtapo(:,:)
    integer(ip),pointer        :: lside(:,:),lsmark(:),lptype(:,:),lline(:)
    integer(ip),pointer        :: lpoin(:),lpsur(:),lcart(:),lptri(:),lfath(:)
    integer(ip),pointer        :: lsnew(:,:),llinenew(:),lsmnew(:),lsconnew(:)
    integer(ip)                :: iside,ip1,ip2,inosi,jside,jp1,jp2
    integer(ip)                :: nsid0,ipoin,ineigh,isconn
    real(rp)                   :: rix,riy,riz,rjx,rjy,rjz,ril,rjl 
    real(rp)                   :: rtx,rty,rtz,rtl
    real(rp)                   :: rnx,rny,rnz,rnl
    real(rp)                   :: rnx2,rny2,rnz2,rnl2
    real(rp)                   :: cscal,c10 
    integer(4)                 :: istat
    !
    !     This sub gets point normals and generates virtual sides
    !
    c10=1.0d+00
    !
    !     Allocate temporary arrays
    !
    allocate(llinenew(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','mptside',llinenew) 
    allocate(lsnew(2,nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSNEW','mptside',lsnew) 
    allocate(lsmnew(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSMNEW','mptside',lsmnew) 
    allocate(lsconnew(nsconn+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LSCONNEW','mptside',lsconnew) 
    !
    !     Loop on side to find virtual sides
    !
    nsid0=nside
    nside=0_ip
    !
    !     Loop on connex component
    !
    do isconn=1,nsconn
       !
       !     First loop on sides of this component
       !
       do iside=lsconn(isconn),lsconn(isconn+1_ip)-2_ip
          !
          !     Copy in new structure
          !
          nside=nside+1  
          call memrea(nside,memor_msh,'LSNEW','mptside',lsnew)
          call memrea(nside,memor_msh,'LSMNEW','mptside',lsmnew)
          call memrea(nside,memor_msh,'LLINENEW','mptside',llinenew)

          lsnew(1,nside)=lside(1,iside)
          lsnew(2,nside)=lside(2,iside)
          lsmnew(nside)=lsmark(iside)
          llinenew(nside)=lline(iside)
          !
          !     Is iside a wetted side?
          ! 
          if(lsmark(iside)==1)then 
             !
             !     Get points of the side
             !
             ip1=lside(1,iside)
             ip2=lside(2,iside)
             !
             !     Get side vector
             !
             rix=coor(1,ip2)-coor(1,ip1)
             riy=coor(2,ip2)-coor(2,ip1)
             riz=coor(3,ip2)-coor(3,ip1)
             ril=sqrt(rix*rix+riy*riy+riz*riz)
             ril=c10/ril
             rix=ril*rix
             riy=ril*riy
             riz=ril*riz
             !
             !     Take next side 
             !
             jside=iside+1
             !
             !     Is jside a wetted side?
             ! 
             if(lsmark(jside)==1)then
                !
                !     Get points of the side
                !
                jp1=lside(1,jside)
                jp2=lside(2,jside)
                !
                !     Get side vector
                !
                rjx=coor(1,jp2)-coor(1,jp1)
                rjy=coor(2,jp2)-coor(2,jp1)
                rjz=coor(3,jp2)-coor(3,jp1)
                rjl=sqrt(rjx*rjx+rjy*rjy+rjz*rjz)
                rjl=c10/rjl
                rjx=rjl*rjx
                rjy=rjl*rjy
                rjz=rjl*rjz
                !
                !     Compute scalar product
                !
                cscal=rix*rjx+riy*rjy+riz*rjz
                !
                !     And the test 
                !
                if(cscal<tolscal)then
                   !
                   !     Check for concave/convex and add normals
                   !           
                   call chkconv(nside,lside,iside,jside,coor,rix,riy,riz,&
                        rjx,rjy,rjz,npoin,rnopo,nnosi,ip2,&
                        ndim,cscal,rtapo,lpoin,lpsur,lcart,&
                        lptri,lptype,lline,lsnew,lsmnew,llinenew,lfath)

                else
                   !
                   !     Get average tangent
                   !
                   rtx=rix+rjx 
                   rty=riy+rjy 
                   rtz=riz+rjz 
                   rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
                   rtl=c10/rtl
                   rtx=rtx*rtl
                   rty=rty*rtl
                   rtz=rtz*rtl
                   !
                   !     Get face normal
                   !
                   rnx=rnopo(1,ip2)
                   rny=rnopo(2,ip2)
                   rnz=rnopo(3,ip2)
                   !
                   !     Compute vector product
                   !
                   rnx2= rny*rtz-rnz*rty
                   rny2=-rnx*rtz+rnz*rtx
                   rnz2= rnx*rty-rny*rtx
                   rnl2=sqrt(rnx2*rnx2+rny2*rny2+rnz2*rnz2)
                   rnl2=c10/rnl2
                   rnx2=rnx2*rnl2
                   rny2=rny2*rnl2
                   rnz2=rnz2*rnl2
                   rtapo(1,ip2)=rnx2
                   rtapo(2,ip2)=rny2
                   rtapo(3,ip2)=rnz2

                endif

             endif

          endif

       enddo
       !
       !     Then handle the last and first side
       !
       iside=lsconn(isconn+1)-1_ip
       jside=lsconn(isconn)
       !
       !     Copy in new structure
       !
       nside=nside+1  
       lsnew(1,nside)=lside(1,iside)
       lsnew(2,nside)=lside(2,iside)
       lsmnew(nside)=lsmark(iside)
       llinenew(nside)=lline(iside)
       !
       !     Is iside a wetted side?
       ! 
       if(lsmark(iside)==1)then 
          !
          !     Get points of the side
          !
          ip1=lside(1,iside)
          ip2=lside(2,iside)
          !
          !     Get side vector
          !
          rix=coor(1,ip2)-coor(1,ip1)
          riy=coor(2,ip2)-coor(2,ip1)
          riz=coor(3,ip2)-coor(3,ip1)
          ril=sqrt(rix*rix+riy*riy+riz*riz)
          ril=c10/ril
          rix=ril*rix
          riy=ril*riy
          riz=ril*riz
          !
          !     Is jside a wetted side?
          ! 
          if(lsmark(jside)==1)then
             !
             !     Get points of the side
             !
             jp1=lside(1,jside)
             jp2=lside(2,jside)
             !
             !     Get side vector
             !
             rjx=coor(1,jp2)-coor(1,jp1)
             rjy=coor(2,jp2)-coor(2,jp1)
             rjz=coor(3,jp2)-coor(3,jp1)
             rjl=sqrt(rjx*rjx+rjy*rjy+rjz*rjz)
             rjl=c10/rjl
             rjx=rjl*rjx
             rjy=rjl*rjy
             rjz=rjl*rjz
             !
             !     Compute scalar product
             !
             cscal=rix*rjx+riy*rjy+riz*rjz
             !
             !     And the test 
             !
             if(cscal<tolscal)then
                !
                !     Check for concave/convex and add normals
                !           
                call chkconv(nside,lside,iside,jside,coor,rix,riy,riz,&
                     rjx,rjy,rjz,npoin,rnopo,nnosi,ip2,&
                     ndim,cscal,rtapo,lpoin,lpsur,lcart,&
                     lptri,lptype,lline,lsnew,lsmnew,llinenew,lfath)
                !
                !     Modify lsnew if lside(1,iside) has been modified
                ! 
                if(lside(1,jside)/=jp1)then
                   lsnew(1,lsconnew(isconn))=lside(1,jside)
                endif
 
             else
                !
                !     Get average tangent
                !
                rtx=rix+rjx 
                rty=riy+rjy 
                rtz=riz+rjz 
                rtl=sqrt(rtx*rtx+rty*rty+rtz*rtz)
                rtl=c10/rtl
                rtx=rtx*rtl
                rty=rty*rtl
                rtz=rtz*rtl
                !
                !     Get face normal
                !
                rnx=rnopo(1,ip2)
                rny=rnopo(2,ip2)
                rnz=rnopo(3,ip2)
                !
                !     Compute vector product
                !
                rnx2= rny*rtz-rnz*rty
                rny2=-rnx*rtz+rnz*rtx
                rnz2= rnx*rty-rny*rtx
                rnl2=sqrt(rnx2*rnx2+rny2*rny2+rnz2*rnz2)
                rnl2=c10/rnl2
                rnx2=rnx2*rnl2
                rny2=rny2*rnl2
                rnz2=rnz2*rnl2
                rtapo(1,ip2)=rnx2
                rtapo(2,ip2)=rny2
                rtapo(3,ip2)=rnz2

             endif

          endif

       endif
       !
       !     Remember for this component
       !
       lsconnew(isconn+1)=nside+1

    enddo

    call memrea(nside,memor_msh,'lside','mptside',lside)
    call memrea(nside,memor_msh,'lline','mptside',lline)
    call memrea(nside,memor_msh,'lsmark','mptside',lsmark)
    !
    !     Transfer in lside
    !
    do iside=1,nside
       lside(1,iside)=lsnew(1,iside)
       lside(2,iside)=lsnew(2,iside)
       lline(iside)=llinenew(iside)
       lsmark(iside)=lsmnew(iside)
    enddo

    do isconn=2,nsconn+1
       lsconn(isconn)=lsconnew(isconn)
    enddo

    call memchk(2_ip,istat,memor_msh,'LSCONNEW','mptside',lsconnew)
    deallocate(lsconnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSCONNEW','mptside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LLINENEW','mptside',llinenew)
    deallocate(llinenew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LLINENEW','mptside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSNEW','mptside',lsnew)
    deallocate(lsnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSNEW','mptside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSMNEW','mptside',lsmnew)
    deallocate(lsmnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSMNEW','mptside',0_ip)

  end subroutine mptside

  subroutine chkconv(nside,lside,iside,jside,coor,rix,riy,riz,rjx,rjy,rjz,&
       npoin,rnopo,nnosi,ipoin,ndim,cscal,rtapo,lpoin,lpsur,&
       lcart,lptri,lptype,lline,lsnew,lsmnew,llinenew,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)    :: iside,jside,nnosi,ipoin,ndim
    integer(ip), intent(inout) :: npoin,nside
    integer(ip), pointer       :: lptype(:,:)
    integer(ip), intent(inout) :: lside(2,nside)
    integer(ip), intent(in)    :: lline(nside)
    integer(ip), pointer       :: lpoin(:),lpsur(:),lcart(:),lptri(:)
    integer(ip), pointer       :: lsnew(:,:),lsmnew(:),llinenew(:),lfath(:)
    real(rp), intent(in)       :: rix,riy,riz,rjx,rjy,rjz,cscal
    real(rp), pointer          :: coor(:,:),rnopo(:,:),rtapo(:,:)
    integer(ip)                :: ipa,ipb,ipc,iside1,iside2,nsubd
    integer(ip)                :: ipbegin,isubd,npnew,ip1,ip2,nsnew
    integer(ip)                :: iline
    real(rp)                   :: rtx1,rty1,rtz1,rtx2,rty2,rtz2,c10,rnl2  
    real(rp)                   :: rnx,rny,rnz,rnx2,rny2,rnz2,rnl,rtl,csca  
    real(rp)                   :: c00,phi,cm1,phito,dphi,cosan,sinan,phian  
    real(rp)                   :: rnx1,rny1,rnz1,rnx4,rny4,rnz4,rtol,rnl1  
    real(rp)                   :: rpx,rpy,rpz  
    real(rp)                   :: p1(3),p2(3)  
    !
    !     This sub checks if the angle formed by iside and jside is convex
    !     We suppose that ri and rj are normalized
    !
    c00=0.0d+00
    c10=1.0d+00
    cm1=-1.0d+00
    rtol=0.2d+00
    phito=35.0*3.1415927/180.  ! A normal every 35 degrees
    !
    !     Get the correct orientation 
    !
    ipa=lside(1,iside)
    ipb=lside(2,iside)
    ipc=lside(2,jside)

    rtx1=rix 
    rty1=riy 
    rtz1=riz 
    rtx2=rjx 
    rty2=rjy 
    rtz2=rjz 

    iside1=iside
    iside2=jside
    !
    !     Get surface normal
    !
    rnx=rnopo(1,ipa)+rnopo(1,ipb)
    rny=rnopo(2,ipa)+rnopo(2,ipb)
    rnz=rnopo(3,ipa)+rnopo(3,ipb)
    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl
    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl
    !
    !     Get binormal to first side
    ! 
    rnx2= rny*rtz1-rnz*rty1
    rny2=-rnx*rtz1+rnz*rtx1
    rnz2= rnx*rty1-rny*rtx1
    rnl2=sqrt(rnx2*rnx2+rny2*rny2+rnz2*rnz2)
    rnl2=c10/rnl2
    rnx2=rnx2*rnl2
    rny2=rny2*rnl2
    rnz2=rnz2*rnl2
    !
    !     Check visibility
    !
    csca=rnx2*rtx2+rny2*rty2+rnz2*rtz2
    !
    !     Do we have a virtual point?
    !
    if(lfath(ipa)==lfath(ipc))then 
       csca=c00
    endif
    !
    !     Do we have a convex corner?
    !
    if(csca<=c00)then
       !
       !     Get the angle
       !
       phi=acos(max(cm1,min(c10,cscal)))
       !
       !      Nr. of inner subdivisions
       !
       nsubd=1+int(phi/phito)
       !
       !     Reallocate coor
       !
       npnew=npoin+nsubd
       call memrea(npnew,memor_msh,'COOR','chkconv',coor)
       call memrea(npnew,memor_msh,'RNOPO','chkconv',rnopo)
       call memrea(npnew,memor_msh,'RTAPO','chkconv',rtapo)
       call memrea(npnew,memor_msh,'LPOIN','chkconv',lpoin)
       call memrea(npnew,memor_msh,'LPSUR','chkconv',lpsur)
       call memrea(npnew,memor_msh,'LTRI','chkconv',lptri)
       call memrea(npnew,memor_msh,'LCART','chkconv',lcart)
       call memrea(npnew,memor_msh,'LPTYPE','chkconv',lptype)
       call memrea(npnew,memor_msh,'LFATH','chkconv',lfath)
       !
       !     Reallocate lsnew
       !
       nsnew=nside+nsubd
       call memrea(nsnew,memor_msh,'LSNEW','chkconv',lsnew)
       call memrea(nsnew,memor_msh,'LSMNEW','chkconv',lsmnew)
       call memrea(nsnew,memor_msh,'LLINENEW','chkconv',llinenew)
       !
       !     Angle for each
       !
       dphi=phi/real(nsubd)

       ipbegin=ipb
       rtapo(1,ipbegin)=rnx2
       rtapo(2,ipbegin)=rny2
       rtapo(3,ipbegin)=rnz2
       iline=lline(iside2)
       !
       !     Initialize phian
       !
       phian=dphi
       !
       !     Add sides and normals 
       !
       do isubd=1,nsubd 

          npoin=npoin+1
          coor(1,npoin)=coor(1,ipb)
          coor(2,npoin)=coor(2,ipb)
          coor(3,npoin)=coor(3,ipb)
          cosan=cos(phian)
          sinan=sin(phian) 
          rtapo(1,npoin)=cosan*rnx2+sinan*rtx1
          rtapo(2,npoin)=cosan*rny2+sinan*rty1
          rtapo(3,npoin)=cosan*rnz2+sinan*rtz1
          lpsur(npoin)=lpsur(ipb)
          lptri(npoin)=lptri(ipb)
          lcart(npoin)=lcart(ipb)
          rnopo(1,npoin)=rnopo(1,ipb)
          rnopo(2,npoin)=rnopo(2,ipb)
          rnopo(3,npoin)=rnopo(3,ipb)
          lptype(1,npoin)=lptype(1,ipb)
          lptype(2,npoin)=lptype(2,ipb)
          lfath(npoin)=ipb
          phian=phian+dphi

          nside=nside+1
          lsnew(1,nside)=ipbegin
          lsnew(2,nside)=npoin
          lsmnew(nside)=2_ip
          llinenew(nside)=iline

          ipbegin=npoin

       enddo
       !
       !     Modify iside2
       !     lside will be copied in lsnew if not last edge
       !
       lside(1,iside2)=npoin

    else
       !
       !     We have a concave corner
       !  
       !
       !     Compute point normal
       ! 
       rnx=rtx2-rtx1
       rny=rty2-rty1
       rnz=rtz2-rtz1
       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       !
       !     Are the sides aligned?
       !
       if(rnl<rtol)then
          !
          !     Take the face normal
          !   
          rnx1=rnopo(1,ipoin) 
          rny1=rnopo(2,ipoin) 
          rnz1=rnopo(3,ipoin) 
          !
          !     Compute vectorial product
          ! 
          rnx= rny1*rtz1-rnz1*rty1
          rny=-rnx1*rtz1+rnz1*rtx1
          rnz= rnx1*rty1-rny1*rtx1
          rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
          rnl=c10/rnl
          rnx=rnx*rnl
          rny=rny*rnl
          rnz=rnz*rnl

       else
          !
          !     Take this normal as bisectrice
          !
          rnl=c10/rnl
          rnx=rnx*rnl
          rny=rny*rnl
          rnz=rnz*rnl

       endif
       !
       !     Compute concavity
       !
       csca=rnx*rtx2+rny*rty2+rnz*rtz2
       csca=c10/csca


       rtapo(1,ipoin)=rnx*csca     
       rtapo(2,ipoin)=rny*csca     
       rtapo(3,ipoin)=rnz*csca

    endif

  end subroutine chkconv

  subroutine gtanisurf(nface,nnofa,lface,lside,lsmark,nside,npoin,lpoin,&
       coor,nnosi,rblay,nblay,rtapo,lplay,lfsid)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)    :: nside,nnofa,nblay,nnosi
    integer(ip), intent(in)    :: lside(2,nside),lsmark(nside)
    integer(ip), intent(inout) :: nface,npoin
    integer(ip),pointer        :: lpoin(:),lface(:,:),lplay(:),lfsid(:)
    real(rp),intent(in)        :: rblay(nblay)
    real(rp),pointer           :: coor(:,:),rtapo(:,:)
    integer(ip)                :: ncont,npnew,iside,ib1,ib2,it1,it2
    integer(ip)                :: ip1,ip2,ilay,inosi,ipoin
    real(rp)                   :: reasont,x0,y0,z0,thickness
    integer(4)                 :: istat
    !
    !     This sub creates the new mesh with approximate coordinates 
    !     of the point
    !   
    !
    !     Clean up lpoin
    !
    do ipoin=1,npoin
       lpoin(ipoin)=0_ip
    enddo
    !
    !     Count the number of new elements and new points
    ! 
    ncont=0_ip
    npnew=npoin
    do iside=1,nside
       if(lsmark(iside)>0)then
          ncont=ncont+1_ip
          do inosi=1,nnosi
             ipoin=lside(inosi,iside)
             if(lpoin(ipoin)==0)then
                lpoin(ipoin)=1_ip
                npnew=npnew+nblay 
             endif
          enddo
       endif
    enddo
    !
    !     Reallocate lpoin && coor
    !
    call memrea(npnew,memor_msh,'LPOIN','gtanisurf',lpoin)
    do ipoin=1,npnew
       lpoin(ipoin)=0
    enddo
    call memrea(npnew,memor_msh,'COOR','gtanisurf',coor)
    call memrea(npnew,memor_msh,'RNOPO','gtanisurf',rtapo)
    call memrea(npnew,memor_msh,'LPLAY','gtanisurf',lplay)
    !
    !     Get the pointer
    ! 
    ncont=0_ip
    !
    !     Loop on sides
    !
    do iside=1,nside
       if(lsmark(iside)>0)then
          ncont=ncont+1_ip
          do inosi=1,nnosi
             ipoin=lside(inosi,iside)
             if(lpoin(ipoin)==0)then
                lpoin(ipoin)=npoin+1
                lplay(ipoin)=0_ip
                reasont=1.0d+00
                !
                !     Build the point position
                !
                x0=coor(1,ipoin)
                y0=coor(2,ipoin)
                z0=coor(3,ipoin)
                !
                !     Loop on layers
                !
                do ilay=1,nblay

                   thickness=rblay(ilay)
                   npoin=npoin+1
                   coor(1,npoin)=x0+rtapo(1,ipoin)*thickness
                   coor(2,npoin)=y0+rtapo(2,ipoin)*thickness
                   coor(3,npoin)=z0+rtapo(3,ipoin)*thickness
                   lplay(npoin)=ilay
                   lpoin(npoin)=ipoin

                   x0=coor(1,npoin)
                   y0=coor(2,npoin)
                   z0=coor(3,npoin)
                enddo
             endif
          enddo
       endif
    enddo
    !
    !     Get new face number
    !
    nface=2*ncont*nblay
    !
    !     Allocate lface
    !
    if(.not.associated(lface))then
       allocate(lface(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LFACE','gtanisurf',lface)
    else
       call memrea(nface,memor_msh,'LFACE','gtanisurf',lface)
    endif
    !
    !     Allocate lfsid
    !    
    allocate(lfsid(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFSID','gtanisurf',lfsid)
    !
    !     Loop on marked sides
    !
    nface=0_ip
    do iside=1,nside

       if(lsmark(iside)>0)then

          ib1=lside(1,iside)
          ib2=lside(2,iside)
          !
          !     Get the pointer
          !
          it1=lpoin(ib1)
          it2=lpoin(ib2)
          !
          !     Loop on layers
          !
          do ilay=1,nblay

             nface=nface+1_ip 
             lface(1,nface)=ib1
             lface(2,nface)=ib2
             lface(3,nface)=it1
             lfsid(nface)=iside

             nface=nface+1_ip 
             lface(1,nface)=it1
             lface(2,nface)=ib2
             lface(3,nface)=it2
             lfsid(nface)=iside

             ib1=it1
             ib2=it2

             it1=it1+1_ip
             it2=it2+1_ip

          enddo
       endif
    enddo

  end subroutine gtanisurf

  subroutine aniadv(nface,nnofa,lface,lside,lsmark,nside,npoin,ndim,&
       nblay,coor,rnopo,rblay,coorold,eltoelold,lfold,nfold,&
       npold,rnofaold,lelemold,lsurfold,isurf,&
       rsize,rsuni,lcell,ncell,lcart,lpsur,lptri,lfmark,&
       rtapo,lpoin,lfsid,nnosi,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip), intent(in)   :: nside,nnofa,ndim,nblay,nfold,npold,isurf
    integer(ip), intent(in)   :: nface,npoin,ncell,nnosi
    real(rp), intent(inout)   :: coor(ndim,npoin)
    real(rp), intent(in)      :: rblay(nblay),rsuni
    real(rp), intent(inout)   :: rnopo(ndim,npoin),rsize(npoin)
    real(rp), intent(inout)   :: rtapo(ndim,npoin)
    integer(ip), intent(in)   :: lside(nnosi,nside),lfsid(nface)
    integer(ip), intent(in)   :: lface(nnofa,nface),lpoin(npoin),lfath(npoin)
    integer(ip),intent(inout) :: lsmark(nside),lcart(npoin),lpsur(npoin)
    integer(ip),intent(inout) :: lptri(npoin)
    integer(ip),intent(in)    :: eltoelold(nnofa,nfold),lfold(nnofa,nfold)
    integer(ip),intent(in)    :: lsurfold(nfold)
    integer(ip),intent(inout) :: lelemold(nfold)
    type(cell),intent(in)     :: lcell(ncell) 
    type(cell),pointer        :: lcellp(:) 
    integer(ip),pointer       :: lcartp(:),lfmark(:),lmarkc(:)
    integer(ip),pointer       :: lctop(:,:),lpmark(:),sitosi(:,:)
    integer(ip),pointer       :: ptoel1(:),ptoel2(:),lstack(:)
    integer(ip),pointer       :: ptosi1(:),ptosi2(:),lmark(:)
    real(rp),intent(in)       :: coorold(ndim,npold),rnofaold(ndim,nfold)
    integer(ip)               :: ilay,iside,ip1,ip2,jlay,ipb1,ipb2,ipt1,ipt2
    integer(ip)               :: ihost,ip1o,ip2o,ip3o,icart,ielem1,npclos,nsclos
    integer(ip)               :: npmax,ncellp,nfac,npmaxm1,nfacl,iface,ierr
    real(rp)                  :: rnx,rny,rnz,rlthick,d1,d2,d3,rsloc,rnofa(3) 
    real(rp)                  :: rnl,qgeo2,q1,c05,qgeotol,rtol,bbox(ndim,2)
    integer(ip),parameter     :: mstackp=500
    integer(ip)               :: lstackp(mstackp),nstackp,mcell    
    integer(4)                :: istat
    !
    !     This sub checks intersection one layer at a time and create the points 
    !
    npmax=9_ip
    npmaxm1=npmax-1_ip
    c05=0.5d+00
    qgeotol=0.5d+00
    nullify(ptoel1,ptoel2,ptosi1,ptosi2,sitosi)
    ierr=0_ip
    !
    !     Allocate lfmark
    !    
    allocate(lfmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFMARK','aniadv',lfmark)
    !
    !     Set all the faces to be kept
    !
    do iface=1,nface
       lfmark(iface)=1_ip
    enddo
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the sides surrounding the points
    !
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2 )
    !
    !     Get the sides surrounding sides
    !
    call sitosa(lside,nside,npoin,nnosi,ptosi1,ptosi2,sitosi)
    !
    !     Allocate lcartp
    !    
    allocate(lcartp(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LCART','aniadv',lcartp)
    allocate(lpmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPMARK','aniadv',lpmark)
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','aniadv',lmark)
    !
    !     First, insert all the points in the cartesian mesh
    ! 
    call cartp(npoin,ndim,lface,nnofa,nface,npmaxm1,lcartp,ptoel1,ptoel2,coor,&
         lcellp,ncellp,rtol,0_ip,bbox)
    !
    !     Allocate arrays for lcellp
    !
    allocate(lctop(npmax,ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LCTOP','aniadv',lctop) 
    allocate(lmarkc(ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARKC','aniadv',lmarkc) 
    allocate(lstack(ncellp),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','aniadv',lstack) 
    mcell=ncell
    !
    !     Build the cell to old point pointers
    !
    call ctopnt2(lcartp,npoin,ncellp,lctop,npmax)
    !
    !     Temporaryyyyyyyyyyyyyyy
    !
    nsclos=0_ip
    !
    !     Outer loop on layer
    !
    do ilay=1,nblay
       !
       !     Set face pointer
       !
       nfac=(ilay-1)*2-2*nblay+1_ip
       !
       !     Inner loop on sides
       !
       do iside=1,nside

          !write(*,*)ilay,iside,lpsur(17)
          !
          !     Is it a wetted side?
          !
          if(lsmark(iside)>0)then
             !
             !     Get points on side  
             !
             ip1=lside(1,iside)
             ip2=lside(2,iside)
             ipb1=ip1
             ipb2=ip2
             ipt1=lpoin(ip1)          
             ipt2=lpoin(ip2) 
             !
             !     Get old side
             !         
             do jlay=2,ilay

                ipb1=ipt1
                ipb2=ipt2
                ipt1=ipt1+1_ip          
                ipt2=ipt2+1_ip  

             enddo
             !
             !     Set rsloc
             !
             rsloc=rblay(ilay)
             !
             !     Get the face number
             !
             nfac=nfac+2*nblay
             !
             !     Have we already inserted ipt1?
             !
             if(lpmark(ipt1)==0)then
                !
                !     Mark ipt1
                !   
                lpmark(ipt1)=1_ip
                !
                !     Project ipt1 on the surface
                !
                call gthost(coor(1,ipt1),lpsur(ipb1),ihost,rnofaold,ndim,nfold,&
                     npold,lfold,nnofa,coorold,d1,d2,d3,eltoelold,lelemold,&
                     lsurfold,isurf,ierr,rsloc)
                !
                !     Get the projected point
                !
                ip1o=lfold(1,ihost)
                ip2o=lfold(2,ihost)
                ip3o=lfold(3,ihost)
                coor(1,ipt1)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
                coor(2,ipt1)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
                coor(3,ipt1)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
                rnopo(1,ipt1)=rnofaold(1,ihost)
                rnopo(2,ipt1)=rnofaold(2,ihost)
                rnopo(3,ipt1)=rnofaold(3,ihost)
                icart=lcart(ipb1)
                call gtelem(ipt1,coor,npoin,ndim,lcell,ncell,icart,rtol)
                lcart(ipt1)=icart
                rsize(ipt1)=rsloc
                lpsur(ipt1)=ihost
                lptri(ipt1)=nfac
                !
                !      Check for closest points
                !
                call getClosIn2(coor(1,ipt1),lcartp(ipt1),ndim,coor,npoin,lcellp,ncellp,&
                     lmarkc,lstack,lctop,rsloc,lstackp,nstackp,mstackp,npmax,mcell,ipt1)
                !
                !     Filter the points of this column
                !
                call filtcol(lstackp,nstackp,npoin,lpoin,ipt1)
                !
                !     Filter the virtual points
                !
                call filtvirt(npoin,nfac,nface,nside,nnosi,lfsid,lsmark,nstackp,lstackp,&
                     lside,sitosi,lmark,lpoin,lfath)

                if(nstackp/=0)then
                   !
                   !     As there are points too close, mark the side and stop
                   !     the progression
                   !
                   lsmark(iside)=-lsmark(iside)
                   !
                   !     Mark the point
                   !
                   lpmark(ip1)=-1_ip
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   !
                   cycle

                endif
                !
                !     Check for closest sides
                !
                if(nsclos/=0)then
                   !
                   !     As there are sides too close, mark the side and stop
                   !     the progression
                   !
                   lsmark(iside)=-lsmark(iside)
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   !
                   cycle 

                endif

             else
                !
                !     Was the point already rejected
                !
                if(lpmark(ipt1)<0)then

                   lsmark(iside)=-lsmark(iside)
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   !
                   cycle

                endif

             endif
             !
             !     Have we already inserted ipt2
             !
             if(lpmark(ipt2)==0)then
                !
                !     Mark ipt2
                !   
                lpmark(ipt2)=1_ip
                !
                !     Project ipt2 on the surface
                !
                call gthost(coor(1,ipt2),lpsur(ip2),ihost,rnofaold,ndim,nfold,&
                     npold,lfold,nnofa,coorold,d1,d2,d3,eltoelold,lelemold,&
                     lsurfold,isurf,ierr,rsloc)
                !
                !     Get the projected point
                !
                ip1o=lfold(1,ihost)
                ip2o=lfold(2,ihost)
                ip3o=lfold(3,ihost)
                coor(1,ipt2)=d1*coorold(1,ip1o)+d2*coorold(1,ip2o)+d3*coorold(1,ip3o)
                coor(2,ipt2)=d1*coorold(2,ip1o)+d2*coorold(2,ip2o)+d3*coorold(2,ip3o)
                coor(3,ipt2)=d1*coorold(3,ip1o)+d2*coorold(3,ip2o)+d3*coorold(3,ip3o)
                rnopo(1,ipt2)=rnofaold(1,ihost)
                rnopo(2,ipt2)=rnofaold(2,ihost)
                rnopo(3,ipt2)=rnofaold(3,ihost)
                icart=lcart(ipb2)
                call gtelem(ipt2,coor,npoin,ndim,lcell,ncell,icart,rtol)
                lcart(ipt2)=icart
                rsize(ipt2)=rsloc
                lpsur(ipt2)=ihost
                lptri(ipt2)=nfac+1_ip
                !
                !      Check for closest points
                !
                call getClosIn2(coor(1,ipt2),lcartp(ipt2),ndim,coor,npoin,lcellp,ncellp,&
                     lmarkc,lstack,lctop,rsloc,lstackp,nstackp,mstackp,npmax,mcell,ipt2)
                !
                !     Filter the points of this column
                !
                call filtcol(lstackp,nstackp,npoin,lpoin,ipt2)
                !
                !     Filter the virtual points
                !
                call filtvirt(npoin,nfac,nface,nside,nnosi,lfsid,lsmark,nstackp,lstackp,&
                     lside,sitosi,lmark,lpoin,lfath)

                if(nstackp/=0)then
                   !
                   !     As there are points too close, mark the side and stop
                   !     the progression
                   !
                   lsmark(iside)=-lsmark(iside)
                   !
                   !     Mark the point
                   !
                   lpmark(ip1)=-1_ip
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   ! 
                   cycle 

                endif
                !
                !     Check for closest sides
                !
                if(nsclos/=0)then
                   !
                   !     As there are sides too close mark the side and stop
                   !     the progression
                   !
                   lsmark(iside)=-lsmark(iside)

                   cycle 

                endif

             else
                !
                !     Was the point already rejected
                !
                if(lpmark(ipt2)<0)then

                   lsmark(iside)=-lsmark(iside)
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   ! 
                   cycle

                endif

             endif
             !
             !     For virtual edges, do not consider geometric quality at lowest level
             !
             if(lsmark(iside)==2 .and. ilay==1)then
                !
                !
                !  
             else
                !
                !     Check quality
                !
                call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa,nfac,rnl)  
                call chkgeo(lface,nface,nnofa,nfac,qgeo2,rnofa,rnopo,npoin,ndim)

                if(qgeo2<qgeotol)then
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   ! 
                   cycle 

                endif

                call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa,nfac+1_ip,rnl)               
                call chkgeo(lface,nface,nnofa,nfac+1_ip,qgeo2,rnofa,rnopo,npoin,ndim)

                if(qgeo2<qgeotol)then
                   !
                   !     Delete the faces in the column 
                   !
                   nfacl=nfac
                   do jlay=ilay,nblay 
                      lfmark(nfacl)=0_ip 
                      lfmark(nfacl+1)=0_ip
                      nfacl=nfacl+2_ip
                   enddo
                   !
                   !     And go to the next side
                   ! 
                   cycle 

                endif
             endif
          endif
       enddo
    enddo
    !
    !     Clean up lsmark
    !
    do iside=1,nside
       lsmark(iside)=abs(lsmark(iside))
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','aniadv',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPMARK','aniadv',lpmark)
    deallocate(lpmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPMARK','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCART','aniadv',lcartp)
    deallocate(lcartp,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCART','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARKC','aniadv',lmarkc)
    deallocate(lmarkc,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARKC','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCTOP','aniadv',lctop)
    deallocate(lctop,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCTOP','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','aniadv',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'SITOSI','aniadv',sitosi)
    deallocate(sitosi,stat=istat)
    if(istat/=0) call memerr(2_ip,'SITOSI','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','aniadv',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','aniadv',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','aniadv',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','aniadv',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','aniadv',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','aniadv',0_ip)

  end subroutine aniadv

  subroutine addcover(lface,nface,nnofa,lfmark,lplay,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(inout) :: nface
    integer(ip),intent(in)    :: npoin,nnofa
    integer(ip),pointer       :: lface(:,:),lfmark(:)   
    integer(ip),pointer       :: lpoin(:),ptoel1(:),ptoel2(:)   
    integer(ip)               :: lplay(npoin)
    integer(ip)               :: ip1,ip2,ip3,iface,nfac0,jface,iflay    
    integer(ip)               :: ip1t,ip2t,ip1h,ip2h,ip1b,jp1b,jflay   
    integer(ip)               :: icont,iph,jp1t,jp2t,ip2b,ipmax,iface1,jface1
    integer(4)                :: istat
    !
    !     This sub adds elements to avoid bad aspect faces
    !
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','addcover',lpoin) 
    !
    !     Mark all the points touched by non deleted faces
    !
    do iface=1,nface
       if(lfmark(iface)/=0)then
          ip1=lface(1,iface) 
          ip2=lface(2,iface) 
          ip3=lface(3,iface) 
          lpoin(ip1)=1_ip 
          lpoin(ip2)=1_ip 
          lpoin(ip3)=1_ip 
       endif
    enddo
    !
    !     Loop on faces
    !
    nfac0=nface 
    iface=-1_ip

    do 
       !
       !     Find the next element marked for deletion
       !
       iface=iface+2_ip
       !
       !     Did we exhaust the list?
       !
       if(iface>nfac0)exit
       !
       !     Is the face still active?
       !
       if(lfmark(iface)==1)cycle 
       !
       !     Bottom side
       !   
       ip1b=lface(1,iface)
       ip2b=lface(2,iface)
       iflay=lplay(ip1b)
       !
       !     Top side
       !
       iface1=iface+1 
       ip1t=lface(1,iface1)
       ip2t=lface(3,iface1)
       !
       !     Initialize the levels
       ! 
       ip1h=0_ip
       ip2h=0_ip

       if(lpoin(ip1t)==1)ip1h=ip1h+1_ip
       if(lpoin(ip2t)==1)ip2h=ip2h+1_ip
       !
       !     Go upwards in the column
       !
       jface=iface

       do 
          !
          !     Update counter
          !
          jface=jface+2 
          !
          !     Did we exhaust the list?
          !
          if(jface>nfac0)exit
          !
          !     Bottom side
          !
          jp1b=lface(1,jface)
          jflay=lplay(jp1b)
          !
          !     Did we exhaust the column?
          !
          if(jflay<=iflay)exit
          !
          !     Top side
          !
          jface1=jface+1_ip
          jp1t=lface(1,jface1)
          jp2t=lface(3,jface1)
          if(lpoin(jp1t)==1)ip1h=ip1h+1_ip
          if(lpoin(jp2t)==1)ip2h=ip2h+1_ip

       enddo
       !
       !     Get the highest point
       !
       if(ip1h>ip2h)then
          ipmax=ip1t
          iph=ip1h
       else
          ipmax=ip2t
          iph=ip2h
       endif
       !
       !     Check that the jump is not higher than 1
       !
       if(iph>1)then
          write(*,*)'Jump higher than 1  --> stop '
          stop
       endif
       !
       !     Then decide the type of element to be added
       !
       icont=ip1h+ip2h

       if(icont==0)then
          !
          !     Nothing to do 
          !
          cycle

       else if(icont==2)then
          !
          !     Reintroduce the elements
          !           
          lfmark(iface)=1_ip
          lfmark(iface+1)=1_ip

       else 
          !
          !     Add the element
          !
          if(ip1h>ip2h)then

             nface=nface+1
             call memrea(nface,memor_msh,'LFACE','addcover',lface)
             lface(1,nface)=ip1b 
             lface(2,nface)=ip2b 
             lface(3,nface)=ip1t 

          else

             nface=nface+1
             call memrea(nface,memor_msh,'LFACE','addcover',lface)
             lface(1,nface)=ip1b 
             lface(2,nface)=ip2b 
             lface(3,nface)=ip2t 

          endif

       endif

    enddo
    !
    !     Resize lfmark
    !
    call memrea(nface,memor_msh,'LFMARK','addcover',lfmark)
    !
    !     Mark the new faces as kept
    !
    do iface=nfac0+1,nface
       lfmark(iface)=1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LPOIN','addcover',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','addcover',0_ip)

  end subroutine addcover

  subroutine blbound(lface,nface,nnofa,npoin,lbsid,nbsid,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)   :: nface,nnofa,npoin 
    integer(ip),intent(inout):: nbsid 
    integer(ip),intent(in)   :: lface(nnofa,nface),eltoel(nnofa,nface) 
    integer(ip),pointer      :: lbsid(:,:)
    integer(ip)              :: iface,inofa
    integer(4)               :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    !
    !     Count boundary sides
    !
    nbsid=0_ip
    do iface=1,nface
       do inofa=1,nnofa
          if(eltoel(inofa,iface)==0)then
             nbsid=nbsid+1_ip
          endif
       enddo
    enddo
    !
    !     Allocate lbsid
    ! 
    allocate(lbsid(2,nbsid),stat=istat)
    call memchk(zero,istat,memor_msh,'LBSID','blbound',lbsid) 
    !
    !     Store boundary sides
    !
    nbsid=0_ip
    do iface=1,nface
       do inofa=1,nnofa
          if(eltoel(inofa,iface)==0)then
             nbsid=nbsid+1_ip
             lbsid(1,nbsid)=lface(ltab(1,inofa),iface)
             lbsid(2,nbsid)=lface(ltab(2,inofa),iface)
          endif
       enddo
    enddo

  end subroutine blbound

  subroutine forcbl(lface,nnofa,nface,lbsid,nbsid,eltoel,lside,nside,npoin,&
       coor,ndim,lptri,lsmark,lpoin,isurf,rnopo,lcell,ncell,lcart,rsize,rsuni)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nnofa,nside,npoin,ndim,isurf,ncell
    integer(ip),intent(inout) :: nbsid
    integer(ip),intent(in)    :: lside(2,nside)
    integer(ip),intent(inout) :: lbsid(2,nbsid)
    integer(ip),intent(in)    :: lsmark(nside),lpoin(npoin),lcart(npoin)
    integer(ip),intent(inout) :: lptri(npoin)
    integer(ip),intent(inout) :: nface
    type(cell),intent(in)     :: lcell(ncell) 
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin),rsuni 
    real(rp),intent(inout)    :: rsize(npoin) 
    integer(ip),pointer       :: lface(:,:),eltoel(:,:)
    real(rp),pointer          :: rnofa(:,:)        
    integer(ip),pointer       :: lmark(:),lelem(:),lsurf(:)
    integer(ip)               :: iside,ip1,ip2,ip3,ipa,ipb,iface,iview
    integer(ip)               :: ipold,ihost,ierr,nbsid0
    real(rp)                  :: d1,d2,d3,p1(3),p2(3),rface(3),dtot
    real(rp)                  :: pproj(3),rx,ry,rz,rscal
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    !
    !     This sub is responsible of forcing the bl boundary in the old mesh 
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','forcbl',lmark) 
    allocate(lelem(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','forcbl',lelem) 
    allocate(lsurf(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURF','forcbl',lsurf) 
    allocate(rnofa(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RNOFA','forcbl',rnofa) 
    ierr=0_ip
    !
    !     Get the face normals
    !
    call gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    !
    !     Initialize lsurf
    !
    do iface=1,nface
       lsurf(iface)=isurf 
    enddo
    !
    !     Mark the points in lside
    !
    do iside=1,nside
       if(lsmark(iside)>0)then
          ip1=lside(1,iside)
          ip2=lside(2,iside)
          lmark(ip1)=1_ip
          lmark(ip2)=1_ip
       endif
    enddo
    !
    !     Compact lbsid
    ! 
    nbsid0=nbsid 
    nbsid=0_ip
    do iside=1,nbsid0
       ip1=lbsid(1,iside)
       ip2=lbsid(2,iside)
       if(lmark(ip1)==0 .or. lmark(ip2)==0)then
          nbsid=nbsid+1_ip
          lbsid(1,nbsid)=lbsid(1,iside)
          lbsid(2,nbsid)=lbsid(2,iside)
       endif
    enddo
    !
    !     Loop on points of lbsid
    !
    do iside=1,nbsid

       !write(*,*)iside
       !write(*,*)lface(1,719),lface(2,719),lface(3,719)
       ipa=lbsid(1,iside)
       !
       !     Has this point been marked or inserted? 
       !
       if(lmark(ipa)==0)then
          lmark(ipa)=-1_ip
          !
          !     Get point in side
          !
          ipold=lpoin(ipa)
          !
          !     Get host face of mesh lface
          !
          call gthost(coor(1,ipa),lptri(ipold),ihost,rnofa,ndim,nface,&
               npoin,lface,nnofa,coor,d1,d2,d3,eltoel,lelem,&
               lsurf,isurf,ierr,rsize(ipa))
          if(ierr==1)then
             write(*,*)'Error in gthost 1 in forcbl'
             stop
          endif 
          !
          !     Remember ihost
          !
          lptri(ipa)=ihost
          !
          !     Insert the point in the mesh
          ! 
          call insertp(nface,ihost,nnofa,ndim,npoin,d1,d2,d3,ipa,&
               lface,eltoel,lptri,lsurf,isurf,rnofa,coor,rnopo,lelem)
          !
          !     Get size of ipa  
          !
          call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipa,ndim,coor,rsuni)

       endif

       ipb=lbsid(2,iside)
       !
       !     Has this point been marked or inserted? 
       !
       if(lmark(ipb)==0)then
          lmark(ipb)=-1_ip
          !
          !     Get point in side
          !
          ipold=lpoin(ipb)
          !
          !     Get host face of mesh lface
          !
          call gthost(coor(1,ipb),lptri(ipold),ihost,rnofa,ndim,nface,&
               npoin,lface,nnofa,coor,d1,d2,d3,eltoel,lelem,&
               lsurf,isurf,ierr,rsize(ipb))
          if(ierr==1)then
             write(*,*)'Error in gthost 2 in forcbl'
             stop
          endif 
          !
          !     Remember ihost
          !
          lptri(ipb)=ihost
          !
          !     Insert the point in the mesh
          ! 
          call insertp(nface,ihost,nnofa,ndim,npoin,d1,d2,d3,ipb,&
               lface,eltoel,lptri,lsurf,isurf,rnofa,coor,rnopo,lelem)
          !
          !     Get size of ipb  
          !
          call gtsiz2(ncell,lcell,npoin,rsize,lcart,ipb,ndim,coor,rsuni)

       endif
       !
       !     Output the mesh
       !
       !call outaniface2(nnofa,nface,npoin,ndim,lface,coor)

    enddo
    !
    !     Output the mesh
    !
    !call outaniface2(nnofa,nface,npoin,ndim,lface,coor)
    !
    !     Now regenerate the sides
    !
    do iside=1,nbsid
       ip1=lbsid(1,iside)
       ip2=lbsid(2,iside)
       !
       !     At least one point must be untouched
       !       
       if(lmark(ip1)/=1_ip .or. lmark(ip2)/=1_ip)then
          write(*,*)iside,ip1,ip2
          !
          !     Regenerate side (ip1,ip2)
          !
          call recoverp(ip1,ip2,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,rnofa,ierr)
          !
          !     Output the mesh
          !
          if(ierr==1)then
             call outaniface2(nnofa,nface,npoin,ndim,lface,coor)
             stop
          endif  

       endif

    enddo
    !
    !     Output the mesh
    !
    !call outaniface2(nnofa,nface,npoin,ndim,lface,coor)

    call memchk(2_ip,istat,memor_msh,'RNOFA','forcbl',rnofa)
    deallocate(rnofa,stat=istat)
    if(istat/=0) call memerr(2_ip,'RNOFA','forcbl',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSURF','forcbl',lsurf)
    deallocate(lsurf,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSURF','forcbl',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','forcbl',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','forcbl',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','forcbl',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','forcbl',0_ip)

  end subroutine forcbl

  subroutine insertp(nface,ihost,nnofa,ndim,npoin,d1,d2,d3,ipnew,&
       lface,eltoel,lptri,lsurf,isurf,rnofa,coor,rnopo,lelem)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(inout) :: nface
    integer(ip),intent(in)    :: ndim,npoin,ihost,nnofa,ipnew,isurf
    real(rp),intent(in)       :: d1,d2,d3,coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lsurf(:),lelem(:)
    real(rp),pointer          :: rnofa(:,:)
    integer(ip),intent(inout) :: lptri(npoin)
    real(rp)                  :: epsil,modul
    integer(ip)               :: ip1,ip2,ip3,neigh1,neigh2,nfac1,nfnew
    integer(ip)               :: ineigh1,ineigh2,ineigh3,ierr
    epsil=1.0d-01
    ierr=0_ip 
    !
    !     Get host face
    !      
    ip1=lface(1,ihost)
    ip2=lface(2,ihost)
    ip3=lface(3,ihost)
    !
    !     Get the neighbors
    !
    ineigh1=eltoel(1,ihost) 
    ineigh2=eltoel(2,ihost) 
    ineigh3=eltoel(3,ihost) 
    !
    !     Check the simplicial coordinates
    !
    if(d1<epsil)then
       if(ineigh1/=0)then
          !
          !     Split on first side
          !  
          call splitp(ihost,1_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
               lface,eltoel,lsurf,isurf,rnofa,coor,lelem,rnopo,ierr)
          if(ierr==0)then
             return
          endif
       endif

    else if(d2<epsil)then
       if(ineigh2/=0)then
          !
          !     Split on second side
          !  
          call splitp(ihost,2_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
               lface,eltoel,lsurf,isurf,rnofa,coor,lelem,rnopo,ierr)
          if(ierr==0)then
             return
          endif
       endif

    else if(d3<epsil)then
       if(ineigh3/=0)then
          !
          !     Split on third side
          !  
          call splitp(ihost,3_ip,nnofa,nface,npoin,ndim,ipnew,lptri,&
               lface,eltoel,lsurf,isurf,rnofa,coor,lelem,rnopo,ierr)
          if(ierr==0)then
             return
          endif
       endif
    endif
    !
    !     It is possible that some of the shape function be below epsil,
    !     but the point is in front of a good element because the old element
    !     in front is very large. Therefore, the element in front must be split
    !     in three elements.  
    !


    !
    !     Insert npoin in the surface mesh 
    !   
    lface(3,ihost)=ipnew
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,ihost),ihost,modul)
    !
    !     Realloc lface 
    !
    nfnew=nface+2
    call memrea(nfnew,memor_msh,'LFACE','insertp',lface)
    call memrea(nfnew,memor_msh,'ELTOEL','insertp',eltoel)
    call memrea(nfnew,memor_msh,'LSURF','insertp',lsurf)
    call memrea(nfnew,memor_msh,'LSURF','insertp',lelem)
    call memrea(nfnew,memor_msh,'RNOFA','insertp',rnofa)

    nface=nface+1
    nfac1=nface
    lface(1,nface)=ip2
    lface(2,nface)=ip3
    lface(3,nface)=ipnew
    lsurf(nface)=isurf
    lelem(nface)=0_ip
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,nface),nface,modul)

    nface=nface+1
    lface(1,nface)=ip3
    lface(2,nface)=ip1
    lface(3,nface)=ipnew
    lsurf(nface)=isurf
    lelem(nface)=0_ip
    call gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa(:,nface),nface,modul)
    !
    !     Update eltoel
    !
    neigh1=eltoel(1,ihost)
    neigh2=eltoel(2,ihost)
    eltoel(1,ihost)=nfac1
    eltoel(2,ihost)=nface

    eltoel(1,nfac1)=nface
    eltoel(2,nfac1)=ihost
    eltoel(3,nfac1)=neigh1

    eltoel(1,nface)=ihost
    eltoel(2,nface)=nfac1
    eltoel(3,nface)=neigh2
    !
    !     Update eltoel outside
    !
    if(neigh1/=0)then
       if(eltoel(1,neigh1)==ihost)then
          eltoel(1,neigh1)=nfac1
       else if(eltoel(2,neigh1)==ihost)then
          eltoel(2,neigh1)=nfac1
       else
          eltoel(3,neigh1)=nfac1
       endif
    endif

    if(neigh2/=0)then
       if(eltoel(1,neigh2)==ihost)then
          eltoel(1,neigh2)=nface
       else if(eltoel(2,neigh2)==ihost)then
          eltoel(2,neigh2)=nface
       else
          eltoel(3,neigh2)=nface
       endif
    endif

    lptri(ipnew)=ihost
    lptri(ip3)=nfac1

  end subroutine insertp

  subroutine splitp(ielem1,iview1,nnofa,nface,npoin,ndim,ipnew,lptri,&
       lface,eltoel,lsurf,isurf,rnofa,coor,lelem,rnopo,ierr)
    use def_kintyp, only       :  ip,rp,lg,cell
    use mod_memchk
    use def_meshin, only     :  memor_msh
    implicit none
    integer(ip),intent(in)    :: ndim,nnofa,ielem1,iview1,ipnew,npoin,isurf
    integer(ip),intent(inout) :: nface,ierr
    integer(ip),intent(inout) :: lptri(npoin)
    integer(ip),pointer       :: lface(:,:),eltoel(:,:),lsurf(:),lelem(:)
    real(rp),pointer          :: rnofa(:,:)
    real(rp),intent(in)       :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip)               :: ip1,ip2,ip3,ip4,iview2,ielem2
    integer(ip)               :: ineigh11,ineigh12,ineigh21,ineigh22
    integer(ip)               :: nfnew,ielem3,ielem4
    integer(ip)               :: lfloc(3,4)
    real(rp)                  :: rnofloc(3,4)
    real(rp)                  :: modul,qgeo2,qgeotol 
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    !
    !     This sub splits the edge separing elements ielem1 and ielem2
    !     by inserting a point at the middle of the edge
    !     Compared to split, point ipnew and its database have already been computed
    !
    qgeotol=0.1d+00

    ielem2=eltoel(iview1,ielem1)  

    ip1=lface(iview1,ielem1)
    ip2=lface(ltab(1,iview1),ielem1)
    ip3=lface(ltab(2,iview1),ielem1)

    if(eltoel(1,ielem2)==ielem1)then
       iview2=1_ip
    else if(eltoel(2,ielem2)==ielem1)then
       iview2=2_ip
    else
       iview2=3_ip
    endif

    ip4=lface(iview2,ielem2)
    !
    !     Create new elements
    ! 
    lfloc(1,1)=ip1  
    lfloc(2,1)=ip2  
    lfloc(3,1)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,1),1_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,1_ip,qgeo2,rnofloc(1,1),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,2)=ip2  
    lfloc(2,2)=ip4  
    lfloc(3,2)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,2),2_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,2_ip,qgeo2,rnofloc(1,2),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,3)=ip4  
    lfloc(2,3)=ip3  
    lfloc(3,3)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,3),3_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,3_ip,qgeo2,rnofloc(1,3),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif

    lfloc(1,4)=ip3  
    lfloc(2,4)=ip1  
    lfloc(3,4)=ipnew
    call gtfnr2(lfloc,4_ip,nnofa,ndim,coor,npoin,rnofloc(1,4),4_ip,modul)
    call chkgeo(lfloc,4_ip,nnofa,4_ip,qgeo2,rnofloc(1,4),rnopo,npoin,ndim)
    if(qgeo2<qgeotol)then
       ierr=1_ip
       return
    endif
    !
    !     We can split
    !
    ineigh11=eltoel(ltab(1,iview1),ielem1)
    ineigh12=eltoel(ltab(2,iview1),ielem1)
    ineigh21=eltoel(ltab(1,iview2),ielem2)
    ineigh22=eltoel(ltab(2,iview2),ielem2)
    !
    !     Allocate the 2 new elements
    !
    nfnew=nface+2
    call memrea(nfnew,memor_msh,'LFACE','split2',lface)
    call memrea(nfnew,memor_msh,'ELTOEL','split2',eltoel)
    call memrea(nfnew,memor_msh,'LSURF','split2',lsurf)
    call memrea(nfnew,memor_msh,'RNOFA','split2',rnofa)
    call memrea(nfnew,memor_msh,'LELEM','split2',lelem)

    nface=nface+1
    ielem3=nface
    nface=nface+1
    ielem4=nface
    !
    !     Create new elements
    ! 
    lface(1,ielem1)=lfloc(1,1)  
    lface(2,ielem1)=lfloc(2,1)
    lface(3,ielem1)=lfloc(3,1)
    lsurf(ielem1)=isurf
    eltoel(1,ielem1)=ielem2
    eltoel(2,ielem1)=ielem3
    eltoel(3,ielem1)=ineigh12
    lelem(ielem1)=0_ip
    rnofa(1,ielem1)=rnofloc(1,1) 
    rnofa(2,ielem1)=rnofloc(2,1) 
    rnofa(3,ielem1)=rnofloc(3,1) 

    lface(1,ielem2)=lfloc(1,2)
    lface(2,ielem2)=lfloc(2,2)
    lface(3,ielem2)=lfloc(3,2)
    lsurf(ielem2)=isurf
    eltoel(1,ielem2)=ielem4
    eltoel(2,ielem2)=ielem1
    eltoel(3,ielem2)=ineigh21
    lelem(ielem2)=0_ip
    rnofa(1,ielem2)=rnofloc(1,2) 
    rnofa(2,ielem2)=rnofloc(2,2) 
    rnofa(3,ielem2)=rnofloc(3,2) 

    lface(1,ielem4)=lfloc(1,3)
    lface(2,ielem4)=lfloc(2,3)
    lface(3,ielem4)=lfloc(3,3)
    lsurf(ielem4)=isurf
    eltoel(1,ielem4)=ielem3
    eltoel(2,ielem4)=ielem2
    eltoel(3,ielem4)=ineigh22
    lelem(ielem4)=0_ip
    rnofa(1,ielem4)=rnofloc(1,3) 
    rnofa(2,ielem4)=rnofloc(2,3) 
    rnofa(3,ielem4)=rnofloc(3,3) 

    lface(1,ielem3)=lfloc(1,4)
    lface(2,ielem3)=lfloc(2,4)  
    lface(3,ielem3)=lfloc(3,4)
    lsurf(ielem3)=isurf
    eltoel(1,ielem3)=ielem1
    eltoel(2,ielem3)=ielem4
    eltoel(3,ielem3)=ineigh11
    lelem(ielem3)=0_ip
    rnofa(1,ielem3)=rnofloc(1,4) 
    rnofa(2,ielem3)=rnofloc(2,4) 
    rnofa(3,ielem3)=rnofloc(3,4) 

    lptri(ip3)=ielem3
    lptri(ipnew)=ielem1

    !
    !     Update outside
    !
    if(ineigh22/=0)then
       if(eltoel(1,ineigh22)==ielem2)then
          eltoel(1,ineigh22)=ielem4
       else if(eltoel(2,ineigh22)==ielem2)then
          eltoel(2,ineigh22)=ielem4
       else
          eltoel(3,ineigh22)=ielem4
       endif
    endif

    if(ineigh11/=0)then 
       if(eltoel(1,ineigh11)==ielem1)then
          eltoel(1,ineigh11)=ielem3
       else if(eltoel(2,ineigh11)==ielem1)then
          eltoel(2,ineigh11)=ielem3
       else
          eltoel(3,ineigh11)=ielem3
       endif
    endif

  end subroutine splitp

  subroutine cleansurf(lface,nnofa,nface,npoin,lside,nside,nnosi,eltoel,&
       coor,ndim,lbsid,nbsid,lptri)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nnofa,npoin,nside,nnosi,ndim,nbsid
    integer(ip),intent(inout) :: nface
    integer(ip),intent(inout) :: lface(nnofa,nface),eltoel(nnofa,nface)
    integer(ip),intent(inout) :: lptri(npoin)
    integer(ip),intent(in)    :: lside(nnosi,nside),lbsid(nnosi,nbsid)
    real(rp),intent(in)       :: coor(ndim,npoin)
    integer(ip),pointer       :: lcolor(:),lstack(:),ptoel1(:),ptoel2(:)  
    integer(ip),pointer       :: ptosi1(:),ptosi2(:),lmark(:)  
    integer(ip)               :: iface,iside,icolor,ncolor,ncolor0,istack,jpoin    
    integer(ip)               :: nstack,jface,inofa,ineigh,ip1,ip2,isto,ipa,ipb,nfac0   
    logical(lg)               :: ifound 
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    !
    !     This sub is responsible of taking out the elements covered by
    !     the bl
    !
    allocate(lmark(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','anisurf',lmark) 
    allocate(lstack(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','anisurf',lstack) 
    !
    !     Get the sides surrounding the points
    !
    call ptoelm(lbsid,nbsid,npoin,nnosi,ptosi1,ptosi2)
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Initialize the stack
    !
    lstack(1)=1_ip
    ncolor=1_ip
    lmark(1)=1_ip
    nstack=1_ip
    !
    !     Outer loop on faces
    !
    istack=0_ip   

    do 

       if(istack==nstack)then
          !
          !     Scan for unmarked face
          !
          do iface=1,nface
             if(lmark(iface)==0)then
                nstack=nstack+1
                lstack(nstack)=iface
                ncolor=ncolor+1_ip 
                lmark(iface)=ncolor
                exit
             endif
          enddo
          !
          !     No more faces?
          !
          if(iface==nface+1)then 
             if(nstack==nface)then 
                exit
             else
                write(*,*)'Error in cleansurf'
                stop
             endif
          endif
       endif

       istack=istack+1_ip
       jface=lstack(istack)
       !
       !     Loop on sides
       !
       do inofa=1,nnofa
          ineigh=eltoel(inofa,jface)
          !
          !     Do we have a neighbor?   
          !
          if(ineigh==0)cycle
          !
          !     Has the neighbor been already marked?
          !
          if(lmark(ineigh)==0)then
             !
             !     Get the points of this side
             !
             ip1=lface(ltab(1,inofa),jface)
             ip2=lface(ltab(2,inofa),jface)
             ifound=.false. 

             do isto=ptosi2(ip1),ptosi2(ip1+1)-1
                iside=ptosi1(isto)
                ipa=lbsid(1,iside) 
                ipb=lbsid(2,iside) 
                if(ipa==ip1)then
                   jpoin=ipb
                else 
                   jpoin=ipa  
                endif
                !
                !     Did we find a boundary side?
                !
                if(jpoin==ip2)then
                   ifound=.true.
                   exit  
                endif
             enddo
             !
             !     We did not find a boundary side
             !
             if(ifound.eqv. .false.)then
                lmark(ineigh)=ncolor 
                nstack=nstack+1 
                lstack(nstack)=ineigh         
             endif
          endif
       enddo
    enddo
    !
    !     Output mesh
    !
    !call outaniface3(nnofa,nface,npoin,ndim,lface,coor,lmark)
    !
    !     Allocate lcolor
    !
    allocate(lcolor(ncolor),stat=istat)
    call memchk(zero,istat,memor_msh,'LCOLOR','anisurf',lcolor) 
    !
    !     Now find the colors touching the sides
    !
    ncolor0=0_ip
    do iside=1,nside
       ip1=lside(1,iside)
       do isto=ptoel2(ip1),ptoel2(ip1+1)-1
          iface=ptoel1(isto)
          icolor=lmark(iface)
          lcolor(icolor)=1_ip
       enddo
    enddo
    !
    !     DBG
    !
    do icolor=1,ncolor 
       write(*,*)icolor,lcolor(icolor)
    enddo
    !
    !     Finally mark the faces to be deleted in lmark
    !
    do iface=1,nface
       icolor=lmark(iface)
       if(lcolor(icolor)==1)then
          lmark(iface)=0
       endif
    enddo
    !
    !     Compress lface
    !
    nfac0=nface
    nface=0_ip 
    do iface=1,nfac0
       if(lmark(iface)/=0)then
          nface=nface+1_ip
          lface(1,nface)=lface(1,iface)
          lface(2,nface)=lface(2,iface)
          lface(3,nface)=lface(3,iface)
          lptri(lface(1,nface))=nface
          lptri(lface(2,nface))=nface
          lptri(lface(3,nface))=nface
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'PTOEL1','cleansurf',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','cleansurf',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','cleansurf',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','cleansurf',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCOLOR','cleansurf',lcolor)
    deallocate(lcolor,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCOLOR','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','cleansurf',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','cleansurf',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','cleansurf',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','cleansurf',0_ip)

  end subroutine cleansurf

  subroutine rmvjmp(lface,nface,npoin,nnofa,lplay,lfmark,nside,nblay)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nnofa,npoin,nside,nblay
    integer(ip),intent(inout) :: nface
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),pointer       :: ptoel2(:),ptoel1(:)
    integer(ip),intent(in)    :: lplay(npoin)  
    integer(ip),intent(inout) :: lfmark(nface)  
    integer(ip)               :: iface,jface,kface,isto,ip1,ip2,ip3,kflay,jflay
    integer(ip)               :: inofa,ipoin,lplist(2),iface2,iflay,iplist
    logical(lg)               :: imarked
    integer(4)                :: istat
    !
    !     This sub guarantees that the jump between neighbors is not higher than 1
    !
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Loop on iteration
    !
    do 
       !
       !     Loop on elements
       !  
       iface=-1_ip
       !
       !     Set flag for exit condition
       ! 
       imarked=.false.

       do 
          !
          !     Update counter
          !
          iface=iface+2 
          !
          !     Did we exhaust the list?
          !
          if(iface>nface)exit
          !
          !     Is the  face marked for deletion?
          !
          if(lfmark(iface)==0)then
             !
             !     Get the level
             !
             iflay=lplay(lface(1,iface))
             !
             !     Create the list of points
             !
             lplist(1)=lface(1,iface+1)
             lplist(2)=lface(3,iface+1)
             !
             !     Mark the faces that are too high in the column
             !
             do iplist=1,2
                ipoin=lplist(iplist) 
                do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
                   kface=ptoel1(isto)  
                   !
                   !     Has this face already been deleted before
                   !
                   if(lfmark(kface)==1)then
                      ip1=lface(1,kface)
                      ip2=lface(2,kface)
                      ip3=lface(3,kface)
                      kflay=min(lplay(ip1),lplay(ip2),lplay(ip3))
                      if(kflay.gt.iflay) then
                         lfmark(kface)=0
                         imarked=.true.
                      endif
                   endif
                enddo
             enddo
          endif
       enddo
       !
       !     Exit condition
       !
       if(imarked.eqv. .false.)exit
       !
       !     Mark the elements not marked in the column
       ! 
       call mrkcol(nside,nblay,lfmark,nface)

    enddo

    call memchk(2_ip,istat,memor_msh,'PTOEL1','rmvjmp',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','rmvjmp',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','rmvjmp',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','rmvjmp',0_ip)

  end subroutine rmvjmp

  subroutine reorside(nside,lside,nnosi,eltoel,lface,nnofa,nface,npoin,&
       lline,ptoel1,ptoel2,lsmark,nsconn,lsconn,coor,rnopo,rtapo,&
       lpoin,lpsur,lptri,lcart,lptype,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :  memor_msh
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: nnosi,nnofa,nface
    integer(ip),intent(inout) :: nside,nsconn,npoin
    integer(ip),intent(in)    :: lface(nnofa,nface),eltoel(nnofa,nface)
    integer(ip),pointer       :: lside(:,:),lline(:),lsmark(:),lsconn(:)
    integer(ip),pointer       :: ptosi1(:),ptosi2(:),lsnew(:,:),lrenu(:)
    integer(ip),pointer       :: llinenew(:),lsmnew(:),lfath(:)
    integer(ip),pointer       :: lpoin(:),lpsur(:),lptri(:),lcart(:),lptype(:,:)
    real(rp),pointer          :: coor(:,:),rnopo(:,:),rtapo(:,:)
    integer(ip)               :: ptoel2(npoin+1),ptoel1(*)
    logical(lg)               :: ifound
    integer(ip)               :: ipa,ipb,ipc,ipd,ip1,ip2,isid,iside,ncont
    integer(ip)               :: iface,inofa,isnew,innerside,ipoin,jsid
    integer(ip)               :: isto,jpmin,jpmax,ipt,ipmin,ipmax
    integer(4)                :: istat
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    !
    !     This sub reorders and orients the sides
    !
    !
    !     On input:  lside(2,nside):  the sides
    !                lline(nside):    the line number
    !                lsmark(nside):   the wetted side marked 
    !
    !     On output: lside(2,nside):  the sides
    !                lline(nside):    the line number
    !                lsmark(nside):   the wetted side marked 
    !
    allocate(llinenew(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','reorside',llinenew) 
    allocate(lrenu(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','reorside',lrenu) 
    allocate(lsnew(2,nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSNEW','reorside',lsnew) 
    allocate(lsmnew(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LSMNEW','reorside',lsmnew) 
    !
    !     Get the sides surrounding the points
    !
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
    !
    !     Get one boundary side 
    !
    do iface=1,nface
       do inofa=1,nnofa
          if(eltoel(inofa,iface)==0)then
             ip1=lface(ltab(1,inofa),iface)
             ip2=lface(ltab(2,inofa),iface)
             goto 100
          endif
       enddo
    enddo
    !
    !     Error 
    ! 
    write(*,*)'Error in reorside, sides not found'
    stop

100 continue
    !
    !     Get first side (Should be on the boundary)
    !
    ifound=.false.
    do isid=ptosi2(ip1),ptosi2(ip1+1)-1
       iside=ptosi1(isid)
       ipa=lside(1,iside) 
       ipb=lside(2,iside) 

       if(ipa==ip1)then
          ipc=ipb
       else
          ipc=ipa
       endif

       if(ipc==ip2)then
          ifound=.true.
          exit
       endif
    enddo
    !
    !     Any error?
    !
    if(ifound.eqv. .false.)then
       write(*,*)'Error in reorside, first point not found'
       stop
    endif
    !
    !     Initialize renumbering
    !  
    allocate(lsconn(2),stat=istat)
    call memchk(zero,istat,memor_msh,'LSCONN','reorside',lsconn)
    lsconn(1)=1_ip 
    nsconn=1_ip
    ncont=1_ip
    lrenu(iside)=ncont
    lsnew(1,ncont)=ip1
    lsnew(2,ncont)=ip2
    llinenew(ncont)=lline(iside)
    lsmnew(ncont)=lsmark(iside)
    !
    !     Form outer contour
    !
    call addside(nnosi,nside,npoin,nnofa,nface,lface,eltoel,lside,lline,&
         ptosi1,ptosi2,lsnew,lrenu,llinenew,ptoel1,ptoel2,ip1,ip2,&
         ncont,iside,lsmark,lsmnew,coor,rnopo,rtapo,lpoin,lpsur,lptri,lcart,&
         lptype,lfath)
    lsconn(2_ip)=ncont+1_ip
    !
    !     Do we have inner sides?
    !
    do iside=1,nside
       if(lrenu(iside)==0)then
          ip1=lside(1,iside)
          ip2=lside(2,iside)
          !
          !     Is it a cusp point?
          !
          if(lptype(1,ip1)==ID_CUSP)then

             ncont=ncont+1_ip
             lrenu(iside)=ncont
             lsnew(1,ncont)=ip1
             lsnew(2,ncont)=ip2
             llinenew(ncont)=lline(iside) 
             lsmnew(ncont)=lsmark(iside)
             !
             !     Form a new contour
             !
             call addside(nnosi,nside,npoin,nnofa,nface,lface,eltoel,lside,lline,&
                  ptosi1,ptosi2,lsnew,lrenu,llinenew,ptoel1,ptoel2,ip1,ip2,&
                  ncont,iside,lsmark,lsmnew,coor,rnopo,rtapo,lpoin,lpsur,lptri,lcart,&
                  lptype,lfath)
             nsconn=nsconn+1_ip
             call memrea(nsconn+1_ip,memor_msh,'lsconn','reorside',lsconn)
             lsconn(nsconn+1)=ncont+1_ip

          endif
          !
          !     Is it a cusp point?
          !
          if(lptype(1,ip2)==ID_CUSP)then
             !
             !     We may have two cusp points on the same side
             !     We do not want to count it twice
             !
             if(lrenu(iside)==0)then  

                ncont=ncont+1_ip
                lrenu(iside)=ncont
                lsnew(1,ncont)=ip2
                lsnew(2,ncont)=ip1
                llinenew(ncont)=lline(iside) 
                lsmnew(ncont)=lsmark(iside)
                !
                !     Form a new contour
                !
                call addside(nnosi,nside,npoin,nnofa,nface,lface,eltoel,lside,lline,&
                     ptosi1,ptosi2,lsnew,lrenu,llinenew,ptoel1,ptoel2,ip2,ip1,&
                     ncont,iside,lsmark,lsmnew,coor,rnopo,rtapo,lpoin,lpsur,lptri,lcart,&
                     lptype,lfath)
                nsconn=nsconn+1_ip
                call memrea(nsconn+1_ip,memor_msh,'lsconn','reorside',lsconn)
                lsconn(nsconn+1)=ncont+1_ip

             endif
          endif
       endif
    enddo
    !
    !     We may still have inner sides without cusp points (damn it)
    !
    do iside=1,nside
       if(lrenu(iside)==0)then
          ip1=lside(1,iside)
          ip2=lside(2,iside)
          !
          !     Initialize new component
          !
          ncont=ncont+1_ip
          lrenu(iside)=ncont
          lsnew(1,ncont)=ip1
          lsnew(2,ncont)=ip2
          llinenew(ncont)=lline(iside) 
          lsmnew(ncont)=lsmark(iside)
          !
          !     Form new contour
          !      
          call addside(nnosi,nside,npoin,nnofa,nface,lface,eltoel,lside,lline,&
               ptosi1,ptosi2,lsnew,lrenu,llinenew,ptoel1,ptoel2,ip1,ip2,&
               ncont,iside,lsmark,lsmnew,coor,rnopo,rtapo,lpoin,lpsur,lptri,lcart,&
               lptype,lfath)
          nsconn=nsconn+1_ip
          call memrea(nsconn+1_ip,memor_msh,'lsconn','reorside',lsconn)
          lsconn(nsconn+1)=ncont+1_ip

       endif
    enddo
    !
    !     Resize
    !
    nside=ncont  
    call memrea(nside,memor_msh,'lside','reorside',lside)
    call memrea(nside,memor_msh,'lline','reorside',lline)
    call memrea(nside,memor_msh,'lsmark','reorside',lsmark)
    !
    !     Transfer in lside
    !
    do iside=1,nside
       lside(1,iside)=lsnew(1,iside)
       lside(2,iside)=lsnew(2,iside)
       lline(iside)=llinenew(iside)
       lsmark(iside)=lsmnew(iside)
    enddo
    !
    !     Deallocate temporal arrays
    !
    call memchk(2_ip,istat,memor_msh,'LRENU','reorside',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','reorside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LLINENEW','reorside',llinenew)
    deallocate(llinenew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LLINENEW','reorside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSNEW','reorside',lsnew)
    deallocate(lsnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSNEW','reorside',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSMNEW','reorside',lsmnew)
    deallocate(lsmnew,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSMNEW','reorside',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','reorside',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','reorside',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','reorside',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','reorside',0_ip)

  end subroutine reorside

  subroutine addside(nnosi,nside,npoin,nnofa,nface,lface,eltoel,lside,lline,&
       ptosi1,ptosi2,lsnew,lrenu,llinenew,ptoel1,ptoel2,ip1,ip2,&
       ncont,iside,lsmark,lsmnew,coor,rnopo,rtapo,lpoin,lpsur,lptri,lcart,&
       lptype,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :  memor_msh
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: nnosi,nnofa,nface
    integer(ip),intent(in)    :: ip1,ip2,iside          
    integer(ip),intent(inout) :: nside,ncont,npoin
    integer(ip),intent(in)    :: lface(nnofa,nface),eltoel(nnofa,nface)
    integer(ip),pointer       :: lside(:,:),lline(:)
    integer(ip),pointer       :: ptosi1(:),ptosi2(:),lsnew(:,:),lrenu(:)
    integer(ip),pointer       :: llinenew(:),lsmark(:),lsmnew(:),lfath(:)
    integer(ip),pointer       :: lpoin(:),lpsur(:),lptri(:),lcart(:),lptype(:,:)
    real(rp),pointer          :: coor(:,:),rnopo(:,:),rtapo(:,:)
    integer(ip)               :: ptoel2(npoin+1),ptoel1(*)
    integer(ip)               :: isid,jside,ipa,ipb,ipmin,ipmax,nsid,ip1t,ip2t
    integer(ip)               :: ipc,ipd,iface,jpmin,jpmax,inofa,isto,ipt,isidet
    integer(ip)               :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  
    real(rp)                  :: c00
    !
    !     This sub adds sides to form the contours of the surface
    !
    c00=0.0d+00
    !
    !     On input:
    ! 
    !     ip1,ip2 the current side iside
    !     
    !
    ip1t=ip1
    ip2t=ip2
    isidet=iside   
    !
    !     Loop on the boundary sides
    !
    do
       !
       !     How many sides has ip2t? (May be globally a corner point and 
       !     locally a ridge point on the patch)
       !
       nsid=ptosi2(ip2t+1)-ptosi2(ip2t)

       if(nsid==2)then
          !
          !     It is a surface ridge point, most common case
          !
          !
          !     Loop on sides of ip2t
          !     Find the side which is not isidet
          !  
          do isid=ptosi2(ip2t),ptosi2(ip2t+1)-1
             jside=ptosi1(isid)
             if(jside/=isidet)then
                !
                !     Did we reach the first side?
                !
                if(jside==iside)then
                   !
                   !     Did we reach the first point?
                   !     One side may appear twice for inner sides 
                   !
                   if(ip2t==ip1)then
                      return
                   endif                   
                endif

                ipa=lside(1,jside) 
                ipb=lside(2,jside) 

                if(ipa==ip2t)then
                   ipc=ipa
                   ipd=ipb
                else
                   ipc=ipb
                   ipd=ipa
                endif

                ncont=ncont+1_ip
                call memrea(ncont,memor_msh,'lsnew','reorside',lsnew)
                call memrea(ncont,memor_msh,'llinenew','reorside',llinenew)
                call memrea(ncont,memor_msh,'lsmnew','reorside',lsmnew)
                !
                !     Is it the first time we reach jside?
                !
                if(lrenu(jside)==0)then       

                   lsnew(1,ncont)=ipc
                   lsnew(2,ncont)=ipd

                else
                   !
                   !     This is the second time we reach iside
                   !     Is ipd a cusp point?
                   !
                   if(lptype(1,ipd)==ID_CUSP)then
                   
                      lsnew(1,ncont)=lsnew(2,ncont-1)
                      lsnew(2,ncont)=ipd

                   else
                      !     
                      !     No, ipd should be duplicated
                      !     Reallocate coor
                      !
                      npoin=npoin+1
                      call memrea(npoin,memor_msh,'COOR','chkconv',coor)
                      call memrea(npoin,memor_msh,'RNOPO','chkconv',rnopo)
                      call memrea(npoin,memor_msh,'RTAPO','chkconv',rtapo)
                      call memrea(npoin,memor_msh,'LPOIN','chkconv',lpoin)
                      call memrea(npoin,memor_msh,'LPSUR','chkconv',lpsur)
                      call memrea(npoin,memor_msh,'LTRI','chkconv',lptri)
                      call memrea(npoin,memor_msh,'LCART','chkconv',lcart)
                      call memrea(npoin,memor_msh,'LPTYPE','chkconv',lptype)
                      call memrea(npoin,memor_msh,'LFATH','chkconv',lfath)
                      
                      coor(1,npoin)=coor(1,ipd)
                      coor(2,npoin)=coor(2,ipd)
                      coor(3,npoin)=coor(3,ipd)
                      rtapo(1,npoin)=c00
                      rtapo(2,npoin)=c00
                      rtapo(3,npoin)=c00
                      lpsur(npoin)=lpsur(ipd)
                      lptri(npoin)=lptri(ipd)
                      lcart(npoin)=lcart(ipd)
                      rnopo(1,npoin)=rnopo(1,ipd)
                      rnopo(2,npoin)=rnopo(2,ipd)
                      rnopo(3,npoin)=rnopo(3,ipd)
                      lptype(1,npoin)=lptype(1,ipd)
                      lptype(2,npoin)=lptype(2,ipd)
                      lfath(npoin)=ipd  
 
                      lsnew(1,ncont)=lsnew(2,ncont-1)
                      lsnew(2,ncont)=npoin
                  
                   endif  
                
                endif

                lrenu(jside)=ncont
                llinenew(ncont)=lline(jside)
                lsmnew(ncont)=lsmark(isidet)
                ip1t=ipc
                ip2t=ipd
                isidet=jside
                exit

             endif

          enddo
       !
       !     This is the surface cusp point
       !
       else if(nsid==1)then
          !
          !     Did we reach the first side?
          !
          if(isidet==iside)then
             !
             !     Did we reach the first point?
             !     One side may appear twice for inner sides 
             !
             if(ip2t==ip1)then
                return
             endif 
          endif
          !
          !     Invert direction
          !  
          ipt=ip1t
          ip1t=ip2t
          ip2t=ipt
          !
          !     Go around the cusp and add a new side
          !
          ncont=ncont+1_ip 
          call memrea(ncont,memor_msh,'lsnew','reorside',lsnew)
          call memrea(ncont,memor_msh,'llinenew','reorside',llinenew)
          call memrea(ncont,memor_msh,'lsmnew','reorside',lsmnew)
          !     
          !     Now add a new point if current ip2t is not a cusp point
          !
          if(lptype(1,ip2t)/=ID_CUSP)then

             npoin=npoin+1
             call memrea(npoin,memor_msh,'COOR','chkconv',coor)
             call memrea(npoin,memor_msh,'RNOPO','chkconv',rnopo)
             call memrea(npoin,memor_msh,'RTAPO','chkconv',rtapo)
             call memrea(npoin,memor_msh,'LPOIN','chkconv',lpoin)
             call memrea(npoin,memor_msh,'LPSUR','chkconv',lpsur)
             call memrea(npoin,memor_msh,'LTRI','chkconv',lptri)
             call memrea(npoin,memor_msh,'LCART','chkconv',lcart)
             call memrea(npoin,memor_msh,'LPTYPE','chkconv',lptype)
             call memrea(npoin,memor_msh,'LFATH','chkconv',lfath)

             coor(1,npoin)=coor(1,ip2t)
             coor(2,npoin)=coor(2,ip2t)
             coor(3,npoin)=coor(3,ip2t)
             rtapo(1,npoin)=c00
             rtapo(2,npoin)=c00
             rtapo(3,npoin)=c00
             lpsur(npoin)=lpsur(ip2t)
             lptri(npoin)=lptri(ip2t)
             lcart(npoin)=lcart(ip2t)
             rnopo(1,npoin)=rnopo(1,ip2t)
             rnopo(2,npoin)=rnopo(2,ip2t)
             rnopo(3,npoin)=rnopo(3,ip2t)
             lptype(1,npoin)=lptype(1,ip2t)
             lptype(2,npoin)=lptype(2,ip2t)
             lfath(npoin)=ip2t  

             lsnew(1,ncont)=ip1t
             lsnew(2,ncont)=npoin

          else
             !
             !     If it is, only add the side 
             ! 
             lsnew(1,ncont)=ip1t
             lsnew(2,ncont)=ip2t

          endif
             
          llinenew(ncont)=lline(isidet)
          lsmnew(ncont)=lsmark(isidet)

       else  
          !
          !     Should find the next side with respect 
          !     to the surface orientation
          !
          !     Must check for last side  
          !
          !     First get the face with side (ip1t,ip2t)
          !  
          ipmin=min(ip1t,ip2t)
          ipmax=max(ip1t,ip2t)

          write(*,*)'Not ready corner' 
          stop

          do isto=ptoel2(ip2t),ptoel2(ip2t+1)-1
             iface=ptoel1(ip2t)
             do inofa=1,nnofa
                ipa=lface(ltab(1,inofa),iface)               
                ipb=lface(ltab(2,inofa),iface)               
                jpmin=min(ipa,ipb)
                jpmax=max(ipa,ipb)
                if(ipmin==jpmin .and. ipmax==jpmax )then 
                   goto 200 
                endif
             enddo
          enddo
          !
          !     Error
          !
          write(*,*)'Error in reorside, iside not found:',ip1t,ip2t
          stop

200       continue
          !
          !     Now roll around ip2t until finding new side
          !


       endif

    enddo


  end subroutine addside

  subroutine outaniface(nnofa,nface,npoin,ndim,lface,coor,lside,nnosi,nside,rnopo)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,nnosi,nside
    integer(ip), intent(in)      :: lface(nnofa,nface),lside(nnosi,nside)
    real(rp), intent(in)         :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh with its surface patch number 
    !

    open(unit=50,file='outaniface.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i)
    enddo

    write(50,6)
    write(50,7)
    write(50,2)
    write(50,3)
    write(50,4)
    write(50,5)
    icont=0_ip
    do i=1,nside
       icont=icont+1
       write(50,400)icont,lside(1,i),lside(2,i)
    enddo
    write(50,6)

    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
400 format(3i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
7   format('MESH dimension 3 ElemType  Linear Nnode 2')

    open(unit=50,file='outaniface.res',status='unknown')
    rewind 50

10  format('GID Post Results File 1.0')
11  format('Result "normals" "Analysis/time" 1 Vector OnNodes')
12  format('ComponentNames "Vx", "Vy" , "Vz" ')
13  format('Values')
14  format('End Values')
15  format('   ')
    write(50,10)
    write(50,11)
    write(50,12)
    write(50,13)
    do  i=1,npoin
       write(50,100)i,rnopo(1,i),rnopo(2,i),rnopo(3,i)
    enddo
    write(50,14)
    write(50,15)

    close(50)

  end subroutine outaniface

  subroutine outaniface2(nnofa,nface,npoin,ndim,lface,coor)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh with its surface patch number 
    !

    open(unit=50,file='outaniface2.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i)
    enddo

    write(50,6)

    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
400 format(3i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')

  end subroutine outaniface2

  subroutine outaniface3(nnofa,nface,npoin,ndim,lface,coor,lmark)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim
    integer(ip), intent(in)      :: lface(nnofa,nface),lmark(nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont
    !
    !     This sub displays a triangular surface mesh with its surface patch number 
    !

    open(unit=50,file='outaniface3.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i),lmark(i)
    enddo

    write(50,6)

    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
400 format(3i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')

  end subroutine outaniface3


  subroutine mrkcol(nside,nblay,lfmark,nface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only        :  memor_msh
    implicit none
    integer(ip),intent(in)    :: nside,nface,nblay
    integer(ip),intent(inout) :: lfmark(nface)
    integer(ip)               :: nfac,iside,iblay,jblay
    !
    !     Initialize the pointer
    !
    nfac=0_ip
    !
    !     Loop on sides to check elements in column no marked
    !   
    do iside=1,nside
       !
       !     Loop on layers
       !  
       do iblay=1,nblay

          nfac=nfac+1_ip
          !
          !     Has this first face been marked for deletion
          ! 
          if(lfmark(nfac)==0)then
             !
             !     Mark next faces in the column
             !
             nfac=nfac+1_ip 
             lfmark(nfac)=0_ip

             do jblay=iblay+1,nblay
                nfac=nfac+1_ip
                lfmark(nfac)=0_ip
                nfac=nfac+1_ip
                lfmark(nfac)=0_ip 
             enddo
             !
             !     And go to the next side
             ! 
             exit

          endif

          nfac=nfac+1_ip
          !
          !     Has this second face been marked for deletion
          ! 
          if(lfmark(nfac)==0)then
             !
             !     Mark next faces in the column
             !
             lfmark(nfac-1)=0_ip

             do jblay=iblay+1,nblay
                nfac=nfac+1_ip
                lfmark(nfac)=0_ip
                nfac=nfac+1_ip
                lfmark(nfac)=0_ip 
             enddo
             !
             !     And go to the next side
             ! 
             exit

          endif

       enddo

    enddo

  end subroutine mrkcol

  subroutine recoverp(ipa,ipb,lface,nnofa,nface,lptri,coor,npoin,ndim,eltoel,rnopo,rnofa,ierr)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: ipa,ipb,nnofa,nface,npoin,ndim
    integer(ip), intent(inout)   :: lface(nnofa,nface),eltoel(nnofa,nface),lptri(npoin),ierr
    real(rp), intent(in)         :: coor(ndim,npoin),rnopo(ndim,npoin)
    real(rp), intent(inout)      :: rnofa(nnofa,nface)
    integer(ip)                  :: ie,j,ip1,ip2,ip3,ienext,iopt,ielem2,iview1
    integer(ip)                  :: lpipe(100),npipe,ieold,ipipe,ienext1,ienext2,iturn,ifirst
    real(rp)                     :: rface(ndim),rx,ry,rz,rscal,pproj(ndim),pnew(ndim),p1(ndim),p2(ndim)
    real(rp)                     :: d1,d2,d3,epsil,dtot  
    integer(ip)                  :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  

    epsil=1.0d-09

    pnew(1)=coor(1,ipb)
    pnew(2)=coor(2,ipb)
    pnew(3)=coor(3,ipb)
    !
    !    Loop until the edges have been regenerated
    !
5   continue   
    !
    !     Find the element in the ball of ipa cut by edge (ipa,ipb)
    !
    ie=lptri(ipa)

    do

       if(lface(1,ie)==ipa)then
          iview1=1_ip
       else if(lface(2,ie)==ipa)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=ipa
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)

       if(ip2==ipb)goto 10
       if(ip3==ipb)goto 10
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)

       rx=pnew(1)-coor(1,ip1)
       ry=pnew(2)-coor(2,ip1)
       rz=pnew(3)-coor(3,ip1)

       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Side (ip1,ip3)
       !
       p1(1)=coor(1,ip3)-pproj(1)
       p1(2)=coor(2,ip3)-pproj(2)
       p1(3)=coor(3,ip3)-pproj(3)
       p2(1)=coor(1,ip1)-pproj(1)
       p2(2)=coor(2,ip1)-pproj(2)
       p2(3)=coor(3,ip1)-pproj(3)

       call orient3D(p1,p2,rface,d1,ndim)
       d1=d1/dtot
       if(d1<-epsil)then
          ie=eltoel(ltab(1,iview1),ie)
          cycle
       endif
       !
       !     Side (ip1,ip2)
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip1)-pproj(1)
       p2(2)=coor(2,ip1)-pproj(2)
       p2(3)=coor(3,ip1)-pproj(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot

       if(d2>epsil)then
          ie=eltoel(ltab(2,iview1),ie)
          cycle
       endif
       !
       !     Element found swap it
       !
       ielem2=eltoel(iview1,ie)

       iopt=0_ip
       call swapf(ie,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt,0_ip,lptri,rnofa)
       if(iopt==0)then
          exit
       endif

    enddo
    !
    !     At this point, we could not swap the first element
    !     Find the pipe
    !
    lpipe(1)=ie
    npipe=1_ip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Go through the pipe
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ieold=ie
    ie=ielem2

    do

       npipe=npipe+1
       lpipe(npipe)=ie

       if(eltoel(1,ie)==ieold)then
          iview1=1_ip
       else if(eltoel(2,ie)==ieold)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=lface(iview1,ie)
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)
       !
       !     Does ipb belong to this element
       !
       if(ipb==ip1)exit
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)

       rx=pnew(1)-coor(1,ip1)
       ry=pnew(2)-coor(2,ip1)
       rz=pnew(3)-coor(3,ip1)

       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Side (ip1,ip2)
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip1)-pproj(1)
       p2(2)=coor(2,ip1)-pproj(2)
       p2(3)=coor(3,ip1)-pproj(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot

       if(d2>epsil)then
          ieold=ie
          ie=eltoel(ltab(2,iview1),ie)
          cycle
       endif
       !
       !     Side (ip1,ip3)
       !
       p1(1)=coor(1,ip1)-pproj(1)
       p1(2)=coor(2,ip1)-pproj(2)
       p1(3)=coor(3,ip1)-pproj(3)
       p2(1)=coor(1,ip3)-pproj(1)
       p2(2)=coor(2,ip3)-pproj(2)
       p2(3)=coor(3,ip3)-pproj(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot

       if(d2>epsil)then
          ieold=ie
          ie=eltoel(ltab(1,iview1),ie)
          cycle
       endif
       !
       !     This is the last element
       !
       exit

    enddo
    !
    !     Swap the hole pipe except the first element
    !

    !
    !     Loop on iterations
    !
    do  
       ipipe=1_ip
       !
       !      Loop on the pipe
       !
       do    

          ipipe=ipipe+1
          if(ipipe>=npipe)exit

          ie=lpipe(ipipe)
          ielem2=lpipe(ipipe+1)

          iopt=0_ip
          call swapf(ie,ielem2,lface,nnofa,nface,eltoel,coor,rnopo,ndim,npoin,iopt,0_ip,lptri,rnofa)
          if(iopt==1)then
             call updatepipep(ipipe,ipipe+1_ip,lpipe,npipe,ipa,ipb,rnofa,nnofa,nface,coor,ndim,npoin,&
                  pnew,eltoel,lface,ierr)
             if(ierr==1)then
                return
             endif
             !
             !     We are almost done
             !
             if(npipe==2)goto 5
          endif

       enddo

    enddo

10  continue

  end subroutine recoverp

  subroutine updatepipep(ibegin,iend,lpipe,npipe,ipbegin,ipb,rnofa,nnofa,nface,coor,ndim,&
       npoin,pnew,eltoel,lface,ierr)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: ibegin,iend,ipbegin,nface,nnofa,ipb
    integer(ip),intent(in)     :: ndim,npoin 
    integer(ip),intent(in)     :: eltoel(nnofa,nface),lface(nnofa,nface)
    integer(ip),intent(inout)  :: npipe,lpipe(npipe),ierr
    real(rp),intent(in)        :: rnofa(nnofa,nface),coor(ndim,npoin),pnew(ndim)
    integer(ip)                :: ieold,ie,lploc(100),nploc,iview1,ip1,ip2,ip3,iploc 
    real(rp)                   :: rface(ndim),pproj(ndim),rscal,p1(ndim),p2(ndim),rx,ry,rz,dtot,d1,d2,epsil
    integer(ip)                :: ltab(2,3)=RESHAPE((/2,3,3,1,1,2/),(/2,3/))  

    epsil=1.0d-09

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Go through the pipe
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    nploc=0_ip

    ieold=lpipe(ibegin-1)
    if(lface(1,ieold)==ipbegin)then
       iview1=1_ip 
    else if(lface(2,ieold)==ipbegin)then
       iview1=2_ip 
    else
       iview1=3_ip 
    endif

    ie=eltoel(iview1,ieold)

    do

       if(ie==iend)exit

       if(eltoel(1,ie)==ieold)then
          iview1=1_ip
       else if(eltoel(2,ie)==ieold)then
          iview1=2_ip
       else
          iview1=3_ip
       endif

       ip1=lface(iview1,ie)
       ip2=lface(ltab(1,iview1),ie)
       ip3=lface(ltab(2,iview1),ie)
       !
       !     Remember for later
       !
       nploc=nploc+1
       lploc(nploc)=ie
       !
       !     Does ipb belong to this element
       !
       if(ipb==ip1)exit
       !
       !     Find which side is crossed
       ! 
       rface(1)=rnofa(1,ie)
       rface(2)=rnofa(2,ie)
       rface(3)=rnofa(3,ie)

       rx=pnew(1)-coor(1,ip1)
       ry=pnew(2)-coor(2,ip1)
       rz=pnew(3)-coor(3,ip1)

       rscal=rface(1)*rx+rface(2)*ry+rface(3)*rz
       pproj(1)= pnew(1)-rscal*rface(1)
       pproj(2)= pnew(2)-rscal*rface(2)
       pproj(3)= pnew(3)-rscal*rface(3)
       !
       !     Compute area to normalize
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       call orient3D(p1,p2,rface,dtot,ndim)
       !
       !     Side (ip1,ip2)
       !
       p1(1)=coor(1,ip2)-pproj(1)
       p1(2)=coor(2,ip2)-pproj(2)
       p1(3)=coor(3,ip2)-pproj(3)
       p2(1)=coor(1,ip1)-pproj(1)
       p2(2)=coor(2,ip1)-pproj(2)
       p2(3)=coor(3,ip1)-pproj(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot

       if(d2>epsil)then
          ieold=ie
          ie=eltoel(ltab(2,iview1),ie)
          cycle
       endif
       !
       !     Side (ip1,ip3)
       !
       p1(1)=coor(1,ip1)-pproj(1)
       p1(2)=coor(2,ip1)-pproj(2)
       p1(3)=coor(3,ip1)-pproj(3)
       p2(1)=coor(1,ip3)-pproj(1)
       p2(2)=coor(2,ip3)-pproj(2)
       p2(3)=coor(3,ip3)-pproj(3)

       call orient3D(p1,p2,rface,d2,ndim)
       d2=d2/dtot

       if(d2>epsil)then
          ieold=ie
          ie=eltoel(ltab(1,iview1),ie)
          cycle
       endif
       !
       !     This is the last element
       !
       exit

    enddo
    !
    !     We must have only one element
    !
    if(nploc/=1)then
       write(*,*)'Error in updatepipep, nploc=',nploc
       ierr=1_ip
       return
       !stop
    endif

    lpipe(ibegin)=lploc(1)

    do iploc=iend,npipe-1
       lpipe(iploc)=lpipe(iploc+1)
    enddo

    npipe=npipe-1

  end subroutine updatepipep

  subroutine compresf(lface,nnofa,nface,lfmark)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nnofa
    integer(ip),intent(inout) :: nface 
    integer(ip),intent(in)    :: lfmark(nface)
    integer(ip),intent(inout) :: lface(nnofa,nface)
    integer(ip)               :: iface,nfac0 
    nfac0=nface
    nface=0_ip

    do iface=1,nfac0
       if(lfmark(iface)==1)then
          nface=nface+1_ip
          lface(1,nface)=lface(1,iface)
          lface(2,nface)=lface(2,iface)
          lface(3,nface)=lface(3,iface)
       endif
    enddo

  end subroutine compresf

  subroutine compresp(lface,nnofa,nface,lfnew,nfnew,npoin,coor,lptri,ndim,lpoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nnofa,nface,nfnew,ndim
    integer(ip),intent(inout) :: npoin
    integer(ip),intent(in)    :: lfnew(nnofa,nfnew)
    integer(ip),intent(in)    :: lface(nnofa,nface)
    integer(ip),intent(inout) :: lptri(npoin),lpoin(npoin)
    real(rp),intent(inout)    :: coor(ndim,npoin)
    integer(ip)               :: iface,npoi0,ipoin,inofa 
    !
    !     This sub compress the point database
    !
    do ipoin=1,npoin 
       lpoin(ipoin)=0_ip
    enddo

    do iface=1,nface
       do inofa=1,nnofa
          lpoin(lface(inofa,iface))=1_ip
       enddo
    enddo

    do iface=1,nfnew
       do inofa=1,nnofa
          lpoin(lfnew(inofa,iface))=1_ip
       enddo
    enddo

    npoi0=npoin
    npoin=0_ip
    do ipoin=1,npoi0
       if(lpoin(ipoin)==1)then
          npoin=npoin+1_ip
          coor(1,npoin)=coor(1,ipoin)
          coor(2,npoin)=coor(2,ipoin)
          coor(3,npoin)=coor(3,ipoin)
          lptri(npoin)=lptri(ipoin)
       endif
    enddo

  end subroutine compresp

  subroutine filtcol(lstackp,nstackp,npoin,lpoin,ipoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: npoin,ipoin
    integer(ip),intent(inout) :: nstackp
    integer(ip),intent(in)    :: lpoin(npoin)
    integer(ip),intent(inout) :: lstackp(nstackp)
    integer(ip)               :: ipbase,jpbase,nstackp0,jpoin,istackp
    !
    !     This sub takes out from lstack the points that belong to the same column
    !   
    !
    !   Get the initial point of the point currently tested 
    !
    ipbase=lpoin(ipoin) 
    !
    !     Take out the points that have the same base point
    !
    nstackp0=nstackp
    nstackp=0
    do istackp=1,nstackp0
       jpoin=lstackp(istackp)
       !
       !     Is this point the base point of the column?
       !
       if(jpoin/=ipbase)then

          jpbase=lpoin(jpoin)
          !
          !     Does this point belong to the same column than ipoin?
          !  
          if(jpbase/=ipbase)then
             nstackp=nstackp+1_ip
             lstackp(nstackp)=jpoin 
          endif
       endif
    enddo

  end subroutine filtcol

  subroutine onemore(nnofa,nface,lface,lfmark,lplay,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)    :: nface,nnofa,npoin
    integer(ip),intent(inout) :: lfmark(nface)
    integer(ip),intent(in)    :: lface(nnofa,nface),lplay(npoin)
    integer(ip)               :: iface,nfac0,jface,ip1b,jp1b,iflay,jflay 
    !
    !     This subroutine marks all the faces for which the next faces 
    !     have been already marked
    !
    !
    !    Loop on faces
    !  
    iface=-1_ip

    do 
       !
       !     Find the next element marked for deletion
       !
       iface=iface+2_ip
       !
       !     Did we exhaust the list?
       !
       if(iface>nface)exit
       !
       !     Is the face still active?
       !
       if(lfmark(iface)==1)cycle 
       !
       !     We have found a marked face
       !     Get iflay
       ! 
       ip1b=lface(1,iface)
       iflay=lplay(ip1b) 
       !
       !     Now go two steps before 
       !
       jface=iface-2_ip
       !
       !     Are we still in the range?
       !
       if(jface<1)cycle
       !
       !     Has jface been deleted?
       !
       if(lfmark(jface)==0)cycle
       !
       !     Get jflay
       ! 
       jp1b=lface(1,jface)
       jflay=lplay(jp1b) 
       !
       !     Is iface above jface?
       !
       if(jflay>=iflay)cycle
       !
       !     Mark jface
       !
       lfmark(jface)=0_ip
       lfmark(jface+1)=0_ip

    enddo

  end subroutine onemore

  subroutine anicol(coor,ndim,npoin,lface,nface,nnofa,lptype,ptoel1,ptoel2,&
       rsize,rnopo,eltoel,nside,&
       nnosi,tolcol,tolcolr,tolref,tolrefr,lside,lsmark,lptri,nboup)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)        :: ndim,nnofa,nnosi,nboup
    integer(ip),intent(inout)     :: nface,npoin
    integer(ip),intent(inout)     :: lptype(2,npoin),lptri(npoin)
    integer(ip),intent(inout)     :: lface(nnofa,nface)
    integer(ip),pointer           :: ptoel1(:),ptoel2(:),ptosi1(:),ptosi2(:),eltoel(:,:)
    integer(ip),pointer           :: lside(:,:),lline(:),lsmark(:),lsurf(:)
    integer(ip),intent(inout)     :: nside
    real(rp),   intent(in)        :: tolcol,tolcolr,tolref,tolrefr
    real(rp),   intent(inout)     :: rsize(npoin)
    real(rp),   intent(inout)     :: coor(ndim,npoin),rnopo(ndim,npoin)
    integer(ip),pointer           :: lstack(:),lmark(:),lelem(:),ledge(:,:),ledg2(:)
    real(rp),pointer              :: redge(:)
    integer(ip)                   :: nstack,ncol,ipoin,iface,isurf,nfloc,ipnt,iboun 
    integer(ip)                   :: ie,ielem,nfold,iside,j,iter,niter,nsurf,nline 
    integer(ip)                   :: nfac0,ip1,ip2 
    integer(4)                    :: istat
    real(rp)                      :: rmaxloc,tolgeo 
    !
    !     This subroutine drives the collapsing by marking and collapsing the surface mesh 
    !
    !
    !     Initialize tolgeo
    !
    tolgeo=0.1d+00
    !
    !     Allocate local arrays
    !
    allocate(lstack(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LSTACK','colglo',lstack) 
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','colglo',lmark) 
    allocate(lelem(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LELEM','colglo',lelem) 
    allocate(lsurf(nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LSURF','colglo',lsurf) 
    allocate(lline(nside),stat=istat)
    call memchk(zero,istat,memor_msh,'LLINE','colglo',lline) 
    allocate(ledge(2,3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGE','reglo',ledge)
    allocate(redge(3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'REDGE','reglo',redge)
    allocate(ledg2(3*nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDG2','reglo',ledg2)
    nullify(ptosi1,ptosi2)
    !
    !     Fake lsurf && lline
    !
    do iface=1,nface
       lsurf(iface)=1_ip
    enddo
    do iside=1,nside
       lline(iside)=1_ip
    enddo
    nsurf=1_ip
    nline=1_ip
    !
    !     Mark ridge points 
    !
    do iside=1,nside
       ip1=lside(1,iside)
       ip2=lside(2,iside)
       lptype(1,ip1)=ID_RIDGE
       lptype(2,ip1)=1_ip         ! fake line number
    enddo
    !
    !     Initialize exit condition
    !   
    nfold=nface
    !
    !     Iter on collapsing iterations
    !     Each iteration is meant to divide the point mesh size by 2
    !     in a multilevel approach
    !
    iter=1_ip
    !
    !     Get the sides surrounding points
    ! 
    call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)

    do 
       write(*,*)'Coarsening Iteration for anisotropic mesh',iter,'nface=',nface
       !
       !     Get the faces surrounding the points for the whole triangulation
       !
       call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2 )
       !
       !     Mark the points to be collapsed on the whole triangulation with lmark
       !
       call mark(lstack,nstack,ptoel1,ptoel2,nnofa,lface,nface,npoin,&
            lptype,lmark,ncol,coor,ndim,rsize,nline,nside,nnosi,nsurf,lsurf,&
            ptosi1,ptosi2,lline,lside)
       !call outmark(lface,nnofa,nface,lelem,coor,npoin,ndim,lmark)
       !
       !     Collapse edges inside each patch
       !
       call colsmo(lface,nnofa,nface,lptype,lmark,npoin,lelem,coor,ndim,ptoel1,ptoel2,&
            rnopo,tolcol,rsize,eltoel,rmaxloc,tolref,ledge,redge,ledg2,lstack,tolgeo)
       !
       !     Compress the faces
       !
       nfac0=nface  
       nface=0_ip
       do iface=1,nfac0
          if(lelem(iface)/=-1)then
             nface=nface+1_ip
             lface(1,nface)=lface(1,iface)
             lface(2,nface)=lface(2,iface)
             lface(3,nface)=lface(3,iface)
             lelem(nface)=0_ip
          endif
       enddo
       !
       !     Did we do something?
       !         
       if(nface==nfold)exit
       nfold=nface
       iter=iter+1
       !
       !     Get the faces surrounding the points
       !
       call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
       !
       !     Get the faces surrounding faces
       !
       call trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
       !
       !     Optimize the whole mesh by swapping
       !
       niter=1_ip
       !
       !     Do we have sides ?
       !
       if(nboup>0)then
          !
          !     Get the sides surrounding the points
          !
          call ptoelm(lside,nside,npoin,nnosi,ptosi1,ptosi2)
       else
          call memrea(npoin+1_ip,memor_msh,'PTOSI2','refsplit',ptosi2)
          do ipoin=1,npoin+1
             ptosi2(ipoin)=0_ip
          enddo

       endif

       call swapsurf(lface,nnofa,nface,coor,npoin,ndim,rnopo,eltoel,niter,ptosi1,&
             ptosi2,lside,nnosi,nside)
       !
       !     Clean up lmark
       !
       do ipoin=1,npoin
          lmark(ipoin)=0_ip
       enddo

    enddo
    !
    !     Update lptri
    !
    do iface=1,nface 
       lptri(lface(1,iface))=iface
       lptri(lface(2,iface))=iface
       lptri(lface(3,iface))=iface
    enddo
    !
    !     Free memory
    !
    call memchk(2_ip,istat,memor_msh,'LSURF','colglo',lsurf)
    deallocate(lsurf,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSURF','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LELEM','colglo',lelem)
    deallocate(lelem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LELEM','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','colglo',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LSTACK','colglo',lstack)
    deallocate(lstack,stat=istat)
    if(istat/=0) call memerr(2_ip,'LSTACK','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI1','colglo',ptosi1)
    deallocate(ptosi1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI1','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOSI2','colglo',ptosi2)
    deallocate(ptosi2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOSI2','colglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'REDGE','refglo',redge)
    deallocate(redge,stat=istat)
    if(istat/=0) call memerr(2_ip,'REDGE','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDG2','refglo',ledg2)
    deallocate(ledg2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDG2','refglo',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','refglo',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','refglo',0_ip)

  end subroutine anicol

  subroutine filtvirt(npoin,iface,nface,nside,nnosi,lfsid,lsmark,nstackp,lstackp,&
       lside,sitosi,lmark,lpoin,lfath)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only     :  memor_msh
    use def_meshin, only        :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    implicit none
    integer(ip),intent(in)    :: iface,nface,nside,nnosi,npoin        
    integer(ip),intent(in)    :: lfsid(nface),lsmark(nside),lpoin(npoin) 
    integer(ip),intent(in)    :: lside(nnosi,nside),sitosi(nnosi,nside),lfath(npoin)
    integer(ip),intent(inout) :: nstackp
    integer(ip),intent(inout) :: lmark(npoin),lstackp(nstackp)
    integer(ip)               :: iside,inosi,ip1,ip2,jside,ineigh,jneigh,nstackp0
    integer(ip)               :: nstackl,istackl,ipbase,istackp,ipoin,ipfath 
    integer(ip),parameter     :: mstackl=500
    integer(ip)               :: lstackl(mstackl)       
    !
    !     This sub filters the close virtual points neighbors of this side 
    !       
    nstackl=0_ip
    !
    !     DBG
    ! 
    do ipoin=1,npoin 
       if(lmark(ipoin)/=0)then
          write(*,*)'Error in filtvirt'
          stop
       endif
    enddo
    !
    !     Original side of this face
    !
    iside=lfsid(iface)
    !
    !     Mark the points 
    ! 
    ip1=lside(1,iside)
    ip2=lside(2,iside)
    lmark(ip1)=1_ip  
    lmark(ip2)=1_ip 
    nstackl=2
    lstackl(1)=ip1   
    lstackl(2)=ip2   
    !
    !     Look at surrounding neighbors
    !
    do inosi=1,nnosi

       ineigh=sitosi(inosi,iside)
       if(ineigh==0)cycle
       if(lsmark(ineigh)/=2)cycle
       !
       !   Mark the points of this side
       !
       ip1=lside(1,ineigh)
       if(lmark(ip1)==0)then
          lmark(ip1)=1_ip  
          nstackl=nstackl+1
          lstackl(nstackl)=ip1
       endif
       ip2=lside(2,ineigh)
       if(lmark(ip2)==0)then
          lmark(ip2)=1_ip 
          nstackl=nstackl+1
          lstackl(nstackl)=ip2
       endif
       jside=ineigh 
       !
       !     Loop on neighbor virtual sides
       !
       do 

          jneigh=sitosi(inosi,jside)
          !
          !     Do we have a neighbor?  
          ! 
          if(jneigh==0)exit
          !
          !     Is this a virtual side?  
          ! 
          if(lsmark(jneigh)/=2)exit
          !
          !   Mark the points of this side
          !
          ip1=lside(1,jneigh)
          if(lmark(ip1)==0)then
             lmark(ip1)=1_ip  
             nstackl=nstackl+1
             lstackl(nstackl)=ip1
          endif
          !
          !   For duplicated sides, mark father 
          !
          ipfath=lfath(ip1)
          if(ipfath/=0)then
             if(lmark(ipfath)==0)then
                lmark(ipfath)=1_ip  
                nstackl=nstackl+1
                lstackl(nstackl)=ipfath
             endif
          endif

          ip2=lside(2,jneigh)
          if(lmark(ip2)==0)then
             lmark(ip2)=1_ip 
             nstackl=nstackl+1
             lstackl(nstackl)=ip2
          endif
          !
          !   For duplicated sides, mark father 
          !
          ipfath=lfath(ip2)
          if(ipfath/=0)then
             if(lmark(ipfath)==0)then
                lmark(ipfath)=1_ip  
                nstackl=nstackl+1
                lstackl(nstackl)=ipfath
             endif
          endif
          jside=jneigh

       enddo

    enddo
    !
    !     Filter the points
    !
    nstackp0=nstackp
    nstackp=0_ip
    do istackp=1,nstackp0
       ipoin=lstackp(istackp)
       ipbase=lpoin(ipoin)
       !
       !     Has the base point of ipoin or ipoin (for the first level)
       !     been marked?
       !     Has the father of ipoin been marked?
       !
       if(lmark(ipbase)/=1 .and. lmark(ipoin)/=1 )then
           !
           !     Check for points of duplicated sides
           !
           ipfath=lfath(ipoin)
           if(ipfath==0)then
              nstackp=nstackp+1
              lstackp(nstackp)=ipoin
           else
              if(lmark(ipfath)/=1)then
                 nstackp=nstackp+1
                 lstackp(nstackp)=ipoin
              endif 
           endif
       endif
    enddo
    !
    !     Clean up 
    !
    do istackl=1,nstackl
       lmark(lstackl(istackl))=0_ip
    enddo

  end subroutine filtvirt

end module mod_anisurf
