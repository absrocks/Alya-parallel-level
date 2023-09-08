subroutine elsest_binpre(&
     ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltype,coord,rpara)
  !
  ! Construct the bin structure
  !
  use def_elsest
  use mod_elsest
  use def_master, only : kfl_paral
  implicit none
  integer(ip), intent(in)  :: ipara(*),imesh,ithre
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode(*)
  integer(ip), intent(in)  :: lnods(mnode,*),ltype(*)
  real(rp),    intent(in)  :: coord(ndime,npoin),rpara(*)
  integer(ip)              :: idime,iboxe,istat,i,j,kk,jj,ii,ll
  integer(ip)              :: imin,imax,jmin,jmax,kmin,kmax,ni,nj,nk
  integer(ip)              :: ielem,box_nr,box_nr1,box_nr2,ninj
  real(rp),    pointer     :: xmima(:,:,:)
  real(rp)                 :: time1,time2,time3,time4
  real(rp)                 :: deltx,delty,deltz,dni,dnj,dnk

  call elsest_cputim(time1)
  bin_struc(imesh)%iallo = 1

  !----------------------------------------------------------------------
  !
  ! Initialize parameters
  !
  !----------------------------------------------------------------------

  !*OMP PARALLEL
  !*OMP SECTIONS
  !*OMP SECTION
  allocate(bin_struc(imesh)%cputi(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'CPUTI','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(bin_struc(imesh)%memor(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'MEMOR','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(bin_struc(imesh)%kstat(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSTAT','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(bin_struc(imesh)%ksear(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSEAR','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(bin_struc(imesh)%kfirs(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KFIRS','elsest_alloc',0_ip)
  !*OMP SECTION
  allocate(bin_struc(imesh)%kseco(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSECO','elsest_alloc',0_ip)
  !*OMP END SECTIONS
  !*OMP END PARALLEL

  !----------------------------------------------------------------------
  !
  ! Point to current mesh (IMESH) structure
  !
  !----------------------------------------------------------------------

  call elsest_binpoi(imesh)
  !*OMP PARALLEL DO PRIVATE(i,j)
  do i=1,nthre
     do j=1,10
        memor(j,i) = 0_ip
        cputi(j,i) = 0.0_rp
        kstat(j,i) = 0_ip
     end do
     kstat(1,i)    = huge(1_ip)
     kstat(3,i)    = huge(1_ip)
     ksear(i)      = 0_ip
     kfirs(i)      = 0_ip
     kseco(i)      = 0_ip
  end do
  !*OMP END PARALLEL DO

  memax      = 0_ip
  nboxx(1)   = ipara(1)
  nboxx(2)   = ipara(2)
  nboxx(3)   = ipara(3)  
  dataf      = ipara(4)
  comin(1)   = 0.0_rp
  comin(2)   = 0.0_rp
  comin(3)   = 0.0_rp
  comax(1)   = 0.0_rp
  comax(2)   = 0.0_rp
  comax(3)   = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! Total number of boxes
  !
  !----------------------------------------------------------------------

  if( ndime == 1 ) then
     nboxe = nboxx(1)
  else if( ndime == 2 ) then
     nboxe = nboxx(1) * nboxx(2)
  else
     nboxe = nboxx(1) * nboxx(2) * nboxx(3)
  end if
  ni   = nboxx(1)
  nj   = nboxx(2)
  nk   = nboxx(3)
  bin_struc(imesh) % nboxx(1) = ni
  bin_struc(imesh) % nboxx(2) = nj
  bin_struc(imesh) % nboxx(3) = nk

  !----------------------------------------------------------------------
  !
  ! Allocate memory for bin structure
  !
  !----------------------------------------------------------------------

  allocate(bin_struc(imesh) % elcod(ndime,mnode,nthre),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_binpre',bin_struc(imesh)%elcod)
  if( dataf == 0 ) then
     allocate(bin_struc(imesh) % tboel(nboxe),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'TBOEL','elsest_binpre',bin_struc(imesh)%tboel)
  else if( dataf == 1 ) then
     allocate(bin_struc(imesh) % pboel(nboxe+1),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'PBOEL','elsest_binpre',bin_struc(imesh)%pboel)
  end if

  elcod => bin_struc(imesh) % elcod
  tboel => bin_struc(imesh) % tboel
  pboel => bin_struc(imesh) % pboel
  lboel => bin_struc(imesh) % lboel


  !----------------------------------------------------------------------
  !
  ! Compute bounding box and bin spacing delta
  !
  !----------------------------------------------------------------------

  call elsest_boubox(ndime,npoin,coord,comin,comax)
  !
  ! Bin spacing DELTA
  !
  do idime = 1,ndime
     delta(idime) = (comax(idime)-comin(idime))/real(nboxx(idime),rp)
  end do
  deltx = 1.0_rp / ( comax(1)-comin(1) )
  delty = 1.0_rp / ( comax(2)-comin(2) )
  dni   = real(ni,rp) * deltx
  dnj   = real(nj,rp) * delty
  if( ndime == 3 ) then
     deltz = 1.0_rp / ( comax(3)-comin(3) )
     dnk   = real(nk,rp) * deltz
  end if
  ninj  = nboxx(1) * nboxx(2)

  call elsest_cputim(time2)
  cputi(1,ithre)=time2-time1

  !----------------------------------------------------------------------
  !
  ! Compute element bounding boxes
  !
  !----------------------------------------------------------------------

  allocate( xmima(3,2,nelem), stat = istat )
  call elsest_elbbox(&
       mnode,ndime,npoin,nelem,nnode,lnods,ltype,coord,xmima)

  !----------------------------------------------------------------------
  !
  ! Number of elements per box
  !
  !----------------------------------------------------------------------

  allocate(nbono(nboxe),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'NBONO','elsest_binpre',nbono)

  if( ndime == 1 ) then

     do ielem = 1,nelem

        imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
        imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1  

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))

        do ii = imin,imax
           box_nr = ii
           nbono(box_nr) = nbono(box_nr) + 1
        end do
     end do

  else if( ndime == 2 ) then

     do ielem = 1,nelem

        imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
        imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1      

        jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1
        jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1      

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))
        jmin = max(jmin,1_ip)
        jmax = min(jmax,nboxx(2))

        do ii = imin,imax
           do jj = jmin,jmax
              box_nr = (jj-1_ip) * ni + ii
              nbono(box_nr) = nbono(box_nr) + 1
           end do
        end do
     end do

  else if( ndime == 3 ) then

     do ielem = 1,nelem

        imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip )  + 1
        imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip )  + 1      
                                                                 
        jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip )  + 1
        jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip )  + 1      
                                                                 
        kmin = int( ( xmima(3,1,ielem) - comin(3)) * dnk , ip )  + 1
        kmax = int( ( xmima(3,2,ielem) - comin(3)) * dnk , ip )  + 1      

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))
        jmin = max(jmin,1_ip)
        jmax = min(jmax,nboxx(2))
        kmin = max(kmin,1_ip)
        kmax = min(kmax,nboxx(3))

        do kk = kmin,kmax
           box_nr2 = ninj * (kk-1)
           do jj = jmin,jmax
              box_nr1 = box_nr2 + nboxx(1) * (jj-1)
              do ii = imin,imax
                 box_nr        = box_nr1 + ii
                 nbono(box_nr) = nbono(box_nr) + 1
              end do
           end do
        end do

     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Fill in box element list
  !
  !----------------------------------------------------------------------

  if( dataf == 0 ) then
     !
     ! Type
     !
     do iboxe = 1,bin_struc(imesh) % nboxe
        allocate( bin_struc(imesh) % tboel(iboxe)%l(nbono(iboxe)) , stat = istat )
        nbono(iboxe) = 0
     end do

     if( ndime == 1 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))

           do ii = imin,imax
              box_nr        = ii
              nbono(box_nr) = nbono(box_nr) + 1
              bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
           end do

        end do

     else if( ndime == 2 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1      

           jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1
           jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1    

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))

           do ii = imin,imax
              do jj = jmin,jmax
                 box_nr        = (jj-1_ip) * ni + ii
                 nbono(box_nr) = nbono(box_nr) + 1
                 bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
              end do

           end do

        end do

     else

        do ielem = 1,nelem

           imin = int( ( xmima(1,1,ielem) - comin(1) ) * dni , ip ) + 1
           imax = int( ( xmima(1,2,ielem) - comin(1) ) * dni , ip ) + 1      

           jmin = int( ( xmima(2,1,ielem) - comin(2) ) * dnj , ip ) + 1
           jmax = int( ( xmima(2,2,ielem) - comin(2) ) * dnj , ip ) + 1      

           kmin = int( ( xmima(3,1,ielem) - comin(3) ) * dnk , ip ) + 1
           kmax = int( ( xmima(3,2,ielem) - comin(3) ) * dnk , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))
           kmin = max(kmin,1_ip)
           kmax = min(kmax,nboxx(3))

           do kk = kmin,kmax
              box_nr2 = ninj * (kk-1)
              do jj = jmin,jmax
                 box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                 do ii = imin,imax
                    box_nr        = box_nr1 + ii
                    nbono(box_nr) = nbono(box_nr) + 1
                    bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
                 end do
              end do

           end do

        end do

     end if

  else
     !
     ! Linked list
     !
     bin_struc(imesh) % pboel(1) = 1

     do iboxe = 1,bin_struc(imesh) % nboxe
        bin_struc(imesh) % pboel(iboxe+1) = bin_struc(imesh) % pboel(iboxe) + nbono(iboxe) 
        nbono(iboxe) = 0
     end do

     allocate(bin_struc(imesh)%lboel(bin_struc(imesh) % pboel(bin_struc(imesh) % nboxe+1)),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'LBOEL','elsest_binlis',bin_struc(imesh)%lboel)
     lboel => bin_struc(imesh) % lboel

     if( ndime == 1 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))

           do ii = imin,imax
              box_nr        = ii
              ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
              nbono(box_nr) = nbono(box_nr) + 1
              lboel(ll)     = ielem
           end do
        end do

     else if( ndime == 2 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * real(ni,rp) , ip ) + 1      

           jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1
           jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * real(nj,rp) , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))

           do ii = imin,imax
              do jj = jmin,jmax
                 box_nr        = (jj-1_ip) * ni + ii
                 ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                 nbono(box_nr) = nbono(box_nr) + 1
                 lboel(ll)     = ielem
              end do
           end do
        end do

     else if( ndime == 3 ) then

        do ielem = 1,nelem

           imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip ) + 1
           imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip ) + 1      

           jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip ) + 1
           jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip ) + 1      

           kmin = int( ( xmima(3,1,ielem) - comin(3)) * dnk , ip ) + 1
           kmax = int( ( xmima(3,2,ielem) - comin(3)) * dnk , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))
           kmin = max(kmin,1_ip)
           kmax = min(kmax,nboxx(3))

           do kk = kmin,kmax
              box_nr2 = ninj * (kk-1)
              do jj = jmin,jmax
                 box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                 do ii = imin,imax
                    box_nr        = box_nr1 + ii
                    ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                    nbono(box_nr) = nbono(box_nr) + 1
                    lboel(ll)     = ielem
                 end do
              end do
           end do
        end do

     end if

  end if

  !
  ! Deallocate memory
  !
  call elsest_cputim(time3)
  deallocate( xmima, stat = istat )
  call elsest_memchk(2_ip,ithre,istat,memor(1,ithre),'NBONO','elsest_binpre',nbono)
  deallocate(nbono,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'NBONO','elsest_binpre',0_ip)
  call elsest_cputim(time4)
  cputi(3,ithre)=time4-time3
  !
  ! Post-process bin
  !
  call elsest_binpos(ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)
  !
  ! Output statistics
  !
  if( ipara(7) /= 0 ) call elsest_statis(1_ip,imesh,ipara,ithre)

  call elsest_binpoi(imesh)

end subroutine elsest_binpre
