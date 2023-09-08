subroutine elsest_binpro(&
     imesh,lmesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltype,ltopo,coord,point_x,rpara,ifoun,shapt,&
     derit,coloc,lchec)
  !
  ! Bin search: look for host element
  !
  use def_elsest
  use mod_elsest
  use mod_elmgeo, only : elmgeo_shapf_deriv_heslo 
  use mod_elmgeo, only : elmgeo_natural_coordinates
  implicit none
  integer(ip), intent(in)    :: imesh,lmesh,ipara(*),ithre
  integer(ip), intent(in)    :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)    :: nnode(*)
  integer(ip), intent(in)    :: lnods(mnode,nelem),ltype(*)
  integer(ip), intent(in)    :: ltopo(*)
  real(rp),    intent(in)    :: coord(ndime,npoin),point_x(*)
  real(rp),    intent(inout) :: rpara(*)
  integer(ip), intent(out)   :: ifoun
  real(rp),    intent(out)   :: shapt(*),derit(*),coloc(*)
  integer(ip), intent(in)    :: lchec(*)
  integer(ip)                :: curr_box_coor(3),box_nr
  integer(ip)                :: array_size,ielem,ii,inode,ipoin,ilook,kk
  integer(ip)                :: pnode,pelty,idime,kelem,ptopo
  real(rp)                   :: time1,time2,time3
  real(rp)                   :: point_y(3),coloc_max
  !
  ! If not allocated, create structure of mesh IMESH
  !
  if( kfl_memor(1) == 0 ) call elsest_alloca(1_ip,ipara)
  if( bin_struc(imesh) % iallo == 0 ) then
     call elsest_binpre(&
          ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
          lnods,ltype,coord,rpara)    
  end if

  if( imesh /= lmesh ) call elsest_binpoi(imesh)
 
  lmini = -rpara(1)
  lmaxi = 1.0_rp+rpara(1)

  ksear(ithre) = ksear(ithre) + 1
  call elsest_cputim(time1)
  !
  ! Check if point is outside the bounding box
  !
  if( ipara(15) == 1 ) then
     point_y(1:ndime) = min(point_x(1:ndime),comax(1:ndime))
     point_y(1:ndime) = max(point_y(1:ndime),comin(1:ndime))
  else
     do idime = 1,ndime
        if( point_x(idime) < comin(idime) .or. point_x(idime) > comax(idime) ) then
           ifoun = -1
           return
        end if
     end do
     point_y(1:ndime) = point_x(1:ndime)
  end if
  !
  ! Determine in which box (i,j,k) the point lies: curr_box_coor
  !
  call elsest_boxijk(ndime,nboxx,point_y,curr_box_coor,comin,comax)
  call elsest_boxnum(ndime,nboxx,curr_box_coor,box_nr)
  !
  ! Geometric CHECK
  !
  if( dataf == 0 ) then 
     array_size = size(bin_struc(imesh) % tboel(box_nr) % l,1)
  else if( dataf == 1 ) then
     array_size = bin_struc(imesh) % pboel(box_nr+1) - bin_struc(imesh) % pboel(box_nr)
     kk         = bin_struc(imesh) % pboel(box_nr) - 1
  end if
  !
  ! First try: Loop over elements in box
  !
  ii       = 0
  ifoun    = 0
  kelem    = 0
  rpara(2) = huge(1.0_rp)
  do while( ifoun == 0 .and. ii < array_size )
     ii = ii + 1
     if( dataf == 0 ) then
        ielem = bin_struc(imesh) % tboel(box_nr) % l(ii)
     else if( dataf == 1 ) then
        ielem = bin_struc(imesh) % lboel(ii+kk)
     end if

     ilook = 1
     if( ipara(14) /= 0 ) then
        if( lchec(ielem) /= ipara(14) ) ilook = 0
     end if
     pelty = ltype(ielem)

     if( ilook == 1 .and. pelty > 0 ) then
        pnode = nnode(pelty)
        do inode = 1,pnode
           ipoin= lnods(inode,ielem)
           do idime = 1,ndime
              bin_struc(imesh) % elcod(idime,inode,ithre) = coord(idime,ipoin)
           end do
        end do

        call elmgeo_natural_coordinates(                       &
             ndime,pelty,pnode,bin_struc(imesh)%elcod(:,:,ithre),&
             shapt,derit,point_y,coloc,ifoun,abs(lmini))

        if( ifoun > 0 ) then
           ifoun        = ielem
           kfirs(ithre) = kfirs(ithre)+1

        else if( ifoun <= 0 .and. ipara(15) == 1 ) then
           ptopo = ltopo(pelty)
           if( ptopo == 0 .or. ptopo == -1 .or. ptopo == 3 ) then 
              !
              ! Lines, Hexa, Pyra
              !
              coloc_max = 0.0_rp
              do idime = 1,ndime
                 coloc_max = coloc_max + coloc(idime)**2
              end do
           else
              !
              ! Other elements
              !
              coloc_max = 0.0_rp
              do idime = 1,ndime
                 coloc_max = coloc_max + (coloc(idime)-0.5_rp)**2
              end do
           end if
           if( coloc_max <= rpara(2) ) then
              kelem = ielem 
              rpara(2) = coloc_max
           end if

        end if
     end if

  end do
  !
  ! If an element is really required, copy information back
  !
  if( ipara(15) == 1 ) then 
     if( ifoun == 0 ) then
        ifoun = kelem
        if( ifoun /= 0 ) then
           pelty = ltype(kelem)
           if( pelty > 0 ) then
              ptopo = ltopo(pelty)
              coloc(1:ndime) = min(coloc(1:ndime),1.0_rp)
              if( ptopo == 0 .or. ptopo == -1 .or. ptopo == 3 ) then 
                 coloc(1:ndime) = max(coloc(1:ndime),-1.0_rp) 
              else
                 coloc(1:ndime) = max(coloc(1:ndime), 0.0_rp)        
                 !call runend('CHKEAR')
              end if
              call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapt,derit)
           end if
        end if
     else
        rpara(2) = -1.0_rp
     end if
  end if

  call elsest_cputim(time2)
  cputi(6,ithre)=cputi(6,ithre)+time2-time1

  call elsest_cputim(time3)
  cputi(7,ithre)=cputi(7,ithre)+time3-time2


end subroutine elsest_binpro
