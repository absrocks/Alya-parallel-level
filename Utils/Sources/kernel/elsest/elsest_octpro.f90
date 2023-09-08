subroutine elsest_octpro(&
     imesh,lmesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltype,ltopo,coord,point_x,rpara,ifoun,shapt,&
     derit,coloc,lchec)
  !
  ! Oct search: look for host element
  !
  use def_elsest, only : ip,rp,poiarr,octbox
  use def_elsest, only : oct_struc
  use def_elsest, only : kfl_memor
  use mod_elmgeo, only : elmgeo_shapf_deriv_heslo 
  use mod_elmgeo, only : elmgeo_natural_coordinates
  implicit none
  integer(ip), intent(in)    :: imesh,lmesh,ipara(*),ithre
  integer(ip), intent(in)    :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)    :: nnode(*)
  integer(ip), intent(in)    :: lnods(mnode,nelem),ltype(*)
  integer(ip), intent(in)    :: ltopo(*)
  integer(ip), intent(out)   :: ifoun
  real(rp),    intent(in)    :: coord(ndime,npoin),point_x(*)
  real(rp),    intent(inout) :: rpara(*)
  real(rp),    intent(out)   :: shapt(*),derit(*),coloc(*)
  integer(ip), intent(in)    :: lchec(*)
  integer(ip)                :: array_size,ielem,ii,inode,ipoin
  integer(ip)                :: pnode,pelty,idime,ilook,kelem
  integer(ip)                :: ptopo,ichild,igene
  real(rp)                   :: time1,time2,time3,elcod(ndime,mnode)
  real(rp)                   :: point_y(3),coloc_max,lmini,rpara2
  real(rp),     pointer      :: comin(:)
  real(rp),     pointer      :: comax(:)
  type(octbox), pointer      :: current_o
  !
  ! If not allocated, create structure of mesh IMESH
  !
  if( kfl_memor(2) == 0 ) call elsest_alloca(2_ip,ipara)

  if( oct_struc(imesh) % iallo == 0 ) then
     call runend('ELSEST_OCTPRO: STRUCTURE HAS NOT BEEN ALLOCATED')
  end if 
  !
  ! This is a new search
  !
  comin     => oct_struc(imesh) % comin
  comax     => oct_struc(imesh) % comax
  lmini     =  -rpara(1)
  !
  !
  !call elsest_octpoi(imesh)
  !!!!ksear(ithre) = ksear(ithre) + 1
  !!!!call elsest_cputim(time1) 
  !
  ! Check if point is outside the bounding box
  !
  if( ipara(15) == 1 ) then
     point_y(1:ndime) = min( point_x(1:ndime) , comax(1:ndime) )
     point_y(1:ndime) = max( point_y(1:ndime) , comin(1:ndime) )
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
  ! Find the quad/oct where the point lies
  !
  !call elsest_octfin(ndime,ithre,point_y,current)

  current_o  => oct_struc(imesh) % tree_root
  igene = 0

  if( ndime == 3 ) then
     
     do while( current_o % whoiam == 0 )    
        igene = igene + 1
        childloop8: do ichild = 1,8           
           if(    point_y(1) >= current_o % children(ichild) % minc(1) .and. &
                & point_y(1) <= current_o % children(ichild) % maxc(1) .and. &
                & point_y(2) >= current_o % children(ichild) % minc(2) .and. &
                & point_y(2) <= current_o % children(ichild) % maxc(2) .and. &
                & point_y(3) >= current_o % children(ichild) % minc(3) .and. &
                & point_y(3) <= current_o % children(ichild) % maxc(3) ) then
              current_o => current_o % children(ichild)
              exit childloop8
           end if
        end do childloop8
     end do

  else if( ndime == 2 ) then

     do while( current_o % whoiam == 0 )    
        igene = igene + 1
        childloop4: do ichild = 1,4
           if(    point_y(1) >= current_o % children(ichild) % minc(1) .and. &
                & point_y(1) <= current_o % children(ichild) % maxc(1) .and. &
                & point_y(2) >= current_o % children(ichild) % minc(2) .and. &
                & point_y(2) <= current_o % children(ichild) % maxc(2)  ) then
              current_o => current_o % children(ichild)
              exit childloop4
           end if
        end do childloop4
     end do

  end if
  ! 
  ! Perform search over elements inside current box
  !
  array_size = current_o % nelembox
  ii         = 0
  ifoun      = 0
  kelem      = 0
  rpara2     = huge(1.0_rp)

  do while( ifoun == 0 .and. ii < array_size )
     ii    = ii+1
     ielem = current_o % elems(ii)

     ilook = 1
     if( ipara(14) /= 0 ) then
        if( lchec(ielem) /= ipara(14) ) ilook = 0
     end if
     
     if( ilook == 1 ) then
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           ptopo = ltopo(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
           call elmgeo_natural_coordinates(          &
                ndime,pelty,pnode,elcod,shapt,derit, &
                point_y,coloc,ifoun,abs(lmini))

           if( ifoun > 0 ) then
              ifoun        = ielem
              !!!!kfirs(ithre) = kfirs(ithre)+1
           else if( ipara(15) == 1 ) then
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
              if( coloc_max <= rpara2 ) then
                 kelem = ielem
                 rpara2 = coloc_max
              end if
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
        rpara2 = -1.0_rp
     end if
  end if

  call elsest_cputim(time2)
  !!!!cputi(6,ithre) = cputi(6,ithre) + time2-time1

  call elsest_cputim(time3)
  !!!!cputi(7,ithre) = cputi(7,ithre) + time3-time2

end subroutine elsest_octpro
