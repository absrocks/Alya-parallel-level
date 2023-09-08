subroutine ibm_skdtib()!(iimbo,nboib)
  !-----------------------------------------------------------------------
  !****f* skdtib/ibm_skdtib
  ! NAME
  !    ibm_skdtib
  ! DESCRIPTION
  !    Skd-tree construction  
  ! USED BY
  !    nepoib
  !***
  !----------------------------------------------------------------------- 
!  use def_kintyp, only     :  ip,rp
!  use def_master, only     :  imbou
!  use def_domain, only     :  ndime
!  implicit none
!  integer(ip), intent(in)  :: iimbo,nboib
!  integer(ip), pointer     :: link(:)
!  integer(ip), pointer     :: vect2(:)
!  integer(ip), pointer     :: perm(:)
!  integer(ip), pointer     :: struc(:,:)
!  real(rp),    pointer     :: fabox(:,:,:),sabox(:,:,:),vect3(:)
!  real(rp),    pointer     :: centr(:,:)
!  integer(ip)              :: idime,iboib,indst,next,curr,imini,imaxi
!  integer(ip)              :: ladim,irang,krang,idim1,jrang
!  integer(ip)              :: jface,nfac2,imedi,elem2
!  integer(4)               :: istat
!  real(rp)                 :: elem1
  
!  allocate( struc(2*nboib-1,3), stat = istat )
!  allocate( perm(nboib),        stat = istat )
!  allocate( centr(ndime,nboib), stat = istat )

!  fabox => imbou(iimbo) % fabox
!  sabox => imbou(iimbo) % sabox
!  link  => imbou(iimbo) % link
  !
  ! Determine the centroid values for each boundig box of a face
  !
!  do iboib = 1,nboib
!     do idime = 1,ndime
!        centr(idime,iboib) = ( fabox(2,idime,iboib) + fabox(1,idime,iboib) ) * 0.5_rp
!     end do
!     perm(iboib) = iboib
!  end do
  !
  ! Store the total bounding box 
  !
!  do idime = 1,ndime
!     sabox(1,idime,1) = imbou(iimbo) % bobox(idime,1)
!     sabox(2,idime,1) = imbou(iimbo) % bobox(idime,2)
!  end do

!  indst      = 1_ip
!  struc(1,1) = 1_ip
!  struc(1,2) = 1_ip
!  struc(1,3) = nboib
!  next       = 2_ip
  !
  ! Tree strucure construction 
  !
!  do while( indst > 0_ip )
!     curr  = struc(indst,1)
!     imini = struc(indst,2)
!     imaxi = struc(indst,3)
!     indst = indst - 1_ip    

!     if( imaxi == imini ) then
!        link(curr) = -perm(imaxi)
!     else
!        allocate( vect3(imaxi-imini+1),stat=istat)        
!        imedi = imaxi+imini
!        imedi = int ( real(imaxi+imini) * 0.5_rp )
        !
        ! Choose the largest dimension of the current bounding box, sabox(curr)
        !
!        ladim = 1_ip
!        do idime = 2,ndime
!           if ( (sabox(2,idime,curr)-sabox(1,idime,curr)) > (sabox(2,idime-1,curr)-sabox(1,idime-1,curr)) ) then
!              ladim = idime
!           end if
!        end do
        !
        ! Reorder perm(imini:imaxi) with the centroid values
        !
!        vect2 => perm(imini:imaxi)
!        krang =  imaxi - imini + 1_ip
!        irang =  0
!        do idim1 = imini,imaxi
!           irang = irang + 1
!           vect3(irang) = centr(ladim,perm(idim1))
!        end do
      
!        do jrang = 2,krang
!           elem1 = vect3(jrang)
!           elem2 = vect2(jrang)
!           do irang = jrang-1,1,-1
!              if( vect3(irang) <= elem1 ) go to 30
!              vect3(irang+1) = vect3(irang)
!              vect2(irang+1) = vect2(irang)
!           end do
!           irang = 0
!30         vect3(irang+1) = elem1
!           vect2(irang+1) = elem2
!        end do
        !
        ! The two children of node curr are locate at link(curr) y link(curr+1)
        ! 
!        link(curr) = next
        !
        ! Determine the total bounding box for each children 
        !
!        do idime = 1,ndime
!           sabox(1,idime,next) =  1.0e6_rp
!           sabox(2,idime,next) = -1.0e6_rp 
!        end do
!        nfac2 = imedi - imini + 1_ip
!        do iboib = 1,nfac2
!           jface = perm(imini+iboib-1)
!           do idime = 1,ndime
!              sabox(1,idime,next) = min( fabox(1,idime,jface) , sabox(1,idime,next) ) 
!              sabox(2,idime,next) = max( fabox(2,idime,jface) , sabox(2,idime,next) ) 
!           end do
!        end do

!        do idime = 1,ndime
!           sabox(1,idime,next+1) =  1.0e6_rp
!           sabox(2,idime,next+1) = -1.0e6_rp
!        end do
!        nfac2 = imaxi-imedi
!        do iboib = 1,nfac2
!           jface = perm(imedi+iboib) 
!           do idime = 1,ndime
!              sabox(1,idime,next+1) = min( fabox(1,idime,jface) , sabox(1,idime,next+1) )               
!              sabox(2,idime,next+1) = max( fabox(2,idime,jface) , sabox(2,idime,next+1) ) 
!           end do
!        end do

        !call unibb(fabox(1:2,1:ndime,perm(imini:imedi)),sabox(1:2,1:ndime,next),imedi-imini+1_ip,ndime)
        !call unibb(fabox(1:2,1:ndime,perm(imedi+1_ip:imaxi)),sabox(1:2,1:ndime,next+1),imaxi-imedi,ndime)

        !
        ! Store the children of current element of the stack (struc)
        !
!        indst          = indst + 1_ip
!        struc(indst,1) = next
!        struc(indst,2) = imini
!        struc(indst,3) = imedi                
!        indst          = indst + 1_ip
!        struc(indst,1) = next  + 1_ip
!        struc(indst,2) = imedi + 1_ip
!        struc(indst,3) = imaxi             
!        next           = next  + 2_ip        
!        deallocate( vect3,stat=istat)
!     end if

!  end do

!  deallocate( centr, stat = istat )
!  deallocate( perm,  stat = istat )
!  deallocate( struc, stat = istat )

end subroutine ibm_skdtib
