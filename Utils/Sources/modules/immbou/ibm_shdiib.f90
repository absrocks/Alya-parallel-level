subroutine ibm_shdiib()! (iimbo,nboib,coord,dista,iboib,proje)
  !-----------------------------------------------------------------------
  ! NAME
  !    nepoib
  ! DESCRIPTION
  !    Shortest signed distance to a given point 
  !    INPUT
  !       iimbo: id of the particle
  !       nboib: number of faces of rthe particle
  !       point: point coordinates
  !    OUTPUT
  !       dista2: is less than 0 is point is inside the particle
  !              is more than 0 is point is outside the particle
  !       
  ! USED BY
  !    inouib
  !----------------------------------------------------------------------- 
  !use def_kintyp, only     :  ip,rp
  !use def_master, only     :  imbou,kfl_paral
  !use def_domain, only     :  ndime,mnoib,nnode
  !implicit none
  !integer(ip), intent(in)  :: iimbo,nboib
  !integer(ip), intent(out) :: iboib
  !real(rp),    intent(in)  :: coord(ndime)
  !real(rp),    intent(out) :: dista,proje(ndime)
  !integer(ip), pointer     :: blink(:),struc(:),canfa(:)
  !real(rp),    pointer     :: sabox(:,:,:),candi(:)
  !integer(ip)              :: idime,indst,ilist,inoib,nlist,curr
  !integer(ip)              :: pblty,pnodb,ipoib
  !integer(4)               :: istat
  !real(rp)                 :: dist,dismi,disma,distu
  !real(rp)                 :: norfa(3),propo(3)
  !real(rp)                 :: bocod(ndime,mnoib)

  !sabox => imbou(iimbo)%sabox
  !blink  => imbou(iimbo)%blink

  !allocate( struc(2*nboib-1),stat=istat)
  !allocate( canfa(nboib),    stat=istat)
  !allocate( candi(nboib),    stat=istat)

  !call ibm_maxdis(coord,sabox(1,1,1),ndime,distu)

  !indst        = 1_ip
  !struc(indst) = 1_ip
  !ilist        = 0_ip
  !
  ! Assemble a list of candidate patches by traversing skd-tree
  !
  !do while( indst > 0_ip )

  !   curr  = struc(indst)
  !   indst = indst - 1_ip
    
  !   call ibm_mindis(coord,sabox(1,1,curr),ndime,dismi)

  !   if( dismi < distu ) then
  !      call ibm_maxdis(coord,sabox(1,1,curr),ndime,disma)
  !      distu = min(distu,disma)
        !
        ! If currnode is a leaf in the tree structure
        !      
  !      if( blink(curr) < 0_ip ) then
  !         ilist = ilist + 1_ip
           !
           ! canli: candidates list
           !
  !         candi(ilist) =  dismi
  !         canfa(ilist) = -blink(curr)
  !      else
  !         indst        = indst + 1_ip
  !         struc(indst) = link(curr)
  !         indst        = indst + 1_ip
  !         struc(indst) = link(curr) + 1_ip
  !      end if
  !   end if
  !end do
  !
  ! Pare down candidate list
  !
  !dist  = 0.0_rp
  !nlist = ilist
  !do ilist = 1,nlist
  !   if( candi(ilist) < distu ) then
        !
        ! Element properties and dimensions
        !
  !      iboib = canfa(ilist)
  !      pblty = imbou(iimbo) % ltyib(iboib)
  !      pnodb = nnode(pblty)        
        !
        ! norfa: Exterior normal           
        !        
  !      call exteib(iimbo,pnodb,iboib,norfa,0_ip)
  !      do inoib = 1,pnodb
  !         ipoib = imbou(iimbo)%lnoib(inoib,iboib)            
  !         do idime = 1,ndime              
  !            bocod(idime,inoib) = imbou(iimbo)%cooib(idime,ipoib)
  !         end do
  !      end do                
        !
        ! Minimun distance from the point to the face
        !
  !      call ibm_pofadi(coord,bocod,norfa,ndime,dist,propo)           
  !      if( abs(dist) < distu ) then
  !         dista = dist
  !         do idime = 1,ndime
  !            proje(idime) = propo(idime)
  !         end do 
  !         iboib = canfa(ilist)
  !         distu = abs(dist)          
  !      end if

  !   end if
  !end do
 
  !deallocate( struc,stat=istat)
  !deallocate( canfa,stat=istat)
  !deallocate( candi,stat=istat)
end subroutine ibm_shdiib

!subroutine ibm_mindis(coord,bobox,ndime,dista)
  !-----------------------------------------------------------------------
  !****f* mindis/mindis
  ! NAME
  !    mindis
  ! DESCRIPTION
  !    Minimun distance between a point and a bounding box
  ! USED BY
  !    ibm_shdiib
  !***
  !----------------------------------------------------------------------- 
!  use def_kintyp, only      :  ip,rp
!  implicit   none
!  integer(ip),intent(in)    :: ndime
!  real(rp),   intent(in)    :: coord(ndime),bobox(2,ndime)
!  real(rp),   intent(out)   :: dista
!  integer(ip)               :: idime
!  real(rp)                  :: temp

!  dista = 0.0_rp

!  do idime = 1,ndime
!     temp  = max(0.0_rp,bobox(1,idime)-coord(idime))  
!     temp  = max(temp,coord(idime)-bobox(2,idime))  
!     dista = dista + temp * temp
!  end do
!  dista = sqrt(dista)

!end subroutine ibm_mindis

!subroutine ibm_maxdis(coord,bobox,ndime,dista)
  !-----------------------------------------------------------------------
  !****f* maxdis/maxdis
  ! NAME
  !    maxdis
  ! DESCRIPTION
  !    Maximum distance between a point and a bounding box
  ! USED BY
  !    ibm_shdiib
  !***
  !----------------------------------------------------------------------- 
!  use def_kintyp, only      :  ip,rp
!  implicit   none
!  integer(ip),intent(in)    :: ndime
!  real(rp),   intent(in)    :: coord(ndime),bobox(2,ndime)
!  real(rp),   intent(out)   :: dista
!  integer(ip)               :: idime
!  real(rp)                  :: temp

!  dista = 0.0_rp
!  do idime = 1,ndime
!     temp  = max(coord(idime)-bobox(1,idime), bobox(2,idime)-coord(idime))
!     dista = dista + temp*temp
!  end do
!  dista = sqrt(dista)

!end subroutine ibm_maxdis
