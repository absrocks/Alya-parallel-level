subroutine meshin(itask)
  !-----------------------------------------------------------------------
  !****f* meshin/meshin
  ! NAME
  !    meshin
  ! DESCRIPTION
  !    Bridge to mesh subroutines
  ! OUTPUT
  ! USED BY
  !    readim
  !    reageo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_parame
  use def_master
  use def_domain, only : nboib,ndime,npoib,mnoib,nimbo,&
       &                 lexis,nelem,npoin,ltype,memor_dom,&
       &                 nnode,ngaib
  use def_inpout
  use def_meshin
  use mod_cart
  use mod_meteo

  ! Pour le case 6 (temporaire)
  use mod_surf
  use mod_extrpr
  use mod_cart
  use mod_vol
  use mod_source
  use mod_memchk
  use mod_optim
  use mod_messages, only : livinf
  ! Pour le case 6 (temporaire)

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem

  ! Pour le case 6 (temporaire)
  integer(ip)             :: ndim,npoint,nx,ny,npsur,nblay,nnofa,nnodet,nfacet,nelemt,nnods
  integer(ip)             :: nline,nnosi,nside,nlinlet,npntr,nlpntr,nface
  integer(ip)             :: iinter,icartout,irefsurf,isizuni
  integer(ip)             :: ipoin,jpoin,jelem,iside,isizcrit,ismoo
  integer(ip)             :: inoib,iboib,ipoib,kpoin,kboun,iimbo,pnodb,pblty,iboun
  integer(ip)             :: pgaub,idime
  real(rp),    pointer    :: coordt(:,:)
  real(rp)                :: bbox(3,2),c00,tolscal,rtol
  integer(ip), pointer    :: lfacet(:,:),llinelet(:),lrenul(:,:),lsurft(:)
  integer(ip),pointer     :: lsmark(:)            
  integer(4)              :: istat
  real(rp)                :: thick,reason,dx,dy,rsize1
  ! Pour le case 6 (temporaire)

  c00=0.0d+00
  ndim=3_ip
  nnods=4_ip 
  nnofa=3_ip
  nnosi=2_ip 


  select case(itask)

  case(1_ip)
     !
     ! Define some mesh arrays and variables (from readim)
     !
     !kfl_camsh    = 1
     lexis(HEX08) = 1
#ifdef NDIMEPAR

#else
     ndime        = 3
#endif      

  case(2_ip)
     !
     ! Read mesh parameters (from reageo)
     !
     call reamsh()

  case(3_ip)

     !-------------------------------------------------------------------
     !
     ! Generate Cartesian mesh (from reageo)
     !
     !-------------------------------------------------------------------

     !
     ! Convert immbou format to carmsh format
     ! - NBOIB, NPOIB, LNOIB, COOIB
     !
     nboib = 0
     npoib = 0
     do iimbo = 1,nimbo 
        nboib = nboib + imbou(iimbo)%nboib
        npoib = npoib + imbou(iimbo)%npoib
     end do
     !allocate( cooib(ndime,npoib) )
     !allocate( lnoib(mnoib,nboib) )

     kpoin = 0
     kboun = 0
     do iimbo = 1,nimbo 
        do iboib = 1,imbou(iimbo)%nboib
           kboun = kboun + 1
           pblty = imbou(iimbo)%ltyib(iboib)
           pnodb = nnode(pblty)
           pgaub = ngaib(pblty)
           do inoib = 1,pnodb
             ! lnoib(inoib,kboun) = imbou(iimbo)%lnoib(inoib,iboib)
           end do
        end do
        do ipoib = 1,imbou(iimbo)%npoib
           kpoin = kpoin + 1
           do idime = 1,ndime
             ! cooib(idime,kpoin) = imbou(iimbo)%cooib(idime,ipoib)
           end do
        end do
     end do
    ! lface => lnoib
    ! coor  => cooib
     !
     ! Generate Cartesian mesh: LCELL, NCELL
     !
     call livinf(60_ip,' ',0_ip)
     isizcrit = 1 ! =1: do nothing. =2: Apply minimum between cell size and uniform size. 3= Apply source size
     iinter   = 0 ! =1: full intersection
     icartout = 0 ! Output cartesian mesh
     irefsurf = 1 ! =1: refine depending on the surface mesh 
     isizuni  = 0 ! =1: Apply uniform size and delete the previous size distribution
     call carmsh(npoib,ndime,nboib,3_ip,nsour,rsuni,rscal,kfl_ifbox,boxin,rtol,&
                 lcell,ncell,coor,lface,rsour,rsgeo,ismoo,isizcrit,iinter,&
                 icartout,irefsurf,isizuni)
     !
     ! Write Cartesian mesh in Alya format:
     ! - Geometrical description
     ! - Hanging nodes treatment
     ! - Boundary conditions
     !
     call cargeo(ncell,lcell)  ! Write
     !call carhan(ncell,lcell)
     call mescek(3_ip)

  case(4_ip)
     !
     ! Mark element with level
     !
     do ielem=1,nelem 
        gisca(ielem)=lcell(ielem)%level 
     end do

  case(-5_ip)
     !
     ! Mark element with mark
     !
     do ielem=1,nelem 
        gesca(ielem)=-real(lcell(ielem)%marked)
     end do

  case(5_ip)
     !
     ! Mark element with mark
     !
     do ielem=1,nelem 
        gisca(ielem)=-lcell(ielem)%marked
     end do

  case(6_ip)
     !
     ! Surface mesh
     !
#ifdef NDIMEPAR

#else
     ndime=3_ip
#endif
     nnofa=3_ip
     nnosi=2_ip
     !
     !     Read the bl distribution
     !
     !call readbl(nblay,rblay)
     !call insurf(nfacet,npoin,nnofa,ndime,nnosi,nside,nsurf,nline,tolscal) ! Read 
     !
     !     DBG
     !
     allocate(lsmark(nside),stat=istat)
     call memchk(zero,istat,memor_msh,'LSMARK','meshin',lsmark)
     do iside=1,nside
        lsmark(iside)=1
     enddo    

     !call readsour(ndime,nsour,rsuni,rscal,rsour,rsgeo,isizcrit,ismoo)
     !call carmsh(npoin,ndime,nfacet,nnofa,nsour,rsuni,rscal,0,bbox,rtol,&
     !               lcell,ncell,coor,lface,rsour,rsgeo,ismoo,isizcrit)
     !call mshsrf(ndime,nnofa,nfacet,npoin,nboup,nsurf,nline,nnosi,nside,&
     !            rsuni,nblay,rblay,tolscal,lsmark,rtol,lcell,ncell,lface,&
     !            rsize,coor,lsurf)
     !call mshvol(ndim,npoin,lface,nnofa,nfacet,rsize,nelem,nnods,coor,&
     !            elem,nnosi)
  case(7_ip)
     !
     ! Meteo
     !
     
     !call meteopre(npoin,nelem,ndim,nnofa,nnods,nface,lface,elem,coor,lsurft)
     !call edlinelet(npoin,ndim,coor,npntr,nlpntr,llinelet,nlinlet,lrenul,elem,&
     !     nnods,nelem,lface,nnofa,nface,lsurft,1,lbous,nsurf)
  
  case(8_ip)

     !call comp() 
  
  case(9_ip)

     !call scaledg() 

  case(10_ip)

     call optimsh() 

     nelem=0_ip
     !
     !     Call openfile Meteo
     !
     call openfi(9_ip)
    
  case(11_ip)
     
     nfacet=0_ip    
  
     bbox(1,1)=c00
     bbox(2,1)=c00
     bbox(3,1)=c00
     bbox(1,2)=3.0d+00
     bbox(2,2)=12.0d+00
     bbox(3,2)=12.0d+00

     !call readsour(ndime,nsour,rsuni,rscal,rsour,rsgeo,isizcrit,ismoo)
     !call carmsh(npoin,ndime,nfacet,nnofa,nsour,rsuni,rscal,1,bbox,rtol,&
     !               lcell,ncell,coor,lface,rsour,rsgeo,ismoo,isizcrit)
 
  end select

end subroutine meshin
