subroutine elsest(&
     itask,imesh,ipara,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltype,ltopo,coord,xcoor,rpara,ifoun,&
     shapt,derit,coloc,lchec)
  !-----------------------------------------------------------------------
  !****f* Elsest
  ! NAME
  !    Elsest
  ! DESCRIPTION
  !    EL      SE     ST
  !      ement   arch   rategies   L I B R A R Y
  ! INPUT
  !    ITASK ....... =0: Allocate memory for the bin and quad/oct
  !                      structures. Should be called only once
  !                  =1: Element search
  !                  =2: Deallocate memory for mesh IMESH
  !                  =3: Deallocate memory for structures
  !    IMESH ....... Working mesh
  !    MNODE ....... Maximum number of nodes per elements
  !    NDIME ....... Dimension
  !    NPOIN ....... Number of nodes
  !    NELEM ....... Number of elements
  !    NNODE(*) .... Number of nodes
  !    IPARA(*) .... Integer input parameters
  !                  IPARA( 1) = # of bins in x                 (BIN)
  !                  IPARA( 2) = # of bins in y                 (BIN)
  !                  IPARA( 3) = # of bins in z                 (BIN)
  !                  IPARA( 4) = Data format (0=type,1=list)    (BIN)
  !                  IPARA( 5) = Max # of background meshes     (BOTH)
  !                  IPARA( 6) = Second try strategy            (BOTH)
  !                            = 1 for neighboring boxes        (BIN)
  !                            = 2 for max segment length       (BOTH)
  !                  IPARA( 7) = output unit (0 dfor no output) (BOTH)
  !                  IPARA( 8) = Search strategy (0=bin,1=oct)
  !                  IPARA( 9) = Max # of nodes per bin         (QUAD/OCT)
  !                  IPARA(10) = Search results freq (0 for no) (QUAD/OCT)
  !                  IPARA(11) = search radius of
  !                            = 0 for no                       (BIN)
  !                            = 1 for yes                      (BIN)
  !                  IPARA(12) = Unit for mesh postprocess      (BOTH)
  !                  IPARA(13) = Unit for result postprocess    (BOTH)
  !                  IPARA(14) = If element flag should be      (BOTH)
  !                              checked during search          
  !                  IPARA(15) = Find an element anyway         (BOTH)
  !    RPARA(*) .... Real input parameters
  !                  RPARA( 1) = Tolerance                      (BOTH)
  !                  RPARA( 2) = Min parametrized Distance
  !                              when using IPARA(15) = 1  
  !                              =-1 when found really inside
  !                              a host element
  ! OUTPUT
  !    IFOUN ....... =-1: Point out of bounding box
  !                  = 0: Element not found
  !                  > 0: Host element
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only     :  ip,rp,iunit,nthre,kfl_memor,lmini,lmaxi
  implicit none
  integer(ip), intent(in)  :: itask,imesh,ipara(*)
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode(*)
  integer(ip), intent(in)  :: lnods(*),ltype(*)
  integer(ip), intent(in)  :: ltopo(*),lchec(*)
  real(rp),    intent(in)  :: coord(*),xcoor(*)
  real(rp),    intent(in)  :: rpara(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: shapt(*),derit(*),coloc(*) 
  integer(ip)              :: ithre,recmeth,ielem,idime,itotn
  integer(ip)              :: ipoin,inode,itote,pnode
  real(rp)                 :: elcod(ndime,mnode)
  integer(ip), save        :: lmesh=0,ipass=0
#ifdef _OPENMPPPPPPPP
  integer                  :: OMP_GET_THREAD_NUM
  integer                  :: OMP_GET_NUM_THREADS
#endif
  !
  ! Errors
  !
  if( imesh > ipara(5) ) call runend('ELSEST: WRONG MESH NUMBER')
  !
  ! Initialization
  !
  !if( ipass == 0 ) then
  !   ipass = 1
  !   kfl_memor(1) = 0
  !   kfl_memor(2) = 0
  !   kfl_memor(3) = 0
  !end if

#ifdef _OPENMPPPPPPPPPPPPPPPPPPPPP
  ithre=OMP_GET_THREAD_NUM()+1
#else
  ithre=1
#endif

  select case(itask)

  case(-1_ip)
     !
     ! Automatically find a method
     !
     call elsest_recomm(&
          mnode,ndime,npoin,nelem,nnode,lnods,ltype,&
          coord,recmeth)

  case(0_ip)
     !
     ! Allocate memory for structures
     !
     call elsest_alloca(0_ip,ipara)

  case(1_ip)
     !
     ! Preprocess
     !
     if(ipara(8)==0) then
        call elsest_binpre(&
             ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&         ! Bin strategy
             lnods,ltype,coord,rpara)
     else if(ipara(8)==1) then
        call elsest_octpre(&                                           ! Quad/Oct strategy
             ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
             lnods,ltype,coord)
     end if

  case(2_ip)
     !
     ! Element search
     !
     if( ipara(8) == 0 ) then

        call elsest_binpro(&                                           ! Bin strategy
             imesh,lmesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
             lnods,ltype,ltopo,coord,xcoor,rpara,ifoun,shapt,&
             derit,coloc,lchec)

     else if( ipara(8) == 1 ) then

        call elsest_octpro(&                                           ! Quad/Oct strategy
             imesh,imesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
             lnods,ltype,ltopo,coord,xcoor,rpara,ifoun,shapt,&
             derit,coloc,lchec)
     end if

  case(3_ip)
     !
     ! Deallocate memory
     !
     call elsest_statis(2_ip,imesh,ipara,1)
     if(ipara(8)==0) then

        call elsest_bindea(ithre,imesh,lmesh)                          ! Bin strategy

     else if(ipara(8)==1) then

        call elsest_octdea(ithre,imesh,lmesh)                          ! Quad/Oct strategy
     end if

  case(4_ip)
     !
     ! Deallocate memory for structures
     !
     call elsest_deallo()

  case(5_ip)
     !
     ! Check if test point is in element IFOUN
     !
     call elsest_inelem(&
          mnode,ndime,nnode,lnods,ltype,ltopo,coord,xcoor,rpara,&
          ifoun,shapt,derit,coloc)

  case(6_ip)
     !
     ! Check if test point is in element IFOUN: stop if not
     !
     ielem = ifoun
     call elsest_inelem(&
          mnode,ndime,nnode,lnods,ltype,ltopo,coord,xcoor,rpara,&
          ifoun,shapt,derit,coloc)
     if( ifoun == 0 ) then
        write(*,*) ielem,coloc(1:ndime)
        call runend('ELSEST: TEST POINT NOT IN ELEMENT')
     end if

  case(7_ip)
     !
     ! Give shape function and deriv in element ielem
     !
     ielem = ifoun     
     itote = (ielem-1) * mnode 
     pnode = nnode(abs(ltype(ielem)))
     do inode = 1,pnode
        itote = itote + 1
        ipoin = lnods(itote)
        itotn = (ipoin-1) * ndime
        do idime = 1,ndime
           itotn = itotn + 1
           elcod(idime,inode) = coord(itotn)
        end do
     end do
     call elsest_chkelm(&
          ndime,ltopo(ielem),pnode,elcod,shapt,derit,&
          xcoor,coloc,ielem,lmini,lmaxi)
     if( ielem <= 0 ) call runend('ELSEST: COULD NOT INTEPROLATE')

  end select
  !
  ! Save last mesh number
  !
  !lmesh=imesh

end subroutine elsest
