!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_kdtree.f90
!> @author  houzeaux
!> @date    2020-09-09
!> @brief   kd-tree
!> @details KD-tree toolbox
!-----------------------------------------------------------------------

module mod_kdtree
  
  use def_kintyp_basic,   only : ip,rp,lg
  use def_kintyp_domain,  only : netyp
  use def_elmtyp,         only : BAR02,TRI03
  use def_domain,         only : ndime,nnode
  use def_domain,         only : memor_dom
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  
  implicit none

  real(rp),    parameter        :: epsil=epsilon(1.0_rp)
  
  real(rp),    pointer, private :: fabox_loc(:,:,:) => null()
  real(rp),    pointer, private :: bobox_loc(:,:)   => null()
  real(rp),    pointer, private :: sabox_loc(:,:,:) => null()
  integer(ip), pointer, private :: blink_loc(:)     => null()
  integer(ip), pointer, private :: stru2_loc(:)     => null()
  real(rp),    pointer, private :: ldist_loc(:)     => null()
  type(netyp), pointer, private :: lnele_loc(:)     => null()
  
  type typ_kdtree
     integer(ip)          :: mnodb
     integer(ip)          :: npoin
     integer(ip)          :: nboun
     real(rp)             :: bobox(3,2)
     real(rp),    pointer :: coord(:,:)
     integer(ip), pointer :: lnodb(:,:)
     integer(ip), pointer :: ltypb(:)
     integer(ip), pointer :: lperm(:)
     real(rp),    pointer :: fabox(:,:,:)
     real(rp),    pointer :: sabox(:,:,:)
     integer(ip), pointer :: blink(:)
     integer(ip), pointer :: stru2(:)
     real(rp),    pointer :: ldist(:)
     type(netyp), pointer :: lnele(:)
  end type typ_kdtree
  
contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/04/2014
  !> @brief   Delete kdtree
  !> @details Delete the kdtree structure
  !>
  !----------------------------------------------------------------------

  subroutine kdtree_deallocate(kdtree_typ)
    
    type(typ_kdtree), intent(inout) :: kdtree_typ
    integer(ip)                     :: ipoib
    
    kdtree_typ % nboun = 0
    kdtree_typ % npoin = 0
    kdtree_typ % mnodb = 0
    kdtree_typ % bobox = 0.0_rp

    call memory_deallo(memor_dom,'COORD','mod_kdtree',kdtree_typ % coord)
    call memory_deallo(memor_dom,'LNODB','mod_kdtree',kdtree_typ % lnodb)
    call memory_deallo(memor_dom,'LTYPB','mod_kdtree',kdtree_typ % ltypb)
    call memory_deallo(memor_dom,'LPERM','mod_kdtree',kdtree_typ % lperm)
    call memory_deallo(memor_dom,'FABOX','mod_kdtree',kdtree_typ % fabox)
    call memory_deallo(memor_dom,'SABOX','mod_kdtree',kdtree_typ % sabox)
    call memory_deallo(memor_dom,'BLINK','mod_kdtree',kdtree_typ % blink)
    call memory_deallo(memor_dom,'STRU2','mod_kdtree',kdtree_typ % stru2)
    call memory_deallo(memor_dom,'LDIST','mod_kdtree',kdtree_typ % ldist)

    if( associated(kdtree_typ % lnele) ) then       
       do ipoib = 1,size(kdtree_typ % lnele,KIND=ip)
          call memory_deallo(memor_dom,'LNELE % ELEME','mod_kdtree',kdtree_typ % lnele(ipoib) % eleme)
          call memory_deallo(memor_dom,'LNELE % LTYPE','mod_kdtree',kdtree_typ % lnele(ipoib) % ltype)
       end do
       deallocate(kdtree_typ % lnele)    
    end if
    
  end subroutine kdtree_deallocate

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/04/2014
  !> @brief   Initialize kdtree
  !> @details Construct the kdtree structure
  !>
  !----------------------------------------------------------------------

  subroutine kdtree_initialize(kdtree_typ)
    type(typ_kdtree), intent(inout) :: kdtree_typ

    kdtree_typ % nboun = 0
    kdtree_typ % npoin = 0
    kdtree_typ % mnodb = 0
    kdtree_typ % bobox = 0.0_rp
    nullify( kdtree_typ % coord )
    nullify( kdtree_typ % lnodb )
    nullify( kdtree_typ % ltypb )
    nullify( kdtree_typ % lperm )
    nullify( kdtree_typ % fabox )
    nullify( kdtree_typ % sabox )
    nullify( kdtree_typ % blink ) 
    nullify( kdtree_typ % stru2 )
    nullify( kdtree_typ % ldist )
    nullify( kdtree_typ % lnele )

  end subroutine kdtree_initialize

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/04/2014
  !> @brief   Initialize kdtree
  !> @details Construct the kdtree structure
  !>
  !----------------------------------------------------------------------

  subroutine kdtree_construct(nboun,npoin,lnodb,ltypb,coord,kdtree_typ,lperm)
    
    integer(ip),                        intent(in)    :: nboun
    integer(ip),                        intent(in)    :: npoin
    integer(ip),      pointer,          intent(in)    :: lnodb(:,:)
    integer(ip),      pointer,          intent(in)    :: ltypb(:)
    real(rp),         pointer,          intent(in)    :: coord(:,:)
    type(typ_kdtree),                   intent(inout) :: kdtree_typ
    integer(ip),      pointer, optional,intent(in)    :: lperm(:)
    integer(ip)                                       :: iboun,ipoin,inodb,ipoin_bou
    integer(ip)                                       :: pblty,mnodb,pnodb,kboun
    integer(ip),      allocatable                     :: gisca(:)
    !
    ! Initialize
    !
    call kdtree_initialize(kdtree_typ)
    !
    ! Dimensions
    !
    mnodb = size(lnodb,1_ip,KIND=ip)

    if( present(lperm) ) then
       if( associated(lperm) ) then          
          kdtree_typ % nboun = size(lperm,KIND=ip)
       else
          kdtree_typ % nboun = 0
       end if
    else
       kdtree_typ % nboun = nboun
    end if
    kdtree_typ % mnodb = mnodb
    !
    ! Permutation
    !
    if( present(lperm) .and. kdtree_typ % nboun > 0 ) then
       call memory_alloca(memor_dom,'LPERM','mod_kdtree',kdtree_typ % lperm,kdtree_typ % nboun)
       do iboun = 1,kdtree_typ % nboun
          kdtree_typ % lperm(iboun) = lperm(iboun)
       end do       
    end if

    if( kdtree_typ % nboun > 0 ) then
       !
       ! Mark nodes
       !
       allocate(gisca(npoin))
       do ipoin = 1,npoin
          gisca(ipoin) = 0
       end do

       if( present(lperm) ) then          
          do kboun = 1,kdtree_typ % nboun
             iboun = lperm(kboun)
             pblty = ltypb(iboun)
             if( pblty > 0 ) then
                pnodb = nnode(pblty)
                do inodb = 1,pnodb
                   ipoin = lnodb(inodb,iboun)
                   if( gisca(ipoin) == 0 ) then
                      kdtree_typ % npoin = kdtree_typ % npoin + 1
                      gisca(ipoin) = kdtree_typ % npoin
                   end if
                end do
             end if
          end do
       else
          do iboun = 1,kdtree_typ % nboun
             pblty = ltypb(iboun)
             if( pblty > 0 ) then
                pnodb = nnode(pblty)
                do inodb = 1,pnodb
                   ipoin = lnodb(inodb,iboun)
                   if( gisca(ipoin) == 0 ) then
                      kdtree_typ % npoin = kdtree_typ % npoin + 1
                      gisca(ipoin) = kdtree_typ % npoin
                   end if
                end do
             end if
          end do
       end if
       !
       ! Allocate memory
       !
       call memory_alloca(memor_dom,'COORD','mod_kdtree',kdtree_typ % coord,ndime,kdtree_typ % npoin)
       call memory_alloca(memor_dom,'LNODB','mod_kdtree',kdtree_typ % lnodb,mnodb,kdtree_typ % nboun)
       call memory_alloca(memor_dom,'LTYPB','mod_kdtree',kdtree_typ % ltypb,      kdtree_typ % nboun)
       !
       ! Bluid kdtree geometry using renumbering of nodes and boundaries
       !
       do ipoin = 1,npoin        
          ipoin_bou = gisca(ipoin)
          if( ipoin_bou > 0 ) then
             kdtree_typ % coord(1:ndime,ipoin_bou) = coord(1:ndime,ipoin)
          end if
       end do

       if( present(lperm) ) then          
          do kboun = 1,kdtree_typ % nboun
             iboun = lperm(kboun)
             pblty = abs(ltypb(iboun))
             kdtree_typ % ltypb(kboun) = pblty
             if( pblty > 0 ) then
                pnodb = nnode(pblty)
                do inodb = 1,pnodb
                   ipoin     = lnodb(inodb,iboun)
                   ipoin_bou = gisca(ipoin)
                   kdtree_typ % lnodb(inodb,kboun) = ipoin_bou
                   !if (ipoin_bou < 0) 
                end do
             end if
          end do
       else
          do iboun = 1,kdtree_typ % nboun
             pblty = abs(ltypb(iboun))
             kdtree_typ % ltypb(iboun) = pblty
             if( pblty > 0 ) then
                pnodb = nnode(pblty)
                do inodb = 1,pnodb
                   ipoin     = lnodb(inodb,iboun)
                   ipoin_bou = gisca(ipoin)
                   kdtree_typ % lnodb(inodb,iboun) = ipoin_bou
                   !if (ipoin_bou < 0) 
                end do
             end if
          end do
       end if

       deallocate( gisca )

       call kdtree(&
            1_ip               , kdtree_typ % mnodb ,                      &
            kdtree_typ % npoin , kdtree_typ % nboun , kdtree_typ % coord , &
            kdtree_typ % lnodb , kdtree_typ % ltypb , kdtree_typ % fabox , &
            kdtree_typ % bobox , kdtree_typ % sabox , kdtree_typ % blink , &
            kdtree_typ % stru2 , kdtree_typ % ldist , kdtree_typ % lnele   )

    end if

  end subroutine kdtree_construct

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/04/2014
  !> @brief   Find nearest boundary
  !> @details Find nearest boundary in the original numbering. If 
  !>          the kdtree was constructed using a permutation array
  !>
  !----------------------------------------------------------------------

  subroutine kdtree_nearest_boundary(coord_test,kdtree_typ,iboun,dista,proje,norma,lperm,lmask,chkdi_)
    implicit none
    real(rp),                        intent(in)    :: coord_test(*)
    type(typ_kdtree),                intent(inout) :: kdtree_typ
    integer(ip),                     intent(out)   :: iboun
    real(rp),                        intent(out)   :: dista
    real(rp),                        intent(out)   :: proje(*)
    real(rp),              optional, intent(out)   :: norma(*)
    integer(ip), pointer,  optional, intent(in)    :: lperm(:)
    logical(lg), pointer,  optional, intent(in)    :: lmask(:)
    real(rp),              optional, intent(in)    :: chkdi_
    real(rp)                                       :: dummy(3),chkdi
    
    if( kdtree_typ % nboun > 0 ) then
       if(present(chkdi_)) then
          chkdi = chkdi_
       else
          chkdi = 1.0e9_rp       
       end if
       call dpopar(&
            1_ip               , coord_test         , &
            kdtree_typ % npoin , kdtree_typ % mnodb , &
            kdtree_typ % nboun , chkdi              , &
            kdtree_typ % ltypb , kdtree_typ % lnodb , &
            kdtree_typ % coord , dista              , &
            dummy              , proje              , &
            iboun              , kdtree_typ % fabox , &
            kdtree_typ % sabox , kdtree_typ % blink , &
            kdtree_typ % stru2 , kdtree_typ % ldist , &
            kdtree_typ % lnele , LMASK=lmask          )   
       if( present(lperm) ) then
          iboun = lperm(iboun)
       else if( associated(kdtree_typ % lperm) ) then
          if( size(kdtree_typ % lperm,KIND=ip) > 0 ) iboun = kdtree_typ % lperm(iboun)
       end if
       if( present(norma) ) norma(1:ndime) = 0.0_rp
    else
       dista          = huge(1.0_rp)
       proje(1:ndime) = huge(1.0_rp)
       if( present(norma) ) norma(1:ndime) = dummy(1:ndime)
       iboun          = 0
    end if

  end subroutine kdtree_nearest_boundary

  recursive function find_height(nod) result(lev) 
  
     implicit none
     integer(ip),      intent(in)   :: nod
     integer(ip)                    :: lev
     integer(ip)                    :: left, right

     if(blink_loc(nod) < 0) then
        lev=0
     else
        left  = find_height(blink_loc(nod))
        right = find_height(blink_loc(nod)+1)
        lev=max(left,right) + 1
     endif

  end function

  subroutine kdtree(                              &
       itask,mnodb,npoib,nboun,coord,lnodb,ltypb, &
       fabox,bobox,sabox,blink,stru2,ldist,lnele, &
       ltest)
    !-----------------------------------------------------------------------
    !****f* domain/kdtree
    ! NAME
    !    kdtree    ! DESCRIPTION
    !    Skd-tree construction. If ITASK=1 memory is allocated for the 
    !    structure at the beginning and for the search at the end
    ! INPUT
    !    ITASK ... If=1 FABOX is computed here
    !    NBOUN ... Number of boundaries
    !    FABOX ... Bounding box of each boundary
    !    BOBOX ... Complete boundaing box
    !    IF ITASK = 1 
    !      COORD ... 
    !      LNODB ... 
    !      LTYPB ... 
    ! OUTPUT
    !    SABOX ... Bounding box kd-tree structure
    !    BLINK ... Index of tree
    ! USED BY
    !    nepoib
    !***
    !----------------------------------------------------------------------- 
    implicit none
    integer(ip), intent(in)                       :: itask
    integer(ip), intent(in)                       :: mnodb
    integer(ip), intent(in)                       :: npoib
    integer(ip), intent(in)                       :: nboun
    real(rp),    intent(in),    pointer           :: coord(:,:)
    integer(ip), intent(in),    pointer           :: lnodb(:,:)
    integer(ip), intent(in),    pointer           :: ltypb(:)

    real(rp),    intent(inout), pointer, optional :: fabox(:,:,:)
    real(rp),    intent(inout), target , optional :: bobox(3,2)
    real(rp),    intent(inout), pointer, optional :: sabox(:,:,:)
    integer(ip), intent(inout), pointer, optional :: blink(:)
    integer(ip),                pointer, optional :: stru2(:)
    real(rp),    intent(inout), pointer, optional :: ldist(:)
    type(netyp), intent(inout), pointer, optional :: lnele(:)
    integer(ip), intent(in),    pointer, optional :: ltest(:)

    integer(ip), pointer                          :: perm(:)    
    integer(ip), pointer                          :: struc(:,:) 
    real(rp),    pointer                          :: vect3(:)   
    real(rp),    pointer                          :: centr(:,:) 
    integer(ip), pointer                          :: blink_lev(:)
    integer(ip)                                   :: idime,ipoib,iboun,inodb
    integer(ip)                                   :: indst,next,curr,imini,imaxi
    integer(ip)                                   :: ladim,irang,krang,idim1
    integer(ip)                                   :: jface,nfac2,imedi
    real(rp)                                      :: toler,dista,time1,time2

    nullify(perm)
    nullify(struc)
    nullify(vect3)
    nullify(centr)
    nullify(blink_lev)
    
    if( itask == 2 ) then
       !
       ! Deallocate kd-tree structure memory
       !
       if( present(fabox) ) then
          if( associated(lnele) ) then
             do ipoib = 1,size(lnele,KIND=ip)
                if( associated(lnele(ipoib) % eleme) ) deallocate( lnele(ipoib) % eleme )
                if( associated(lnele(ipoib) % ltype) ) deallocate( lnele(ipoib) % ltype )
             end do
             deallocate( lnele )
          end if
          if( associated(ldist) ) deallocate( ldist )
          if( associated(stru2) ) deallocate( stru2 )
          if( associated(blink) ) deallocate( blink )
          if( associated(fabox) ) deallocate( fabox )
          if( associated(sabox) ) deallocate( sabox )
       else
          if( associated(lnele_loc) ) then
             do ipoib = 1,size(lnele_loc,KIND=ip)
                if( associated(lnele_loc(ipoib) % eleme) ) deallocate( lnele_loc(ipoib) % eleme )
                if( associated(lnele_loc(ipoib) % ltype) ) deallocate( lnele_loc(ipoib) % ltype )
             end do
             deallocate( lnele_loc )
          end if
          if( associated(ldist_loc) ) deallocate( ldist_loc )
          if( associated(stru2_loc) ) deallocate( stru2_loc )
          if( associated(blink_loc) ) deallocate( blink_loc )
          if( associated(fabox_loc) ) deallocate( fabox_loc )
          if( associated(sabox_loc) ) deallocate( sabox_loc )
          if( associated(bobox_loc) ) deallocate( bobox_loc )
       end if

    else if( itask == 1 ) then
       !
       ! Check
       !
       if( ndime == 2 ) then
          do iboun = 1,nboun
             if( ltypb(iboun) /= BAR02 ) then
                call runend('MOD_KDTREE: WRONG BOUNDARY TYPE')
             end if
          end do
       else          
          do iboun = 1,nboun
             !
             ! OJO DESDOCUEMNTAR ESTO:
             !
             !if( ltypb(iboun) /= TRI03 ) then
             !   call runend('MOD_KDTREE: WRONG BOUNDARY TYPE')
             !end if
          end do
       end if
       !
       ! Construct FABOX_LOC and BOBOX_LOC
       !       
       allocate( blink_lev(max(1_ip,2_ip*nboun-1_ip)) )
       if( present(fabox) ) then
          call memory_alloca(memor_dom,'FABOX','mod_kdtree',fabox,2_ip,ndime,max(1_ip,nboun)           )
          call memory_alloca(memor_dom,'SABOX','mod_kdtree',sabox,2_ip,ndime,max(1_ip,2_ip*nboun-1_ip) )
          call memory_alloca(memor_dom,'BLINK','mod_kdtree',blink,max(1_ip,2_ip*nboun-1_ip)            )
          call memory_alloca(memor_dom,'STRU2','mod_kdtree',stru2,2_ip*nboun-1_ip                      )
          call memory_alloca(memor_dom,'LDIST','mod_kdtree',ldist,2_ip*nboun-1_ip                      )
          allocate( lnele(npoib) )
          fabox_loc => fabox 
          bobox_loc => bobox 
          sabox_loc => sabox 
          blink_loc => blink 
          stru2_loc => stru2 
          ldist_loc => ldist 
          lnele_loc => lnele        
       else
          allocate( fabox_loc(2,ndime,max(1_ip,nboun)) )
          allocate( sabox_loc(2,ndime,max(1_ip,2_ip*nboun-1_ip)) )
          allocate( blink_loc(max(1_ip,2_ip*nboun-1_ip)) )
          allocate( lnele_loc(npoib) )
          allocate( stru2_loc(2*nboun-1))
          allocate( ldist_loc(2*nboun-1))
          allocate( bobox_loc(3,2))
       end if
       !
       ! Nullify(lnele)
       !
       do ipoib = 1,npoib
          lnele_loc(ipoib) % nelem = 0_ip
          nullify(lnele_loc(ipoib) % eleme)
          nullify(lnele_loc(ipoib) % ltype)
       end do
       !
       ! Construct LNELE_LOC. List of connected elements for each node 
       !       
       if( present(ltest) ) then

          do iboun = 1,nboun
             if( ltest(iboun) /= 0 ) then 
               do inodb = 1,nnode(ltypb(iboun))
                   ipoib = lnodb(inodb,iboun)            
                   lnele_loc(ipoib)%nelem = lnele_loc(ipoib)%nelem + 1_ip
                end do
             end if
          end do
          do ipoib = 1,npoib
             call memory_alloca(memor_dom,'LNELE % ELEME','mod_kdtree',lnele_loc(ipoib) % eleme,lnele_loc(ipoib) % nelem)
             call memory_alloca(memor_dom,'LNELE % LTYPE','mod_kdtree',lnele_loc(ipoib) % ltype,lnele_loc(ipoib) % nelem)
          end do
          do ipoib = 1,npoib
             lnele_loc(ipoib)%nelem = 0_ip
          end do
          do iboun = 1,nboun
             if( ltest(iboun) /= 0 ) then
                do inodb = 1,nnode(ltypb(iboun))
                   ipoib = lnodb(inodb,iboun)       
                   lnele_loc(ipoib)%nelem = lnele_loc(ipoib)%nelem + 1_ip     
                   lnele_loc(ipoib)%eleme(lnele_loc(ipoib)%nelem) = iboun
                   lnele_loc(ipoib)%ltype(lnele_loc(ipoib)%nelem) = ltypb(iboun)
                end do
             end if
          end do

       else

          do iboun = 1,nboun
             do inodb = 1,nnode(ltypb(iboun))
                ipoib = lnodb(inodb,iboun)   
                lnele_loc(ipoib)%nelem = lnele_loc(ipoib)%nelem + 1_ip
             end do
          end do
          do ipoib = 1,npoib
             call memory_alloca(memor_dom,'LNELE % ELEME','mod_kdtree',lnele_loc(ipoib) % eleme,lnele_loc(ipoib) % nelem)
             call memory_alloca(memor_dom,'LNELE % LTYPE','mod_kdtree',lnele_loc(ipoib) % ltype,lnele_loc(ipoib) % nelem)
          end do
          do ipoib = 1,npoib
             lnele_loc(ipoib)%nelem = 0_ip
          end do
          do iboun = 1,nboun
             do inodb = 1,nnode(ltypb(iboun))
                ipoib = lnodb(inodb,iboun)       
                lnele_loc(ipoib)%nelem = lnele_loc(ipoib)%nelem + 1_ip     
                lnele_loc(ipoib)%eleme(lnele_loc(ipoib)%nelem) = iboun
                lnele_loc(ipoib)%ltype(lnele_loc(ipoib)%nelem) = ltypb(iboun)
             end do
          end do

       end if

    end if

    if( itask /= 2 ) then

       if( present(ltest) ) then

          toler = 1.0e-6_rp
          do iboun = 1,nboun
             if( ltest(iboun) /= 0 ) then
                do idime = 1,ndime
                   fabox_loc(1,idime,iboun) =  1.0e12_rp
                   fabox_loc(2,idime,iboun) = -1.0e12_rp
                   do inodb = 1,nnode(ltypb(iboun))
                      ipoib = lnodb(inodb,iboun)
                      fabox_loc(1,idime,iboun) = min( fabox_loc(1,idime,iboun) , coord(idime,ipoib) )
                      fabox_loc(2,idime,iboun) = max( fabox_loc(2,idime,iboun) , coord(idime,ipoib) )
                   end do
                   dista = fabox_loc(2,idime,iboun) - fabox_loc(1,idime,iboun)
                   if (dista < 2*toler) then
                      fabox_loc(1,idime,iboun) = fabox_loc(1,idime,iboun) - toler
                      fabox_loc(2,idime,iboun) = fabox_loc(2,idime,iboun) + toler
                   else
                      fabox_loc(1,idime,iboun) = fabox_loc(1,idime,iboun) - abs(dista*toler)
                      fabox_loc(2,idime,iboun) = fabox_loc(2,idime,iboun) + abs(dista*toler)
                   end if
                end do
             end if
          end do
          do idime = 1,ndime
             bobox_loc(idime,1) =  1.0e12_rp
             bobox_loc(idime,2) = -1.0e12_rp
          end do
          do iboun = 1,nboun
             if( ltest(iboun) /= 0 ) then
                do idime = 1,ndime
                   bobox_loc(idime,1) = min( bobox_loc(idime,1) , fabox_loc(1,idime,iboun) )
                   bobox_loc(idime,2) = max( bobox_loc(idime,2) , fabox_loc(2,idime,iboun) )
                end do
             end if
          end do

       else

          toler = 1.0e-6_rp
          do iboun = 1,nboun
             do idime = 1,ndime
                fabox_loc(1,idime,iboun) =  1.0e12_rp
                fabox_loc(2,idime,iboun) = -1.0e12_rp
                do inodb = 1,nnode(ltypb(iboun))
                   ipoib = lnodb(inodb,iboun)
                   fabox_loc(1,idime,iboun) = min( fabox_loc(1,idime,iboun) , coord(idime,ipoib) )
                   fabox_loc(2,idime,iboun) = max( fabox_loc(2,idime,iboun) , coord(idime,ipoib) )
                end do
                dista = fabox_loc(2,idime,iboun) - fabox_loc(1,idime,iboun)
                if( dista < 2.0_rp*toler ) then
                   fabox_loc(1,idime,iboun) = fabox_loc(1,idime,iboun) - toler
                   fabox_loc(2,idime,iboun) = fabox_loc(2,idime,iboun) + toler
                else
                   fabox_loc(1,idime,iboun) = fabox_loc(1,idime,iboun) - abs(dista*toler)
                   fabox_loc(2,idime,iboun) = fabox_loc(2,idime,iboun) + abs(dista*toler)
                end if
             end do
          end do
          do idime = 1,ndime
             bobox_loc(idime,1) =  1.0e12_rp
             bobox_loc(idime,2) = -1.0e12_rp
          end do
          do iboun = 1,nboun
             do idime = 1,ndime
                bobox_loc(idime,1) = min( bobox_loc(idime,1) , fabox_loc(1,idime,iboun) )
                bobox_loc(idime,2) = max( bobox_loc(idime,2) , fabox_loc(2,idime,iboun) )
             end do
          end do

       end if
       !
       ! Allocate memory
       !
       allocate( perm(nboun)        )
       allocate( centr(ndime,nboun) )
       allocate( struc(2*nboun-1,3) )
       !
       ! Determine the centroid values for each boundig box of a face
       !
       if( present(ltest) ) then
          do iboun = 1,nboun
             if( ltest(iboun) /= 0 ) then
                do idime = 1,ndime
                   centr(idime,iboun) = ( fabox_loc(2,idime,iboun) + fabox_loc(1,idime,iboun) ) * 0.5_rp
                end do
                perm(iboun) = iboun
             end if
          end do
       else
          do iboun = 1,nboun
             do idime = 1,ndime
                centr(idime,iboun) = ( fabox_loc(2,idime,iboun) + fabox_loc(1,idime,iboun) ) * 0.5_rp
             end do
             perm(iboun) = iboun
          end do
       end if
       !
       ! Store the total bounding box 
       !
       do idime = 1,ndime
          sabox_loc(1,idime,1) = bobox_loc(idime,1)
          sabox_loc(2,idime,1) = bobox_loc(idime,2)
       end do

       indst      = 1_ip
       struc(1,1) = 1_ip
       struc(1,2) = 1_ip
       struc(1,3) = nboun
       next       = 2_ip
       blink_loc    = 0
       blink_lev(1) = 1
       call cputim(time1)
       !
       ! Tree strucure construction 
       !
       do while( indst > 0_ip )
          
          curr  = struc(indst,1)
          imini = struc(indst,2)
          imaxi = struc(indst,3)
          indst = indst - 1_ip    

          if( imaxi == imini ) then
             blink_loc(curr) = -perm(imaxi)
          else
             allocate( vect3(imaxi-imini+1) )
             imedi = int ( real(imaxi+imini,rp) * 0.5_rp , KIND=ip)
             !
             ! Choose the largest dimension of the current bounding box, sabox_loc(curr)
             !
             ladim = 1_ip
             if ( (sabox_loc(2,2,curr)-sabox_loc(1,2,curr)) > (sabox_loc(2,1,curr)-sabox_loc(1,1,curr)) ) ladim = 2
             if ( (sabox_loc(2,ndime,curr)-sabox_loc(1,ndime,curr)) > (sabox_loc(2,2,curr)-sabox_loc(1,2,curr)) ) ladim = ndime
             !
             ! Reorder perm(imini:imaxi) with the centroid values
             !
             krang =  imaxi - imini + 1_ip
             irang =  0
             do idim1 = imini,imaxi
                irang = irang + 1
                vect3(irang) = centr(ladim,perm(idim1))
             end do
             !
             ! Sort VECT2 using VECT3 in increasing order
             !
             call heapsortri(2_ip,krang,vect3,perm(imini))
             !
             ! The two children of node curr are locate at blink_loc(curr) y blink_loc(curr+1)
             ! 
             blink_loc(curr) = next
             !
             ! Determine the total bounding box for each children 
             !
             sabox_loc(1,1,next)     =  1.0e6_rp
             sabox_loc(2,1,next)     = -1.0e6_rp
             sabox_loc(1,2,next)     =  1.0e6_rp
             sabox_loc(2,2,next)     = -1.0e6_rp
             sabox_loc(1,ndime,next) =  1.0e6_rp
             sabox_loc(2,ndime,next) = -1.0e6_rp 
             nfac2 = imedi - imini + 1_ip

             if( present(ltest) ) then
                do iboun = 1,nfac2
                   if( ltest(iboun) /= 0 ) then
                      jface = perm(imini+iboun-1)
                      sabox_loc(1,1,next)     = min( fabox_loc(1,1,jface) , sabox_loc(1,1,next) ) 
                      sabox_loc(2,1,next)     = max( fabox_loc(2,1,jface) , sabox_loc(2,1,next) ) 
                      sabox_loc(1,2,next)     = min( fabox_loc(1,2,jface) , sabox_loc(1,2,next) ) 
                      sabox_loc(2,2,next)     = max( fabox_loc(2,2,jface) , sabox_loc(2,2,next) ) 
                      sabox_loc(1,ndime,next) = min( fabox_loc(1,ndime,jface) , sabox_loc(1,ndime,next) ) 
                      sabox_loc(2,ndime,next) = max( fabox_loc(2,ndime,jface) , sabox_loc(2,ndime,next) ) 
                   end if
                end do
             else
                do iboun = 1,nfac2
                   jface = perm(imini+iboun-1)
                   sabox_loc(1,1,next)     = min( fabox_loc(1,1,jface) , sabox_loc(1,1,next) ) 
                   sabox_loc(2,1,next)     = max( fabox_loc(2,1,jface) , sabox_loc(2,1,next) ) 
                   sabox_loc(1,2,next)     = min( fabox_loc(1,2,jface) , sabox_loc(1,2,next) ) 
                   sabox_loc(2,2,next)     = max( fabox_loc(2,2,jface) , sabox_loc(2,2,next) ) 
                   sabox_loc(1,ndime,next) = min( fabox_loc(1,ndime,jface) , sabox_loc(1,ndime,next) ) 
                   sabox_loc(2,ndime,next) = max( fabox_loc(2,ndime,jface) , sabox_loc(2,ndime,next) ) 
                end do
             end if

             sabox_loc(1,1,next+1)     =  huge(1.0_rp)*0.1_rp
             sabox_loc(2,1,next+1)     = -huge(1.0_rp)*0.1_rp
             sabox_loc(1,2,next+1)     =  huge(1.0_rp)*0.1_rp
             sabox_loc(2,2,next+1)     = -huge(1.0_rp)*0.1_rp
             sabox_loc(1,ndime,next+1) =  huge(1.0_rp)*0.1_rp
             sabox_loc(2,ndime,next+1) = -huge(1.0_rp)*0.1_rp

             nfac2 = imaxi-imedi
             if( present(ltest) ) then
                do iboun = 1,nfac2
                   if( ltest(iboun) /= 0 ) then
                      jface = perm(imedi+iboun) 
                      sabox_loc(1,1,next+1)     = min( fabox_loc(1,1,jface) , sabox_loc(1,1,next+1) )               
                      sabox_loc(2,1,next+1)     = max( fabox_loc(2,1,jface) , sabox_loc(2,1,next+1) ) 
                      sabox_loc(1,2,next+1)     = min( fabox_loc(1,2,jface) , sabox_loc(1,2,next+1) )               
                      sabox_loc(2,2,next+1)     = max( fabox_loc(2,2,jface) , sabox_loc(2,2,next+1) ) 
                      sabox_loc(1,ndime,next+1) = min( fabox_loc(1,ndime,jface) , sabox_loc(1,ndime,next+1) )               
                      sabox_loc(2,ndime,next+1) = max( fabox_loc(2,ndime,jface) , sabox_loc(2,ndime,next+1) ) 
                   end if
                end do
             else
                do iboun = 1,nfac2
                   jface = perm(imedi+iboun) 
                   sabox_loc(1,1,next+1)     = min( fabox_loc(1,1,jface) , sabox_loc(1,1,next+1) )               
                   sabox_loc(2,1,next+1)     = max( fabox_loc(2,1,jface) , sabox_loc(2,1,next+1) ) 
                   sabox_loc(1,2,next+1)     = min( fabox_loc(1,2,jface) , sabox_loc(1,2,next+1) )               
                   sabox_loc(2,2,next+1)     = max( fabox_loc(2,2,jface) , sabox_loc(2,2,next+1) ) 
                   sabox_loc(1,ndime,next+1) = min( fabox_loc(1,ndime,jface) , sabox_loc(1,ndime,next+1) )               
                   sabox_loc(2,ndime,next+1) = max( fabox_loc(2,ndime,jface) , sabox_loc(2,ndime,next+1) ) 
                end do
             end if
             !
             ! Store the children of current element of the stack (struc)
             !
             indst          = indst + 1_ip
             struc(indst,1) = next
             struc(indst,2) = imini
             struc(indst,3) = imedi                
             indst          = indst + 1_ip
             struc(indst,1) = next  + 1_ip
             struc(indst,2) = imedi + 1_ip
             struc(indst,3) = imaxi             
             blink_lev(next)    = blink_lev(curr) + 1
             blink_lev(next+1)  = blink_lev(curr) + 1
             next           = next  + 2_ip       

             deallocate( vect3 )
          end if

       end do

       call cputim(time2)

       deallocate( centr )
       deallocate( perm  )
       deallocate( struc )
       deallocate( blink_lev)

    end if
  end subroutine kdtree

  subroutine dpopar(                              &
       itask,xcoor,npoib,mnodb,nboun,chkdi,ltypb, &
       lnodb,coord,dista,norma,proje,bound,fabox, &
       sabox,blink,stru2,ldist,lnele,lmask)

    !-----------------------------------------------------------------------
    ! NAME
    !    nepoib
    ! DESCRIPTION
    !    Shortest signed distance from a given point to the particle
    !    INPUT
    !       XCOOR ... Test point coordinates
    !       SABOX_LOC ... Tree stru2_locture
    !       BLINK_LOC ... Tree stru2_locture
    !       NBOUN ... Number of boundaries
    !       LTYPB ... Types of boundaries
    !       LNODB ... Boundary connectivity
    !       COORD ... Boundary node coordinates
    !       CHKDI ... Distance use to check the boundary box of the faces in the particle 
    !                 for find the shortest distance (use chkdi = 1.0e10_rp in general)
    !    OUTPUT
    !       DISTA ... is less than 0 is point is inside the particle
    !                 is more than 0 is point is outside the particle
    !       NORMA ... Normal to the plane that contains the projection of the point
    !       PROJE ... Nearest projection point to the particle surface
    !       
    ! USED BY
    !    inouib
    !----------------------------------------------------------------------- 
    
    integer(ip), intent(in)                       :: itask
    integer(ip), intent(in)                       :: npoib
    integer(ip), intent(in)                       :: mnodb
    integer(ip), intent(in)                       :: nboun
    integer(ip), intent(in)                       :: ltypb(nboun)
    integer(ip), intent(in)                       :: lnodb(mnodb,nboun)
    real(rp),    intent(in)                       :: xcoor(ndime)
    real(rp),    intent(in)                       :: chkdi
    real(rp),    intent(in)                       :: coord(ndime,*)
    real(rp),    intent(out)                      :: dista
    real(rp),    intent(out)                      :: norma(ndime)
    real(rp),    intent(out)                      :: proje(ndime)
    integer(ip), intent(out)                      :: bound
    real(rp),    intent(in),    target, optional  :: sabox(2,ndime,2*nboun-1)
    integer(ip), intent(in),    target, optional  :: blink(2*nboun-1)
    real(rp),    intent(in),    target, optional  :: fabox(2,ndime,nboun)
    type(netyp), intent(inout), target, optional  :: lnele(npoib)
    integer(ip), intent(inout), target, optional  :: stru2(2*nboun-1)
    real(rp),    intent(inout), target, optional  :: ldist(2*nboun-1)
    logical(lg), intent(in),    pointer,optional  :: lmask(:)
    integer(ip)                                   :: iboun
    integer(ip)                                   :: idime,ipoib,indst,inodb,curr,ifoun
    integer(ip)                                   :: pblty,pnodb,qnodb, jboun,kboun,ilist
    integer(ip)                                   :: jlist,node1,node2,inde1,ntria
    real(rp)                                      :: temdi,dismi1,dismi2,temp,tmpfa,facto
    real(rp)                                      :: tenor(3),propo(3),toler,pladi,numer,dummi
    real(rp)                                      :: bocod(ndime,mnodb),bari1,bari2,denom,t
    real(rp)                                      :: x01(3),x02(3),x12(3),x012(3)
    logical(lg)                                   :: if_mask
    logical(lg)                                   :: consider_iboun

    if_mask = .false.
    if( present(lmask) ) then
       if( associated(lmask) ) if_mask = .true.
    end if

    if( present(fabox) ) then
       fabox_loc => fabox
       sabox_loc => sabox
       blink_loc => blink
       stru2_loc => stru2
       ldist_loc => ldist
       lnele_loc => lnele
    end if

    toler = 1.0e-6_rp
    dista = 0.0_rp

    do idime = 1,ndime
       temp  = max(xcoor(idime)-sabox_loc(1,idime,1), sabox_loc(2,idime,1)-xcoor(idime))
       dista = dista + temp*temp
    end do

    bound = 0_ip
    if ( chkdi*chkdi < dista ) then       
       dista = chkdi*chkdi       
    end if

    indst            = 1_ip
    stru2_loc(indst) = 1_ip

    dismi1 = 0.0_rp
    do idime = 1,ndime
       temp   = max(0.0_rp,sabox_loc(1,idime,1)-xcoor(idime))  
       temp   = max(temp,xcoor(idime)-sabox_loc(2,idime,1))  
       dismi1 = dismi1 + temp * temp
    end do
    ldist_loc(1) = dismi1
    !
    ! Assemble a list of candidate patches by traversing skd-tree
    !
    do while( indst > 0_ip )
       curr  = stru2_loc(indst)      
       indst = indst - 1_ip
       !
       ! Minimum distance
       !       
       if ( ldist_loc(curr) < dista ) then
          !
          ! If currnode is a leaf in the tree stru2_locture
          !      
          if( blink_loc(curr) < 0_ip ) then
             !
             ! Find the closest projection to the particle 
             !
             iboun = -blink_loc(curr)
             if( if_mask ) then
                consider_iboun = lmask(iboun)
             else
                consider_iboun = .true.
             end if
             if( consider_iboun ) then

                pblty =  ltypb(iboun)
                pnodb =  nnode(pblty)                  

                do inodb = 1,pnodb
                   ipoib              = lnodb(inodb,iboun)     
                   bocod(1    ,inodb) = coord(1    ,ipoib)
                   bocod(2    ,inodb) = coord(2    ,ipoib)
                   bocod(ndime,inodb) = coord(ndime,ipoib)                
                end do

                if( ndime == 3 .and. pnodb == 2 ) then
                   !
                   ! Distance to a BAR02 element (segment) in 3D
                   !
                   ! http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
                   !
                   ! Find the value of t that minimizes the distance from the point to the line. If t is between 0.0 and 1.0, 
                   ! then the closest point lies on the segment, otherwise the closest point is one of the segment's end points.
                   ! 
                   do idime = 1,3
                      x01(idime) = bocod(idime,1) - xcoor(idime)
                      x12(idime) = bocod(idime,2) - bocod(idime,1)
                      x02(idime) = bocod(idime,2) - xcoor(idime)
                   end do
                   dismi1 = dot_product(x01(1:3),x12(1:3))
                   dismi2 = dot_product(x12(1:3),x12(1:3))
                   t      = - dismi1 / (dismi2+epsil)
                   if( t < 0.0_rp ) then
                      temdi = sqrt(dot_product(x01(1:3),x01(1:3)))
                      t     = 0.0_rp 
                   else if( t > 1.0_rp ) then
                      temdi = sqrt(dot_product(x02(1:3),x02(1:3)))
                      t     = 1.0_rp 
                   else
                      call kdtree_vecpro(x01,x02,x012,3_ip)
                      dismi1 = dot_product(x012(1:3),x012(1:3))
                      temdi  = sqrt(dismi1) / (sqrt(dismi2)+epsil)
                   end if
                   temdi      = temdi * temdi
                   propo(1:3) = bocod(1:3,1) + x12(1:3) * t
                   
                else
                   !
                   ! Distance to plane in 3D and segment in 2D
                   ! TENOR= Exterior normal           
                   !
                   call extbou(itask,pnodb,lnodb(1,iboun),coord,tenor)

                   denom = dot_product(tenor(1:ndime),tenor(1:ndime))
                   if( denom > epsil ) then
                      call dpofac(pnodb,xcoor,bocod,fabox_loc(1,1,iboun),tenor,temdi,tmpfa,propo)
                   else
                      temdi = huge(1.0_rp)
                   end if

                end if

                if( temdi < dista ) then               
                   dista          = temdi
                   bound          = iboun
                   facto          = tmpfa
                   norma(1:ndime) = tenor(1:ndime)
                   proje(1:ndime) = propo(1:ndime)         
                end if

             end if

          else                         

             dismi1 = 0.0_rp
             dismi2 = 0.0_rp
             do idime = 1,ndime
                temp   = 0.0_rp
                temp   = max(0.0_rp,sabox_loc(1,idime,blink_loc(curr))-xcoor(idime))  
                temp   = max(temp,xcoor(idime)-sabox_loc(2,idime,blink_loc(curr)))  
                dismi1 = dismi1 + temp * temp
                temp   = 0.0_rp
                temp   = max(0.0_rp,sabox_loc(1,idime,blink_loc(curr)+1_ip)-xcoor(idime))  
                temp   = max(temp,xcoor(idime)-sabox_loc(2,idime,blink_loc(curr)+1_ip) )  
                dismi2 = dismi2 + temp * temp
             end do
             ldist_loc(blink_loc(curr))    = dismi1             
             ldist_loc(blink_loc(curr)+1)  = dismi2

             indst          = indst + 2_ip
             if (dismi1 > dismi2) then               
                stru2_loc(indst-1) = blink_loc(curr)
                stru2_loc(indst)   = blink_loc(curr) + 1_ip
             else
                stru2_loc(indst-1) = blink_loc(curr) + 1_ip
                stru2_loc(indst)   = blink_loc(curr)                              
             end if

          end if
       end if
    end do

    if( bound /= 0_ip ) then

       if(  ndime == 3 .and. pnodb == 2 ) then
          !
          ! Particular case of BAR02 elements in 3D
          ! The normal does not make sense for BAR02 elements in 3D!
          ! Neither the signed distance, so distance is positive
          !
          dista = sqrt(dista) 
          norma = 0.0_rp

       else
          !
          ! Determine if the projection point is near to a node or a segment of the face
          !
          iboun = bound
          pblty = ltypb(iboun)
          pnodb = nnode(pblty) 
          do inodb = 1,pnodb
             ipoib = lnodb(inodb,iboun)            
             bocod(1    ,inodb) = coord(1    ,ipoib)
             bocod(2    ,inodb) = coord(2    ,ipoib)
             bocod(ndime,inodb) = coord(ndime,ipoib)                
          end do

          if( ndime == 3 ) then
             !
             ! 3D
             !
             call instr2(pnodb,proje,bocod,ifoun,bari1,bari2,ntria)
             node1 = 0_ip
             node2 = 0_ip
             !
             ! Check if the projection point is near to a node
             !
             if (pnodb == 3) then
                if ( bari1 > 1.0_rp-toler .and. mnodb > 2 ) then
                   node1 = lnodb(3,iboun)       
                   inde1 = 3
                else if ( bari2 > 1.0_rp-toler ) then
                   node1 = lnodb(2,iboun)       
                   inde1 = 2
                else if ( 1.0_rp-bari1-bari2  > 1.0_rp-toler ) then
                   node1 = lnodb(1,iboun)       
                   inde1 = 1
                else if ( bari1 < toler ) then
                   node1 = lnodb(1,iboun)     
                   node2 = lnodb(2,iboun)     
                   inde1 = 1
                else if ( bari2 < toler .and. mnodb > 2 ) then
                   node1 = lnodb(1,iboun)     
                   node2 = lnodb(3,iboun)     
                   inde1 = 1
                else if ( 1.0_rp-bari1-bari2 < toler .and. mnodb > 2 ) then
                   node1 = lnodb(2,iboun)     
                   node2 = lnodb(3,iboun)     
                   inde1 = 2
                end if
             else if (pnodb == 4) then
                if (ntria == 1) then
                   if ( bari1 > 1.0_rp-toler ) then
                      node1 = lnodb(3,iboun)       
                      inde1 = 3
                   else if ( bari2 > 1.0_rp-toler ) then
                      node1 = lnodb(2,iboun)       
                      inde1 = 4
                   else if ( 1.0_rp-bari1-bari2  > 1.0_rp-toler ) then
                      node1 = lnodb(1,iboun)       
                      inde1 = 1
                   else if ( bari1 < toler ) then
                      node1 = lnodb(1,iboun)     
                      node2 = lnodb(2,iboun)     
                      inde1 = 1
                   else if ( 1.0_rp-bari1-bari2 < toler ) then
                      node1 = lnodb(2,iboun)     
                      node2 = lnodb(3,iboun)     
                      inde1 = 4
                   end if
                else if (ntria == 2) then
                   if ( bari1 > 1.0_rp-toler ) then
                      node1 = lnodb(1,iboun)       
                      inde1 = 1
                   else if ( bari2 > 1.0_rp-toler ) then
                      node1 = lnodb(4,iboun)       
                      inde1 = 4
                   else if ( 1.0_rp-bari1-bari2  > 1.0_rp-toler ) then
                      node1 = lnodb(1,iboun)       
                      inde1 = 1
                   else if ( bari1 < toler ) then
                      node1 = lnodb(3,iboun)     
                      node2 = lnodb(4,iboun)     
                      inde1 = 3
                   else if ( 1.0_rp-bari1-bari2 < toler ) then
                      node1 = lnodb(4,iboun)     
                      node2 = lnodb(1,iboun)     
                      inde1 = 4
                   end if
                end if
             end if
             !
             ! It is near to a segmnent
             !
             if ( node2 /= 0_ip ) then
                norma = 0.0_rp
                !
                ! Sum all the exterior normals that use the segment
                !                    
                !
                ! Determine if the point are really inside or outside the particle
                !        
                facto = -1.0_rp          
                do ilist = 1,lnele_loc(node1)%nelem
                   do jlist = 1,lnele_loc(node2)%nelem
                      jboun = lnele_loc(node1)%eleme(ilist)
                      kboun = lnele_loc(node2)%eleme(jlist)
                      if (jboun == kboun) then
                         qnodb = nnode(lnele_loc(node1)%ltype(ilist))
                         call extbou(itask,qnodb,lnodb(1,jboun),coord,tenor)
                         do idime = 1,ndime  
                            norma(idime) = norma(idime) + tenor(idime)
                         end  do
                         call vecuni(ndime,norma,dummi)

                         pladi = dot_product(tenor(1:ndime),(xcoor(1:ndime) - bocod(1:ndime,inde1)))
                         if ( pladi >= 0.0_rp ) facto = 1.0_rp                
                      end if
                   end  do
                end  do
             else if ( node1 /= 0_ip ) then
                !
                ! It is near a node
                !
                norma = 0.0_rp
                !
                ! Sum all the exterior normals that use the node 
                !       
                !
                ! Determine if the point are really inside or outside the particle
                !       
                facto = -1.0_rp   
                do ilist = 1,lnele_loc(node1)%nelem
                   jboun = lnele_loc(node1)%eleme(ilist)
                   qnodb = nnode(lnele_loc(node1)%ltype(ilist))
                   call extbou(itask,qnodb,lnodb(1,jboun),coord,tenor)
                   do idime = 1,ndime  
                      norma(idime) = norma(idime) + tenor(idime)
                   end  do
                   call vecuni(ndime,norma,dummi)
                   pladi = dot_product(tenor(1:ndime),(xcoor(1:ndime) - bocod(1:ndime,inde1)))
                   if ( pladi >= 0.0_rp ) facto = 1.0_rp
                end  do

             end if

          else
             !
             ! 2D
             !
             node1 = 0_ip
             !
             ! Check if the projection point is near to a node
             !
             numer = 0.0_rp
             denom = 0.0_rp
             do idime = 1,ndime
                numer = numer + (bocod(idime,2) - bocod(idime,1)) * (xcoor(idime)   - bocod(idime,1))
                denom = denom + (bocod(idime,2) - bocod(idime,1)) * (bocod(idime,2) - bocod(idime,1))
             end do
             bari1 = numer / (denom+epsil)
             if( bari1 < toler ) then
                node1 = lnodb(1,iboun)
             else if( bari1 > 1.0_rp - toler ) then
                node1 = lnodb(2,iboun)
             end if
             !
             ! It is near to a node
             !
             if ( node1 /= 0_ip ) then
                norma(1)     = 0.0_rp
                norma(2)     = 0.0_rp
                !
                ! Sum all the exterior normals that use the node 
                !       
                !
                ! Determine if the point are really inside or outside the particle
                !       
                facto = -1.0_rp   
                do ilist = 1,lnele_loc(node1)%nelem
                   jboun = lnele_loc(node1)%eleme(ilist)
                   qnodb = nnode(lnele_loc(node1)%ltype(ilist))
                   call extbou(itask,qnodb,lnodb(1,jboun),coord,tenor)
                   do idime = 1,ndime  
                      norma(idime) = norma(idime) + tenor(idime)
                      call vecuni(ndime,norma,dummi)
                   end  do

                   pladi = 0.0_rp
                   do idime = 1,ndime
                      pladi = pladi + tenor(idime) * ( xcoor(idime) - bocod(idime,1) )
                   end do
                   if ( sign(1.0_rp,pladi) > 0.0_rp ) facto = 1.0_rp
                end  do
             end if

          end if
          !
          ! Signed distance
          !           
          dista = facto * sqrt(dista)
       end if

    else
       dista = 1.0e10_rp
    end if

  end subroutine dpopar

  pure subroutine kdtree_vecpro(v1,v2,v3,n)

    !-----------------------------------------------------------------------
    !
    ! Two and three-dimensional vectorial product of two vectors  v3 = v1 x v2.
    ! The same pointer as for v1 or v2 may be used for v3. If N = 2, it is
    !  assumed that v1 = (0,0,v1_3) and v2 = (v2_1,v2_2,0).      
    !
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v2(n)
    real(rp),    intent(in)  :: v1(3)
    real(rp),    intent(out) :: v3(n)
    real(rp)                 :: c1,c2,c3

    if(n==2) then
       c1=-v1(3)*v2(2)
       c2= v1(3)*v2(1)
       v3(1)=c1
       v3(2)=c2
    else if(n==3) then
       c1=v1(2)*v2(3)-v1(3)*v2(2)
       c2=v1(3)*v2(1)-v1(1)*v2(3)
       c3=v1(1)*v2(2)-v1(2)*v2(1)
       v3(1)=c1
       v3(2)=c2
       v3(3)=c3
    end if

  end subroutine kdtree_vecpro

end module mod_kdtree
