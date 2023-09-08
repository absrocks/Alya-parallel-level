module mod_maillage

!!$  use, intrinsic :: ISO_C_BINDING
!!$  use def_kintyp,    only : ip,rp
!!$  use mod_memory,    only : memory_alloca
!!$  use mod_memory,    only : memory_deallo
!!$  use def_domain,    only : memor_dom
!!$  use def_elmtyp,    only : TRI03,TET04,BAR02
!!$  use def_domain,    only : ndime,npoin,nelem,lnods,coord,mnode
!!$  use def_domain,    only : ltype,lelch,lesub,lmast,lnoch
!!$  use def_domain,    only : nboun,lnodb,mnodb,lmate,lnodb,lboel,lboch,ltypb
!!$  use mod_mesh_type
!!$  implicit none
!!$  !
!!$  ! https://docs.oracle.com/cd/E19059-01/stud.9/817-6694/11_cfort.html
!!$  !  
!!$  type, BIND (C) :: typ_adr 
!!$     integer(C_INT) :: X,T       
!!$  end type typ_adr
!!$
!!$  type, BIND (C) :: typ_lng 
!!$     integer(C_INT) :: m
!!$     integer(C_INT) :: x,t,sx,st  
!!$  end type typ_lng
!!$
!!$  type, BIND (C) :: typ_maillage 
!!$     integer(C_INT) :: d,dim
!!$     integer(C_INT) :: nX,nT
!!$     type(C_PTR)    :: T
!!$     type(C_PTR)    :: X
!!$     integer(C_INT) :: nF
!!$     type(C_PTR)    :: I,T2N,T2T
!!$     type(C_PTR)    :: Ir,T2Nr
!!$     type(typ_adr)  :: libre
!!$     type(C_PTR)    :: Nom
!!$     type(typ_lng)  :: L
!!$     type(C_PTR)    :: H
!!$     type(C_PTR)    :: metrique 
!!$     type(C_PTR)    :: S
!!$  end type typ_maillage
!!$
!!$  !interface
!!$  !   !extern maillage * procedure_mtc_metrique(const maillage *M,char *,int,int,const double*,double **);
!!$  !   function procedure_mtc_metrique(M,c1,i1,i2,metrique,d2) BIND (C,NAME='procedure_mtc_metrique')
!!$  !     type(C_PTR)        :: M
!!$  !     type(C_PTR)        :: c1
!!$  !     integer(ip), value :: i1
!!$  !     integer(ip), value :: i2
!!$  !     type(C_PTR)        :: metrique
!!$  !     type(C_PTR)        :: d2
!!$  !     type(C_PTR)        :: procedure_mtc_metrique
!!$  !   end function procedure_mtc_metrique
!!$  !end interface
!!$
!!$  interface
!!$     subroutine maillage_proc(Mesh_in,Mesh_out,metrique) BIND (C,NAME='maillage_proc')
!!$       use iso_c_binding
!!$       type(C_PTR), VALUE  :: Mesh_in
!!$       type(C_PTR)         :: Mesh_out
!!$       type(C_PTR), VALUE  :: metrique
!!$     end subroutine maillage_proc
!!$  end interface 
!!$
!!$  interface
!!$     subroutine maillage_detruire(Mesh_in) BIND (C,NAME='maillage_detruire')
!!$       use iso_c_binding
!!$       type(C_PTR),VALUE        :: Mesh_in
!!$     end subroutine maillage_detruire
!!$  end interface
!!$
!!$  character(100), PARAMETER :: vacal = "mod_maillage"
!!$
!!$contains
!!$
!!$  subroutine maillage_invocation(metrique)
!!$
!!$    implicit none 
!!$
!!$    real(rp),    pointer, intent(in) :: metrique(:,:,:)
!!$
!!$    type(typ_maillage)               :: Mesh_in
!!$    type(typ_maillage), pointer      :: Mesh_out
!!$
!!$
!!$    integer(ip)                      :: d
!!$    integer(ip)                      :: nenti
!!$    real(rp),    pointer             :: X(:)
!!$    integer(ip), pointer             :: T(:)
!!$    real(rp),    pointer             :: M(:)
!!$
!!$    integer(ip)                      :: nenti_new
!!$    integer(ip)                      :: npoin_new
!!$    real(rp),    pointer             :: X_new(:)
!!$    integer(ip), pointer             :: T_new(:)
!!$
!!$    type(c_ptr)                      :: px_new,pnx_new
!!$    type(c_ptr)                      :: pt_new,pnt_new
!!$
!!$    integer(ip)                      :: ipoin,inode,Tdim,iboun
!!$    integer(ip)                      :: idime,jdime,kpoin,kelem,ielem
!!$    integer(ip)                      :: npoin_
!!$
!!$        
!!$    nullify(X,T,M)
!!$    ! nullify(X_new,T_new)
!!$
!!$    d     = ndime + 1
!!$    nenti = nelem + nboun
!!$    
!!$    allocate(T(0:nenti*d-1),X(0:(npoin+1)*ndime-1),M(0:(npoin+1)*ndime*ndime-1))
!!$    !
!!$    ! LNODS
!!$    !
!!$    kelem = 0
!!$    do ielem = 1,nelem
!!$       do inode = 1,d
!!$          T(kelem) = lnods(inode,ielem)
!!$          kelem = kelem + 1
!!$       end do
!!$    end do
!!$    do iboun = 1,nboun
!!$       do inode = 1,mnodb
!!$          T(kelem) = lnodb(inode,iboun)
!!$          kelem = kelem + 1
!!$       end do
!!$       T(kelem) = 0_ip
!!$       kelem = kelem + 1
!!$    end do
!!$    !
!!$    ! COORD
!!$    !
!!$    X     = 0.0_rp
!!$    kpoin = ndime
!!$    do ipoin = 1,npoin
!!$       do idime = 1,ndime
!!$          X(kpoin) = coord(idime,ipoin)
!!$          kpoin = kpoin + 1
!!$       end do
!!$    end do
!!$    !
!!$    ! METRIC
!!$    !
!!$    allocate(M(0:(npoin+1)*ndime*ndime-1))
!!$    M     = 0.0_rp
!!$    kpoin = ndime*ndime
!!$    do ipoin = 1,npoin       
!!$       do idime = 1,ndime
!!$          do jdime = 1,ndime
!!$             M(kpoin) = metrique(idime,jdime,ipoin)
!!$             kpoin = kpoin + 1 
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$    print*,'GENERATE MESH'
!!$    call maillage_bridge(Mesh_in,Mesh_out,T,X,M)
!!$    deallocate(T,X,M)
!!$    !
!!$    ! Remeshing
!!$    ! 
!!$    !call maillage_proc(d,ndime,npoin,nenti,X,T,M,pnx_new,pnt_new,px_new,pt_new)
!!$    print*,'END GENERATE MESH'
!!$    !
!!$    ! Reconstruct mesh
!!$    !
!!$    !call c_f_pointer(pnx_new,npoin_new)
!!$    !print*,'nx,nt=',npoin_new,nenti_new
!!$    !print*,'alloc=',associated(X_new),size(X_new,KIND=ip)
!!$    !call c_f_pointer(px_new,X_new,shape=[npoin_new*ndime])
!!$
!!$    npoin = npoin_new
!!$
!!$
!!$    call outdom(1_ip)
!!$
!!$    deallocate(X,T) 
!!$    call runend("O.K.") 
!!$
!!$  end subroutine maillage_invocation
!!$
!!$  subroutine maillage_bridge(Mesh_in,Mesh_out,T,X,M)
!!$    
!!$    type(typ_maillage), target   :: Mesh_in
!!$    type(typ_maillage), pointer  :: Mesh_out
!!$    type(c_ptr)                  :: Mesh_out_ptr
!!$    
!!$    integer(ip),    intent(inout), pointer :: T(:)
!!$    real(rp),       intent(inout), pointer :: X(:)
!!$    real(rp),       intent(inout), pointer :: M(:)
!!$    
!!$    integer(C_INT),                target  :: T_c(size(T,DIM=1,KIND=ip))
!!$    real(C_DOUBLE),                target  :: X_c(size(X,DIM=1,KIND=ip))
!!$    real(C_DOUBLE),                target  :: M_c(size(M,DIM=1,KIND=ip))
!!$    
!!$    integer(ip),                   pointer :: T_out(:)
!!$    real(C_DOUBLE),                pointer :: X_out(:)
!!$    real(rp),                      pointer :: M_out(:)
!!$    integer(ip)                            :: Tdim,Xdim,npoin_new,nenti_new,nboun_new
!!$
!!$    T_c = T
!!$    X_c = X 
!!$    M_c = M
!!$
!!$    Mesh_in % d   = ndime+1
!!$    Mesh_in % dim = ndime
!!$    Mesh_in % nX  = npoin
!!$    Mesh_in % nT  = nelem+nboun    
!!$    Mesh_in % T   = C_LOC(T_c)
!!$    Mesh_in % X   = C_LOC(X_c)
!!$
!!$    call maillage_proc(C_LOC(Mesh_in),Mesh_out_ptr,C_LOC(M_c))
!!$    CALL C_F_POINTER(Mesh_out_ptr, Mesh_out)
!!$    
!!$    Xdim = (Mesh_out % nX + 1)* Mesh_out % dim
!!$    call c_f_pointer(Mesh_out % X, X_out, [Xdim]) 
!!$    Tdim = (Mesh_out % nT )* Mesh_out % d
!!$    call c_f_pointer(Mesh_out % T, T_out, [Tdim])
!!$    
!!$    nboun_new = Mesh_out % nF 
!!$    npoin_new = Mesh_out % nX
!!$    nenti_new = Mesh_out % nT
!!$    
!!$    call maillage_newmesh(npoin_new,nenti_new,nboun_new,T_out,X_out)
!!$    !
!!$    ! Destroy mesh
!!$    !
!!$    print*,'a=',associated(T_out),size(T_out,KIND=ip)
!!$    call maillage_detruire(Mesh_out_ptr)
!!$
!!$    print*,'b=',associated(T_out),size(T_out,KIND=ip)
!!$    
!!$  end subroutine maillage_bridge
!!$
!!$  subroutine maillage_newmesh(npoin_new,nenti_new,nboun_new,T_new,X_new)
!!$
!!$    integer(ip)                 :: npoin_new
!!$    integer(ip)                 :: nenti_new
!!$    integer(ip)                 :: nboun_new
!!$    integer(ip),     intent(inout), pointer :: T_new(:)
!!$    real(C_DOUBLE),                 pointer :: X_new(:)
!!$    integer(ip)        :: d,kelem,inode,ielem,idime,iboun,kpoin,ipoin,nbou1
!!$    !
!!$    !
!!$    !
!!$    npoin = npoin_new
!!$    nboun = nboun_new
!!$    nelem = nenti_new-nboun
!!$    nbou1 = max(1_ip,nboun)
!!$    d     = ndime+1
!!$    
!!$    !kelem = 0
!!$    !inode = 0
!!$    !loop_nenti: do ielem = 1,nenti_new
!!$    !   if( T_new(ielem*d) == 0 ) then
!!$    !      nelem = ielem-1
!!$    !      exit loop_nenti
!!$    !   end if
!!$    !end do loop_nenti
!!$
!!$    print*,'A=',npoin,nelem,nboun
!!$    
!!$    call memory_deallo(memor_dom,'LTYPE'    ,'memgeo' , ltype    )
!!$    call memory_deallo(memor_dom,'LELCH'    ,'memgeo' , lelch    )
!!$    call memory_deallo(memor_dom,'LNODS'    ,'memgeo' , lnods    )
!!$    call memory_deallo(memor_dom,'LESUB'    ,'memgeo' , lesub    )
!!$    call memory_deallo(memor_dom,'LMATE'    ,'memgeo' , lmate    )
!!$
!!$    call memory_deallo(memor_dom,'COORD'    ,'memgeo' , coord    )
!!$    call memory_deallo(memor_dom,'LNOCH'    ,'memgeo' , lnoch    )
!!$    call memory_deallo(memor_dom,'LMAST'    ,'memgeo' , lmast    )
!!$
!!$    call memory_deallo(memor_dom,'LNODB'    ,'memgeo' , lnodb    )
!!$    call memory_deallo(memor_dom,'LBOEL'    ,'memgeo' , lboel    )
!!$    call memory_deallo(memor_dom,'LTYPB'    ,'memgeo' , ltypb    )
!!$    call memory_deallo(memor_dom,'LBOCH'    ,'memgeo' , lboch    )
!!$    
!!$    call memory_alloca(memor_dom,'LTYPE'    ,'memgeo' , ltype    , nelem   )
!!$    call memory_alloca(memor_dom,'LELCH'    ,'memgeo' , lelch    , nelem   )
!!$    call memory_alloca(memor_dom,'LNODS'    ,'memgeo' , lnods    , mnode   , nelem )
!!$    call memory_alloca(memor_dom,'LESUB'    ,'memgeo' , lesub    , nelem   )
!!$    call memory_alloca(memor_dom,'LMATE'    ,'memgeo' , lmate    , nelem )
!!$
!!$    call memory_alloca(memor_dom,'COORD'    ,'memgeo' , coord    , ndime   , npoin )
!!$    call memory_alloca(memor_dom,'LNOCH'    ,'memgeo' , lnoch    , npoin   )
!!$    call memory_alloca(memor_dom,'LMAST'    ,'memgeo' , lmast    , npoin   )
!!$
!!$     call memory_alloca(memor_dom,'LNODB'   ,'memgeo' , lnodb   , mnodb   , nbou1 )
!!$     call memory_alloca(memor_dom,'LBOEL'   ,'memgeo' , lboel   , mnodb   , nbou1 )
!!$     call memory_alloca(memor_dom,'LTYPB'   ,'memgeo' , ltypb   , nbou1   )
!!$     call memory_alloca(memor_dom,'LBOCH'   ,'memgeo' , lboch   , nbou1   )    
!!$    !
!!$    ! LNODS
!!$    !
!!$    kelem = 1
!!$    do ielem = 1,nelem
!!$       do inode = 1,d
!!$          lnods(inode,ielem) = T_new(kelem) 
!!$          kelem = kelem + 1
!!$       end do
!!$       ltype(ielem) = TRI03
!!$    end do
!!$    do iboun = 1,nboun
!!$       do inode = 1,mnodb
!!$          lnodb(inode,iboun) = T_new(kelem) 
!!$          kelem = kelem + 1
!!$       end do
!!$       ltypb(iboun) = BAR02
!!$       kelem = kelem + 1
!!$    end do
!!$    !
!!$    ! COORD
!!$    !
!!$    kpoin = ndime+1
!!$    do ipoin = 1,npoin
!!$       do idime = 1,ndime
!!$          coord(idime,ipoin) = X_new(kpoin) 
!!$          kpoin = kpoin + 1
!!$       end do
!!$    end do
!!$
!!$    write(90,*) ' '
!!$    write(90,*) ' '
!!$    write(90,*) ' '
!!$    write(90,*) 'coordinates'
!!$    do ipoin=1,npoin
!!$       write(90,*) ipoin,coord(1:ndime,ipoin)
!!$    end do
!!$    write(90,*) 'end coordinates'
!!$    write(90,*) 'elements'
!!$    do ielem=1,nelem
!!$       write(90,*) ielem,lnods(:,ielem)
!!$    end do
!!$    write(90,*) 'end elements'
!!$
!!$    
!!$    call mesh_type_update_last_mesh()
!!$    call outdom(1_ip)
!!$
!!$    call runend('O.K.!')
!!$    
!!$  end subroutine maillage_newmesh
  
end module mod_maillage
 
