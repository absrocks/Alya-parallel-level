!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_domain.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_domain

  use def_kintyp_basic

  integer(ip), parameter  :: &
       nelty =  60, &                          ! # of element types
       mfiel = 500                             ! Maximum number of fields

  !----------------------------------------------------------------------
  !
  ! Element data base type
  !
  !----------------------------------------------------------------------

  type elm
     ! User integration rule
     integer(ip)          :: pgaus            ! Number of Gauss points
     real(rp),    pointer :: shape(:,:)       ! pnode,pgaus
     real(rp),    pointer :: deriv(:,:,:)     ! ndime,pnode,pgaus
     real(rp),    pointer :: heslo(:,:,:)     ! ntens,pnode,pgaus
     real(rp),    pointer :: weigp(:)
     real(rp),    pointer :: shaga(:,:)
     real(rp),    pointer :: posgp(:,:)
     ! Bubble
     real(rp),    pointer :: shape_bub(:)     ! pgaus
     real(rp),    pointer :: deriv_bub(:,:)   ! ndime,pgaus
     real(rp),    pointer :: heslo_bub(:,:)   ! ntens,pgaus
     ! Center of gravity
     real(rp),    pointer :: shacg(:)
     real(rp),    pointer :: dercg(:,:)
     real(rp),    pointer :: hescg(:,:)
     real(rp)             :: weicg
     ! Close rule
     real(rp),    pointer :: shapc(:,:)
     real(rp),    pointer :: deric(:,:,:)
     real(rp),    pointer :: heslc(:,:,:)
     real(rp),    pointer :: weigc(:)
     ! IB integration rule
     real(rp),    pointer :: shaib(:,:)
     real(rp),    pointer :: derib(:,:,:)
     real(rp),    pointer :: weiib(:)
  end type elm

  type elmgp
     real(rp),    pointer :: gpvol(:)
     real(rp),    pointer :: gpcar(:,:,:)
     real(rp),    pointer :: gphes(:,:,:)
     real(rp),    pointer :: hleng(:)
     real(rp),    pointer :: tragl(:,:)
  end type elmgp
  
  !----------------------------------------------------------------------
  !
  ! Mesh
  !
  !----------------------------------------------------------------------

  type cell

     real(rp)      :: coor(3,2)
     integer(ip)   :: neigh(6)
     integer(ip)   :: level
     integer(ip)   :: marked
     real(rp)      :: rsize

  end type cell
  
   type typ_lobas
     real(rp),   pointer :: local_basis(:,:,:)
     integer(ip)         :: type
     real(rp)            :: param(9)
  end type typ_lobas

  !----------------------------------------------------------------------
  !
  ! Element bin
  !
  !----------------------------------------------------------------------

  type typ_element_bin
     real(rp)             :: comin(3)
     real(rp)             :: comax(3)
     integer(ip)          :: boxes(3)
     integer(ip), pointer :: bin_size(:,:,:)
     type(i1p),   pointer :: list_elements(:,:,:)
  end type typ_element_bin

  !----------------------------------------------------------------------
  !
  ! OMPSS_DOMAIN
  !
  !----------------------------------------------------------------------

  type ompss_domain
     integer(ip), pointer :: neighbours(:)
     integer(ip), pointer :: elements(:)
     integer(ip)          :: neighIdx
     integer(ip)          :: elemIdx
  end type ompss_domain

  !------------------------------------------------------------------------
  !
  ! KD-TREE
  !
  !------------------------------------------------------------------------

  type netyp
     integer(ip)          :: nelem             ! Number of elements for a node
     integer(ip), pointer :: eleme(:)          ! Ids of elements
     integer(ip), pointer :: ltype(:)          ! Types of elements
  end type netyp
  
end module def_kintyp_domain
!> @}
