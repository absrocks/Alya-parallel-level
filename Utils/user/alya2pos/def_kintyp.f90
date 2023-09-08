module def_kintyp
  !-----------------------------------------------------------------------
  !****f* defmod/def_kintyp
  ! NAME
  !   def_kintyp
  ! DESCRIPTION
  !   Definition of kinds and types. 
  !   "The range of the default integers is not specified in the language
  !   but on a computer with a word size of n bits, is often from 
  !   -2^{n-1} to +2^{n-1}-1. Thus on a 32-bit computer the range is
  !   often -2.14*10^9 to +2.14*10^9."
  !   M. Metclaf and J. Reid, FORTRAN 90/95 explained, 2nd edition.
  !
  !   Defaults are:
  !   Integers: 4-bytes 
  !   Reals:    8-bytes
  !
  !***
  !-----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !
  ! Symbolc names for integers, reals and logicals
  !
  !----------------------------------------------------------------------

  !
  ! Symbolic names for integers
  !
# ifdef I8
  integer, parameter  :: ip = 8             ! 8-byte integer
# else
  integer, parameter  :: ip = 4             ! 4-byte integer
# endif
  integer, parameter  :: rp = 8             ! Double precision 
  !
  ! Symbolic name for kind type of default logical
  !
  integer, parameter  :: lg = kind(.true.)

  !----------------------------------------------------------------------
  !
  ! General types
  !
  !----------------------------------------------------------------------

  type i1p
     integer(ip), pointer :: l(:)
  end type i1p
  type i2p
     integer(ip), pointer :: l(:,:)
  end type i2p
  type r1p
     real(rp),    pointer :: a(:)
  end type r1p
  type r2p
     real(rp),    pointer :: a(:,:)
  end type r2p
  type r3p
     real(rp),    pointer :: a(:,:,:)
  end type r3p

  integer(ip), parameter   :: &
       nelty = 100 
  integer(ip)              :: &
       ldime(nelty),          &      ! List of element dimensions
       ltopo(nelty),          &      ! List of element topology
       ngaus(nelty),          &      ! List of element Laplacian
       llapl(nelty),          &      ! List of element Laplacian
       lrule(nelty),          &      ! List of element integration rules
       lruib(nelty),          &      ! List of IB integration rules
       lorde(nelty),          &      ! List of element order
       nnode(-nelty:nelty),   &      ! List of element # of nodes
       nface(nelty),          &      ! List of element # of faces
       needg(nelty),          &      ! List of element # of edges
       lexis(nelty),          &      ! List of existing element
       lbxis(nelty),          &      ! List of existing boundary
       lnuty(nelty),          &      ! Number of element types                  
       jttim,                 &
       ncoun_pos,             &
       nllll,                 &
       kttim
  integer(ip), pointer     :: &
       lllll(:),              &
       ltyp2(:),              &
       llll2(:) 
  character(7)             :: &
       cenal(nelty)                  ! List of element names (lower case)
  character(13)            :: &
       cetop(nelty)                  ! List of element topology name
  character(5)             :: &
       cenam(nelty)                  ! List of element names
  character(15)            :: &
       cepos(nelty)                  ! List of element names (upper case)
  character(150)           :: &
       fil_pdata
  !
  ! Graph
  !
  integer(ip)              :: &
       mepoi,                 &      ! Max. # of elements by node
       mpopo                         ! Max. # of point-point connectivity
  integer(ip), pointer :: &
       nepoi(:),              &      ! # of neighbor elements
       pelpo(:),              &      ! Pointer node/element connectivity 
       lelpo(:),              &      ! List node/element connectivities  
       pelel(:),              &      ! Pointer element/element connectivity 
       lelel(:),              &      ! List element/element connectivities 
!!!2 new arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       peled(:),              &      ! Pointer node/edge connectivity 
       leled(:)                      ! List node/edge connectivities 
  !
  ! ENSIGHT special variables
  !
  integer(ip)              :: &
       npart_pos                     ! # of geometry part
  character(5)             :: &
       varna_pos(4,100)              ! Variables names to postprocess
  character(5)             :: &
       wopow(2)                      ! Variables names to postprocess
  integer(ip)              :: &      
       varnu_pos                     ! Number of variables to postprocess
  character(20)            :: &
       nunam_pos                     ! Postprocess name
  type partt_pos 
     character(50) :: name
     integer(ip)   :: numepart
     integer(ip)   :: npoin
  end type partt_pos
  type(partt_pos)  :: parts_pos(10)
  integer(ip)              :: &
       nppti_ens,             &      ! Ensight postprocess time counter
       nppva_ens,             &      ! Ensight postprocess time counter
       kfl_statu_ens                 ! Ensight postprocess status           
  real(rp)                 :: &
       tipoe_ens(10000)              ! Postprocessing times (Ensight)
  character(30)            :: &      ! Ensight variable type
       ensty_ens(50)
  character(15)            :: &      ! Ensight variable name
       ensva_ens(50)
  character(150)           :: &      ! Ensight variable file name
       ensfi_ens(2,50)

contains

  function intost(integ)
    !-------------------------------------
    !
    !  Convert an integer(ip) to a string
    !
    !-------------------------------------
    implicit none
    integer(ip)   :: integ
    integer(4)    :: integ4
    character(20) :: intost
    character(20) :: intaux

    integ4=int(integ,4)
    write(intaux,*) integ4
    intost=adjustl(intaux)

  end function intost

end module def_kintyp
