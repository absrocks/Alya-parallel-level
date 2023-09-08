module mod_redistribution

  use def_kintyp,          only : ip, rp, r1p, comm_data_par,lg, i1p
  use mod_parall,          only : PAR_INITIALIZE_COMMUNICATION_ARRAY
  use mod_parall,          only : PAR_DEALLOCATE_COMMUNICATION_ARRAY
  use mod_parall,          only : par_memor,PAR_COMM_WORLD
  use def_parall,          only : method_redistribution_par
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  use mod_memory,          only : memory_size
  use mod_communications,  only : PAR_MAX
  use mod_communications,  only : PAR_SEND_RECEIVE 
  use mod_communications,  only : PAR_ALLTOALL
  use mod_communications,  only : PAR_ALLTOALLV
  use mod_communications,  only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications,  only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,  only : PAR_START_NON_BLOCKING_COMM
  use mod_communications,  only : PAR_END_NON_BLOCKING_COMM
  use mod_communications,  only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use def_domain,          only : nnode,memor_dom
  use def_domain
  use def_master
  use mod_messages

  use mod_communications

  implicit none

  private

  interface redistribution_array
     module procedure &
          &           redistribution_array_IP_1, &
          &           redistribution_array_IP_2, &
          &           redistribution_array_RP_1, &
          &           redistribution_array_RP_2, &
          &           redistribution_array_RP_3           
  end interface redistribution_array

  integer(ip)         :: PAR_COMM_RED
  integer(4)          :: PAR_COMM_RED4
  integer(ip)         :: PAR_COMM_RED_RANK
  integer(ip)         :: PAR_COMM_RED_SIZE
  type(comm_data_par) :: commd_nelem      
  type(comm_data_par) :: commd_npoin
  type(comm_data_par) :: commd_nboun     

  type :: hash_chunkdist
     integer(ip)             :: chunk
     integer(ip)             :: nchunk
     integer(ip),  pointer   :: offset(:)
     type(i1p),    pointer   :: perm(:)
  end type hash_chunkdist


#ifndef MPI_OFF
  include 'mpif.h'
#endif

  character(100), PARAMETER :: vacal = "mod_redistribution"
  interface generate_comm_chunkdist
     module procedure generate_comm_chunkdist_a,&
          &           generate_comm_chunkdist_b
  end interface generate_comm_chunkdist
  !
  ! Public functions
  !
  public :: generate_comms
  public :: generate_comm_elem 
  public :: generate_comm_poin 
  public :: generate_comm_chunkdist
  public :: generate_comm_sizes
  public :: redistribution_domain
  public :: redistribution_array
  public :: redistribution_array_RP_3
  public :: redistribution_initialization
  public :: redistribution_lid_chunkdist
  public :: redistribution_generate_hash_chunkdist
  public :: redistribution_deallocate_hash
  public :: hash_chunkdist

contains

  subroutine redistribution_initialization(PAR_COMM)

    integer(ip), intent(in) :: PAR_COMM

    PAR_COMM_RED  = PAR_COMM
    PAR_COMM_RED4 = int(PAR_COMM,4_4)
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_RED,PAR_COMM_RED_RANK,PAR_COMM_RED_SIZE)

  end subroutine redistribution_initialization

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell 
  !> @date    28/05/2018 
  !> @brief   Generate communicators
  !> @details Generate elements, nodes and boundary nodes communicators
  !           for mesh redistribution
  !
  !----------------------------------------------------------------------

  subroutine generate_comms(new_part,comm_elem,comm_poin,comm_boun)

    implicit none

    integer(ip), pointer, intent(in)    :: new_part(:)    
    type(comm_data_par),  intent(inout) :: comm_elem      
    type(comm_data_par),  intent(inout) :: comm_poin 
    type(comm_data_par),  intent(inout) :: comm_boun      

    call generate_comm_elem(new_part,comm_elem)
    call generate_comm_poin(new_part,comm_poin)
    call generate_comm_boun(new_part,comm_boun)

  end subroutine generate_comms

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate elements communicator
  !> @details
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_elem(new_part,comm_elem, nenti_)

    use def_domain,  only : nelem

    implicit none

    integer(ip), pointer, intent(in)    :: new_part(:)   
    type(comm_data_par),  intent(inout) :: comm_elem      
    integer(ip), optional,intent(in)    :: nenti_

    integer(ip)           :: ielem, ipart
    integer(ip)           :: isend, nsendp, ibuff, iperm
    integer(ip), pointer  :: all2send(:)
    integer(ip), pointer  :: nperm(:)
    integer(4),  pointer  :: lnsend(:)
    type(i1p),   pointer  :: send_perm(:)
    integer(ip)           :: nenti
    integer(4)            :: PAR_COMM_RED_SIZE4

    call PAR_DEALLOCATE_COMMUNICATION_ARRAY(comm_elem,par_memor)
    call PAR_INITIALIZE_COMMUNICATION_ARRAY(comm_elem)
    !
    ! Nullify
    !
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)      
    !
    ! Check arguments and allocate pointers
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti) call runend("mod_redistribution: wrong size new_part")           
    !
    ! Evaluate lnsend
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)
    lnsend = 0_4
    do ielem = 1,nenti
       if( new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip ) then
          call runend("mod_redistribution: wrong values partition ranks in new_part")
       else
          lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4
       endif
    enddo
    !
    ! Eval comm_elem sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_elem)      
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'all2send', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0_ip,PAR_COMM_RED_SIZE-1
       if(lnsend(ipart) > 0_4) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do ielem = 1,nenti
       isend = all2send(new_part(ielem))
       nperm(isend) = nperm(isend) + 1_ip
       send_perm(isend) % l( nperm(isend) ) = ielem
    enddo

    iperm=1_ip
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_elem % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,' lnsend ',   vacal,lnsend)
    call memory_deallo(par_memor,' send_perm ',vacal,send_perm)
    call memory_deallo(par_memor,' all2send ', vacal,all2send)
    call memory_deallo(par_memor,' nperm ',    vacal,nperm)

  end subroutine generate_comm_elem

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate nodes communicator
  !> @details nenti_ refers to number of elements to consider (outer loop is
  !           through elements).
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_poin(new_part,comm_poin,nenti_)

    use def_domain,  only : npoin, nnode, mnode
    use def_domain,  only : lnods, nelem, ltype

    implicit none

    integer(ip),         pointer,  intent(in)    :: new_part(:)   
    type(comm_data_par),           intent(inout) :: comm_poin      
    integer(ip),         optional, intent(in)    :: nenti_

    integer(ip)                                  :: ipoin, ielem, inode
    integer(ip)                                  :: isend, ipart, nsendp
    integer(ip)                                  :: melem, ibuff, iperm
    integer(ip)                                  :: nenti, isin
    integer(ip),         pointer                 :: lelems(:,:)
    integer(ip),         pointer                 :: lnelem(:)
    integer(ip),         pointer                 :: all2send(:)
    integer(ip),         pointer                 :: nperm(:)
    integer(4),          pointer                 :: lnsend(:)
    type(i1p),           pointer                 :: send_perm(:)
    integer(4)                                   :: PAR_COMM_RED_SIZE4

    call PAR_DEALLOCATE_COMMUNICATION_ARRAY(comm_poin,par_memor)
    call PAR_INITIALIZE_COMMUNICATION_ARRAY(comm_poin)
    !
    ! Nullify
    !
    nullify(lelems)
    nullify(lnelem)
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)
    !
    ! Check arguments
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti)   call runend("mod_redistribution: wrong size new_part")
    if(nenti /= 0_ip .and. nenti /= nelem) call runend("mod_redistribution: option nenti not encoded")
    !
    ! Generate lnelem and lelems 
    !
    call memory_alloca(par_memor,vacal,' lnelem',lnelem,max(1_ip,npoin))

    lnelem = 0_ip
    melem = 0_ip
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          lnelem(inode) = lnelem(inode) + 1_ip 
          if(lnelem(inode) > melem ) melem = lnelem(inode)
       enddo
    enddo

    lnelem = 0_ip
    call memory_alloca(par_memor,vacal,' lelems',lelems,melem,max(1_ip,npoin))
    do ielem = 1,nenti
       do ipoin = 1,nnode(ltype(ielem))
          inode = lnods(ipoin,ielem)
          isin=0_ip;
          do ibuff=1_ip,lnelem(inode)
             if(new_part(lelems(ibuff,inode))==new_part(ielem)) isin = isin + 1_ip;
          enddo
          if(isin == 0_ip) then
             lnelem(inode) = lnelem(inode) + 1_ip
             lelems(lnelem(inode),inode) = ielem 
          else if(isin > 1_ip) then
             call runend("mod_redistribution: something broken isin > 1")
          endif
       enddo
    enddo
    !
    ! Eval lnsend
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)
    do ipoin = 1,npoin
       do ibuff = 1,lnelem(ipoin)
          ielem = lelems(ibuff,ipoin)
          if(new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip) then
             call runend("mod_redistribution: wrong values partition ranks in new_part")
          else
             lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4 
          endif
       enddo
    enddo
    !
    ! Eval comm_poin sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_poin)
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'all2send', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0,PAR_COMM_RED_SIZE-1
       if( lnsend(ipart) > 0_4 ) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do ipoin = 1,npoin
       do ibuff = 1,lnelem(ipoin)
          ielem = lelems(ibuff,ipoin)
          isend = all2send(new_part(ielem))
          nperm(isend) = nperm(isend) + 1_ip
          send_perm(isend) % l( nperm(isend) ) = ipoin
       enddo
    enddo

    iperm=1
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_poin % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,' lnelem '   ,vacal,lnelem)
    call memory_deallo(par_memor,' lelems '   ,vacal,lelems)
    call memory_deallo(par_memor,' lnsend '   ,vacal,lnsend)
    call memory_deallo(par_memor,' send_perm ',vacal,send_perm)
    call memory_deallo(par_memor,' all2send ' ,vacal,all2send)
    call memory_deallo(par_memor,' nperm '    ,vacal,nperm)

  end subroutine generate_comm_poin

  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    28/05/2018 
  !> @brief   Generate elements communicator
  !> @details
  !
  !----------------------------------------------------------------------

  subroutine generate_comm_boun(new_part,comm_boun,nenti_)

    use def_domain,  only : lelbo, nboun, ltypb
    use def_domain,  only : npoin, nnode, mnode
    use def_domain,  only : lnods, nelem, ltype
    use mod_memory      
    implicit none

    integer(ip), pointer,          intent(in)    :: new_part(:)  
    type(comm_data_par),           intent(out)   :: comm_boun
    integer(ip),         optional, intent(in)    :: nenti_

    integer(ip)                                  :: ipoin, ielem
    integer(ip)                                  :: isend, ipart, nsendp
    integer(ip)                                  :: melem, ibuff, iperm
    integer(ip)                                  :: iboun, nenti
    integer(ip),         pointer                 :: all2send(:)
    integer(ip),         pointer                 :: nperm(:)
    integer(4),          pointer                 :: lnsend(:)
    type(i1p),           pointer                 :: send_perm(:)
    integer(4)                                   :: PAR_COMM_RED_SIZE4

    call PAR_DEALLOCATE_COMMUNICATION_ARRAY(comm_boun,par_memor)
    call PAR_INITIALIZE_COMMUNICATION_ARRAY(comm_boun)
    !
    ! Nullify pointers and communicator 
    !
    nullify(all2send)
    nullify(nperm)
    nullify(lnsend)
    nullify(send_perm)
    !
    ! Check arguments
    !
    nenti=nelem;
    if(present(nenti_)) nenti=nenti_
    if( int(memory_size(new_part),ip) <  nenti)   call runend("mod_redistribution: wrong size new_part")
    if(nenti /= 0_ip .and. nenti /= nelem) call runend("mod_redistribution: option nenti not encoded")
    !
    ! Generate lnelem and lelems
    !
    PAR_COMM_RED_SIZE4 = int(PAR_COMM_RED_SIZE,4)
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,PAR_COMM_RED_SIZE4,lboun=0_4)

    do iboun = 1,nboun
       ielem = lelbo(iboun)
       if(new_part(ielem) > PAR_COMM_RED_SIZE-1 .or. new_part(ielem) < 0_ip) then
          call runend("mod_redistribution: wrong values partition ranks in new_part")
       else
          lnsend(new_part(ielem)) = lnsend(new_part(ielem)) + 1_4 
       endif
    end do
    !
    ! Eval comm_poin sizes (lsend_dim, lsens_size, lrecv_dim ...)
    !
    call generate_comm_sizes(lnsend,nsendp,comm_boun)      
    !
    ! Evaluate send permutation (lsend_perm)
    !
    call memory_alloca(par_memor,'all2send', vacal,all2send, PAR_COMM_RED_SIZE,lboun=0_ip)
    call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
    call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))

    isend = 1_ip
    all2send = 0_ip
    do ipart = 0,PAR_COMM_RED_SIZE-1
       if( lnsend(ipart) > 0_4 ) then
          all2send(ipart) = isend
          call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(ipart),ip)) 
          isend = isend + 1_ip
       endif
    end do

    nperm = 0_ip
    do iboun = 1,nboun
       ielem        = lelbo(iboun)
       isend        = all2send(new_part(ielem))
       nperm(isend) = nperm(isend) + 1_ip
       send_perm(isend) % l( nperm(isend) ) = iboun
    enddo

    iperm=1_ip
    do isend=1,nsendp
       if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
          call runend("mod_redistribution: wrong evaluation of cont_perm")
       endif
       do ibuff=1,nperm(isend)
          comm_boun % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    enddo
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,' lnsend ',   vacal,lnsend)
    call memory_deallo(par_memor,' send_perm ',vacal,send_perm)
    call memory_deallo(par_memor,' all2send ', vacal,all2send)
    call memory_deallo(par_memor,' nperm ',    vacal,nperm)

  end subroutine generate_comm_boun
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell
  !> @date    29/05/2018 
  !> @brief   Generate communicator sizes (private subroutine)
  !> @details Generate communicator size structures from the send 
  !           the send requirements   
  ! 
  !----------------------------------------------------------------------

  subroutine generate_comm_sizes(lnsend,nsendp,comm,PAR_COMM_IN)

    implicit none

    integer(4), pointer,  intent(in)    :: lnsend(:)    
    integer(ip),          intent(out)   :: nsendp 
    type(comm_data_par),  intent(out)   :: comm 
    integer(ip),          intent(in),    optional :: PAR_COMM_IN

    integer(ip)          :: ipart, ineig
    integer(4), pointer  :: lnrecv(:)
    integer(4)           :: istat
    integer(4)           :: PAR_MYCOMM_SIZE4
    integer(4)           :: PAR_MYCOMM_RANK4    
    integer(4)           :: PAR_MYCOMM4    
    !
    ! Allocate pointers
    !
    nullify(lnrecv)
    if(.not. present(PAR_COMM_IN)) then
       PAR_MYCOMM_SIZE4 = int(PAR_COMM_RED_SIZE,4)
       PAR_MYCOMM4      = PAR_COMM_RED4
    else
       PAR_MYCOMM4         = int(PAR_COMM_IN,4) 
       call PAR_COMM_RANK_AND_SIZE(PAR_MYCOMM4,PAR_MYCOMM_RANK4,PAR_MYCOMM_SIZE4)
    endif
    call memory_alloca(par_memor,'lnrecv',vacal,lnrecv,PAR_MYCOMM_SIZE4,lboun=0_4)
    lnrecv = 0_ip
    nsendp = 0_ip 
    ! 
    ! Communication of communication requirements
    !
#ifndef MPI_OFF      
    call MPI_Alltoall(lnsend,1_ip,MPI_INTEGER,lnrecv,1_ip,MPI_INTEGER,PAR_MYCOMM4,istat)
#endif
    ! 
    ! Allocate and fill send structures (size, dime, neights)
    ! 
    comm  % PAR_COMM_WORLD = int(PAR_MYCOMM4,ip)
    comm  % nneig = 0_ip
    do ipart = 0_ip,PAR_MYCOMM_SIZE4-1
       if(lnrecv(ipart) > 0_ip .or. lnsend(ipart) > 0_ip) then
          comm  % nneig = comm  % nneig + 1_ip
       endif
       if( lnsend(ipart) > 0_ip ) nsendp = nsendp + 1_ip
    end do
    nullify(comm % neights,comm % lsend_size,comm % lrecv_size)
    call memory_alloca(par_memor,'neights'   ,vacal,comm % neights,   max(1_ip, comm % nneig   ))
    call memory_alloca(par_memor,'lsend_size',vacal,comm % lsend_size,max(1_ip, comm % nneig+1 ))
    call memory_alloca(par_memor,'lrecv_size',vacal,comm % lrecv_size,max(1_ip, comm % nneig+1 ))

    comm % lsend_dim = 0_ip
    comm % lrecv_dim = 0_ip
    comm % lsend_size(1_ip) = 1_ip
    comm % lrecv_size(1_ip) = 1_ip

    ineig = 1_ip
    do ipart = 0_ip,PAR_MYCOMM_SIZE4-1
       if( lnrecv(ipart) > 0_ip .or. lnsend(ipart) > 0_ip ) then
          comm % neights(ineig)      = ipart
          comm % lsend_dim           = comm % lsend_dim + lnsend(ipart)
          comm % lsend_size(ineig+1) = comm % lsend_size(ineig) + lnsend(ipart)
          comm % lrecv_dim           = comm % lrecv_dim + lnrecv(ipart)
          comm % lrecv_size(ineig+1) = comm % lrecv_size(ineig) + lnrecv(ipart)
          ineig                      = ineig + 1_ip
       endif
    end do

    nullify( comm % lsend_perm )
    call memory_alloca(par_memor,'lsend_perm',vacal,comm % lsend_perm, comm % lsend_dim)
    !
    ! Deallocate pointers
    !
    call memory_deallo(par_memor,' lnrecv ',vacal,lnrecv)

  end subroutine generate_comm_sizes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-30
  !> @brief   ???
  !> @details NENTI = NELEM, NPOIN, NBOUN
  !>          XX(NENTI)
  !>          XX(:,NENTI), XX(NENTI,:)
  !>          XX(:,:,NENTI), XX(:,NENTI,:), XX(NENTI,:,:)
  !>          XX(:,:,:,NENTI), XX(:,:,NENTI,:), XX(:,NENTI,:,:), XX(NENTI,:,:,:)
  !> 
  !-----------------------------------------------------------------------

!!$   subroutine redistribution_array_ip_2(xx,WENTI,posit_opt)
!!$  
!!$     integer(ip),  intent(inout)        :: xx(:,:)
!!$     character(*), intent(in)           :: WENTI
!!$     integer(ip),  intent(in), optional :: posit_opt
!!$     integer(ip)                        :: ndim1,ndim2,posit,ndofn
!!$     integer(ip)                        :: ndim,ldim(10),ii
!!$     integer(ip)                        :: ldim_perm(10)
!!$     !
!!$     ! Find dimensions of XX
!!$     !
!!$     if( associated(xx) ) then
!!$        ndim = size(shape(xx))
!!$        ldim = shape(xx)
!!$     else
!!$        ndim = 0
!!$        ldim = 0
!!$     end if
!!$     call PAR_MAX(ndim)
!!$     call PAR_MAX(ndim,ldim)
!!$     !
!!$     ! Position of entity
!!$     !
!!$     if( present(posit_opt) ) then
!!$        posit = posit_opt
!!$     else
!!$        posit = ndim
!!$     end if
!!$     !
!!$     ! Check dimensions
!!$     !
!!$     if( associated(xx) ) then
!!$        do ii = 1,ndim
!!$           if( ii /= posit .and. size(xx,ii) /= ldim(ii) ) call runend('WRONG DIMENSIONS')
!!$        end do
!!$     end if
!!$
!!$     call redistribution_array(xx,WENTI,posit,ndim,ldim)
!!$     
!!$   end subroutine redistribution_array_ip_2
!!$
!!$   subroutine redistribution_array(xx,WENTI,posit,ndim,ldim)
!!$
!!$     integer(ip),  intent(inout)           :: xx(*)
!!$     character(*), intent(in)              :: WENTI
!!$     integer(ip)                           :: posit
!!$     integer(ip)                           :: ndim
!!$     integer(ip)                           :: ldim(10)
!!$     type(comm_data_par),         pointer  :: COMMD 
!!$     integer(ip)                           :: ndim1,ndim2,posit,ndofn
!!$     integer(ip)                           :: ldim_perm(10),ii
!!$     integer(ip),                 pointer  :: bsend(:,:),brecv(:,:)
!!$     !
!!$     ! => Independent of XX
!!$     !
!!$     nullify(bsend,brecv)
!!$     !
!!$     ! Dimensions
!!$     ! LDIM      = n1,n2,nelem,n3
!!$     ! LDIM_PERM = n1,n2,n3,nelem
!!$     ! NODFN     = n1*n2*n3
!!$     !
!!$     ndofn = 1
!!$     jj    = 0
!!$     do ii = 1,ndim
!!$        if( ii /= posit ) then
!!$           ndofn         = ndofn * ldim(ii)
!!$           jj            = jj + 1
!!$           ldim_perm(jj) = ldim(ii)
!!$        end if
!!$     end do
!!$     ldim_perm(ndim) = ndim
!!$     !
!!$     !
!!$     !
!!$     if(      WENTI == 'NELEM' ) then
!!$        COMMD => commd_nelem
!!$        nenti =  nelem
!!$     else if( WENTI == 'NPOIN' ) then
!!$        COMMD => commd_npoin
!!$        nenti =  npoin
!!$     else if( WENTI == 'NBOUN' ) then
!!$        COMMD => commd_nboun
!!$        nenti =  nboun
!!$     end if
!!$     !
!!$     ! Prepare buffers
!!$     !
!!$     call memory_alloca(par_memor,vacal,'bsend',bsend,ndofn,COMMD % lsend_dim)
!!$     call memory_alloca(par_memor,vacal,'brecv',brecv,ndofn,COMMD % lrecv_dim)
!!$     !
!!$     ! Reshape
!!$     !
!!$     if( posit == ndim ) then
!!$        
!!$     else
!!$     do ii = 1_ip,nenti
!!$        idofn = 0
!!$        do idim = 1,ndim-1
!!$           jdim = ldim_perm(idim) ! Dimension in XX
!!$           do idofi = 1,ldim(idim)
!!$              idofn = idofn + 1       ! Location in BSEND
!!$              !jdofn = 
!!$              bsend(idim,ii) = xx(jj)
!!$           end do
!!$        end do
!!$     end do
!!$     !
!!$     ! Perform redistribution
!!$     !
!!$     call PAR_SEND_RECEIVE_TO_ALL(ndofn,bsend,brecv,COMMD,'SYNCHRONOUS',COMMD % lsend_perm)
!!$     !
!!$     ! Reallocate depends on XX
!!$     !        
!!$     !call memory_deallo(par_memor,vacal,'xx',xx)
!!$     !if(      posit == ndim ) then
!!$     !   call memory_alloca(par_memor,vacal,'xx',xx,ldim(1),COMMD % lrecv_dim)
!!$     !else if( posit == 1    ) then           
!!$     !   call memory_alloca(par_memor,vacal,'xx',xx,COMMD % lrecv_dim,ldim(2))
!!$     !end if
!!$     
!!$   end subroutine redistribution_array_ip
!!$

  subroutine redistribution_domain(lpart,PAR_COMM_IN,VERBOSE)

    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaini
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades

    integer(ip),        intent(in),           pointer     :: lpart(:)
    integer(ip),        intent(in), optional              :: PAR_COMM_IN
    logical(lg),        intent(in), optional              :: VERBOSE
    integer(ip)                                           :: ipoin,jpoin,indx,pnode
    integer(ip)                                           :: inode,ii,icont,ielem,ineig
    integer(ip)                                           :: nenti,kpoin,offset,iboun
    integer(ip)                                           :: ipoin_loc,dom_i,jneig,pnodb,kelem,jelem
    integer(ip)                                           :: PAR_COMM,ifiel,istep,knode,knodb
    integer(4)                                            :: PAR_COMM4
    logical(lg)                                           :: lfoun
    integer(ip),                              pointer     :: bsend(:,:)
    integer(ip),                              pointer     :: irecv(:,:)
    integer(ip),                              pointer     :: irecv1(:)
    real(rp),                                 pointer     :: rrecv(:,:)
    integer(ip),                              allocatable :: permn(:)
    integer(ip),                              allocatable :: perme(:)
    integer(ip),                              pointer     :: lninv_tmp(:)
    real(rp),                                 pointer     :: dumm2(:,:)
    logical(lg)                                           :: if_verbose
    type(hash_t)                                          :: ht
    integer(ip)                                           :: lid

    if_verbose = .true.
    if( present(VERBOSE) ) if_verbose = VERBOSE

    if( if_verbose ) call messages_live('REDISTRIBUTION OF DOMAIN ARRAYS','START SECTION')
    !
    ! Compute communicator
    !
    if( present(PAR_COMM_IN) ) then
       PAR_COMM  = int(PAR_COMM_IN)
    else
       PAR_COMM  = int(PAR_COMM_WORLD)
    end if
    PAR_COMM4 = int(PAR_COMM,4)
    call redistribution_initialization(PAR_COMM)
    !
    ! Nullify
    !
    nullify(bsend)
    nullify(irecv)
    nullify(irecv1)
    nullify(rrecv)
    nullify(lninv_tmp)
    !
    ! Domain dimensions
    !
    knode = memory_size(lnods,1_ip)
    knodb = memory_size(lnodb,1_ip)
    call PAR_MAX(knode,PAR_COMM_IN4=PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    call PAR_MAX(knodb,PAR_COMM_IN4=PAR_COMM_RED4,INCLUDE_ROOT=.true.)
    !
    ! Generate communicator
    !
    nenti = memory_size(lpart)

    if( nenti /= nelem .and. nenti /= 0 ) call runend('WRONG LPART')

    if( if_verbose ) call messages_live('GENERATE ELEMENT COMMUNICATOR')
    call generate_comm_elem(lpart,commd_nelem,nenti)
    if( if_verbose ) call messages_live('GENERATE NODE COMMUNICATOR')
    call generate_comm_poin(lpart,commd_npoin,nenti)
    if( if_verbose ) call messages_live('GENERATE BOUNDARY COMMUNICATOR')
    call generate_comm_boun(lpart,commd_nboun,nenti)     
    !
    ! Renumber LNODS/LNODB/LBOEL/LELBO in fragmented local numbering
    !
    if( if_verbose ) call messages_live('RENUMBER ELEMENT/BOUNDARY CONNECTIVITY')
    allocate( permn(max(1_ip,npoin)) )
    permn = -1_ip
    do ineig = 1,commd_npoin % nneig
       icont = 0_ip
       do ii = commd_npoin % lsend_size(ineig),commd_npoin % lsend_size(ineig+1)-1
          ipoin        = commd_npoin % lsend_perm(ii)
          icont        = icont + 1_ip
          permn(ipoin) = icont
       end do
       jneig = 1
       lfoun = .false.
       do while( jneig <= commd_nelem % nneig )
          if( commd_npoin % neights(ineig) == commd_nelem % neights(jneig) ) then
             lfoun = .true.
             exit
          end if
          jneig = jneig + 1
       end do
       if( lfoun ) then
          do ii = commd_nelem  % lsend_size(jneig),commd_nelem % lsend_size(jneig+1)-1
             ielem = commd_nelem % lsend_perm(ii)
             do inode = 1,nnode(abs(ltype(ielem)))
                ipoin              = lnods(inode,ielem)
                kpoin              = permn(ipoin)
                if( kpoin == -1 ) call runend('WE ARE IN BIG TROUBLES 1')
                lnods(inode,ielem) = kpoin
             end do
          end do
       end if
       if( commd_nboun % nneig /= 0 ) then
          jneig = 1
          lfoun = .false.
          do while( jneig <= commd_nboun % nneig )
             if( commd_npoin % neights(ineig) == commd_nboun % neights(jneig) ) then
                lfoun = .true.
                exit
             end if
             jneig = jneig + 1
          end do
          if( lfoun ) then
             do ii = commd_nboun  % lsend_size(jneig),commd_nboun % lsend_size(jneig+1)-1
                iboun = commd_nboun % lsend_perm(ii)
                pnodb = nnode(abs(ltypb(iboun)))
                do inode = 1,pnodb
                   ipoin              = lnodb(inode,iboun)
                   kpoin              = permn(ipoin)
                   if( kpoin == -1 ) call runend('WE ARE IN BIG TROUBLES 2')
                   lnodb(inode,iboun) = kpoin
                end do
             end do
          end if
       end if
       permn = -1_ip
    end do
    !
    ! Renumber LELBO
    !
    allocate(perme(nelem))
    do ineig = 1,commd_nelem % nneig        
       kelem = 0
       do ii = commd_nelem  % lsend_size(ineig),commd_nelem % lsend_size(ineig+1)-1
          ielem = commd_nelem % lsend_perm(ii)
          kelem = kelem + 1
          perme(ielem) = kelem
       end do
    end do
    do iboun = 1,nboun
       ielem = lelbo(iboun)
       kelem = perme(ielem)
       lelbo(iboun) = kelem
    end do
    deallocate( perme )
    deallocate( permn )
    !
    ! New dimensions
    !
    nelem = commd_nelem % lrecv_dim
    nboun = commd_nboun % lrecv_dim
    !
    ! LNODS
    !
    call redistribution_array(lnods,'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LNODS')
    !
    ! LNODB
    !
    call redistribution_array(lnodb,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LNODB')
    !
    ! LBOEL and LELBO
    !
    call redistribution_array(lboel,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBOEL') 
    call redistribution_array(lelbo,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LELBO')
    !
    ! LNINV_TMP
    !
    call messages_live('REDISTRIBUTE GLOBAL NUMBERING')
    call memory_alloca(par_memor,'IRECV1',vacal,irecv1,commd_npoin % lrecv_dim)
    call PAR_SEND_RECEIVE_TO_ALL(lninv_loc,irecv1,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    !
    ! NPOIN and node permutation: unique node numbering from fragmented geometries
    !
    call memory_alloca(par_memor,'commd_npoin % lrecv_perm',vacal,commd_npoin % lrecv_perm,commd_npoin % lrecv_dim)         
    npoin = 0_ip
    call htaini( ht, commd_npoin % lrecv_dim+1, lidson=.true., AUTOMATIC_SIZE=.true. )
    do ipoin = 1,commd_npoin % lrecv_dim
       call htaadd( ht, irecv1(ipoin), lid) 
       commd_npoin % lrecv_perm(ipoin) = lid
    end do
    npoin = ht % nelem 
    call memory_deallo(par_memor,'IRECV1',vacal,irecv1)
    !
    ! Renumber LNODS
    !
    ielem = 0
    do ineig = 1,commd_nelem % nneig
       offset = commd_npoin % lrecv_size(ineig)-1

       do ii = commd_nelem % lrecv_size(ineig),commd_nelem % lrecv_size(ineig+1)-1
          ielem = ielem + 1
          inode = 0
          loop_inode: do while( inode < knode )
             inode              = inode + 1
             ipoin              = lnods(inode,ielem)
             if( ipoin == 0 ) exit loop_inode
             lnods(inode,ielem) = commd_npoin % lrecv_perm(ipoin+offset)
          end do loop_inode
       end do

    end do
    !
    ! Renumber LNODB
    !
    iboun = 0
    do ineig = 1,commd_nboun % nneig

       dom_i = commd_nboun % neights(ineig)
       jneig = 1
       do while( commd_npoin % neights(jneig) /= dom_i )
          jneig = jneig + 1
       end do
       offset = commd_npoin % lrecv_size(jneig)-1

       do ii = commd_nboun % lrecv_size(ineig),commd_nboun % lrecv_size(ineig+1)-1
          iboun = iboun + 1
          inode = 0
          loop_inodb: do while( inode < knodb )
             inode              = inode + 1
             ipoin              = lnodb(inode,iboun)
             if( ipoin == 0 ) exit loop_inodb
             lnodb(inode,iboun) = commd_npoin % lrecv_perm(ipoin+offset)
          end do loop_inodb
       end do

    end do
    !
    ! Deallocate hash table
    !
    call htades( ht )

    !-------------------------------------------------------------------
    !
    ! NELEM arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(leinv_loc,'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LEINV_LOC')
    call redistribution_array(ltype,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LTYPE')
    call redistribution_array(lelch,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LELCH')
    call redistribution_array(lesub,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LESUB')
    call redistribution_array(lmate,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LMATE')
    call redistribution_array(leset,    'NELEM',MEMOR=memor_dom,VARIABLE_NAME='LESET')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NELEM_TYPE ) then
          call redistribution_array_RP_3(xfiel(ifiel) % a,'NELEM',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='FIELD')
       end if
    end do

    !-------------------------------------------------------------------
    !
    ! NPOIN arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(coord,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='COORD')
    call redistribution_array(lninv_loc,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNINV_LOC')
    call redistribution_array(lnoch,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNOCH')
    call redistribution_array(lmast,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LMAST')
    call redistribution_array(kfl_codno,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODNO')
    call redistribution_array(lnset,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNSET')
    call redistribution_array(lgrou_dom,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LGROU_DOM')

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
          call redistribution_array_RP_3(xfiel(ifiel) % a,'NPOIN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='FIELD')
       end if
    end do

    !-------------------------------------------------------------------
    !
    ! NBOUN arrays
    !
    !-------------------------------------------------------------------

    call redistribution_array(ltypb,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LTYPB')
    call redistribution_array(lboch,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBOCH')
    call redistribution_array(lbset,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBSET')
    call redistribution_array(kfl_codbo,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODBO')
    call redistribution_array(lbinv_loc,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBINV_LOC')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
          call redistribution_array_RP_3(xfiel(ifiel) % a,'NBOUN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='FIELD')
       end if
    end do
    !
    ! Renumber LELBO
    !
    iboun = 0
    do ineig = 1,commd_nboun % nneig
       offset = 0_ip
       jneig_loop: do jneig = 1,commd_nelem % nneig
          if( commd_nboun % neights(ineig) == commd_nelem % neights(jneig) ) then
             offset = commd_nelem % lrecv_size(jneig)-1
             exit jneig_loop
          endif
       end do jneig_loop
       do ii = commd_nboun % lrecv_size(ineig),commd_nboun % lrecv_size(ineig+1)-1
          iboun = iboun + 1
          kelem = lelbo(iboun)
          lelbo(iboun) = kelem + offset
       end do
    end do
    if(iboun/=nboun) call runend("Something wrong here")

    if( if_verbose ) call messages_live('REDISTRIBUTION OF DOMAIN ARRAYS','END SECTION')

  end subroutine redistribution_domain

  subroutine redistribution_array_IP_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    integer(ip),                   pointer,  intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV
    integer(ip),         optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    integer(ip),                   pointer                 :: irecv(:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ipoin,jpoin,nsize
    integer(ip)                                            :: ienti,jenti,nenti
    logical(lg)                                            :: if_permute_recv

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(irecv)
    !
    ! Check if array should be redistributed
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))
    else
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.)
    end if

    if( nsize == 0 ) return     

    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if

    if( present(COMM) ) then

       nenti = COMM % lrecv_dim
       call memory_alloca(memor_loc,'IRECV',vacal,irecv,COMM % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   OUTPUT_ARRAY(jenti) = irecv(ienti)
                end do
             else
                OUTPUT_ARRAY(1:nenti) = irecv(1:nenti)
             end if
          end if
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   xx(jenti) = irecv(ienti)
                end do
             else
                xx(1:nenti) = irecv(1:nenti)
             end if
          end if
       end if

    else if(      trim(wtype) == 'NELEM' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_nelem % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nelem)
          if( commd_nelem % lrecv_dim > 0 ) OUTPUT_ARRAY(:) = irecv(:)
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nelem)
          if( commd_nelem % lrecv_dim > 0 ) xx(:) = irecv(:)
       end if

    else if( trim(wtype) == 'NPOIN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_npoin % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             OUTPUT_ARRAY(jpoin) = irecv(ipoin)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             xx(jpoin) = irecv(ipoin)
          end do
       end if

    else if( trim(wtype) == 'NBOUN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_nboun % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nboun)
          if( commd_nboun % lrecv_dim > 0 ) OUTPUT_ARRAY(:) = irecv(:)
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nboun)
          if( commd_nboun % lrecv_dim > 0 ) xx(:) = irecv(:)
       end if

    end if

    call memory_deallo(memor_loc,'IRECV',vacal,irecv)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_IP_1

  subroutine redistribution_array_IP_2(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,MERGE_RECV,OUTPUT_ARRAY)

    integer(ip),                   pointer,  intent(inout) :: xx(:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    logical(lg),         optional,           intent(in)    :: MERGE_RECV     
    integer(ip),         optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:)
    integer(ip),  pointer                                  :: irecv(:,:)
    integer(ip),  pointer                                  :: icont(:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,ipoin,jpoin,nsize
    integer(ip)                                            :: nenti,ienti,jenti,idime
    integer(ip)                                            :: jdime
    logical(lg)                                            :: if_permute_recv
    logical(lg)                                            :: if_merge_recv

    nullify(irecv)
    nullify(icont)

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if

    if_permute_recv = .true.
    if_merge_recv   = .false.

    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    if( present(MERGE_RECV) )   if_merge_recv   = MERGE_RECV


    nsize = memory_size(xx)
    ndim1 = 0
    if( present(posit) ) then
       if(      posit == 1 ) then
          ndim1 = memory_size(xx,2_ip)
       else if( posit == 2 ) then
          ndim1 = memory_size(xx,1_ip)
       else
          call runend('REDISTRIBUITION: WRONG POSITION')
       end if
    else
       ndim1 = memory_size(xx,1_ip)
    end if
    if( present(COMM) ) then
       call PAR_MAX(ndim1,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))     
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))
    else
       call PAR_MAX(ndim1,INCLUDE_ROOT=.true.)     
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    if( present(COMM) ) then
       if( COMM % lrecv_dim        == 0 .and. COMM % lsend_dim        == 0 ) return
    else
       if( commd_nelem % lrecv_dim == 0 .and. commd_nelem % lsend_dim == 0 ) return
    end if

    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if


    if( present(COMM) ) then

       nenti = COMM % lrecv_dim
       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,COMM % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti)
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nenti)
       end if
       if( COMM % lrecv_dim > 0 ) then
          if( associated(COMM % lrecv_perm) ) then
             if( if_permute_recv ) then
                !
                ! Permute
                !
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,jenti) = irecv(1:ndim1,ienti)
                   else
                      xx(1:ndim1,jenti) = irecv(1:ndim1,ienti)
                   end if
                end do
             else if( if_merge_recv .and. if_permute_recv ) then
                !
                ! Permute and merge
                !
                call memory_alloca(memor_loc,'ICONT',vacal,icont,nenti)
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   do idime = 1,ndim1                       
                      if( irecv(idime,ienti) /= 0 ) then
                         icont(jenti)    = icont(jenti) + 1
                         jdime           = icont(jenti)
                         if( present(OUTPUT_ARRAY) ) then
                            OUTPUT_ARRAY(jdime,jenti) = irecv(idime,ienti)
                         else
                            xx(jdime,jenti) = irecv(idime,ienti)
                         end if
                      end if
                   end do
                end do
                call memory_deallo(memor_loc,'ICONT',vacal,icont)                 
             else
                if( present(OUTPUT_ARRAY) ) then
                   do ienti = 1,nenti
                      OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                   end do
                else
                   do ienti = 1,nenti
                      xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                   end do
                end if
             end if
          else
             if( present(OUTPUT_ARRAY) ) then
                do ienti = 1,nenti
                   OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                end do
             else
                do ienti = 1,nenti
                   xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                end do
             end if
          end if
       end if

    else if( trim(wtype) == 'NELEM' ) then
       
       if( commd_nelem % lrecv_dim > 0 .or. commd_nelem % lsend_dim > 0 ) then 
          call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_nelem % lrecv_dim)
          call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
          if( present(OUTPUT_ARRAY) ) then
             call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nelem)
             do ienti = 1,commd_nelem % lrecv_dim 
                OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
             end do
          else
             call memory_deallo(memor_loc,'XX',vacal,xx)
             call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nelem)
             do ienti = 1,commd_nelem % lrecv_dim
                xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
             end do
          end if
       end if

    else if( trim(wtype) == 'NPOIN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_npoin % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             OUTPUT_ARRAY(1:ndim1,jpoin) = irecv(1:ndim1,ipoin)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             xx(1:ndim1,jpoin) = irecv(1:ndim1,ipoin)
          end do
       end if

    else if( trim(wtype) == 'NBOUN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_nboun % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       end if

    end if

    call memory_deallo(memor_loc,'IRECV',vacal,irecv)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_IP_2

  subroutine redistribution_array_RP_1(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    real(rp),                      pointer,  intent(inout) :: xx(:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:)
    real(rp),                      pointer                 :: irecv(:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ipoin,jpoin,nsize
    integer(ip)                                            :: ienti,jenti,nenti
    logical(lg)                                            :: if_permute_recv

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(irecv)
    !
    ! Check if array should be redistributed
    !
    nsize = memory_size(xx)
    if( present(COMM) ) then
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))
    else
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.)
    end if

    if( nsize == 0 ) return     

    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if

    if( present(COMM) ) then

       nenti = COMM % lrecv_dim
       call memory_alloca(memor_loc,'IRECV',vacal,irecv,COMM % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   OUTPUT_ARRAY(jenti) = irecv(ienti)
                end do
             else
                do ienti = 1,COMM % lrecv_dim
                   OUTPUT_ARRAY(ienti) = irecv(ienti)
                end do
             end if
          end if
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   xx(jenti) = irecv(ienti)
                end do
             else
                do ienti = 1,nenti
                   xx(ienti) = irecv(ienti)
                end do
             end if
          end if
       end if

    else if(      trim(wtype) == 'NELEM' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_nelem % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nelem)
          do ienti = 1,commd_nelem % lrecv_dim 
             OUTPUT_ARRAY(ienti) = irecv(ienti)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nelem)
          do ienti = 1,commd_nelem % lrecv_dim 
             xx(ienti) = irecv(ienti)
          end do
       end if

    else if( trim(wtype) == 'NPOIN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_npoin % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             OUTPUT_ARRAY(jpoin) = irecv(ipoin)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             xx(jpoin) = irecv(ipoin)
          end do
       end if

    else if( trim(wtype) == 'NBOUN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,commd_nboun % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             OUTPUT_ARRAY(ienti) = irecv(ienti)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             xx(ienti) = irecv(ienti)
          end do
      end if

    end if

    call memory_deallo(memor_loc,'IRECV',vacal,irecv)

    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_1

  subroutine redistribution_array_RP_2(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    real(rp),                      pointer,  intent(inout) :: xx(:,:)
    character(*),                            intent(in)    :: wtype
    integer(ip),         optional,           intent(in)    :: posit
    integer(8),          optional,           intent(inout) :: memor(2)
    character(*),        optional,           intent(in)    :: variable_name
    type(comm_data_par), optional,           intent(in)    :: COMM
    logical(lg),         optional,           intent(in)    :: PERMUTE_RECV     
    real(rp),            optional, pointer,  intent(inout) :: OUTPUT_ARRAY(:,:)
    real(rp),                      pointer                 :: irecv(:,:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: ndim1,ipoin,jpoin,nsize
    integer(ip)                                            :: ienti,jenti,nenti
    logical(lg)                                            :: if_permute_recv

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(irecv)
    !
    ! Check dimensions
    !
    nsize = memory_size(xx)
    ndim1 = 0
    if( present(posit) ) then
       if(      posit == 1 ) then
          ndim1 = memory_size(xx,2_ip)
       else if( posit == 2 ) then
          ndim1 = memory_size(xx,1_ip)
       else
          call runend('REDISTRIBUITION: WRONG POSITION')
       end if
    else
       ndim1 = memory_size(xx,1_ip)
    end if
    if( present(COMM) ) then
       call PAR_MAX(ndim1,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))     
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))
    else
       call PAR_MAX(ndim1,INCLUDE_ROOT=.true.)     
       call PAR_MAX(nsize,INCLUDE_ROOT=.true.)
    end if
    if( nsize == 0 ) return
    !
    ! Redistribute and permute if required
    !
    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if

    if( present(COMM) ) then

       nenti = COMM % lrecv_dim
       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,COMM % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,COMM,trim(method_redistribution_par),PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   OUTPUT_ARRAY(1:ndim1,jenti) = irecv(1:ndim1,ienti)
                end do
             else
                do ienti = 1,COMM % lrecv_dim
                   OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                end do
             end if
          end if
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nenti)
          if( COMM % lrecv_dim > 0 ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   xx(1:ndim1,jenti) = irecv(1:ndim1,ienti)
                end do
             else
                do ienti = 1,COMM % lrecv_dim
                   xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
                end do
             end if
          end if
       end if

    else if( trim(wtype) == 'NELEM' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_nelem % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nelem)
          do ienti = 1,commd_nelem % lrecv_dim 
             OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nelem)
          do ienti = 1,commd_nelem % lrecv_dim 
             xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       end if

    else if( trim(wtype) == 'NPOIN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_npoin % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             OUTPUT_ARRAY(1:ndim1,jpoin) = irecv(1:ndim1,ipoin)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,npoin)
          do ipoin = 1,commd_npoin % lrecv_dim
             jpoin = commd_npoin % lrecv_perm(ipoin)
             xx(1:ndim1,jpoin) = irecv(1:ndim1,ipoin)
          end do
       end if

    else if( trim(wtype) == 'NBOUN' ) then

       call memory_alloca(memor_loc,'IRECV',vacal,irecv,ndim1,commd_nboun % lrecv_dim)
       call PAR_SEND_RECEIVE_TO_ALL(xx,irecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
       if( present(OUTPUT_ARRAY) ) then
          call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
          call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             OUTPUT_ARRAY(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       else
          call memory_deallo(memor_loc,'XX',vacal,xx)
          call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nboun)
          do ienti = 1,commd_nboun % lrecv_dim
             xx(1:ndim1,ienti) = irecv(1:ndim1,ienti)
          end do
       end if

    else

       call runend('NOT POSSIBLE')

    end if

    call memory_deallo(memor_loc,'IRECV',vacal,irecv)     
    if( present(memor) ) then
       memor     = memor_loc
    else
       par_memor = memor_loc 
    end if

  end subroutine redistribution_array_RP_2

  subroutine redistribution_array_RP_3(xx,wtype,posit,memor,variable_name,COMM,PERMUTE_RECV,OUTPUT_ARRAY)

    real(rp),                      pointer,   intent(inout) :: xx(:,:,:)
    character(*),                             intent(in)    :: wtype
    integer(ip),         optional,            intent(in)    :: posit
    integer(8),          optional,            intent(inout) :: memor(2)
    character(*),        optional,            intent(in)    :: variable_name
    type(comm_data_par), optional,            intent(in)    :: COMM
    logical(lg),         optional,            intent(in)    :: PERMUTE_RECV
    real(rp),            optional, pointer,   intent(inout) :: OUTPUT_ARRAY(:,:,:)
    integer(8)                                              :: memor_loc(2)
    integer(ip)                                             :: ndim1,ndim2,ipoin,jpoin,nenti,idim1,idim2
    integer(ip)                                             :: ndime_recv,ndime_send,nsize,ienti,jenti
    integer(ip)                                             :: my_posit
    real(rp),                       pointer                 :: xsend(:,:,:)
    real(rp),                       pointer                 :: xrecv(:,:,:)
    logical(lg)                                             :: if_permute_recv

    if( present(memor) ) then
       memor_loc = memor
    else
       memor_loc = par_memor
    end if
    if_permute_recv = .true.
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV

    nullify(xsend)
    nullify(xrecv)

    ndim1      = 0
    ndim2      = 0
    ndime_recv = memory_size(xx)
    ndime_send = 0

    if( present(posit) ) then

       my_posit = posit

       if(      posit == 1 ) then
          ndime_send = memory_size(xx,1_ip)
          ndim1      = memory_size(xx,2_ip)
          ndim2      = memory_size(xx,3_ip)
       else if( posit == 2 ) then
          ndime_send = memory_size(xx,2_ip)
          ndim1      = memory_size(xx,1_ip)
          ndim2      = memory_size(xx,3_ip)
       else if( posit == 3 ) then
          ndime_send = memory_size(xx,3_ip)
          ndim1      = memory_size(xx,1_ip)
          ndim2      = memory_size(xx,2_ip)
       else
          call runend('REDISTRIBUITION: WRONG POSITION')
       end if
    else
       my_posit   = 3
       ndime_send = memory_size(xx,3_ip)
       ndim1      = memory_size(xx,1_ip)
       ndim2      = memory_size(xx,2_ip)
    end if
    if( present(COMM) ) then
       call PAR_MAX(ndim1,     INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))     
       call PAR_MAX(ndim2,     INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))     
       call PAR_MAX(ndime_recv,INCLUDE_ROOT=.true.,PAR_COMM_IN4=int(COMM % PAR_COMM_WORLD,4))
    else
       call PAR_MAX(ndim1,     INCLUDE_ROOT=.true.)     
       call PAR_MAX(ndim2,     INCLUDE_ROOT=.true.)     
       call PAR_MAX(ndime_recv,INCLUDE_ROOT=.true.)
    end if
    if( ndime_recv == 0 ) return

    if( present(variable_name) ) then
       call messages_live('REDISTRIBUTE '//trim(variable_name))
    end if

    if( present(COMM) ) then
       ndime_recv = COMM % lrecv_dim
       nenti      = COMM % lrecv_dim
    else if( trim(wtype) == 'NELEM' ) then
       ndime_recv = commd_nelem % lrecv_dim
       nenti      = nelem
    else if( trim(wtype) == 'NPOIN' ) then
       ndime_recv = commd_npoin % lrecv_dim
       nenti      = npoin
    else if( trim(wtype) == 'NBOUN' ) then
       ndime_recv = commd_nboun % lrecv_dim
       nenti      = nboun
    end if

    call memory_alloca(memor_loc,'XSEND',vacal,xsend,ndim1,ndim2,ndime_send)
    call memory_alloca(memor_loc,'XRECV',vacal,xrecv,ndim1,ndim2,ndime_recv)

    if(      my_posit == 1 ) then
       do ienti = 1,ndime_send
          xsend(1:ndim1,1:ndim2,ienti) = xx(ienti,1:ndim1,1:ndim2)
       end do
    else if( my_posit == 2 ) then
       do ienti = 1,ndime_send
          xsend(1:ndim1,1:ndim2,ienti) = xx(1:ndim1,ienti,1:ndim2)
       end do
    else if( my_posit == 3 ) then
       do ienti = 1,ndime_send
          xsend(1:ndim1,1:ndim2,ienti) = xx(1:ndim1,1:ndim2,ienti)
       end do
    end if

    if( present(COMM) ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,COMM,       'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NELEM' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nelem,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NPOIN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_npoin,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    else if( trim(wtype) == 'NBOUN' ) then
       call PAR_SEND_RECEIVE_TO_ALL(xsend,xrecv,commd_nboun,'SYNCHRONOUS',PERMUTE_SEND=.true.)
    end if

    if( present(OUTPUT_ARRAY) ) then
       call memory_deallo(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY)
    else
       call memory_deallo(memor_loc,'XX',vacal,xx)
    end if

    if( ndime_recv > 0 ) then

       if(      my_posit == 1 ) then

          if( present(OUTPUT_ARRAY) ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,nenti,ndim1,ndim2)
          else
             call memory_alloca(memor_loc,'XX',vacal,xx,nenti,ndim1,ndim2)
          end if

          if( present(COMM) ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ipoin = 1,COMM % lrecv_dim
                   jpoin = COMM % lrecv_perm(ipoin)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(jpoin,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   else
                      xx(jpoin,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   end if
                end do
             else
                if( present(OUTPUT_ARRAY) ) then
                   OUTPUT_ARRAY(1:ndime_recv,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,1:ndime_recv)
                else
                   xx(1:ndime_recv,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,1:ndime_recv)
                end if
             end if
          else
             if( trim(wtype) == 'NPOIN' ) then
                do ipoin = 1,commd_npoin % lrecv_dim
                   jpoin = commd_npoin % lrecv_perm(ipoin)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(jpoin,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   else
                      xx(jpoin,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   end if
                end do
             else
                if( present(OUTPUT_ARRAY) ) then
                   OUTPUT_ARRAY(1:ndime_recv,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,1:ndime_recv)
                else
                   xx(1:ndime_recv,1:ndim1,1:ndim2) = xrecv(1:ndim1,1:ndim2,1:ndime_recv)
                end if
             end if
          end if

       else if( my_posit == 2 ) then

          if( present(OUTPUT_ARRAY) ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,nenti,ndim2)
          else
             call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,nenti,ndim2)
          end if

          if( present(COMM) ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ienti = 1,COMM % lrecv_dim
                   jenti = COMM % lrecv_perm(ienti)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,jenti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   else
                      xx(1:ndim1,jenti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   end if
                end do
             else
                do ienti = 1,ndime_recv
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,ienti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   else
                      xx(1:ndim1,ienti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   end if
                end do
             end if
          else
             if( trim(wtype) == 'NPOIN' ) then
                do ipoin = 1,commd_npoin % lrecv_dim
                   jpoin = commd_npoin % lrecv_perm(ipoin)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,jpoin,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   else
                      xx(1:ndim1,jpoin,1:ndim2) = xrecv(1:ndim1,1:ndim2,ipoin)
                   end if
                end do
             else
                do ienti = 1,ndime_recv
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,ienti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   else
                      xx(1:ndim1,ienti,1:ndim2) = xrecv(1:ndim1,1:ndim2,ienti)
                   end if
                end do
             end if
          end if

       else if( my_posit == 3 ) then

          if( present(OUTPUT_ARRAY) ) then
             call memory_alloca(memor_loc,'OUTPUT_ARRAY',vacal,OUTPUT_ARRAY,ndim1,ndim2,nenti)
          else
             call memory_alloca(memor_loc,'XX',vacal,xx,ndim1,ndim2,nenti)
          end if

          if( present(COMM) ) then
             if( associated(COMM % lrecv_perm) .and. if_permute_recv ) then
                do ipoin = 1,COMM % lrecv_dim
                   jpoin = COMM % lrecv_perm(ipoin)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,1:ndim2,jpoin) = xrecv(1:ndim1,1:ndim2,ipoin)
                   else
                      xx(1:ndim1,1:ndim2,jpoin) = xrecv(1:ndim1,1:ndim2,ipoin)
                   end if
                end do
             else
                if( present(OUTPUT_ARRAY) ) then
                   do ienti = 1,ndime_recv
                      OUTPUT_ARRAY(1:ndim1,1:ndim2,ienti) = xrecv(1:ndim1,1:ndim2,ienti)
                   end do
                else
                   do ienti = 1,ndime_recv
                      xx(1:ndim1,1:ndim2,ienti) = xrecv(1:ndim1,1:ndim2,ienti)
                   end do
                end if
             end if
          else
             if( trim(wtype) == 'NPOIN' ) then
                do ipoin = 1,commd_npoin % lrecv_dim
                   jpoin = commd_npoin % lrecv_perm(ipoin)
                   if( present(OUTPUT_ARRAY) ) then
                      OUTPUT_ARRAY(1:ndim1,1:ndim2,jpoin) = xrecv(1:ndim1,1:ndim2,ipoin)
                   else
                      xx(1:ndim1,1:ndim2,jpoin) = xrecv(1:ndim1,1:ndim2,ipoin)
                   end if
                end do
             else
                if( present(OUTPUT_ARRAY) ) then
                   do ienti = 1,ndime_recv
                      OUTPUT_ARRAY(1:ndim1,1:ndim2,ienti) = xrecv(1:ndim1,1:ndim2,ienti)
                   end do
                else
                   do ienti = 1,ndime_recv
                      xx(1:ndim1,1:ndim2,ienti) = xrecv(1:ndim1,1:ndim2,ienti)
                   end do
                end if
             end if
          end if

       end if
    end if

    call memory_deallo(memor_loc,'XSEND',vacal,xsend)
    call memory_deallo(memor_loc,'XRECV',vacal,xrecv)


  end subroutine redistribution_array_RP_3

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    8/11/2018 
  !> @brief   Generate comm_data_par to transfer requests from/to an 
  !>          array distributed in chunks
  !> @details
  !>
  !>        Inputs:
  !>            - reqs:     global indices of the requested components
  !>            - chunk:    chunk size
  !>            - dist_dom: SOURCE (defalult) / TARGET
  !>
  !>        Note: when the distributed domain (dist_dom) is the  TARGET,    
  !>              the requests (reqs) are in fact the inputs that will 
  !>              be tranfered to the target domain partitioned into 
  !>              chunks of size "chunk"
  !>
  !>        Outputs:
  !>             - comm:   comm_data_par filled in
  !>
  !----------------------------------------------------------------------------

  subroutine generate_comm_chunkdist_a(reqs,chunk,PAR_COMM4,comm,dist_dom)

    implicit none

    integer(ip),         pointer,  intent(in)    :: reqs(:)   
    integer(ip),                   intent(in)    :: chunk
    integer(4),                    intent(in)    :: PAR_COMM4
    type(comm_data_par),           intent(inout) :: comm
    character(*),        optional, intent(in)    :: dist_dom
    integer(ip)                                  :: nreqs
    integer(ip)                                  :: reqs_tmp(1)

    nreqs = memory_size(reqs)
    if( nreqs == 0 ) then
       call generate_comm_chunkdist_b(reqs_tmp,nreqs,chunk,par_comm4,comm,dist_dom)
    else
       call generate_comm_chunkdist_b(reqs,nreqs,chunk,par_comm4,comm,dist_dom)
    end if
    
  end subroutine generate_comm_chunkdist_a

  subroutine generate_comm_chunkdist_b(reqs,nreqs,chunk,par_comm4,comm,dist_dom)

    implicit none

    integer(ip),                   intent(in)    :: reqs(*)   
    integer(ip),                   intent(in)    :: nreqs   
    integer(ip),                   intent(in)    :: chunk
    integer(4),                    intent(in)    :: PAR_COMM4
    type(comm_data_par),           intent(inout) :: comm
    character(*),        optional, intent(in)    :: dist_dom

    logical(lg)                                  :: dist_src
    integer(ip)                                  :: nrank, rank, nneig    
    integer(ip)                                  :: nsend, nrecv,ineig, isize  
    integer(ip)                                  :: ireq,irank,idloc
    integer(ip)                                  :: ibuff, iperm
    integer(ip),                   pointer       :: dummi(:)
    integer(ip),                   pointer       :: lnreq_send(:)
    integer(ip),                   pointer       :: lnreq_recv(:)
    type(i1p),                     pointer       :: buff_send(:)
    type(i1p),                     pointer       :: buff_recv(:)
    integer(ip),                   pointer       :: libuf(:)
    integer(ip),                   pointer       :: rank2neig(:)
    integer(ip),                   pointer       :: lsend_perm(:)
    integer(ip),                   pointer       :: lrecv_perm(:)
    integer(ip),                   pointer       :: lrecv_size(:)
    integer(4)                                   :: istat
    integer(ip)                                  :: mpi_sumsend
    integer(ip)                                  :: mpi_sumrecv
    integer(ip),                   pointer       :: mpi_sendbuf(:)
    integer(ip),                   pointer       :: mpi_recvbuf(:)
    integer(4),                    pointer       :: mpi_sendcounts(:)
    integer(4),                    pointer       :: mpi_recvcounts(:)

    character(100), PARAMETER :: vacal = "generate_comm_chunkdist"

    !
    ! Process inputs
    !
    dist_src = .true.
    if(present(dist_dom))then
       if(dist_dom == 'TARGET') then
          dist_src = .false.
       elseif(dist_dom /= 'SOURCE') then
          call runend('generate_comm_chunkdist: wrong dist_dom argument')
       endif
    endif

    !
    ! Eval rank and nrank and allocate memory
    !
    call PAR_COMM_RANK_AND_SIZE(int(PAR_COMM4,ip),rank,nrank)
    call PAR_INITIALIZE_COMMUNICATION_ARRAY(comm)

    nullify(lnreq_send,lnreq_recv,libuf,buff_send,buff_recv,rank2neig,dummi)
    call memory_alloca(par_memor,'lnreq_send',vacal,lnreq_send,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'lnreq_recv',vacal,lnreq_recv,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'buff_send' ,vacal,buff_send ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'buff_recv' ,vacal,buff_recv ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'libuf'     ,vacal,libuf     ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'rank2neig' ,vacal,rank2neig ,nrank,lboun=0_ip)
    call memory_alloca(par_memor,'dummi'     ,vacal,dummi,2_ip,lboun=0_ip)

    !
    !  Generate lists with #requests to send/recv
    !
    do ireq = 1_ip,nreqs
       irank             = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       lnreq_send(irank) = lnreq_send(irank) + 1_ip
    enddo

    !
    ! Send/recv #requests   
    !
    call  PAR_ALLTOALL(lnreq_send,lnreq_recv,PAR_COMM_IN4 = PAR_COMM4)

    !
    ! Send/recv requests
    !

    !1) Allocate buffers to receive and send
    do irank = 0,nrank-1
       call memory_alloca(par_memor,"buff_send % l",vacal,buff_send(irank) % l,max(lnreq_send(irank),1_ip))
       call memory_alloca(par_memor,"buff_recv % l",vacal,buff_recv(irank) % l,max(lnreq_recv(irank),1_ip))
    enddo

    !2) Fill buffers with to send (with receiver local order)
    do ireq = 1_ip,nreqs
       irank                              = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       idloc                              = mod((reqs(ireq)-1_ip),chunk) + 1_ip
       libuf(irank)                       = libuf(irank) + 1_ip
       buff_send(irank) % l(libuf(irank)) = idloc 
    end do

    !3) Send and recv
    ! isize = 2_ip*nrank
    ! call PAR_START_NON_BLOCKING_COMM(1_ip,isize)
    ! call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    nneig = 0_ip
    do irank = 0,nrank-1
       nsend = lnreq_send(irank)
       nrecv = lnreq_recv(irank)       
       !    if(nrecv > 0_ip) then 
       !       call PAR_SEND_RECEIVE(0_ip,nrecv, dummi,buff_recv(irank) % l,&
       !          & dom_i=irank, wsynch='NON BLOCKING', PAR_COMM_IN4=PAR_COMM4)
       !    endif
       if(nsend > 0_ip .or. nrecv > 0_ip) then 
          nneig = nneig+1_ip
          rank2neig(irank) = nneig
       endif
       ! enddo
       ! do irank = 0,nrank-1
       !    nsend = lnreq_send(irank)
       !    if(nsend > 0_ip) then 
       !       call PAR_SEND_RECEIVE(nsend,0_ip, buff_send(irank) % l,dummi,&
       !          & dom_i=irank, wsynch='NON BLOCKING', PAR_COMM_IN4=PAR_COMM4)
       !    endif
    end do
    ! call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Fill mpi arrays
    !
    nullify(mpi_sendbuf,mpi_recvbuf,mpi_sendcounts,mpi_recvcounts)
    call memory_alloca(par_memor,'mpi_sendcounts',vacal,mpi_sendcounts,int(nrank,4),lboun=0_4)
    call memory_alloca(par_memor,'mpi_recvcounts',vacal,mpi_recvcounts,int(nrank,4),lboun=0_4)

    mpi_sumsend = 0
    mpi_sumrecv = 0
    do irank = 0,nrank-1
       mpi_sendcounts(irank) = int(lnreq_send(irank),4)
       mpi_recvcounts(irank) = int(lnreq_recv(irank),4)
       mpi_sumsend           = mpi_sumsend + lnreq_send(irank)
       mpi_sumrecv           = mpi_sumrecv + lnreq_recv(irank)
    end do
    call memory_alloca(par_memor,'mpi_sendbuf',vacal,mpi_sendbuf,max(mpi_sumsend,1_ip),lboun=0_ip)
    call memory_alloca(par_memor,'mpi_recvbuf',vacal,mpi_recvbuf,max(mpi_sumrecv,1_ip),lboun=0_ip)
    mpi_sumsend = 0_4
    do irank = 0,nrank-1
       do ibuff = 1,mpi_sendcounts(irank) 
          mpi_sendbuf(mpi_sumsend + ibuff -1) = buff_send(irank) % l(ibuff)
       enddo
       mpi_sumsend = mpi_sumsend + lnreq_send(irank)
    enddo
    !
    ! Perform mpi communications
    !
    call PAR_ALLTOALLV(mpi_sendbuf,mpi_recvbuf,mpi_sendcounts,mpi_recvcounts,PAR_COMM_IN4=PAR_COMM4)     
    ! 
    ! Deallocate buffers
    !
    mpi_sumrecv = 0_4
    do irank = 0,nrank-1
       do ibuff = 1,  mpi_recvcounts(irank) 
          buff_recv(irank) % l(ibuff) = mpi_recvbuf(mpi_sumrecv + ibuff -1)
       end do
       mpi_sumrecv = mpi_sumrecv + lnreq_recv(irank)
    end do
    call memory_deallo(par_memor,'mpi_sendcounts',vacal,mpi_sendcounts)
    call memory_deallo(par_memor,'mpi_recvcounts',vacal,mpi_recvcounts)
    call memory_deallo(par_memor,'mpi_sendbuf'   ,vacal,mpi_sendbuf)
    call memory_deallo(par_memor,'mpi_recvbuf'   ,vacal,mpi_recvbuf)
    !
    ! Generate comm_data_par 
    !
    call memory_alloca(par_memor,'neights'   ,vacal,comm % neights,   max(1_ip, nneig   ))
    call memory_alloca(par_memor,'lsend_size',vacal,comm % lsend_size,max(1_ip, nneig+1 ))
    call memory_alloca(par_memor,'lrecv_size',vacal,comm % lrecv_size,max(1_ip, nneig+1 ))

    comm % PAR_COMM_WORLD   = PAR_COMM4
    comm % RANK4            = int(rank,4)
    comm % lsend_dim        = 0_ip
    comm % lrecv_dim        = 0_ip
    comm % lsend_size(1_ip) = 1_ip
    comm % lrecv_size(1_ip) = 1_ip
    comm % nneig            = nneig

    ineig = 1_ip
    do irank = 0_ip,nrank-1
       if(dist_src) then
          nrecv = lnreq_send(irank) ! you receive the requests sent
          nsend = lnreq_recv(irank) ! you send the requests received
       else
          nsend = lnreq_send(irank) ! inputs sent
          nrecv = lnreq_recv(irank) ! inputs received
       endif

       if( nsend > 0_ip .or. nrecv > 0_ip ) then
          comm % neights(ineig)      = irank
          comm % lsend_dim           = comm % lsend_dim + nsend
          comm % lsend_size(ineig+1) = comm % lsend_size(ineig) + nsend
          comm % lrecv_dim           = comm % lrecv_dim + nrecv
          comm % lrecv_size(ineig+1) = comm % lrecv_size(ineig) + nrecv
          ineig                      = ineig + 1_ip
       endif
    end do

    call memory_alloca(par_memor,'lsend_perm',vacal,comm % lsend_perm, comm % lsend_dim)
    call memory_alloca(par_memor,'lrecv_perm',vacal,comm % lrecv_perm, comm % lrecv_dim)

    nullify(lsend_perm,lrecv_perm,lrecv_size)

    if(dist_src) then
       lsend_perm => comm % lsend_perm
       lrecv_perm => comm % lrecv_perm
       lrecv_size => comm % lrecv_size
    else
       lsend_perm => comm % lrecv_perm
       lrecv_perm => comm % lsend_perm
       lrecv_size => comm % lsend_size
    endif

    !lsend_perm
    iperm = 1_ip
    do irank = 0_ip,nrank-1
       do ibuff = 1_ip, lnreq_recv(irank)
          lsend_perm(iperm)=buff_recv(irank) % l(ibuff)
          iperm = iperm + 1_ip
       enddo
    end do

    !lrecv_perm
    libuf(:)=0_ip
    do ireq = 1_ip,nreqs
       irank = min((reqs(ireq)-1_ip)/chunk,nrank-1_ip) 
       ineig = rank2neig(irank)
       lrecv_perm(lrecv_size(ineig)+libuf(irank)) = ireq
       libuf(irank) = libuf(irank) + 1_ip 
    enddo

    nullify(lsend_perm,lrecv_perm,lrecv_size) 
    !
    ! Deallocate memory
    !
    call memory_deallo(par_memor,'lnreq_send',vacal,lnreq_send)     
    call memory_deallo(par_memor,'lnreq_recv',vacal,lnreq_recv)     
    call memory_deallo(par_memor,'buff_send', vacal,buff_send)
    call memory_deallo(par_memor,'buff_recv', vacal,buff_recv)
    call memory_deallo(par_memor,'libuf',     vacal,libuf)
    call memory_deallo(par_memor,'rank2neig', vacal,rank2neig)
    call memory_deallo(par_memor,'dummi' ,vacal,dummi)

  end subroutine generate_comm_chunkdist_b


  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Generate a hash_chunkdist in order to perform permuations from
  !>          gloval to local id of the nodes. The global id (gid) is the 
  !>          global identification of the node, the local id (lid) is its
  !>          position in the list of requests (reqs) used in the subroutione: 
  !>          generate_comm_chunkdist and sent here also as an argument.
  !>  
  !> @details
  !>
  !>        Inputs:
  !>            - reqs:     global indices of the requested components
  !>            - chunk:    chunk size
  !>            - comm:     comm_data_par evaluated in generate_comm_chunkdist
  !>
  !>        Outputs:
  !>             - htable: hash_chunkdist for permutation
  !>
  !----------------------------------------------------------------------------
  subroutine redistribution_generate_hash_chunkdist(reqs,chunk,comm,htable)

    implicit none

    integer(ip),          pointer, intent(in)    :: reqs(:)   
    integer(ip),                   intent(in)    :: chunk
    type(comm_data_par),           intent(inout) :: comm
    type(hash_chunkdist),          intent(inout) :: htable

    integer(ip)                                  :: nrank, rank
    integer(ip)                                  :: ineig,nneig,neigh
    integer(ip)                                  :: isize,irecv,nrecv
    integer(ip)                                  :: iperm,ilid
    integer(ip)                                  :: ireq
    integer(ip)                                  :: mingid,maxgid

    character(100), PARAMETER :: vacal = "redistribution_generate_hash_chunkdist"

    call PAR_COMM_RANK_AND_SIZE( comm % PAR_COMM_WORLD,rank,nrank)
    nullify(htable % offset,htable % perm)
    call memory_alloca(par_memor,'htable % offset',vacal,htable % offset, nrank,lboun=0_ip)
    call memory_alloca(par_memor,'htable % perm'  ,vacal,htable % perm,   nrank,lboun=0_ip)

    htable % chunk  = chunk
    htable % nchunk = nrank

    ilid = 1_ip
    do ineig = 1_ip, comm % nneig

       neigh  = comm % neights(ineig)

       mingid = huge(1_ip)
       maxgid = -1_ip
       nrecv  =  comm % lrecv_size(ineig+1_ip)-comm % lrecv_size(ineig)  

       if(nrecv > 0_ip) then

          !
          ! Eval range within each block
          !
          do irecv = comm % lrecv_size(ineig) , comm % lrecv_size(ineig+1_ip) - 1_ip
             ireq = reqs(comm % lrecv_perm(irecv))
             if(ireq < mingid) mingid = ireq
             if(ireq > maxgid) maxgid = ireq
          enddo

          !
          ! Store offset and allocte local permutation array
          !
          htable % offset(neigh) = mingid - 1_ip
          call memory_alloca(par_memor,'hatable % perm % l',vacal,htable % perm(neigh) % l,maxgid-mingid+1_ip)

          !
          ! Eval local permutaiton array
          !
          do irecv = comm % lrecv_size(ineig) , comm % lrecv_size(ineig+1_ip) - 1_ip
             ireq = reqs(comm % lrecv_perm(irecv))
             htable % perm(neigh) % l(ireq-mingid+1) = 1_ip
          enddo
          do iperm = 1_ip,maxgid-mingid+1_ip
             if(htable % perm(neigh) % l(iperm) == 1_ip) then
                htable % perm(neigh) % l(iperm) = ilid
                ilid = ilid + 1_ip 
             endif
          enddo
       endif

    enddo

  end subroutine redistribution_generate_hash_chunkdist

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Deallocate a hash_chunkdist (htable)
  !> @details Deallocate pointers within the type hash_chunkdist
  !>
  !----------------------------------------------------------------------------
  subroutine redistribution_deallocate_hash(htable)

    type(hash_chunkdist), intent(inout) :: htable

    character(100), PARAMETER :: vacal = "redistribution_deallocate_hash"

    call memory_deallo(par_memor,'htable % offset',vacal,htable % offset)
    call memory_deallo(par_memor,'htable % perm'  ,vacal,htable % perm)

  end subroutine redistribution_deallocate_hash

  !----------------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    13/11/2018 
  !> @brief   Permute from global id (gid) to local id (lid) using a  
  !>          hash_chunkdist
  !> @details
  !>
  !>        Inputs:
  !>            - gid: global id of the node
  !>
  !>        Outputs:
  !>             - local id (position of the gid in reqs)
  !>
  !----------------------------------------------------------------------------
  integer(ip) function redistribution_lid_chunkdist(gid,htable)

    integer(ip),                   intent(in)   :: gid 
    type(hash_chunkdist),          intent(in)   :: htable

    integer(ip)                                 :: ichunk

    ichunk = min((gid-1_ip)/htable % chunk,htable % nchunk-1_ip) 
    redistribution_lid_chunkdist = htable % perm(ichunk) % l(gid-htable % offset(ichunk))

  end function redistribution_lid_chunkdist


end module mod_redistribution
