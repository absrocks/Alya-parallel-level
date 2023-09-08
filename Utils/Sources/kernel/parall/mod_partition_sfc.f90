module mod_partition_sfc
 
  use mod_memory,          only : memory_alloca, memory_deallo
  use def_kintyp,          only : ip, rp, r1p, i1p, comm_data_par,lg
  use mod_parall,          only : par_memor
  use mod_communications,  only : PAR_SUM, PAR_MAX,PAR_BARRIER
  use mod_communications,  only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications,  only : PAR_COMM_RANK_AND_SIZE

  implicit none

  integer(ip),PARAMETER   :: DIM_BIN_CORE        = 128_ip    ! Number of boxes per direction per partitioning process  
  integer(ip),PARAMETER   :: DEFAULT_COARSE_BIN  = 1_ip      ! <1> The minimum larger than the number of procs.
  ! <0> The maximum lower than the number of procs. 
  private

#ifndef MPI_OFF
  include 'mpif.h'
#endif

  integer(ip)             :: dboxf(3)        ! #bin boxes per direction in the fine bin
  integer(ip)             :: dboxc(3)        ! #bin boxes per direction in the coarse bin
  integer(ip)             :: dboxl(3)        ! #bin boxes per direction in the local bin
  integer(ip)             :: nboxc           ! #boxes of the coarse bin
  integer(ip)             :: nboxl           ! #boxes of the local bin

  real(rp)                :: min_coord(3)    ! bounding box min
  real(rp)                :: max_coord(3)    ! bounding bos max

  real(rp),    pointer    :: weigc(:)        ! Weight of the boxes in the coarse bin
  real(rp),    pointer    :: weigc_loc(:)    ! Weight of the boxes in the coarse bin within one phase
  real(rp),    pointer    :: weigl(:)        ! Weight of the boxes in the local bin
  integer(ip), pointer    :: partl(:)        ! Partition of the local bin

  type(comm_data_par)     :: comm1           ! First point-to-point communication structure
  type(comm_data_par)     :: comm2           ! Second point-to-point communication structure
  real(rp), pointer       :: send_buff(:)    ! Send commuincation buffer
  real(rp), pointer       :: recv_buff(:)    ! Receive communication buffer

  integer(ip), pointer    :: lboxes(:)       ! local box assigned to each element/node (if negative is in the owned coarse box) !rick:better explain
  real(rp), pointer       :: lcorr(:)        ! dystribution conrrection

  integer(ip), pointer    :: lepar(:)     
  real(rp), pointer       :: lenti(:,:)
  real(rp), pointer       :: lweig(:)
  integer(ip)             :: nenti
  integer(ip)             :: sdim
  integer(ip)             :: npart_sfc
  integer(4)              :: PAR_COMM
  integer(4)              :: PAR_MY_RANK
  integer(4)              :: PAR_MY_SIZE

  integer(ip)             :: boxes_coarse(3)
  integer(ip)             :: boxes_fine(3)

  real(rp)                :: time_ini
  real(rp)                :: time_partition

  integer(ip)             :: nphase   
  integer(ip)             :: iphase   
  integer(ip)             :: iboxc0,iboxc1,iboxc1_
  integer(ip)             :: gbuff0

  real(rp)                :: totalw

  integer(ip), pointer   :: lboxc(:)
  integer(ip), pointer   :: lboxl(:)
  integer(ip), pointer   :: neighs(:)     ! 1: I have elements/points in the corresponding coarse box, 0: none
  integer(4), pointer    :: lnsend(:)     !number of "local" boxes send/received
  integer(4), pointer    :: lnsend_loc(:) 
  integer(4), pointer    :: lnrecv_loc(:) 
  real(rp),    pointer   :: bufwei(:)     ! Local buffer to accumulate weights
  type(r1p),   pointer   :: bufse2(:)     ! Buffer to send "local" boxes info to the partitioning processes
  integer(ip)            :: nneigs
  integer(ip)            :: ineis

  integer(ip), pointer   :: reord(:), invor(:) 
  !
  ! Output statistics
  !
  real(rp)                :: number_touched_bin
  real(rp)                :: percentage_touched_bin
  real(rp)                :: max_weight_bin
  real(rp)                :: ave_weight_bin
  !
  ! Public functions
  !
  public :: partition_sfc
  public :: partition_sfc_statistics

contains
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell & Juan Cajas
  !> @date    29/11/2016
  !> @brief   Partition of the mesh using Space Filling Curves (SFC)
  !> @details
  !
  !----------------------------------------------------------------------
  subroutine partition_sfc(lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,PAR_COMM_)

    implicit none

    integer(ip), pointer,           intent(inout) :: lepar_(:)         ! partition rank of each entity 
    integer(ip),                    intent(inout) :: npart_sfc_            ! number of partitions
    real(rp),     target,           intent(in)    :: lenti_(:,:)       ! coordinates of the entities 
    integer(ip),                    intent(in)    :: nenti_            ! number of entities
    integer(ip),                    intent(in)    :: sdim_             ! space dimension (2D/3D)
    real(rp),     target, optional, intent(in)    :: lweig_(:)         ! weight of each entitiy 
    integer(ip),          optional, intent(in)    :: boxes_coarse_(3)  ! coarse boxed (to distribute bounding box)
    integer(ip),          optional, intent(in)    :: boxes_fine_(3)    ! find boxes (to define SFC)            
    integer(4),           optional, intent(in)    :: PAR_COMM_         ! communicator (empty if sequential)

    call constructor(lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,PAR_COMM_)

    call setup()
    call define_distribution()

    do iphase=1_ip,nphase


       call redistribute_weights()
       call local_partition()
       call redistribute_result()

    enddo

    call destructor()

  end subroutine partition_sfc
  !
  ! Module constructor
  !
  subroutine constructor(lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,PAR_COMM_)

    implicit none

    integer(ip),  pointer,          intent(inout) :: lepar_(:)      ! List of elements partition indices !rick: pq ha de ser in tambe?
    integer(ip),                    intent(inout) :: npart_sfc_         ! #partitions
    real(rp),     target,           intent(in)    :: lenti_(:,:) 
    integer(ip),                    intent(in)    :: nenti_
    integer(ip),                    intent(in)    :: sdim_
    real(rp),     target, optional, intent(in)    :: lweig_(:) 
    integer(ip),          optional, intent(in)    :: boxes_coarse_(3)
    integer(ip),          optional, intent(in)    :: boxes_fine_(3)                
    integer(4),           optional, intent(in)    :: PAR_COMM_
    integer(ip)                                   :: ii, ipart
    logical(lg)                                   :: auxl

    character(100), PARAMETER :: vacal = "constructor"
    !
    ! Take module arguments
    !
    if( .not. associated(lepar_) ) then
       nullify(lepar_);
       call memory_alloca(par_memor,' lepar_ ',vacal,lepar_,nenti_)
    end if
    lepar => lepar_;
    npart_sfc =  npart_sfc_
    lenti => lenti_
    nenti =  nenti_
    sdim =   sdim_
    if(present(lweig_)) then 
       lweig => lweig_
    else
       nullify(lweig)
       call memory_alloca(par_memor,' lweig ',vacal,lweig,max(1_ip,nenti_))
       lweig = 1.0_rp
    endif
    boxes_coarse(1:3)=0_ip
    boxes_fine(1:3)=0_ip
    if(present(boxes_coarse_)) boxes_coarse(1:3)=boxes_coarse_(1:3)
    if(present(boxes_fine_))   boxes_fine(1:3)=boxes_fine_(1:3)
    if(present(PAR_COMM_)) then
       PAR_COMM=PAR_COMM_
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM,PAR_MY_RANK,PAR_MY_SIZE)
    else   
       PAR_MY_SIZE=1_4
    endif
    !
    !  Check boxes arguments
    !
    if(        present(boxes_coarse_)  .and. ( .not. present(boxes_fine_)) .or. &
         (.not. present(boxes_coarse_)) .and.         present(boxes_fine_)) then
       call runend("ERROR mod_partition_sfc: both coarse boxes and fine boxes must be present or absent in the arguments")
    endif
    if( present(boxes_coarse_) .and. present(boxes_fine_) ) then
       if(      (boxes_coarse_(1)/=0_ip .and. boxes_fine_(1)==0_ip) .or. &
            &   (boxes_coarse_(1)==0_ip .and. boxes_fine_(1)/=0_ip) )then
          call runend("ERROR mod_partition_sfc: both coarse boxes and fine boxes must be defined or not defined")
       endif
    end if
    do ii=2,sdim
       if(boxes_coarse(1) /= boxes_coarse(ii)) then
          call runend("ERROR mod_partition_sfc: the boxes grid must be isotropic")
       endif
       if(boxes_fine(1) /= boxes_fine(ii)) then
          call runend("ERROR mod_partition_sfc: the boxes grid must be isotropic")
       endif
    enddo
    !
    ! Get correction coeficients from file
    !
    nullify(lcorr)  !rick: These should be obtained as an argument of the module
    call memory_alloca(par_memor,' lcorr ','kk',lcorr,max(1_ip,npart_sfc))

    inquire( file="rank-elements.dat", exist=auxl )
    if(auxl) then
       open(unit=1405873,file='rank-elements.dat')
       do ii=1, npart_sfc
          read(1405873,*) ipart, lcorr(ii)
       end do
       lcorr(npart_sfc) = real(npart_sfc,rp) - sum(lcorr(1:npart_sfc-1))
    else
       lcorr(1:npart_sfc) = 1.0_rp
    endif

    ineis = 1_ip

    call PAR_BARRIER(PAR_COMM_IN4=PAR_COMM_)
    call cputim(time_ini)

  end subroutine constructor
  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-02-05
  !> @brief   SFC statistics
  !> @details x1 = Ratio of touched bin (give an idea on the 
  !>               anisotropy of the mesh)
  !>          x2 = Max weight
  !>          x3 = Average weight
  !>          x4 = Granularity = max/sum (realtive weight of the 
  !>               heaviest bin)
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_statistics(x1,x2,x3,x4,x5,i1,i2,COMM4)

    use def_master
    real(rp),    intent(out)          :: x1 
    real(rp),    intent(out)          :: x2
    real(rp),    intent(out)          :: x3
    real(rp),    intent(out)          :: x4
    real(rp),    intent(out)          :: x5
    integer(ip), intent(out)          :: i1(:)
    integer(ip), intent(out)          :: i2(:)
    integer(4),  intent(in), optional :: COMM4
    integer(ip)                       :: ii,boxes(2)

    x1          = number_touched_bin 
    x2          = max_weight_bin
    x3          = ave_weight_bin
    x5          = real(nboxc,rp)*real(nboxl,rp)

    if( present(COMM4) ) then
       call PAR_SUM(x1,PAR_COMM_IN4=COMM4)
       call PAR_MAX(x2,PAR_COMM_IN4=COMM4)
       call PAR_SUM(x3,PAR_COMM_IN4=COMM4)
       call PAR_MAX(x5,PAR_COMM_IN4=COMM4)
    else
       call PAR_SUM(x1,PAR_COMM_IN4=PAR_COMM)
       call PAR_MAX(x2,PAR_COMM_IN4=PAR_COMM)
       call PAR_SUM(x3,PAR_COMM_IN4=PAR_COMM)
       call PAR_MAX(x5,PAR_COMM_IN4=PAR_COMM)
    end if

    if( x1 /= 0.0_rp ) then
       x4 = x2 / x3
       x3 = x3 / x1
       x1 = x1 / x5
    else
       x1 = 0.0_rp
       x2 = 0.0_rp
       x3 = 0.0_rp
       x4 = 0.0_rp
       x5 = 0.0_rp
    end if

    boxes(1) = dboxf(1)
    boxes(2) = dboxc(1)
    if( present(COMM4) ) then
       call PAR_MAX(2_ip,boxes,PAR_COMM_IN4=COMM4)
    else
       call PAR_MAX(2_ip,boxes,PAR_COMM_IN4=PAR_COMM)
    end if
    i1(:) = boxes(1)
    i2(:) = boxes(2)
   
  end subroutine partition_sfc_statistics
  !
  ! Define bins
  !
  subroutine setup()

    use mod_maths,   only : maths_sfc_1d_to2d3d_tab

    implicit none
    real(rp), pointer     :: aux_coord(:)
    integer(ip)           :: idime, lev,ii
    real(rp)              :: localw, r1
    integer(ip)           :: iproc,x,y,z,coo1d
    character(100), PARAMETER :: vacal = "setup"
    !
    ! Set number of slaves according to hilbert restriction
    !

    if(DEFAULT_COARSE_BIN==0_ip) then
       lev = int(log(real(PAR_MY_SIZE,rp))/log((2.0_rp)**(sdim)),ip)
    elseif(DEFAULT_COARSE_BIN==1_ip) then
       r1 = log(real(PAR_MY_SIZE,rp))/log((2.0_rp)**(sdim))
       lev = int(r1,ip)
       if(r1 > real(lev,rp)) then
          lev = lev + 1_ip
       endif
    else
       call runend("ERROR mod_partition_sfc: unexpected DEFAULT_COARSE_BIN")
    endif
    dboxc(1:3) = 2_ip**lev
    dboxf(1:3) = DIM_BIN_CORE*dboxc(1:3)
    if( boxes_fine(1) /= 0 .and. boxes_coarse(1) /= 0 ) then
       dboxc(1:3) = boxes_coarse(1)
       dboxf(1:3) = boxes_fine(1)               
    endif
    dboxl(1:3) = dboxf(1:3) / dboxc(1:3)
    !
    ! Bounding box
    !
    nullify(aux_coord)
    allocate(aux_coord(6))
    aux_coord(1:3)=-huge(1.0_rp)*0.1_rp
    aux_coord(4:6)=-huge(1.0_rp)*0.1_rp
    if(nenti>0) then
       do idime = 1_ip,sdim
          aux_coord(idime)        =  maxval(lenti(idime,1:nenti))
          aux_coord(3_ip+idime)   = -minval(lenti(idime,1:nenti))
       end do
    endif
    call PAR_MAX(aux_coord,PAR_COMM_IN4=PAR_COMM)
    !max_coord(1:3) =  aux_coord(1:3)*1.001_rp
    !min_coord(1:3) = -aux_coord(4:6)*1.001_rp
    max_coord(1:3) =  aux_coord(1:3)
    min_coord(1:3) = -aux_coord(4:6)
    max_coord(1:3) = max_coord(1:3) + (max_coord(1:3)-min_coord(1:3))*0.001_rp
    min_coord(1:3) = min_coord(1:3) - (max_coord(1:3)-min_coord(1:3))*0.001_rp
    deallocate(aux_coord)
    !
    ! Eval bin sizes
    !
    nboxc = dboxc(1)**sdim
    nboxl = dboxl(1)**sdim
    !
    ! Phases
    !
    nphase = int(nboxc/PAR_MY_SIZE,ip)
    if(mod(nboxc,int(PAR_MY_SIZE,ip)) > 0_ip) nphase = nphase + 1_ip;
    !
    ! totalw
    !
    totalw=0.0_rp
    do ii = 1_ip,nenti
       totalw=totalw+lweig(ii) 
    enddo
    call PAR_SUM(totalw, PAR_COMM_IN4=PAR_COMM)
    !
    ! Ordering of coarse bins
    !
    gbuff0=0_ip
    nullify(reord,invor)
    call memory_alloca(par_memor,'reord'    ,vacal,reord,nboxc)
    call memory_alloca(par_memor,'invor'    ,vacal,invor,nboxc)
    do iproc=1,nboxc
       if(sdim==2_ip) then
          call maths_sfc_1d_to2d3d_tab(dboxc(1),iproc,x,y)
          coo1d=(y-1)*dboxc(1)+x
       else
          call maths_sfc_1d_to2d3d_tab(dboxc(1),iproc,x,y,z)
          coo1d=(z-1)*dboxc(1)*dboxc(1)+(y-1)*dboxc(1)+x
       endif
       reord(iproc) = coo1d
       invor(coo1d) = iproc
    enddo
    !
    ! SFC statistics
    !
    number_touched_bin     =  0.0_rp
    percentage_touched_bin =  0.0_rp
    max_weight_bin         = -1.0_rp
    ave_weight_bin         =  0.0_rp

  end subroutine setup
  !
  ! Define distribution
  !
  subroutine define_distribution()

    use mod_maths,               only : maths_mapping_coord_to_3d

    implicit none
    integer(ip)           :: bcoof(3),bcooc(3),bcool(3) 
    integer(ip)           :: idime, ielem, gbuff, ibuff
    integer(ip)           :: iboxc, iboxl, ineig, kbuff
    type(i1p),   pointer  :: lenbc(:)
    integer(ip), pointer  :: auxbf(:)

    character(100), PARAMETER :: vacal = "define_distribution"
    !
    ! Allocate pointers and initialize
    !
    nullify(weigc,weigc_loc,lboxes,weigl,lnsend,lnsend_loc,lnrecv_loc)
    nullify(lboxc,lboxl,neighs,bufwei,bufse2,lenbc,auxbf)

    call memory_alloca(par_memor,' weigc_loc ',vacal,weigc_loc, max(1_ip,int(PAR_MY_SIZE,ip)))
    call memory_alloca(par_memor,' weigc ',vacal,weigc, max(1_ip,nboxc))
    call memory_alloca(par_memor,' lboxes ',vacal,lboxes,max(1_ip,nenti))
    call memory_alloca(par_memor,' weigl  ',vacal,weigl, max(1_ip,nboxl))
    call memory_alloca(par_memor,' lnsend ',vacal,lnsend,int(max(1_ip,nboxc),4),lboun=0_4)
    call memory_alloca(par_memor,' lnsend_loc ',vacal,lnsend_loc,int(max(1_ip,int(PAR_MY_SIZE,ip)),4),lboun=0_4)
    call memory_alloca(par_memor,' lnrecv_loc ',vacal,lnrecv_loc,int(max(1_ip,int(PAR_MY_SIZE,ip)),4),lboun=0_4)
    call memory_alloca(par_memor,' neighs ',vacal,neighs,max(1_ip,nboxc))
    call memory_alloca(par_memor,' lenbc  ',vacal,lenbc,max(1_ip,nboxc))
    call memory_alloca(par_memor,' lboxc  ',vacal,lboxc, max(1_ip,nenti))
    call memory_alloca(par_memor,' lboxl  ',vacal,lboxl, max(1_ip,nenti))
    call memory_alloca(par_memor,' bufwei ',vacal,bufwei,max(1_ip,2*nboxl))
    call memory_alloca(par_memor,' auxbf ',vacal,auxbf,max(1_ip,nenti))

    lnsend(:) = 0_4
    neighs(:) = 0_ip
    nneigs    = 0_ip
    weigc(:)  = 0_ip
    bcool     = 1_ip

    do ielem = 1_ip, nenti
       !
       ! Eval. index of the mc box (box coordinates) at the different bins: bcoof,bcooc,bcool
       !
       call maths_mapping_coord_to_3d(sdim,dboxf,min_coord,max_coord,lenti(1:sdim,ielem),bcoof(1),bcoof(2),bcoof(3))
       call maths_mapping_coord_to_3d(sdim,dboxc,min_coord,max_coord,lenti(1:sdim,ielem),bcooc(1),bcooc(2),bcooc(3))

       bcoof(1:sdim) = min(bcoof(1:sdim),dboxf(1:sdim)) ! just in case
       bcooc(1:sdim) = min(bcooc(1:sdim),dboxc(1:sdim)) ! just in case
       bcool(1:sdim) = bcoof(1:sdim)-(bcooc(1:sdim)-1)*dboxl(1:sdim)
       !
       ! Eval. position of the element in the boxes grid according to lexicografic order
       !
       lboxc(ielem) = invor(dboxc(2_ip)*dboxc(1_ip)*(bcooc(3)-1_ip) + dboxc(1_ip)*(bcooc(2)-1_ip) + bcooc(1))
       lboxl(ielem) = dboxl(2_ip)*dboxl(1_ip)*(bcool(3)-1_ip) + dboxl(1_ip)*(bcool(2)-1_ip) + bcool(1)

       if(neighs(lboxc(ielem)) == 0_ip) then
          nneigs = nneigs + 1_ip
       endif
       neighs(lboxc(ielem)) = neighs(lboxc(ielem)) + 1_ip
    end do
    !
    !  Allocate and fill lenbc
    !
    do iboxc = 1_ip, nboxc
       if(neighs(iboxc) > 0_ip) then
          call memory_alloca(par_memor,' lenbc % l ',vacal,lenbc(iboxc) % l,neighs(iboxc))
       endif
    enddo
    neighs(:)=0_ip
    do ielem = 1_ip, nenti
       iboxc = lboxc(ielem)
       neighs(iboxc)  = neighs(iboxc) + 1_ip
       lenbc(iboxc) % l(neighs(iboxc)) = ielem 
    end do

    !
    ! Store information being sent in temporal buffers
    !
    ineig = 0_ip
    gbuff = 1_ip !accumulates total of "local" boxes sent
    call memory_alloca(par_memor,' bufse2 ',vacal,bufse2,max(1_ip,nneigs))

    do iboxc = 1_ip, nboxc

       if( neighs(iboxc) > 0_ip) then
          ineig = ineig + 1_ip
          !
          ! Count boxes to be sent and its weight
          !
          do ibuff = 1_ip, neighs(iboxc) 
             ielem = lenbc(iboxc) % l(ibuff)
             iboxl = lboxl(ielem)
             !
             ! In bufwei(iboxl*2) we store store iboxf
             ! In bufwei((iboxl-1)*2+1) we accumulate the weight of the bi
             !
             if(bufwei(iboxl*2) == 0.0_rp) then
                lnsend(iboxc-1_ip) = lnsend(iboxc-1_ip) + 1_4
                auxbf(lnsend(iboxc-1)) = iboxl
                bufwei( iboxl*2 ) =  1.0_rp
             endif
             bufwei((iboxl-1)*2+1) = bufwei((iboxl-1)*2+1) + lweig(ielem)
             lboxes(ielem)    = iboxl
          end do
          !
          ! Save info to be sent in temporal buffer
          !
          call memory_alloca(par_memor,' bufse2 ',vacal,bufse2(ineig) % a,max(1_ip,2_ip*lnsend(iboxc-1)))
          ibuff = 1_ip
          !
          ! we send:
          !  1. local index of the box
          !  2. asssociated weight
          !
          do kbuff=1,lnsend(iboxc-1_ip)
             iboxl                            = auxbf(kbuff)
             bufse2(ineig) % a((ibuff-1)*2+1) = real(iboxl,rp)         !iboxl
             bufse2(ineig) % a(ibuff*2)       = bufwei((iboxl-1)*2+1)  !weight
             bufwei(iboxl*2)                  = real(gbuff*2,rp)
             ibuff                            = ibuff + 1_ip
             gbuff                            = gbuff + 1_ip
          end do
          !
          ! Store the place where the result will be obtained in the buffer
          !
          do ibuff = 1_ip, neighs(iboxc) 
             ielem = lenbc(iboxc) % l(ibuff)
             lboxes(ielem) = int(bufwei(lboxes(ielem)*2),ip)
          end do
          !
          ! Initialize to 0.0 the components of bufwei used
          !
          do kbuff=1,lnsend(iboxc-1_ip)
             iboxl=auxbf(kbuff)
             bufwei((iboxl-1)*2+1)=0.0_rp
             bufwei(iboxl*2)=0.0_rp
          end do
       end if
    end do
    call memory_deallo(par_memor,' lenbc ',vacal,lenbc)
    call memory_deallo(par_memor,' auxbf ',vacal,auxbf)

  end subroutine define_distribution
  !
  ! This subroutine creates a bin structure of the geometrical domain storing the elements in the boxes
  !
  subroutine redistribute_weights()

    integer(ip)        :: ineig, neigh
    integer(ip)        :: ielem, idime, icont
    integer(ip)        :: iboxf, iboxc,iboxl, ibuff
    integer(4)         :: istat

    character(100), PARAMETER :: vacal = "redistribute_weights"

    nullify(send_buff,recv_buff)
    !
    ! iboxc0 and iboxc1 limits of the coarse bins treated in this phase
    !
    iboxc0=(iphase-1_ip)*PAR_MY_SIZE
    iboxc1=iboxc0+PAR_MY_SIZE-1_ip
    iboxc1_=min(iboxc1,nboxc-1)
    weigl(:)  = 0.0_rp
#ifndef MPI_OFF   
    lnsend_loc(:)=0_4
    lnrecv_loc(:)=0_4
    do iboxc=iboxc0,iboxc1_
       lnsend_loc(iboxc-iboxc0)=lnsend(iboxc)
    end do
    !
    ! Share communication requierements with others slaves
    !
    call MPI_Alltoall(lnsend_loc,1,MPI_INTEGER,lnrecv_loc,1,MPI_INTEGER,PAR_COMM,istat)
    !
    ! Allocate memory for communication structure
    !
    comm1 % PAR_COMM_WORLD = PAR_COMM
    comm2 % PAR_COMM_WORLD = PAR_COMM
    comm1 % nneig = 0_ip
    do iboxc = iboxc0, iboxc1 
       if(lnrecv_loc(iboxc-iboxc0) > 0_ip .or. lnsend_loc(iboxc-iboxc0) > 0_ip) then
          comm1 % nneig = comm1 % nneig + 1
       endif
    end do
    comm2 % nneig = comm1 % nneig
    nullify(comm1 % neights,comm1 % lsend_size,comm1 % lrecv_size)
    nullify(comm2 % neights,comm2 % lsend_size,comm2 % lrecv_size)
    call memory_alloca(par_memor,'NEIGHTS'   ,vacal,comm1 % neights,   max(1_ip, comm1 % nneig   ))
    call memory_alloca(par_memor,'LSEND_SIZE',vacal,comm1 % lsend_size,max(1_ip, comm1 % nneig+1 ))
    call memory_alloca(par_memor,'LRECV_SIZE',vacal,comm1 % lrecv_size,max(1_ip, comm1 % nneig+1 ))
    call memory_alloca(par_memor,'NEIGHTS'   ,vacal,comm2 % neights,   max(1_ip, comm2 % nneig   ))
    call memory_alloca(par_memor,'LSEND_SIZE',vacal,comm2 % lsend_size,max(1_ip, comm2 % nneig+1 ))
    call memory_alloca(par_memor,'LRECV_SIZE',vacal,comm2 % lrecv_size,max(1_ip, comm2 % nneig+1 ))
    comm1 % lsend_dim = 0_ip
    comm1 % lrecv_dim = 0_ip
    comm1 % lsend_size(1_ip) = 1_ip
    comm1 % lrecv_size(1_ip) = 1_ip
    ineig = 1_ip

    do iboxc = iboxc0, iboxc1
       if(lnrecv_loc(iboxc-iboxc0) > 0_ip .or. lnsend_loc(iboxc-iboxc0) > 0_ip) then
          comm1 % neights(ineig) = iboxc-iboxc0
          comm1 % lsend_dim = comm1 % lsend_dim + 2*lnsend_loc(iboxc-iboxc0)
          comm1 % lsend_size(ineig+1) = comm1 % lsend_size(ineig) + 2_ip*lnsend_loc(iboxc-iboxc0)
          comm1 % lrecv_dim = comm1 % lrecv_dim + 2*lnrecv_loc(iboxc-iboxc0)
          comm1 % lrecv_size(ineig+1) = comm1 % lrecv_size(ineig) + 2_ip*lnrecv_loc(iboxc-iboxc0)
          ineig = ineig + 1_ip
       endif
    end do
    comm2 % lsend_dim = comm1 % lrecv_dim
    comm2 % lrecv_dim = comm1 % lsend_dim
    comm2 % neights(:) = comm1 % neights(:)
    comm2 % lsend_size(:) = comm1 % lrecv_size(:)
    comm2 % lrecv_size(:) = comm1 % lsend_size(:)
    !
    ! Store data in unified buffer, perform boxes transfer, save local boxes infor in weiglh/sub_bins
    !
    call memory_alloca( par_memor,' SEND_BUFF ',vacal,send_buff,max(1_ip,2*comm1 % lsend_dim))
    call memory_alloca( par_memor,' RECV_BUFF ',vacal,recv_buff,max(1_ip,2*comm1 % lrecv_dim))
    icont = 1_ip
    do ineig = 1_ip, comm1 % nneig
       neigh = comm1 % neights(ineig)
       if(lnsend_loc(neigh) > 0) then
          do ibuff = 1_ip, lnsend_loc(neigh)
             send_buff(((icont-1)*2)+1)  = bufse2(ineis) % a((ibuff-1)*2+1) !iboxl
             send_buff(icont*2) = bufse2(ineis) % a(ibuff*2)                !weight
             icont = icont + 1_ip
          end do
          ineis = ineis + 1_ip
       end if
    enddo
    call PAR_SEND_RECEIVE_TO_ALL(send_buff,recv_buff,comm1,'ASYNCHRONOUS')
    icont = 1_ip
    do ineig = 1_ip, comm1 % nneig
       neigh = comm1 % neights(ineig)
       do ibuff = 1_ip, lnrecv_loc(neigh)
          iboxl          = int(recv_buff((icont-1)*2+1),ip)
          if( weigl(iboxl) == 0.0_rp ) then
             number_touched_bin = number_touched_bin + 1.0_rp
          end if
          ave_weight_bin = ave_weight_bin + recv_buff(icont*2)
          weigl(iboxl)   = weigl(iboxl) + recv_buff(icont*2)
          icont          = icont + 1             
          max_weight_bin = max(max_weight_bin,weigl(iboxl))
       end do
    end do
    !
    ! Eval weigh distribution in coarse bin
    !
    weigc_loc(:)=0_ip
    do icont = 1_ip, nboxl
       weigc_loc(PAR_MY_RANK + 1_ip) = weigc_loc(PAR_MY_RANK + 1_ip) + weigl(icont)
    end do
    call PAR_SUM(weigc_loc, PAR_COMM_IN4=PAR_COMM)

#endif

  end subroutine redistribute_weights
  subroutine local_partition

    use mod_maths,   only : maths_sfc_par_part_2
    !
    ! Partition of the local bins by means of SFC
    !
    weigc(iboxc0+1:iboxc1_+1)=weigc_loc(1:iboxc1_-iboxc0+1)
    if(weigc_loc(PAR_MY_RANK+1_ip)>0) then
       call maths_sfc_par_part_2(sdim,dboxl,weigl,npart_sfc,dboxc,iboxc0+PAR_MY_RANK+1_ip,weigc,partl,lcorr,totalw)
    endif

  end subroutine local_partition
  subroutine redistribute_result()

    implicit none
    integer(ip)                         :: ielem, ibuff,gbuffb,iaux
    character(100), PARAMETER :: vacal = "par_redistribute"
    !
    ! Here we use the receive buffer to send and the sent buffer to receive...
    !
    gbuffb=gbuff0
    do ibuff = 1_ip,comm1 % lrecv_dim/2_ip
       recv_buff(ibuff*2) = real(partl(int(recv_buff((ibuff-1)*2+1),ip)),rp)
    end do
    call PAR_SEND_RECEIVE_TO_ALL(recv_buff,send_buff,comm2,'ASYNCHRONOUS')
    !
    ! Assign elments to slave
    !
    do ielem = 1_ip, nenti
       if(lboxc(ielem)>= iboxc0+1_ip .and. lboxc(ielem)<=iboxc1_+1_ip) then
          lepar(ielem) = int(send_buff(lboxes(ielem)-gbuffb))
          gbuff0=max(gbuff0,lboxes(ielem))
       endif
    end do
    !
    ! Deallocate communication pointers
    !
    call memory_deallo(par_memor,'SEND_BUFF'         ,vacal,send_buff         )
    call memory_deallo(par_memor,'RECV_BUFF'         ,vacal,recv_buff         )
    call memory_deallo(par_memor,'PARTL'             ,vacal,partl             )
    call memory_deallo(par_memor,'COMM1 % NEIGHTS'   ,vacal,comm1 % neights   )
    call memory_deallo(par_memor,'COMM2 % NEIGHTS'   ,vacal,comm2 % neights   )
    call memory_deallo(par_memor,'COMM1 % LSEND_SIZE',vacal,comm1 % lsend_size)
    call memory_deallo(par_memor,'COMM2 % LSEND_SIZE',vacal,comm2 % lsend_size)
    call memory_deallo(par_memor,'COMM1 % LRECV_SIZE',vacal,comm1 % lrecv_size)
    call memory_deallo(par_memor,'COMM2 % LRECV_SIZE',vacal,comm2 % lrecv_size)

  end subroutine redistribute_result
  subroutine destructor

    character(100), PARAMETER :: vacal = "par_dealocate_sfc"

    if(associated(lboxes)) call memory_deallo(par_memor,' lboxes ',vacal,lboxes)
    if(associated(weigl))  call memory_deallo(par_memor,' weigl ',vacal,weigl)
    if(associated(weigc))  call memory_deallo(par_memor,' weigc ',vacal,weigc)
    if(associated(weigc_loc))  call memory_deallo(par_memor,' weigc_loc ',vacal,weigc_loc)
    if(associated(bufwei)) call memory_deallo(par_memor,' bufwei ',vacal,bufwei)
    if(associated(lnsend)) call memory_deallo(par_memor,' lnsend ',vacal,lnsend)
    if(associated(lnsend_loc)) call memory_deallo(par_memor,' lnsend_loc ',vacal,lnsend_loc)
    if(associated(lnrecv_loc)) call memory_deallo(par_memor,' lnrecv_loc ',vacal,lnrecv_loc)
    if(associated(lboxc))  call memory_deallo(par_memor,' lboxc  ',vacal,lboxc)
    if(associated(lboxl))  call memory_deallo(par_memor,' lboxl  ',vacal,lboxl)
    if(associated(neighs)) call memory_deallo(par_memor,' neighs ',vacal,neighs)
    if(associated(bufse2)) call memory_deallo(par_memor,' bufse2 ',vacal,bufse2)
    if(associated(lcorr)) call memory_deallo(par_memor,' lcorr ',' kk ',lcorr)

    call PAR_BARRIER(PAR_COMM_IN4=PAR_COMM)
    call cputim(time_partition)
    time_partition=time_partition-time_ini

  end subroutine destructor
end module mod_partition_sfc
!> @}
