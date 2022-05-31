!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    Alya sparse direct solver
!> @file    mod_alya_direct_solver.f90
!> @author  Guillaume Houzeaux
!> @date    13/06/2016
!> @brief   Alya sparse direct solver module 
!> @details Alya sparse direct solver module 
!>          Input: iA_in, jA_in, An_in
!>          1. Initialization
!>          2. Ordering
!>             Allocate Permr and Invpr if necessary
!>             Compute Invpr, Permr, iA, jA
!>          3. Symbolic factorization:  alya_Symbolical_CSR_LU
!>             Compute iL,jL,iU,jU
!>          +---> 
!>          |   
!>          4. Numerical factorization: alya_Numerical_CSR_LU
!>             Permute An_in into An
!>             Compute Ln,Un
!>          5. Solve An x = b
!>          6. Partial clean
!>             Deallocate LN,UN 
!>              |
!>          <---+
!>          7. Clean
!>             Deallocate iL,jL,iU,jU
!>             Deallocate iA,jA,An,Invpr,Permr
!------------------------------------------------------------------------

module mod_alya_direct_solver

  use def_kintyp, only :  ip,rp,i1p,lg
  use mod_memory, only :  memory_alloca
  use mod_memory, only :  memory_deallo
  implicit none 
  integer(8)           :: alya_CSR_memor(2) = 0_8
  private

  public :: alya_Symbolical_CSR_LU            ! Symbolical factorization: iL,jL,iU,jU
  public :: alya_Numerical_CSR_LU             ! Numerical factorization; Ln,Un
  public :: alya_CSR_LUSol                    ! Solution: x
  public :: alya_Numerical_CSR_LU_Deallocate  ! Deallocate Un,Ln
  public :: alya_Symbolical_CSR_LU_Deallocate ! Deallocate 
  public :: alya_CSR_LUfin                    ! Deallocate everything

  public :: alya_cholesky_factorization       ! Cholesky factorization
  public :: alya_cholesky_solution            ! Cholesky solution

  public :: alya_direct_solver_initialization ! Initializaiton of the module

  public :: alya_CSR_memor                    ! Memory counter

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Direct solver initialization
  !> @details Initialize some variables of this module
  !> 
  !-----------------------------------------------------------------------

  subroutine alya_direct_solver_initialization()

    alya_CSR_memor = 0_8

  end subroutine alya_direct_solver_initialization

  !---------------------------------------------------------------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-11
  !> @brief   Direct solver initialization
  !> @details This routine performs Symbolical CSR LU factorization of a matrix A stored in the CSR format.
  !>          INPUT ARGUMENTS
  !>            NBROWS .... Number of rows of matrix A (dimension of matrix A)
  !>            IA ........ CSR format: index vector for beginning of a row block for matrix A
  !>            JA ........ CSR format: index vector for column numbers for matrix A
  !>          OUTPUT ARGUMENTS 
  !>            IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
  !>            JL ........ CSR format: index vector for column numbers for lower triangular matrix L
  !>            IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
  !>            JU ........ CSR format: index vector for column numbers for upper triangular matrix U
  !>
  !---------------------------------------------------------------------------------------------------------------------------

  subroutine alya_Symbolical_CSR_LU(nbrows,iA,jA,iL,jL,iU,jU,rfillin,MEMORY_COUNTER)

    implicit none
    integer(ip), intent(in)              :: nbrows
    integer(ip), intent(in)              :: iA(*),jA(*)
    integer(ip), pointer                 :: iL(:),jL(:)
    integer(ip), pointer                 :: iU(:),jU(:)
    real(rp),    intent(out),   optional :: rfillin
    integer(8),  intent(inout), optional :: MEMORY_COUNTER(2)
    integer(ip)                          :: i,j,k,col,nxt,totalL,totalU
    integer(ip)                          :: rowIA_ini,rowIA_diag,rowIA_fin  
    integer(ip)                          :: minColL,maxColL,minColU,maxColU
    integer(ip), pointer                 :: firstL(:),firstU(:),seen(:),nextU(:),iwa(:)   
    integer(ip), pointer                 :: nzL(:),nzU(:)       !Local arrays for L/U matrices 
    type(i1p),   pointer                 :: ptL(:),ptU(:)
    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = alya_CSR_memor
    end if

    if( associated(iL) .or. associated(jL) .or. associated(iU) .or. associated(jU) ) then
       write(*,*) 'Error!! Pointers already associated.'
    else
       !
       ! Alloc local work arrays
       !
       nullify(firstL)
       nullify(firstU)
       nullify(seen)
       nullify(nextU)
       nullify(iwa)
       nullify(nzL)
       nullify(nzU)
       nullify(ptL)
       nullify(ptU)

       call memory_alloca(memor,'FIRSTL','Symbolical_CSR_LU',firstL,nbrows)
       call memory_alloca(memor,'FIRSTU','Symbolical_CSR_LU',firstU,nbrows)
       call memory_alloca(memor,'SEEN'  ,'Symbolical_CSR_LU',seen,nbrows)
       call memory_alloca(memor,'NEXTU' ,'Symbolical_CSR_LU',nextU,nbrows)
       call memory_alloca(memor,'IWA'   ,'Symbolical_CSR_LU',iwa,nbrows)
       call memory_alloca(memor,'NZL'   ,'Symbolical_CSR_LU',nzL,nbrows)
       call memory_alloca(memor,'NZU'   ,'Symbolical_CSR_LU',nzU,nbrows)
       call memory_alloca(memor,'PTL'   ,'Symbolical_CSR_LU',ptL,nbrows)
       call memory_alloca(memor,'PTU'   ,'Symbolical_CSR_LU',ptU,nbrows)
       !
       ! Initialize the fill-in links and local L/U arrays
       !
       do i= 1, nbrows
          firstL(i) = -1
          firstU(i) = -1
          seen(i)   = -1
          nzL(i)    =  0
          nzU(i)    =  0
          nullify( ptL(i) % l , ptU(i) % l )   
       end do
       !
       ! Initialize total number of non-zero elements in L/U
       !
       totalL = 0
       totalU = 0
       !
       ! Main loop in rows
       !
       do i= 1, nbrows
          minColL = nbrows
          maxColL = -1
          minColU = nbrows
          maxColU = -1
          !
          !For all the elements in the i-th row of the matrix A
          !
          rowIA_ini  = iA(i)
          rowIA_fin  = iA(i+1)-1
          rowIA_diag = iA(i+1)
          !
          ! Find a diagonal element in the i-th row of the matrix A
          !
          loop1: do k = rowIA_ini, rowIA_fin
             if( jA(k) == i ) then
                rowIA_diag = k
                exit loop1
             end if
          end do loop1
          !
          ! For all the elements in the i-th row of the matrix A before diagonal
          !
          do j = rowIA_ini,rowIA_diag-1
             col = jA(j)
             if( seen(col) /= i ) then
                seen(col) = i
                minColL   = min(minColL,col)
                maxColL   = max(maxColL,col)
                !
                ! Compute Reachable Set from L
                !
                col = firstL(col)
                loop2:do while (col /= -1)
                   if( seen(col) /= i ) then
                      minColL   = min(minColL,col)
                      maxColL   = max(maxColL,col)
                      seen(col) = i
                      col = firstL(col)
                   else
                      exit loop2
                   end if
                end do loop2
             end if
          end do
          !
          ! For all the elements in the i-th row of the matrix A after diagonal
          !
          do j = rowIA_diag+1, rowIA_fin
             col       = jA(j)
             seen(col) = i
             minColU   = min(minColU,col)
             maxColU   = max(maxColU,col)
          end do
          !
          ! Compute Reachable Set from U
          ! 
          k = firstU(i)
          do while (k /= -1)
             !
             ! For all the elements in the k-th row of the matrix U without the diagonal element
             !
             do j= 2, nzU(k)
                col = ptU(k) % l(j)
                if( col > i ) then
                   minColU = min(minColU,col)
                   maxColU = max(maxColU,col)
                   seen(col) = i
                end if
             end do
             k = nextU(k)
          end do
          !
          ! For all the non-zero elements of the matrix L. L matrix is stored without the diagonal element
          !
          nxt = 0
          do j= minColL, maxColL
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
                if( firstL(j) == -1) firstL(j) = i
             end if
          end do
          !
          ! Allocate space for L values and copy iwa to ptL 
          !
          totalL = totalL + nxt
          nzL(i) = nxt
          call memory_alloca(memor,'PTL % L','Symbolical_CSR_LU',ptL(i) % l,nxt,'DO_NOT_INITIALIZE')
          !if( istat /= 0 ) call memerr(0_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
          do k= 1, nxt
             ptL(i) % l(k) = iwa(k)
          end do
          !
          ! For all the non-zero elements of the matrix U. The diagonal element must be included
          !
          iwa(1) = i       !The diagonal element 
          nxt    = 1
          do j= minColU, maxColU
             if( seen(j) == i ) then
                iwa(nxt+1) = j
                nxt = nxt + 1
             end if
          end do
          !
          ! Put the proper link for future fill-in generation 
          !
          if( nxt > 2 ) then
             col         = iwa(2)
             nextU(i)    = firstU(col)
             firstU(col) = i
          end if
          !
          ! Allocate space for U values and copy iwa to ptU 
          !
          totalU = totalU + nxt
          nzU(i) = nxt 
          call memory_alloca(memor,'PTU % L','Symbolical_CSR_LU',ptU(i) % l,nxt,'DO_NOT_INITIALIZE')
          !if( istat /= 0 ) call memerr(0_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
          do k= 1, nxt
             ptU(i) % l(k) = iwa(k)
          end do
       end do
       !
       ! Put L and U in CSR format 
       !
       call memory_alloca(memor,'IL','Symbolical_CSR_LU',iL,nbrows+1_ip)
       call memory_alloca(memor,'IU','Symbolical_CSR_LU',iU,nbrows+1_ip)
       call memory_alloca(memor,'JL','Symbolical_CSR_LU',jL,totalL)
       call memory_alloca(memor,'JU','Symbolical_CSR_LU',jU,totalU)
       !
       ! CSR format for the L matrix
       !
       iL(1) = 1
       do i= 1, nbrows
          iL(i+1) = iL(i) + nzL(i)
          nxt = iL(i) - 1
          do k= 1, nzL(i)
             jL(nxt+k) = ptL(i) % l(k)
          end do
       end do
       !
       ! CSR format for the U matrix
       !
       iU(1) = 1
       do i= 1, nbrows
          iU(i+1) = iU(i) + nzU(i)

          nxt = iU(i) - 1
          do k= 1, nzU(i)
             jU(nxt+k) =  ptU(i) % l(k)
          end do
       end do
       !
       ! Free local arrays
       !
       !do i= 1, nbrows
       !   deallocate(ptL(i) % l) 
       !   if( istat /= 0 ) call memerr(2_ip,'ptL(i)%l','Symbolical_CSR_LU',0_ip)
       !   deallocate(ptU(i) % l)
       !   if( istat /= 0 ) call memerr(2_ip,'ptU(i)%l','Symbolical_CSR_LU',0_ip)
       !end do

       call memory_deallo(memor,'PTU'   ,'Symbolical_CSR_LU',ptU)
       call memory_deallo(memor,'PTL'   ,'Symbolical_CSR_LU',ptL)
       call memory_deallo(memor,'NZU'   ,'Symbolical_CSR_LU',nzU)
       call memory_deallo(memor,'NZL'   ,'Symbolical_CSR_LU',nzL)
       call memory_deallo(memor,'IWA'   ,'Symbolical_CSR_LU',iwa)
       call memory_deallo(memor,'NEXTU' ,'Symbolical_CSR_LU',nextU)
       call memory_deallo(memor,'SEEN'  ,'Symbolical_CSR_LU',seen)
       call memory_deallo(memor,'FIRSTU','Symbolical_CSR_LU',firstU)
       call memory_deallo(memor,'FIRSTL','Symbolical_CSR_LU',firstL)

       if( present(rfillin) ) then
          rfillin = (real(iL(nbrows+1)-1,rp)+real(iU(nbrows+1)-1,rp))/real(iA(nbrows+1)-1,rp)
       end if

       if( present(MEMORY_COUNTER) ) then
          MEMORY_COUNTER = memor
       else
          alya_CSR_memor = memor 
       end if

    end if

  end subroutine Alya_Symbolical_CSR_LU

  subroutine alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numerical_CSR_LU1
    ! DESCRIPTION
    !    This routine performs Numerical CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
    !    This routine deals with matrices whose number of degrees of freedom is 1.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    ! OUTPUT ARGUMENTS 
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !    SING ...... Singularity marker
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                     :: nbrows
    integer(ip), intent(in)                     :: iA(*), jA(*)
    real(rp),    intent(in)                     :: An(*)
    integer(ip), intent(in)                     :: iL(*), jL(*)
    real(rp),    intent(inout), pointer         :: Ln(:)
    integer(ip), intent(in)                     :: iU(*), jU(*)
    real(rp),    intent(inout), pointer         :: Un(:)
    integer(ip), intent(out)                    :: sing
    integer(ip), intent(in),  pointer, optional :: invpr(:)
    integer(ip), intent(in),  pointer, optional :: permr(:)
    integer(ip), intent(in),  pointer, optional :: iAold(:)
    integer(ip), intent(in),  pointer, optional :: jAold(:)
    !
    ! Local variables
    !
    integer(ip)                                 :: i,j,k,col,coli,dia,nz
    integer(ip)                                 :: inew,jnew,iold,jold
    real(rp)                                    :: pivot
    real(rp),    pointer                        :: wa(:) 
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(alya_CSR_memor,'WA','Numerical_CSR_LU1',wa,nbrows,'DO_NOT_INITIALIZE')
    sing = 0
    !
    ! Main loop in rows
    !
    do i = 1, nbrows
       !
       ! Initialize wa with zeros and values of the matrix A 
       !
       do j = iL(i), iL(i+1)-1
          col     = jL(j)
          wa(col) = 0.0_rp
       end do
       do j = iU(i), iU(i+1)-1
          col     = jU(j)
          wa(col) = 0.0_rp
       end do
       if( present(invpr) ) then
          iold = invpr(i) 
          do jold = iAold(iold), iAold(iold+1)-1
             col     = permr(jAold(jold))
             wa(col) = An(jold) 
          end do
       else
          do j = iA(i), iA(i+1)-1
             col     = jA(j)
             wa(col) = An(j)
          end do
       end if
       !
       ! Factorize the i-th row of the matrix A
       ! For all the elements in the i-th row before diagonal
       !
       do j = iL(i), iL(i+1)-1
          col     = jL(j) 
          dia     = iU(col)
          pivot   = wa(col) * Un(dia)       !U diagonal is stored inverted 
          wa(col) = pivot
          !For all the elements in the col-th row after diagonal
          do k = iU(col)+1, iU(col+1)-1
             coli     = jU(k)
             wa(coli) = wa(coli) - pivot * Un(k)
          end do
       end do
       !
       ! Store factorized i-th row in L/U 
       !
       pivot = wa(i)
       if( pivot /= 0.0_rp ) then
          !Store the elements before diagonal in L
          do j     = iL(i), iL(i+1)-1
             col   = jL(j)
             Ln(j) = wa(col)
          end do
          !Store the diagonal element in U, at the first place in the i-th row of the matrix U
          dia     = iU(i) 
          Un(dia) = 1.0_rp / pivot       !Diagonal element is stored inverted
          !Store the elements after diagonal in U
          do j = iU(i)+1, iU(i+1)-1
             col   = jU(j)
             Un(j) = wa(col)
          end do
       else
          write(*,*) 'Singular pivot in row ', i
          sing = -1
       end if
    end do

    !Deallocate working array
    call memory_deallo(alya_CSR_memor,'WA','Numerical_CSR_LU1',wa)

  end subroutine Alya_Numerical_CSR_LU1

  subroutine alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    Numerical_CSR_LU
    ! DESCRIPTION
    !    This routine performs Numerical CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    NDOF ...... Number of degrees of freedom in each node
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    ! OUTPUT ARGUMENTS 
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !    SING ...... Singularity marker
    !---------------------------------------------------------------------------------------------------------------------------
    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                     :: nbrows, ndof
    integer(ip), intent(in)                     :: iA(*)
    integer(ip), intent(in)                     :: jA(*)
    real(rp),    intent(in)                     :: An(ndof,ndof,*)
    integer(ip), intent(in)                     :: iL(*)
    integer(ip), intent(in)                     :: jL(*)
    real(rp),    intent(out)                    :: Ln(ndof,ndof,*)
    integer(ip), intent(in)                     :: iU(*)
    integer(ip), intent(in)                     :: jU(*)
    real(rp),    intent(out)                    :: Un(ndof,ndof,*)
    integer(ip), intent(out)                    :: sing
    integer(ip), intent(in),  pointer, optional :: invpr(:)
    integer(ip), intent(in),  pointer, optional :: permr(:)
    integer(ip), intent(in),  pointer, optional :: iAold(:)
    integer(ip), intent(in),  pointer, optional :: jAold(:)
    !
    ! Local variables
    !
    integer(ip)                                 :: i,j,k,l,s,t,col,dia,coli
    integer(ip)                                 :: iold,jold
    real(rp)                                    :: pivot
    real(rp), pointer                           :: wa(:,:,:)     
    !
    ! Allocate working array, wa
    !
    nullify(wa)
    call memory_alloca(alya_CSR_memor,'WA','Numerical_CSR_LUM',wa,ndof,ndof,nbrows,'DO_NOT_INITIALIZE')

    sing = 0
    !
    ! Main loop in rows
    !
    do i= 1, nbrows
       !
       ! Initialize wa with zeros and values of the matrix A 
       !
       do j= iL(i), iL(i+1)-1
          col = jL(j)
          do l= 1, ndof
             do k= 1, ndof
                wa(k,l,col) = 0.0_rp
             end do
          end do
       end do
       do j= iU(i), iU(i+1)-1
          col = jU(j)
          do l= 1, ndof
             do k= 1, ndof
                wa(k,l,col) = 0.0_rp
             end do
          end do
       end do
       if( present(invpr) ) then
          iold = invpr(i) 
          do jold = iAold(iold), iAold(iold+1)-1
             col = permr(jAold(jold))
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = An(l,k,jold)
                end do
             end do
          end do
       else
          do j= iA(i), iA(i+1)-1
             col = jA(j)
             do l= 1, ndof
                do k= 1, ndof
                   wa(k,l,col) = An(l,k,j)
                end do
             end do
          end do
       end if
       !
       ! Factorize all the blocks in the i-th row of the matrix A
       ! For all the elements (blocks) in the i-th row before diagonal
       !
       do j= iL(i), iL(i+1)-1
          col = jL(j)
          dia = iU(col)
          !Factorize the block in column col
          do l= 1, ndof
             do k= 1, ndof      
                pivot = wa(k,l,col) * Un(l,l,dia)       !U diagonal is stored inverted
                wa(k,l,col) = pivot
                !For all the elements in the k-th row inside the current block
                do s= l+1, ndof
                   wa(k,s,col) = wa(k,s,col) - pivot * Un(l,s,dia)
                end do
             end do
          end do
          !For all the blocks after diagonal
          do s= iU(col)+1, iU(col+1)-1
             coli = jU(s)
             do l= 1, ndof
                do k= 1, ndof
                   do t= 1, ndof 
                      wa(k,t,coli) = wa(k,t,coli) - Un(l,t,s) * wa(k,l,col)
                   end do
                end do
             end do
          end do
       end do
       !
       ! Factorize the lower triangle of the diagonal block
       !
       do k= 1, ndof
          do l= 1, k-1
             pivot = wa(k,l,i) * wa(l,l,i)
             wa(k,l,i) = pivot
             !For all the elements in the k-th row inside the diagonal block
             do s= l+1, ndof
                wa(k,s,i) = wa(k,s,i) - pivot * wa(l,s,i)
             end do
          end do
          !Check if some of the diagonal elements of the diagonal block is zero
          if( wa(k,k,i) /= 0.0_rp ) then
             wa(k,k,i) = 1.0_rp / wa(k,k,i)
          else
             sing = -1
             print*,'singular matrix'
             exit
          end if
       end do
       !
       ! For all the blocks after diagonal
       !
       do s= iU(i)+1, iU(i+1)-1
          coli = jU(s)
          do k= 1, ndof
             do l= 1, k-1
                do t= 1, ndof 
                   wa(k,t,coli) = wa(k,t,coli) - wa(l,t,coli) * wa(k,l,i)
                end do
             end do
          end do
       end do
       !
       ! Store factorized i-th row of blocks in L/U matrices
       !
       if( sing == 0 ) then
          !Store the elements of the blocks before diagonal in L
          do j= iL(i), iL(i+1)-1
             col= jL(j)
             do l= 1, ndof
                do k= 1, ndof   
                   Ln(k,l,j) = wa(k,l,col)
                end do
             end do
          end do
          !Store the elements of the diagonal block in U
          !Note that the WHOLE diagonal block is stored in the matrix U
          col = iU(i)
          do l= 1, ndof
             do k= 1, ndof   
                Un(k,l,col) = wa(k,l,i)
             end do
          end do
          !Store the elements of the blocks after diagonal in U
          do j= iU(i)+1, iU(i+1)-1
             col= jU(j)
             do l= 1, ndof
                do k= 1, ndof
                   Un(k,l,j) = wa(k,l,col)    
                end do
             end do
          end do
       end if
    end do
    !
    ! Deallocate working array   
    !
    call memory_deallo(alya_CSR_memor,'WA','Numerical_CSR_LUM',wa)

  end subroutine Alya_Numerical_CSR_LUM

  subroutine alya_Numerical_CSR_LU(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LU_Factorization
    ! DESCRIPTION
    !    This routine performs CSR LU factorization of a matrix A stored in the CSR Row (CSR) format.
    ! INPUT ARGUMENTS
    !    NBROWS .... Number of rows of matrix A (dimension of matrix A)
    !    NDOF ...... Number of degrees of freedom in each node
    !    IA ........ CSR format: index vector for beginning of a row block for matrix A
    !    JA ........ CSR format: index vector for column numbers for matrix A
    !    AN ........ Sparse complex matrix A in BCSR (Blocked CSR Row) format
    ! OUTPUT ARGUMENTS 
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !    SING ...... Singularity marker
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)                    :: nbrows  !< Number of rows of matrix A (dimension of matrix A)
    integer(ip), intent(in)                    :: ndof    !< Number of degrees of freedom in each node
    integer(ip), intent(in)                    :: iA(*)   !< CSR format: index vector for beginning of a row block for matrix A
    integer(ip), intent(in)                    :: jA(*)   !< CSR format: index vector for column numbers for matrix A
    real(rp),    intent(in)                    :: An(*)   !< Sparse complex matrix A in BCSR (Blocked CSR Row) format
    integer(ip), pointer                       :: iL(:)   !< CSR format: index vector for beginning of a row block for lower triangular matrix L
    integer(ip), pointer                       :: jL(:)   !< CSR format: index vector for column numbers for lower triangular matrix L
    real(rp),    pointer                       :: Ln(:)   !< Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    integer(ip), pointer                       :: iU(:)   !< CSR format: index vector for beginning of a row block for upper triangular matrix U
    integer(ip), pointer                       :: jU(:)   !< CSR format: index vector for column numbers for upper triangular matrix U
    real(rp),    pointer                       :: Un(:)   !< Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    integer(ip), intent(out)                   :: sing    !< Singularity marker
    integer(ip), intent(in), optional, pointer :: invpr(:)
    integer(ip), intent(in), optional, pointer :: permr(:)
    integer(ip), intent(in), optional, pointer :: iAold(:)
    integer(ip), intent(in), optional, pointer :: jAold(:)
    logical(lg)                                :: lpermute
    !
    ! Allocate Un and Ln if necessary
    !
    call alya_Numerical_CSR_LU_Allocate(nbrows,ndof,iL,iU,Ln,Un)
    !
    ! Permute or not
    !
    lpermute = present(invpr) .and. present(permr) .and. present(iAold) .and. present(jAold)
    !
    ! Perform Numerical factorization dependig on number of degrees of freedom
    !
    if( ndof == 1 ) then
       if( lpermute ) then
          call alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)
       else
          call alya_Numerical_CSR_LU1(nbrows,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)
       end if
    else
       if( lpermute ) then 
          call alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing,invpr,permr,iAold,jAold)  
       else
          call alya_Numerical_CSR_LUM(nbrows,ndof,iA,jA,An,iL,jL,Ln,iU,jU,Un,sing)  
       end if
    end if

  end subroutine Alya_Numerical_CSR_LU

  subroutine alya_CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Lsol
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is L - lower triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: OLD = INVPR(NEW)
    !    IL ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JL ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    LN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !    IU ........ CSR format: index vector for beginning of a row block for upper triangular matrix U
    !    JU ........ CSR format: index vector for column numbers for upper triangular matrix U
    !    UN ........ Sparse complex upper triangular matrix U in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows,ndof
    integer(ip), pointer      :: invpR(:)
    integer(ip), intent(in)   :: iL(*),jL(*)
    real(rp),    intent(in)   :: Ln(ndof,ndof,*)
    integer(ip), intent(in)   :: iU(*),jU(*)
    real(rp),    intent(in)   :: Un(ndof,ndof,*)
    real(rp),    intent(in)   :: b(ndof,*)
    real(rp),    intent(out)  :: x(ndof,*)
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,l,col

    if( .not. associated(invpR) ) then
       !
       ! x <= b
       !
       do j = 1,nbrows
          x(1:ndof,j) = b(1:ndof,j)
       end do
    else 
       !
       ! x (new) <= b (old)
       !
       if( ndof == 1 ) then
          do j = 1,nbrows
             k = invpR(j)
             x(1,j) = b(1,k)
          end do
       else
          do j = 1,nbrows
             k = invpR(j)
             x(1:ndof,j) = b(1:ndof,k)
          end do
       end if
    end if

    if( ndof == 1 ) then
       do i = 1, nbrows
          do j = iL(i),iL(i+1)-1
             col = jL(j)
             x(1,i) = x(1,i) - Ln(1,1,j) * x(1,col)
          end do
       end do
    else
       do i = 1, nbrows
          do j = iL(i),iL(i+1)-1
             col = jL(j)
             do l = 1,ndof  
                do k = 1,ndof
                   x(k,i) = x(k,i) - Ln(k,l,j) * x(l,col)
                end do
             end do
          end do
          !Update using the lower triangle of the diagonal block stored in the U matrix
          j = iU(i)
          do l = 1,ndof
             do k = l+1,ndof
                x(k,i) = x(k,i) - Un(k,l,j) * x(l,i)
             end do
          end do
       end do
    end if

  end subroutine Alya_CSR_Lsol

  subroutine Alya_CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_Usol
    ! DESCRIPTION
    !    This routine solves a triangular linear system of equations where the system matrix is U - upper triangular matrix
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPC ..... 
    !    IU ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JU ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    UN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows, ndof
    integer(ip), pointer      :: invpC(:)
    integer(ip), intent(in)   :: iU(*), jU(*)
    real(rp),    intent(in)   :: Un(ndof,ndof,*)
    real(rp),    intent(inout):: b(ndof,*)
    real(rp),    intent(out)  :: x(ndof,*)
    !
    ! Local variables
    !
    integer(ip)               :: i,j,k,l,col    

    if( ndof == 1 ) then
       DO i= nbrows, 1, -1
          !Update rhs with the block after diagonal block
          do j= iU(i)+1, iU(i+1)-1
             col = jU(j)
             b(1,i) = b(1,i) - Un(1,1,j) * b(1,col)
          end do
          !Solve the upper triangle in the diagonal of the U matrix
          j = iU(i)
          b(1,i) = b(1,i) * Un(1,1,j)
       END DO
    else
       DO i= nbrows, 1, -1
          !Update rhs with the block after diagonal block
          do j= iU(i)+1, iU(i+1)-1
             col = jU(j)
             do l= 1, ndof
                do k= 1, ndof
                   b(k,i) = b(k,i) - Un(k,l,j) * b(l,col)
                end do
             end do
          end do
          !Solve the upper triangle in the diagonal of the U matrix
          j = iU(i)
          do l= ndof, 1, -1
             b(l,i) = b(l,i) * Un(l,l,j)
             do k= 1, l-1
                b(k,i) = b(k,i) - Un(k,l,j) * b(l,i)
             end do
          end do
       END DO
    end if

    if( .not. associated(invpC) ) then
       !
       ! x <= b
       !
       do j= 1, nbrows
          do i= 1, ndof
             x(i,j) = b(i,j)
          end do
       end do
    else
       !
       ! x (old) <= b (new)
       !
       if( ndof == 1 ) then
          do j= 1, nbrows
             col = invpC(j)
             x(1,col) = b(1,j)
          end do
       else
          do j= 1, nbrows
             col = invpC(j)
             do i= 1, ndof
                x(i,col) = b(i,j)
             end do
          end do
       end if
    end if

  end subroutine Alya_CSR_Usol

  subroutine alya_CSR_LUsol(nbrows,ndof,invpR,invpC,iL,jL,Ln,iU,jU,Un,b,x)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUsol
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.
    !
    !          +---------+
    !    b --> | LUx = b | --> x
    !          +---------+
    !   (old)    (new)        (old)
    !
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... 
    !    INVPC ..... 
    !    IU ........ CSR format: index vector for beginning of a row block for lower triangular matrix L
    !    JU ........ CSR format: index vector for column numbers for lower triangular matrix L
    !    UN ........ Sparse complex lower triangular matrix L in BCSR (Blocked CSR Row) format
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)  :: nbrows,ndof
    integer(ip), pointer     :: invpR(:),invpC(:)
    integer(ip), intent(in)  :: iL(*),jL(*)
    real(rp),    intent(in)  :: Ln(ndof,ndof,*)
    integer(ip), intent(in)  :: iU(*),jU(*)
    real(rp),    intent(in)  :: Un(ndof,ndof,*)
    real(rp),    intent(in)  :: b(ndof,*)
    real(rp),    intent(out) :: x(ndof,*)
    real(rp),    pointer     :: b_tmp(:,:)
    !
    ! Local variables
    !
    integer(ip)               :: i,j

    nullify(b_tmp)
    call memory_alloca(alya_CSR_memor,'B_TMP','CSR_LUsol',b_tmp,ndof,nbrows,'DO_NOT_INITIALIZE')

    b_tmp(1:ndof,1:nbrows) = b(1:ndof,1:nbrows)

    call alya_CSR_Lsol(nbrows,ndof,invpR,iL,jL,Ln,iU,jU,Un,b_tmp,x)
    call alya_CSR_Usol(nbrows,ndof,invpC,iU,jU,Un,x,b_tmp)

    do j = 1,nbrows
       x(1:ndof,j) = b_tmp(1:ndof,j)
    end do

    call memory_deallo(alya_CSR_memor,'B_TMP','CSR_LUsol',b_tmp)

  end subroutine Alya_CSR_LUsol

  subroutine alya_CSR_LUfin(iL,jL,Ln,iU,jU,Un,invpR,invpC)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none

    !Dummy arguments
    integer(ip), pointer           :: iL(:),jL(:)
    real(rp),    pointer           :: Ln(:)
    integer(ip), pointer           :: iU(:),jU(:)
    real(rp),    pointer           :: Un(:)
    integer(ip), pointer, optional :: invpR(:)
    integer(ip), pointer, optional :: invpC(:)

    call memory_deallo(alya_CSR_memor,'IL','CSR_LUfin',iL)
    call memory_deallo(alya_CSR_memor,'JL','CSR_LUfin',jL)
    call memory_deallo(alya_CSR_memor,'LN','CSR_LUfin',Ln)
    call memory_deallo(alya_CSR_memor,'IU','CSR_LUfin',iU)
    call memory_deallo(alya_CSR_memor,'JU','CSR_LUfin',jU)
    call memory_deallo(alya_CSR_memor,'UN','CSR_LUfin',Un)
    if( present(invpR) ) call memory_deallo(alya_CSR_memor,'invpR','CSR_LUfin',invpR)
    if( present(invpC) ) call memory_deallo(alya_CSR_memor,'invpC','CSR_LUfin',invpC)

  end subroutine Alya_CSR_LUfin

  subroutine alya_Symbolical_CSR_LU_Deallocate(iL,jL,iU,jU,MEMORY_COUNTER)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none

    !Dummy arguments
    integer(ip),           pointer       :: iL(:),jL(:)
    integer(ip),           pointer       :: iU(:),jU(:)
    integer(8),  optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                           :: memor(2)
    
    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = alya_CSR_memor
    end if

    call memory_deallo(memor,'iL','CSR_LUfin',iL)
    call memory_deallo(memor,'jL','CSR_LUfin',jL)
    call memory_deallo(memor,'iU','CSR_LUfin',iU)
    call memory_deallo(memor,'jU','CSR_LUfin',jU)

    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor
    else
       alya_CSR_memor = memor 
    end if

  end subroutine Alya_Symbolical_CSR_LU_Deallocate

  subroutine alya_Numerical_CSR_LU_Deallocate(Ln,Un)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_LUfin
    ! DESCRIPTION
    !    This routine solves LU triangular linear systems of equations.eallocates memoryd
    ! INPUT ARGUMENTS
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none

    !Dummy arguments
    real(rp), pointer :: Ln(:),Un(:)

    call memory_deallo(alya_CSR_memor,'LN','CSR_LUfin',Ln)
    call memory_deallo(alya_CSR_memor,'UN','CSR_LUfin',Un)

  end subroutine Alya_Numerical_CSR_LU_Deallocate

  subroutine alya_Numerical_CSR_LU_Allocate(nbrows,ndof,iL,iU,Ln,Un)

    integer(ip), intent(in)             :: nbrows
    integer(ip), intent(in)             :: ndof
    integer(ip), intent(in),    pointer :: iL(:),iU(:)
    real(rp),    intent(inout), pointer :: Ln(:),Un(:)
    integer(ip)                         :: nnzL,nnzU,ndof2
    !
    ! Allocate memory for L/U
    !
    ndof2 = ndof*ndof
    if( .not. associated(Ln) ) then
       nnzL = iL(nbrows+1)-1
       call memory_alloca(alya_CSR_memor,'LN','alya_Numerical_CSR_LU_Allocate',Ln,nnzL*ndof2)
    end if
    if( .not. associated(Un) ) then
       nnzU = iU(nbrows+1)-1
       call memory_alloca(alya_CSR_memor,'UN','alya_Numerical_CSR_LU_Allocate',Un,nnzU*ndof2) 
    end if

  end subroutine alya_Numerical_CSR_LU_Allocate

  subroutine alya_CSR_Permute_and_Copy_matrix(nbrows,ndof,iAin,jAin,iAout,jAout,invpR,Ain,Aout)

    !---------------------------------------------------------------------------------------------------------------------------
    ! NAME
    !    CSR_SMVP
    ! DESCRIPTION
    !    This routine copy a matrix to another
    ! INPUT ARGUMENTS
    !    NBROWS .... Dimension of the sysytem
    !    NDOF ...... Number of degrees of freedom in each node
    !    INVPR ..... Permutation array: OLD = INVPR(NEW)
    !     B ........ RHS of the system
    ! OUTPUT ARGUMENTS 
    !     X ........ Solution vector of the system
    !---------------------------------------------------------------------------------------------------------------------------

    implicit none
    !
    ! Dummy arguments
    !
    integer(ip), intent(in)   :: nbrows            !< Number of rows
    integer(ip), intent(in)   :: ndof              !< Number of dof per row
    integer(ip), intent(in)   :: iAin(*)           !< Matrix graph
    integer(ip), intent(in)   :: jAin(*)           !< Matrix graph
    integer(ip), intent(in)   :: iAout(*)          !< Matrix graph
    integer(ip), intent(in)   :: jAout(*)          !< Matrix graph
    integer(ip), pointer      :: invpR(:)          !< Inverse permutation
    real(rp),    intent(in)   :: Ain(ndof,ndof,*)  !< Input matrix
    real(rp),    intent(out)  :: Aout(ndof,ndof,*) !< Output matrix
    integer(ip)               :: ii,jj,iz,ki,kj
    integer(ip)               :: jjold,iiold,kz

    if( .not. associated(invpR) ) then
       do ii = 1,nbrows
          do iz = iAin(ii),iAin(ii+1)-1
             Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,iz)
          end do
       end do
    else
       do ii = 1,nbrows
          iiold = invpR(ii)
          do iz = iAout(ii),iAout(ii+1)-1
             jj    = jAout(iz)
             jjold = invpR(jj)
             kz    = iAin(iiold)
             do while( jAin(kz) /= jjold )
                kz = kz + 1
             end do
             Aout(1:ndof,1:ndof,iz) = Ain(1:ndof,1:ndof,kz)
          end do
       end do
    end if

  end subroutine Alya_CSR_Permute_and_Copy_matrix

  !-----------------------------------------------------------------------
  !****f* domain/skygro
  ! NAME
  !    skygro
  ! DESCRIPTION
  !    Set up the skyline structure
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------

  subroutine alya_cholesky_initialization(nbrows,ndof,kfl_symmetric,ia,ja,nskyl,iskyl,idiag)

    integer(ip), intent(in)             :: nbrows
    integer(ip), intent(in)             :: ndof
    integer(ip), intent(in)             :: kfl_symmetric
    integer(ip), intent(out)            :: nskyl
    integer(ip), intent(in),    pointer :: ia(:)
    integer(ip), intent(in),    pointer :: ja(:)
    integer(ip), intent(inout), pointer :: iskyl(:)
    integer(ip), intent(inout), pointer :: idiag(:)
    integer(ip)                         :: kskyl,idof,jdof
    integer(ip)                         :: izdom,ii,jj,kk,ll
    !
    ! Allocate memory
    !
    if( .not. associated(iskyl) ) then
       call memory_alloca(alya_CSR_memor,'ISKYL','alya_cholesky_initialization',iskyl,ndof*nbrows+1)
    end if
    !
    ! Skyline format
    !
    do ii = 1,nbrows*ndof+1
       iskyl(ii) = nbrows*ndof
    end do
    do ii = 1,nbrows
       do izdom = ia(ii),ia(ii+1)-1
          jj = ja(izdom)  
          if( jj > 0 .and. ii >= jj ) then
             do idof = 1,ndof 
                kk = (ii-1)*ndof+idof+1  
                do jdof = 1,ndof
                   ll = (jj-1)*ndof+jdof
                   if( ll < iskyl(kk) ) iskyl(kk) = ll
                end do
             end do
          end if
       end do
    end do

    !call PAR_MIN(ndof*nbrows+1_ip,iskyl)

    nskyl    = 1
    iskyl(1) = 1

    if( kfl_symmetric == 1 ) then
       ! 
       ! For the symmetric case, do not need idiag
       !
       do kk = 1,nbrows*ndof
          kskyl       = kk - iskyl(kk+1) + 1
          nskyl       = nskyl + kskyl
          iskyl(kk+1) = nskyl
       end do

    else
       !
       ! For the nonsymmetric case, set idiag 
       !
       if( .not. associated(idiag) ) then
          call memory_alloca(alya_CSR_memor,'IDIAG','alya_cholesky_initialization',idiag,ndof*nbrows)
       end if
       do kk = 1,nbrows*ndof
          kskyl       = kk - iskyl(kk+1)
          idiag(kk)   = nskyl + kskyl
          kskyl       = 2 * kskyl + 1  
          nskyl       = nskyl + kskyl
          iskyl(kk+1) = nskyl
       end do

    end if
    nskyl = nskyl - 1 

  end subroutine alya_cholesky_initialization

  subroutine alya_cholesky_factorization( N, NNZ, NZ, A, INFO )

    !---------------------------------------------------------------------------------------------------------------------------
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  chol_fact computes the Cholesky factorization of a real symmetric
    !  positive definite matrix A stored in skyline format.
    !
    !  The factorization has the form
    !     A = LL^T
    !  where L is a matrix stored in skyline format.
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.
    !
    !  NNZ     (input) INTEGER
    !          The number of non-zeros of the matrix A.  NNZ >= 0
    !
    !  NZ      (input) INTEGER array, dimension (N+1)
    !          The non-zeros of each row/column.
    !
    !  A       (input/output) DOUBLE PRECISION array, dimension (NNZ)
    !          On entry, the symmetri! matrix A stored in skyline format.  
    !
    !          On exit, if INFO = 0, the factor C from the Cholesky
    !          factorization A = C*C' stored in skyline format.
    !
    !  INFO    (output) INTEGER
    !          = 0: successful exit
    !          < 0: if INFO = -k, the k-th argument had an illegal value
    !          > 0: if INFO = k, the leading minor of order k is not
    !               positive definite, and the factorization could not be
    !               completed.
    !
    !---------------------------------------------------------------------------------------------------------------------------
    !
    integer(ip), intent(in)    :: N
    integer(ip), intent(in)    :: NNZ
    integer(ip), intent(in)    :: NZ( N+1 )
    real(rp),    intent(inout) :: A( NNZ )
    integer(ip), intent(out)   :: INFO
    INTEGER(ip)                :: I, J,  K
    INTEGER(ip)                :: I0, J0
    INTEGER(ip)                :: IPOS,IDIF,IROW,JCOL
    real(rp)                   :: TEMP, RDIAG
    !
    ! Test the input parameters
    !
    INFO = 0
    if( N<0 ) then
       write(*,*)'N < 0' 
       INFO = -1
    else if( NNZ<0 ) then
       write(*,*)'NNZ < 0' 
       INFO = -2
    end if
    if( INFO.NE.0 ) then
       return
    end if
    !
    !     Quick return if possible.
    !
    if( N ==0  .or. NNZ == 0 )  return
    !
    !     Compute the Cholesky factorization A = LL^T.
    !
    !
    !     Initialize pointer
    !
    IPOS = 2
    !
    !     Compute first pivot
    !
    if(A(1)<0.0_rp)then
       write(*,*) 'first pivot negative'
       INFO=1
       return
    endif
    ! 
    A(1) = SQRT ( A(1))
    !    
    do  I = 2, N
       !
       !     -----SET TO 0 DIAGONAL TERM
       !
       RDIAG = 0.0_rp
       !
       !     -----COLUMN NUMBER OF THE FIRST NON ZERO IN ROW I
       !
       JCOL=I-(NZ(I+1)-NZ(I))+1
       !     
       !     -----Compute elements JCOL:(I-1) of line I.
       !     
       do J = JCOL , I-1
          !
          !     -----ROW NUMBER OF FIRST NON ZERO IN COLUMN J
          !
          IROW=J-(NZ(J+1)-NZ(J))+1
          !
          !     -----Check for non zero bounds to setup the pointers for the scalar product
          !
          if(JCOL>IROW) then
             IDIF = J-JCOL
             I0 = NZ(I)
             J0 = NZ(J+1)-IDIF-1
          else
             IDIF = J-IROW
             I0 = IPOS-IDIF
             J0 = NZ(J)
          endif
          !
          !     -----compute scalar product
          !
          TEMP = 0.0_rp
          do K= 1,IDIF
             TEMP = TEMP + A( I0 )*A( J0 )
             I0 = I0 + 1
             J0 = J0 + 1
          end do
          !
          A(IPOS)=(A(IPOS)-TEMP)/A(NZ(J+1)-1)
          !
          !     -----ACCUMULATE IN RDIAG
          !
          RDIAG=RDIAG+A(IPOS)*A(IPOS)
          !
          !     -----MOVE POINTER
          !
          IPOS=IPOS+1
          !
       end do
       !
       !     -----COMPUTE THE DIAGONAL
       !

       A(IPOS) = A(IPOS)-RDIAG
       if(A(IPOS)<0.0_rp)then
          write(*,*)'Pivot negative in line ',I
          INFO=I
          return
       endif
       !
       A(IPOS)=SQRT(A(IPOS)) 
       IPOS=IPOS+1
       !
    end do

  end subroutine alya_cholesky_factorization

  subroutine alya_cholesky_solution( N, NNZ, NZ, NRHS, A, B, LDB, INFO )
    !     ..
    !
    !  Purpose
    !  =======
    !!  cho_solve solves a system of linear equations A*X = B with a symmetric
    !  positive definite matrix A stored in skyline format using the Cholesky
    !  factorization A = LL^T computed by chol_fact.
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.
    !
    !  NNZ     (input) INTEGER
    !          The number of non-zeros of the matrix A.  NNZ >= 0
    !
    !  NZ      (input) INTEGER array, dimension (N)
    !          The non-zeros of each row/column.
    !
    !  A       (input) DOUBLE PRECISION array, dimension (NNZ)
    !          The triangular factor C stored in skyline format from the 
    !          Cholesky factorization A = L*L**T, as computed by chol_fact.
    !
    !  NRHS    (input) INTEGER
    !          The number of right hand sides, i.e., the number of columns
    !          of the matrix B.  NRHS >= 0.
    !
    !  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
    !          On entry, the right hand side matrix B.
    !          On exit, the solution matrix X.
    !
    !  LDB     (input) INTEGER
    !          The leading dimension of the array B.  LDB >= max(1,N).
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !
    !  =====================================================================
    !
    INTEGER(ip)          :: INFO, LDB, N, NNZ, NRHS
    INTEGER(ip)          :: NZ( N+1 )
    real(rp)             :: A( NNZ ), B( LDB, * )
    INTEGER(ip)          :: J, K
    INTEGER(ip)          :: K0
    INTEGER(ip)          :: JCNT, I
    real(rp)             :: TEMP
    !     ..
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    if( N<0 ) then
       write(*,*)'N <0' 
       INFO = -1
    else if( NNZ<0 ) then
       write(*,*)'NNZ <0' 
       INFO = -2
    else if( NRHS<0 ) then
       write(*,*)'NRHS <0' 
       INFO = -4
    else if( LDB<MAX( 1_ip, N ) ) then
       write(*,*)'LDB < MAX(1,N)' 
       INFO = -7
    end if
    if( INFO.NE.0 ) then
       return
    end if
    !
    !     Quick return if possible.
    !
    if( N==0 .OR. NNZ==0 .OR. NRHS==0 )  return
    !
    !     Loop over NRHS
    !
    do  I = 1, NRHS
       !     
       !     Forward substitution.
       !    
       B(1,I)=B(1,I)/A(1)
       JCNT = 2 
       !
       do  J = 2, N
          !     
          K0 = J - (NZ(J+1)-JCNT) + 1
          !     
          TEMP = B(J,I)
          do  K = K0, J-1
             !     
             TEMP = TEMP - A(JCNT)*B(K,I)
             JCNT = JCNT + 1
             !     
          enddo
          !     
          B(J,I) = TEMP /A(JCNT)
          JCNT=JCNT+1 
          !     
       enddo
       !     
       !     Backward substitution.
       !     
       do  J = N, 1, -1
          !     
          JCNT=JCNT-1
          TEMP = B(J,I) / A(JCNT)
          B(J,I) = TEMP
          !     
          K0 = J - (NZ(J+1)-NZ(J)) + 1
          !
          do  K =  J-1, K0, -1
             !     
             JCNT = JCNT - 1
             B(K,I) = B(K,I) - A(JCNT)*TEMP
             !     
          enddo
          !     
       enddo
       !     
    enddo
    !
    return
  end subroutine alya_cholesky_solution

end module mod_alya_direct_solver
!> @} 
