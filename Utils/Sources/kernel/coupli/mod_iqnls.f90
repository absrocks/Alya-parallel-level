  module mod_iqnls
    use def_kintyp

    implicit none

    private

    public :: compute_alpha

    contains
    !----------------------------------------------------
    !>
    !> @author  Alfonso Santiago
    !> @date    15/06/2016
    !> @brief   Parallel implementation of a portion of the interface quasi newton algorithm
    !> @details 
    !>
    ! This subrotuine obtains the alpha vector for the interface
    ! quasi newton algoritm.
    ! With this alpha:
    !
    !     x^k+1 = x^k + W*alpha
    !
    ! As alpha =f(V,r) tricks can be made in the inside
    ! to avoid buling some matrices and hereafter generate
    ! a more efficient parallel code.
    !
    ! This algoritm is based in the SEQUENTIAL
    ! economy QR decomposition.
    !
    ! The overall algorithm in this funciton reads:
    !
    !     1. Decompose matrix A in a set of vectors
    !        v_i where:
    !         Q=(I - 2*v_1*v_1^T)*(I - 2*v_2*v_2^T)...
    !
    !     2. Obtain upper triangular matrix U where:
    !         Q*U=A
    ! 
    !     3. Compute Q^T*r, where r is the residue and
    !        Q=f(v1,v2,v3...) as a matrix vector loop
    !
    !     4. Backsubstitute R*alpha=Q^T*r (rhs previously
    !        computed)
    !
    !
    !                 n
    !            +--+--+--+
    !            |  |  |  |
    !            +--+--+--+
    !            |  |  |  |                 +-+-+-+
    !  A(m,n)=   +--+--+--+ m     alpha(n)= | | | |
    !            |  |  |  |                 +-+-+-+
    !            +--+--+--+
    !            |  |  |  |
    !            +--+--+--+
    !
    !
    !                n
    !            +--+--+--+       
    !            |  |  |  |                   n
    !            +--+--+--+               +--+--+--+  
    !            |  |  |  |               |  |  |  |
    !  Q(m,n)=   +--+--+--+ m             +--+--+--+  
    !            |  |  |  |      U(m,n)=  |  |  |  | n
    !            +--+--+--+               +--+--+--+  
    !            |  |  |  |               |  |  |  |
    !            +--+--+--+               +--+--+--+ 
    !             
    !
    !
    !----------------------------------------------------
    subroutine compute_alpha(A, residue, mmax, nmax, alpha)


      use mod_communications, only : PAR_SUM
      use mod_communications, only : PAR_MAX
      use mod_parall,         only : PAR_MY_CODE_RANK

      use mod_maths,          only : maths_backsu

      implicit none

      real(rp), intent (in)               :: A(:,:)
      real(rp), intent (in)               :: residue(:)
      real(rp), intent (out), target      :: alpha(:)
      integer(ip), intent(in)             :: mmax, nmax !! maximum values of the matrix to decompose

      real(rp), pointer                   :: palpha(:)

      real(rp)                            :: Q(mmax,nmax)
      real(rp), allocatable, target       :: R(:,:)
      real(rp), pointer                   :: pR(:,:)
      integer(ip)                         :: i_col, i, j, k
      real(rp)                            :: norm_column, re_aux
      real(rp)                            :: v(mmax,nmax)
      real(rp)                            :: A_aux(mmax,nmax)

      real(rp)                            :: vecaux(mmax)
      real(rp)                            :: Qt_times_r(mmax)

      integer(ip)                         :: leader_rank
      real(rp)                            :: maxdofs ! declarated as real to enforce compatibility with PAR_MAX

      allocate(R(nmax,nmax))
      pR => R
      palpha => alpha
      Q=0.0_rp

      A_aux=A 

      !! IDEAS PARALL
      !!
      !! Elejir el subdominio con mas nodos mojados como el 
      !! "subdominio lider", el subdominio que tiene la cabeza
      !! de la matriz a descomponer.
      !!

      !! TODO Al agregar el PAR_MAX hace que el loop se me vaya out of bounds

      ! Find leader rank: the subdomain that heads the decomposition
      !
      maxdofs=real(mmax,rp)
      call PAR_MAX(maxdofs,'IN CURRENT TARGET COLOR','EXCLUDE MASTER',leader_rank)
      !!
      !! END IDEAS PARALL

      columns: do i_col=1, nmax 
        ! a_i=A(i_col,i_col:mmax)
        
        v(:,i_col)=0.0_rp
        !
        ! If I'm the leader subdomain, I crop the vector
        ! If Im not the leader, I use the full vector
        !
        if(PAR_MY_CODE_RANK .eq. leader_rank) then
          v(i_col:mmax,i_col) = A_aux(i_col:mmax,i_col)
        else
          v(:,i_col) = A_aux(:,i_col)
        endif

        ! ||alpha|| = sqrt(sum(A(i)^2))
        !
        norm_column = 0.0_rp
        !
        ! If Im the leader use the cropped vector, If'm not
        ! use the full vector
        !
        if(PAR_MY_CODE_RANK .eq. leader_rank) then
          do i=i_col,mmax
              norm_column =norm_column + v(i,i_col) * v(i,i_col)
          enddo
        else
          do i=1,mmax
              norm_column =norm_column + v(i,i_col) * v(i,i_col)
          enddo
        endif
        
        call PAR_SUM(norm_column,'IN CURRENT TARGET COLOR')
        
        norm_column = sqrt(norm_column)

        ! u=a_i-alpha*e_1
        !
        !! IMPORTANT THING HERE. Some algoritmhs (i.e. matlab)
        !! uses a_i + alpha* e_i. This algorithm also converges
        !! but slightly different. The minus sign has been
        !! left here trivially.

        ! If I'm the leader, for sure I will have the i_col
        ! element. So If Im the leader do the substraction
        !
        if(PAR_MY_CODE_RANK .eq. leader_rank) then
          v(i_col,i_col) = v(i_col,i_col) - norm_column
        endif

        ! ||u|| = sqrt(sum(u(i)^2))
        !
        norm_column = 0.0_rp
        !
        ! If Im the leder use the cropped vector, 
        ! if not, use the full vector.
        !
        if(PAR_MY_CODE_RANK .eq. leader_rank) then
          do i=i_col,mmax
            norm_column =norm_column + v(i,i_col) * v(i,i_col)
          enddo
        else
          do i=1,mmax
            norm_column =norm_column + v(i,i_col) * v(i,i_col)
          enddo
        endif

        call PAR_SUM(norm_column,'IN CURRENT TARGET COLOR')
        
        norm_column = sqrt(norm_column)
        
        ! v=u/||u||
        !
        !
        ! If Im the leder use the cropped vector, 
        ! if not, use the full vector.
        !
        if(PAR_MY_CODE_RANK .eq. leader_rank) then
          do i=i_col,mmax
            v(i,i_col)=v(i,i_col)/norm_column   
          enddo
        else
          do i=1,mmax
            v(i,i_col)=v(i,i_col)/norm_column   
          enddo
        endif

        ! Q_i=I-2*v*v^T
        !
        ! Step not necesary, because Q_i is obtained in function Q_mat(v,i_ini,i,j)

        ! A_i+1=Q_i * A
        !
        if(.not.(i_col.eq.nmax))then
           call vQ_times_Qaux(v(:,i_col),i_col,A_aux)
        endif
        
      enddo columns


     
      ! Now we need the upper triangular matrix R
      !
      !    A=Q*R  -> R=Q^T*A
      !
      ! This means
      ! 
      ! Q^T*r=(Q_1*Q_2*Q_3*....Q_n)^T * A
      !
      !     = Q_n^T * ... Q_3^T * Q_2^T * Q_1^T * A
      !
      ! In that way we can multiply from right to left.

      ! First we do Q_1^T *A
      A_aux=A
      call vQ_times_Qaux(v(:,1), 1_ip, A_aux(:,:))

      ! Now we keep multipling 
      do i_col=2,nmax
        call vQ_times_Qaux(v(:,i_col), i_col, A_aux(:,:))
      enddo
      
      !! Copy the first nxn square matrix different from zero
      !
      R=0.0_rp
      !! IDEAS PARALL
      !! 
      !! este bucle solo lo realizara el lider porque nmax
      !! (el numero mas grande de rank de iteraciones)
      !! es mas chico que su numero de degrees of freedom.
      !! 
      !! END IDEAS PARALL

      !
      ! If I'm the master, I'll have the head of
      ! the matrix, and only I should have the R
      !
      if(PAR_MY_CODE_RANK .eq. leader_rank) then
        do i=1,nmax
          do j=1,nmax
            R(i,j)=A_aux(i,j)
          enddo
        enddo
      endif

      ! Now sum all the results of every subdomain
      !
      call PAR_SUM(pR,'IN CURRENT TARGET COLOR')





      ! Finally in the backsubstitution we are going to
      ! compute:
      !
      !      R*alpha = - Q^T * r
      !
      ! So first we are going to compute
      ! the RHS as:
      !
      ! Q^T*r=(Q_1*Q_2*Q_3*....Q_n)^T * r
      !
      !     = Q_n^T * ... Q_3^T * Q_2^T * Q_1^T * r
      !
      ! In that way we can multiply from right to left.

      ! First we do Q_1^T *r
      Qt_times_r =0.0_rp
      Qt_times_r = vQ_times_vector(v(:,1),1_ip,.true.,residue(:))


      ! Now we keep multipling 
      do i_col=2,nmax
          Qt_times_r = vQ_times_vector(v(:,i_col),i_col,.true.,Qt_times_r(:))
      enddo

      Qt_times_r = - Qt_times_r


      !
      ! Everyvody should have R and Qt_times_r in this point

      ! Now we call backsubstitution to obtain alpha
      !
      ! R*alpha = - Q^T * r

      alpha=0.0_rp
      !
      ! If I'm the master, Ive the good R and I've
      ! the Qt_times_r that is going to be backsu
      ! bstituted, I should do the backsu
      !
      if(PAR_MY_CODE_RANK .eq. leader_rank) then
        call maths_backsu(R(:,:),alpha(:), Qt_times_r,nmax)
      endif

      call PAR_SUM(palpha, 'IN CURRENT TARGET COLOR')

      contains
        !--------------------------------------
        ! Function that gives me the result of
        ! the matrix without having the matrix
        !--------------------------------------
        function Q_mat(v,i_ini,i,j) result(q_ij)   
            implicit none

            real(rp), intent(in)        :: v(:)
            integer(ip), intent(in)     :: i_ini,i,j
            real(rp)                    :: q_ij
           
            if (i .lt. i_ini) then
                if(i .eq. j) then
                    q_ij = 1.0_rp
                else
                    q_ij = 0.0_rp 
                endif
            else
                if(i .eq. j) then
                    q_ij = 1.0_rp - 2.0_rp*v(i)*v(j)
                else  
                    q_ij = -2.0_rp*v(i)*v(j)
                endif
            endif

        end function Q_mat
        !--------------------------------------
        ! Function that gives me the result of
        ! the matrix without having the matrix
        !--------------------------------------
        function obtain_Q_mat(v,i_ini,max_col) result(Q)   
            implicit none

            real(rp), intent(in)        :: v(:)
            integer(ip), intent(in)     :: i_ini, max_col
            integer(ip)                 :: i, j, m
            real(rp)                    :: Q(size(v,1),size(v,1))
          
            m=size(v)

            
            if(m.eq.0_ip) then
              Q=0.0_rp
              return
            endif

            if ( (max_col .gt. m)     .or. &
                 (max_col .lt. i_ini) ) then
                  call runend('mod_maths: problem with obtain_Q_mat. Wrong max column number')
            endif
            do i=1,m
                do j=1,max_col
                    Q(i,j) = Q_mat(v,i_ini,i,j) 
                enddo
            enddo

        end function obtain_Q_mat

        !--------------------------------------
        ! Function that gives me the result of
        ! the matrix without having the matrix
        !--------------------------------------
        subroutine vQ_times_Qaux(v,i_ini,Qaux)   
            implicit none

            real(rp), intent(in)          :: v(:)
            integer(ip), intent(in)       :: i_ini
            real(rp), intent(inout)       :: Qaux(:,:)
            real(rp), allocatable         :: mat_aux(:,:)
            real(rp), allocatable, target :: vec_aux(:)
            real(rp), pointer             :: pvec_aux(:)
            real(rp)                      :: re_aux
            integer(ip)                   :: mv,mQ,nQ
            integer(ip)                   :: i, j
          
            mv=size(v)
            mQ=size(Qaux,1)
            nQ=size(Qaux,2)
        
            if ( (mv .ne. mQ ) ) call runend('mod_maths: vQ_times_Qaux: wrong dimensions')

            allocate(mat_aux(mQ,nQ))

            allocate(vec_aux(nQ))
            pvec_aux => vec_aux

            ! vQ_times_Qaux is :
            ! 
            !  (I-2*v*v^T) * A
            !
            !  A - 2*v*v^T * A
            !
            ! So I can start with v^T * A
            !
            do i=1,nQ
              re_aux=0.0_rp
              do j=1,mQ
                !! IDEAS PARALL
                !!
                !! en este if habria que agregar " Y si soy el subdominio lider"
                !! ya que i_ini siempre sera menor que el rank de iteraciones
                !! guardadas por lo tanto unicamente correspondera a el lider.
                !! algo asi:
                !!
                !!  if ((j .lt. i_ini) .and. IAMLEADER) then
                !!
                !! END IDEAS PARALL

                if (j .lt. i_ini) then
                  re_aux = re_aux + 0.0_rp
                else
                  re_aux = re_aux + v(j)*Qaux(j,i) !TODO aca creo que es el segudno indice
                endif
              enddo
              vec_aux(i)=re_aux ! here v^T * A is stored
            enddo

            call PAR_SUM(pvec_aux,'IN CURRENT TARGET COLOR')
            
            ! Now we can multiply v * (v^T * A)
            ! and add the A - 2*(...) in the same
            ! line
            !
            do i=1,mv
              do j=1,nQ
                !! IDEAS PARALL
                !!
                !! En este if pasa lo mismo que en el anterior
                !! solo se hace lo de i_ini si soy el subdomino
                !! lider
                !! 
                !!  if ((j .lt. i_ini) .and. IAMLEADER) then
                !!
                !! END IDEAS PARALL
                if (j .lt. i_ini) then
                  mat_aux(i,j) =  Qaux(i,j) - 0.0_rp
                else
                  mat_aux(i,j) =  Qaux(i,j) - 2.0_rp * v(i) * vec_aux(j)
                endif
              enddo
            enddo


            Qaux=mat_aux

            deallocate(mat_aux)

        end subroutine vQ_times_Qaux


        !--------------------------------------
        ! Function that multiplies two Q
        ! special matrices
        !--------------------------------------
        function vQ_times_vQ(v1, i_ini_1, v2, i_ini_2) result(Q)
            implicit none

            real(rp), intent(in)        :: v1(:), v2(:)
            integer(ip), intent(in)     :: i_ini_1, i_ini_2
            real(rp)                    :: Q(size(v1,1),size(v1,1))
            real(rp)                    :: rea_aux,part_V1, part_v2
            integer(ip)                 :: m,i,j,k
            
            m=size(v1,1)
            Q=0.0_rp
            


        do i=1,m !swipe in rows
            do j=1,m
                rea_aux=0.0_rp
                do k=1,m !swipe in columns
                    rea_aux = rea_aux + Q_mat(v1, i_ini_1, i, k) * Q_mat(v2, i_ini_2, k, j)!contract columns
                enddo
                Q(i,j)=rea_aux
            enddo
        enddo


        end function vQ_times_vQ
        !--------------------------------------
        ! Function that multiplies a special
        ! matrix Q with a vector
        !--------------------------------------
        function vQ_times_vector(v1, i_ini_1, transp, vector) result(vec_res)
            implicit none

            real(rp), intent(in)          :: v1(:)
            integer(ip), intent(in)       :: i_ini_1
            real(rp), intent(in)          :: vector(:)
            logical, intent(in)           :: transp

            real(rp), allocatable, target :: vec_res(:)
            real(rp), pointer             :: pvec_res(:)

            real(rp)                      :: rea_aux
            integer(ip)                   :: m,n
            integer(ip)                   :: i,j

            m=size(v1,1)
            n=size(vector,1)

            allocate(vec_res(n))
            pvec_res => vec_res

            vec_res=0.0_rp
          

            if((m .ne. n) .and. (m.ne.0_ip) ) then
              call runend('vq_times_vector: vector and matrix with different dimension')
            endif



            !
            ! This subroutine does:
            !
            !   (I - 2*v*v^T)*r
            !
            ! by computing:
            !
            !    r - 2*v*v^T*r
            !
            ! THe solution for the transpose vQ matrix is trivial as:
            !
            !   (I - 2*v*v^T)^T
            !
            !   I - 2 * (v^T)^T * v^T
            !
            !   I - 2 * v*v^T
            !
            ! So finally:
            !
            !   (I - 2*v*v^T)^T = I - 2*v*v^T
            !
            ! wich is understandable because it is symetric
            !

            !
            ! so we start with v^T * r
            !
            rea_aux=0.0_rp
            do i=1,m
              !! IDEAS PARALL
              !!
              !! En este if pasa lo mismo que en el anterior
              !! solo se hace lo de i_ini si soy el subdomino
              !! lider
              !! 
              !!  if ((j .lt. i_ini) .and. IAMLEADER) then
              !!
              !! END IDEAS PARALL
              if(i .lt. i_ini_1) then
                  rea_aux=rea_aux + 0.0_rp
              else
                  rea_aux=rea_aux + v1(i)*vector(i)
              endif
            enddo

            call PAR_SUM(rea_aux,'IN CURRENT TARGET COLOR')


            !
            ! Now we can mulyiply v * (v^T * r)
            ! and we can take profit of the loop
            ! and add r - 2*(...)
            !
            do i=1,m
              !! IDEAS PARALL
              !!
              !! En este if pasa lo mismo que en el anterior
              !! solo se hace lo de i_ini si soy el subdomino
              !! lider
              !! 
              !!  if ((j .lt. i_ini) .and. IAMLEADER) then
              !!
              !! END IDEAS PARALL
              if(i .lt. i_ini_1) then
                vec_res(i) = vector(i) - 0.0_rp
              else
                vec_res(i) = vector(i) - 2.0_rp * v1(i) * rea_aux 
              endif
            enddo


        end function vQ_times_vector
    end subroutine compute_alpha
  end module mod_iqnls
