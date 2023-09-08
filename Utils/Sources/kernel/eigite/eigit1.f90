subroutine eigit1(amatr,eigva,eigen,bmatr,iter)

  !-----------------------------------------------------------------------
  !
  ! This routine solves eigenvalues systems using the iterative method of 
  ! conjugated gradients. used sparce system but (allways there are a BUTT!)
  ! you need to shift the pencyl or there are not convergence
  ! It is assumed that the mesh graph has been properly renumbered.      
  !
  ! NOTE: Pencil --> A x = lambda B x  (ehh!!! do you know that? ;o) )
  !-----------------------------------------------------------------------
  use def_solver
  use def_master, only      :  IMASTER,INOTSLAVE,INOTMASTER,ISLAVE
  use def_domain, only      :  c_sym,r_sym,npoin,r_dom,c_dom
  use mod_memchk
  implicit none
  real(rp),   intent(inout) :: amatr(*)
  real(rp),   intent(in)    :: bmatr(*)
  real(rp),   intent(out)   :: eigen(npoin,*)
  real(rp),   intent(out)   :: eigva(*)
  integer(ip),intent(in)    :: iter
  integer(ip)               :: indic,keig,NPASO,kk,jj,ii,nbnodes
  integer(4)                :: istat
  real(rp)                  :: funmaxi,FUNMAXMIN,dummr,fmaxi
  real(rp)                  :: ABSTOL,ALANDA,sigma,YMAX
  real(rp),   allocatable   :: X(:),Y(:),R1(:),D(:),R(:)

  !----------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !----------------------------------------------------------------------

  if( IMASTER ) then
     nbnodes = 1
  else
     nbnodes = npoin
  end if

  allocate(X(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'X','eigit1',x)
  allocate(Y(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'Y','eigit1',y)
  allocate(R(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'R','eigit1',r)
  allocate(R1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'R1','eigit1',r1)
  allocate(D(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'D','eigit1',d)

  if( IMASTER ) nbnodes  = 0
  !
  ! Tell solver RHS for algebraic system is already global
  !
  solve_sol(1)%kfl_recov = 2 

  !----------------------------------------------------------------------
  !
  ! Initialization: A <= A - sig * B
  !
  !----------------------------------------------------------------------

  sigma = eigen_sol(1)%shift ! - 20.0_rp

  if( INOTMASTER ) then

     if( eigen_sol(1)%kfl_massm == 0 ) then
        do indic = 1,solve_sol(1)%nzmat
           amatr(indic) = amatr(indic) 
        end do 
        if( solve_sol(1) % kfl_symme == 1 ) then
           do ii = 1,nbnodes
              indic = r_sym(ii+1)-1
              amatr(indic) = amatr(indic) - sigma * bmatr(ii) ! bmatr should be local
              if( ISLAVE ) call runend('POPO')
           end do
        else
           do ii = 1,nbnodes
              indic = r_dom(ii)
              do while( c_dom(indic) /= II )
                 indic = indic + 1
              end do
              amatr(indic) = amatr(indic) - sigma * bmatr(ii) ! bmatr should be local
              if( ISLAVE ) call runend('POPO')
           end do           
        end if
     else
        do indic = 1,solve_sol(1)%nzmat
           amatr(indic) = amatr(indic) - sigma * bmatr(indic)
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Loop over eigenvalues
  !
  !----------------------------------------------------------------------

  do keig = 1,eigen_sol(1)%neiva

     ALANDA = 1.0_rp
     ABSTOL = 1.0_rp
     NPASO  = 0
     do KK = 1,nbnodes
        Y(KK)            = 0.0_rp
        X(KK)            = 1.0_rp
        D(KK)            = 1.0_rp
        if(iter==1)  eigen(KK,KEIG) = 1.0_rp
     end do

     if( keig > 1 ) then
        !
        ! Orthonormalize eigen vectors
        !
        call ORTOGGS_BIS(nbnodes,KEIG,eigen_sol(1)%kfl_massm,eigen,bmatr,c_sym,r_sym) 

        do JJ = 1,nbnodes
           Y(JJ) =  eigen(JJ,KEIG)
        end do

        fmaxi = funmaxi(nbnodes,eigen,KEIG)
        if( fmaxi /= 0.0_rp ) then
           fmaxi = 1.0_rp/fmaxi
           do JJ = 1,nbnodes
              X(JJ) = eigen(JJ,KEIG)*fmaxi
           end do
        end if

     end if

     !solve_sol(1)%solco = 1.0e-2_rp

     do while( ABSTOL > eigen_sol(1)%solco .and. NPASO < eigen_sol(1)%miter ) 

        NPASO = NPASO + 1
        !
        ! R = B.x -  ( A - sig*B ).y
        !
        if( solve_sol(1)%kfl_symme == 1 ) then
           call bsymax(1_ip,nbnodes,1_ip,amatr,c_sym,r_sym,Y,R)
        else
           call bcsrax(1_ip,nbnodes,1_ip,amatr,c_dom,r_dom,Y,R) 
        end if

        if( eigen_sol(1)%kfl_massm == 0 ) then
           do KK = 1,nbnodes
              R(KK) = bmatr(KK)*X(KK) - R(KK)  
           end do
        else
           if( solve_sol(1)%kfl_symme == 1 ) then
              call bsymax(1_ip,nbnodes,1_ip,bmatr,c_sym,r_sym,X,R1)
           else
              call bcsrax(1_ip,nbnodes,1_ip,bmatr,c_dom,r_dom,X,R1)              
           end if
           do KK = 1,nbnodes
              R(KK) = R1(KK) - R(KK) 
           end do
        end if
        !
        ! Solve ( A - sig*B ) . D = R
        !
        call solver(R,D,amatr,dummr)
        !
        ! y = y + d
        !
        do KK = 1,nbnodes
           Y(KK) = Y(KK) + D(KK)
        end do
        !
        ! Orthonormalize eigenvectors
        !
        if( keig > 1 ) then
           call orto_AM(nbnodes,KEIG,eigen_sol(1)%kfl_massm,eigen,Y,bmatr,c_sym,r_sym)
        endif       
        YMAX   = funmaxmin(npoin,Y)
        ALANDA = 1.0_rp/YMAX
        !
        ! Compute residual ABSTOL and update eigenvector
        !
        call eigup1(nbnodes,keig,ALANDA,y,x,abstol,eigen)

        !solve_sol(1)%solco = 0.1_rp*abstol

     end do

     call NORMALIZA_AM(nbnodes,KEIG,eigen_sol(1)%kfl_massm,eigen,bmatr,c_sym,r_sym)

     if( INOTMASTER ) eigva(keig) = ALANDA + sigma

     if( eigen_sol(1)%lun_solei /= 0 .and. INOTSLAVE ) then
        write(6,*) 'autovalor: ', keig,ALANDA + sigma,NPASO,ABSTOL
        write(eigen_sol(1)%lun_solei,10) keig,npaso,ABSTOL,ALANDA + sigma
     end if

  end do

  !----------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !----------------------------------------------------------------------

  call memchk(two,istat,memit,'D','eigit1',d)
  deallocate(d,stat=istat)
  if(istat/=0) call memerr(two,'D','',0_ip)
  call memchk(two,istat,memit,'R1','eigit1',r1)
  deallocate(r1,stat=istat)
  if(istat/=0) call memerr(two,'R1','',0_ip)
  call memchk(two,istat,memit,'R','eigit1',r)
  deallocate(r,stat=istat)
  if(istat/=0) call memerr(two,'R','',0_ip)
  call memchk(two,istat,memit,'Y','eigit1',y)
  deallocate(y,stat=istat)
  if(istat/=0) call memerr(two,'Y','',0_ip)
  call memchk(two,istat,memit,'X','eigit1',x)
  deallocate(x,stat=istat)
  if(istat/=0) call memerr(two,'X','',0_ip)

10 format(i10,i10,4(1x,e12.6))

end subroutine eigit1
