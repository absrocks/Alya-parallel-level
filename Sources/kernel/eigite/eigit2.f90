subroutine eigit2(amatr,eigva,eigen,bmatr,iter)

  !-----------------------------------------------------------------------
  !
  ! This routine solves eigenvalues systems using the iterative method of 
  ! CG. It used sparce system 
  ! (It is assumed that the mesh graph has been properly renumbered.I repeat 
  ! this Fa.. phrase but really I don't know what means)      
  !
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
  integer(ip)               :: indic,keig,NPASO,kk,jj,ii,nbnodes,in1,in2,in3,in4,IM
  integer(4)                :: istat
  real(rp)                  :: eigdot
  real(rp)                  :: ABSTOL,ALANDA,sigma,YMAX,DENOM,NUMER,NORM
  real(rp)                  :: ERROR,DUMP,BETA,MU,CTHETA,TITA,STHETA,BETA0

  real(rp),   allocatable   :: X(:),Y(:),RR(:),DD(:),DX(:),AX(:),DDP(:),PRR(:),AD(:)
  real(rp),   allocatable   :: DT(:), X1(:),AX1(:),DDP1(:),dx1(:)

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
  call memchk(0_ip,istat,memit,'X','eigit2',x)
  allocate(AD(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'AD','eigit2',AD)
  allocate(Y(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'Y','eigit2',y)
  allocate(RR(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'RR','eigit2',rr)
  allocate(PRR(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'PRR','eigit2',prr)
  allocate(DD(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DD','eigit2',dd)
  allocate(DX(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DX','eigit2',dx)


  allocate(AX(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'AX','eigit2',ax)

  allocate(X1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'X1','eigit2',x1)
  allocate(AX1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'AX1','eigit2',ax1)
  allocate(DDP1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DDP1','eigit2',ddp1)
  allocate(DX1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DX1','eigit2',dx1)


  allocate(DDP(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DDP','eigit2',ddp)

  allocate(DT(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'DT','eigit2',DT)

  if( IMASTER ) nbnodes  = 0

  if( eigen_sol(1)%kfl_massm == 0 ) then
     ! matriz diagonal de masa
     ! calculo inv(sqrt(D)) 
     if( solve_sol(1) % kfl_symme == 1 ) then
        DO II=1,nbnodes
           DT(II)=sqrt(1.0_rp/bmatr(ii))
           indic = r_sym(ii+1)-1
           AD(II)=amatr(indic)/bmatr(ii)
        ENDDO
     else
       DO II=1,nbnodes
           DT(II)=sqrt(1.0_rp/bmatr(ii))
           indic = r_dom(ii)
           do while( c_dom(indic) /= II )
              indic = indic + 1
           end do
           AD(II)=amatr(indic)/bmatr(ii)
        ENDDO        
     end if
  else

     call runend('BAD EIGENSOLVER')

  end if



  !    OPEN(UNIT=111, FILE='a_mat_alya.dat',STATUS='UNKNOWN')
  !       write(111,*) nbnodes       
  !         do kk=1,nbnodes
  !	  indic = r_sym(kk+1)-1
  !      write(111,*) amatr(indic)
  !	   enddo

  !      do kk=1,2000
  !	      write(111,*) kk,amatr(kk)
  !	   enddo
  !	 close(111)
  !	 OPEN(UNIT=111, FILE='b_mat_alya.dat',STATUS='UNKNOWN')
  !      write(111,*) nbnodes    
  !       do kk=1,nbnodes
  !	      write(111,*) bmatr(kk)
  !	   enddo
  !	 close(111)


  DO IM = 1,eigen_sol(1)%neiva

     ! VALOR INICIAL AL AUTOVALOR QUE QUEIRO
     DO KK = 1,nbnodes
        if(iter==1) then  
           eigen(KK,IM) = 1.0_rp   ! LO HAGO ORTONORMAL 
        endif
        DD(KK)       = 0.0_rp
     ENDDO

     ! ORTOGONALIZO EL AUTOVECTOR "NUEVO"
     !IF( IM > 1 )  THEN
     !   CALL ORTOGONALIZA(nbnodes,eigen,IM)
     !ENDIF

     call norm2x(1_ip,eigen(1,IM),denom)
     if( denom == 0.0_rp ) then
        call runend('ZERO NORM')
     else
        denom = 1.0_rp/denom
     end if

     DO KK = 1,nbnodes
        X(KK) = eigen(KK,IM)*denom
     ENDDO

     do kk = 1,nbnodes
        x1(kk) = dt(kk)*x(kk)
     enddo

     if( solve_sol(1)%kfl_symme == 1 ) then
        call bsymax(1_ip,nbnodes,1_ip,amatr,c_sym,r_sym,x1,ax1)
     else
        call bcsrax(1_ip,nbnodes,1_ip,amatr,c_dom,r_dom,x1,ax1)
     end if

     do kk=1,nbnodes
        ax(kk)=dt(kk)*ax1(kk)
     enddo

     ERROR = 1.0_rp
     MU    = 1.0_rp
     NPASO = 0

     DO WHILE( ERROR > eigen_sol(1)%solco .AND. NPASO < eigen_sol(1)%miter )

        NPASO = NPASO+1
        call prodxy(1_ip,nbnodes,x,ax,alanda)

        DO KK = 1,nbnodes
           RR(KK) = ALANDA * X(KK) - AX(KK)
        END DO

        DO KK=1,nbnodes
           PRR(KK)=RR(KK)/AD(KK)
        ENDDO

        ! 4.3. Get steepest descent vector. EN QDOT USA LA ORTOGONALIZACION!!
        !       sd = lambda*psi - phi   SD == RR

        DO KK = 1, IM-1
           DO JJ = 1,nbnodes
              Y(JJ) = eigen(JJ,KK)
           END DO
           !           call prodxy(1_ip,nbnodes,y,rr,dump)
           call prodxy(1_ip,nbnodes,y,prr,dump)
           DO JJ = 1,nbnodes
              !              RR(JJ) = RR(JJ) - dump*eigen(JJ,KK)
              PRR(JJ) = PRR(JJ) - dump*eigen(JJ,KK)
           END DO
        enddo

        !        call prodxy(1_ip,nbnodes,rr,rr,beta0)
        call prodxy(1_ip,nbnodes,rr,prr,beta0)
        BETA = BETA0/MU
        MU   = BETA0

        DO KK=1,nbnodes
           !           DD(KK) = RR(KK) + BETA*DD(KK)
           DD(KK) = PRR(KK) + BETA*DD(KK)
        END DO

        call prodxy(1_ip,nbnodes,x,dd,dump)
        DO KK = 1,nbnodes
           DDP(KK) = DD(KK) - DUMP*X(KK)
        END DO
        call norm2x(1_ip,ddp,dump)
        dump = 1.0_rp/dump

        DO KK = 1,nbnodes
           DDP(KK)  = DDP(KK) * DUMP
           ddp1(kk) = dt(kk)  * ddp(kk)
        END DO

        if( solve_sol(1)%kfl_symme == 1 ) then        
           call bsymax(1_ip,nbnodes,1_ip,amatr,c_sym,r_sym,DDP1,DX1)
        else
           call bcsrax(1_ip,nbnodes,1_ip,amatr,c_dom,r_dom,DDP1,DX1)
        end if

        do kk=1,nbnodes
           dx(kk)=dt(kk)*dx1(kk)
        end do

        call prxyxz(1_ip,ddp,ax,dx,numer,denom)
        numer = 2.0_rp*numer
        denom = -alanda + denom

        TITA = 0.5_rp * ABS(ATAN(-NUMER/DENOM))

        CTHETA = COS(TITA)
        STHETA = SIN(TITA)

        DO KK = 1,nbnodes
           X(KK)  = CTHETA * X(KK)  + STHETA * DDP(KK)
           AX(KK) = CTHETA * AX(KK) + STHETA * DX(KK)
        END DO

        call proaxy(1_ip,ALANDA,AX,X,ERROR)

     END DO

     DO KK = 1,nbnodes
        eigen(KK,IM) = X(KK)
     END DO
     eigva(IM) = ALANDA

     if( eigen_sol(1)%lun_solei /= 0 .and. INOTSLAVE ) then
        write(6,*) 'autovalor: ', IM,ALANDA,ERROR,NPASO 
        write(eigen_sol(1)%lun_solei,*) im,ALANDA
     end if


  END DO

  ! re-scaleo los autovectores

  do IM = 1,eigen_sol(1)%neiva
     DO KK = 1,nbnodes
        eigen(KK,IM) = eigen(KK,IM)*dt(KK)
     ENDDO
  enddo

  deALLOCATE(X,Y,RR,DD,DX,AX,DDP,dt,x1,ax1,ddp1,dx1,PRR,AD)
end subroutine eigit2
