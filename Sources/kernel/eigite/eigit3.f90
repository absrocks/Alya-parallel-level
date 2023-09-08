subroutine eigit3(an,eigva,eigen,bmatr,iter)

  !-----------------------------------------------------------------------
  !
  ! This routine solves eigenvalues systems using the iterative method of 
  ! CG. It used sparce system 
  ! (It is assumed that the mesh graph has been properly renumbered.I repeat 
  ! this Fa.. phrase but really I don't know what means)      
  !
  !-----------------------------------------------------------------------
  use def_solver
  use def_master, only       :  IMASTER,INOTSLAVE,INOTMASTER,ISLAVE
  use def_master, only       :  npari,parin,npoi1,npoi2,npoi3,kfl_paral
  use def_domain, only       :  c_sym,r_sym,npoin,r_dom,c_dom,nzdom
  use mod_memchk
  implicit none
  real(rp),    intent(inout) :: an(*)
  real(rp),    intent(in)    :: bmatr(*)
  real(rp),    intent(out)   :: eigen(npoin,*)
  real(rp),    intent(out)   :: eigva(*)
  integer(ip), intent(in)    :: iter
  integer(ip)                :: indic,keig,npaso,kk,jj,ii,nbnodes,in1
  integer(ip)                :: in2,in3,in4,im,npos,nestates,Ninterna
  integer(ip)                :: interna,nsuma,PASAJE,im2,ipoin
  integer(4)                 :: istat
  real(rp)                   :: eigdot
  real(rp)                   :: ABSTOL,ALANDA,sigma,YMAx,denom,numer
  real(rp)                   :: ErrOR,DUMP,beta,mu,CTHETA,TITA,STHETA
  real(rp)                   :: beta0,dot_p,NORM

  real(rp),    allocatable   :: x(:),Y(:),rr(:),dd(:),dx(:),Ax(:)
  real(rp),    allocatable   :: ddp(:),prr(:),ad(:)
  real(rp),    allocatable   :: dt(:), x1(:),Ax1(:),ddp1(:),dx1(:)
  real(rp),    allocatable   :: rrp(:), px(:)
  real(rp),    allocatable   :: H_esp(:),lan_esp(:),eig_esp(:,:),w_esp(:,:)

  integer(ip), allocatable   :: napunte(:)
  integer(ip), target        :: npoint(1)

  !----------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !----------------------------------------------------------------------

  if( imASTER ) then
     nbnodes = 1
  else
     nbnodes = npoin
  end if

  nestates = eigen_sol(1) % neiva

  allocate(x(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'x','eigit3',x)
  allocate(ad(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'ad','eigit3',ad)
  allocate(Y(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'Y','eigit3',y)
  allocate(rr(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'rr','eigit3',rr)
  allocate(prr(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'prr','eigit3',prr)
  allocate(dd(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'dd','eigit3',dd)
  allocate(dx(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'dx','eigit3',dx)

  allocate(Ax(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'Ax','eigit3',ax)

  allocate(x1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'x1','eigit3',x1)
  allocate(Ax1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'Ax1','eigit3',ax1)
  allocate(ddp1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'ddp1','eigit3',ddp1)
  allocate(dx1(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'dx1','eigit3',dx1)

  allocate(ddp(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'ddp','eigit3',ddp)

  allocate(dt(nbnodes),stat=istat)
  call memchk(0_ip,istat,memit,'dt','eigit3',dt)

  allocate(Napunte(nestates),stat=istat)
  call memchk(0_ip,istat,memit,'Napunte','eigit3',Napunte)

  allocate(H_esp(nestates*(nestates+1)/2),stat=istat)
  call memchk(0_ip,istat,memit,'H_esp','eigit3',H_esp)
  allocate(lan_esp(nestates),stat=istat)
  call memchk(0_ip,istat,memit,'lan_esp','eigit3',lan_esp)
  allocate(eig_esp(nestates,nestates))
  call memchk(0_ip,istat,memit,'eig_esp','eigit3',eig_esp)
  allocate(w_esp(nbnodes,nestates),stat=istat)
  call memchk(0_ip,istat,memit,'w_esp','eigit3',w_esp)

  if( imASTER ) nbnodes  = 0

  if( eigen_sol(1) % kfl_massm /= 0 ) then
     call runend('Bad EIGENSOLVER')
  endif
  !
  ! DT = 1 / sqrt( M )
  ! AD = diag(A) / M
  !
  if( solve_sol(1) % kfl_symme == 1 ) then
     call diagon(nbnodes,1_ip,1_ip,0_ip,solve_sol(1) % kfl_full_rows,r_sym,c_sym,an,ad)
  else
     call diagon(nbnodes,1_ip,solve_sol(1) % kfl_symme,solve_sol(1) % kfl_full_rows,r_dom,c_dom,an,ad)
  end if

  do ii = 1,nbnodes
     !dt(ii) = 1.0_rp / sqrt( bmatr(ii) )
     dt(ii) = sqrt( 1.0_rp / bmatr(ii) )
     ad(ii) = ad(ii) / bmatr(ii)
  end do

  Ninterna = 100
  ErrOR    = 1.0_rp
  mu       = 1.0_rp 
  npaso    = 0
  pasaje   = 0

  do kk = 1,nestates
     Napunte(kk) = 0  
  end do 
  !
  ! Total number of nodes
  !
  npoint = npoin
  if( imASTER ) then
     npoint = 0
  else if( ISLAVE ) then
     npoint = npoi1 + npoi3 - npoi2 + 1
  end if
  call parari('SUM',0_ip,1_ip,npoint)

  do while( pasaje == 0 .and. npaso < eigen_sol(1) % miter )  
 
    npaso = npaso + 1

     do im = 1,nestates
        !
        ! VALOR INICIAL AL AUTOVALOR QUE QUIERO
        !
        if( npaso == 1 ) then
           do kk = 1,nbnodes
              eigen(kk,im) = 1.0_rp / real(npoint(1),rp)  ! LO HAGO ORTONORMAL 
              dd(kk)       = 0.0_rp
           end do
           eigva(im)       = 1.0_rp
        else
           do kk = 1,nbnodes
              eigen(kk,im) = w_esp(kk,im)
              dd(kk)       = 0.0_rp
           end do
        endif
        !
        ! DOT_P = 1 / sqrt( EIGEN(:,IM)^2 )
        !
        call norm2x(1_ip,eigen(1,im),dot_p)
        if( dot_p == 0.0_rp ) then
           call runend('ZERO NORM')
        else
           dot_p = 1.0_rp / dot_p
        end if
        !
        ! X  = DOT_P * EIGEN(:,IM)
        ! X1 = X . M^{-1/2}
        !
        do kk=1,nbnodes
           x(kk) = eigen(kk,im) * dot_p
        end do

        do kk = 1,nbnodes
           x1(kk) = dt(kk) * x(kk)
        end do
        !
        ! AX1 = A * X1
        !
        if( solve_sol(1) % kfl_symme == 1 ) then        
           call bsymax(1_ip,nbnodes,1_ip,an,c_sym,r_sym,x1,Ax1)
        else
           call bcsrax(1_ip,nbnodes,1_ip,an,c_dom,r_dom,x1,Ax1)
        end if

        do kk = 1,nbnodes
           ax(kk) = dt(kk) * ax1(kk)
        end do
        !
        ! ERROR = ( AX - lambda * X )^2 / XX^2
        !
        call proaxy(1_ip,eigva(im),Ax,x,ErrOR)

        if( error > eigen_sol(1) % solco ) then

           mu  = 1.0_rp

           do interna = 1,Ninterna
              !
              ! ALANDA = X . AX
              !
              call prodxy(1_ip,nbnodes,x,ax,alanda)
              !
              ! RR  = ALANDA * X - AX
              ! PRR = RR / ( diag(A) / M )
              !
              do kk = 1,nbnodes
                 rr(kk)  = alanda * x(kk)-Ax(kk)
                 prr(kk) = rr(kk) / ad(kk)
              end do
              !
              ! Get steepest descent vector. EN QdoT USA LA ORTOGONALIZACION!!
              ! sd = lambda*psi - phi   SD == rr
              !
              do kk = 1,im-1
                 call prodxy(1_ip,nbnodes,eigen(1,kk),prr,dump)
                 !
                 ! PRR(:) = PRR(:) - ( EIGEN(:,KK) . PRR^t ) * EIGEN(:,KK)
                 !
                 do jj = 1,nbnodes
                    prr(jj) = prr(jj) - dump * eigen(JJ,kk)
                 end do
              end do
              !
              ! BETA0= RR . PRR^t
              !
              call prodxy(1_ip,nbnodes,rr,prr,beta0)

              beta = beta0 / mu
              mu   = beta0

              do kk = 1,nbnodes
                 dd(kk) = prr(kk) + beta * dd(kk)
              end do
              !
              ! DUMP    = X . DD^t
              ! DDP(:)  = DD(:) - DUMP * X(:)
              !
              call prodxy(1_ip,nbnodes,x,dd,dump)
              do kk = 1,nbnodes
                 ddp(kk) = dd(kk) - dump * x(kk)
              end do
              !
              ! DUMP = 1 / sqrt( DPP^2 )
              !
              call norm2x(1_ip,ddp,dump)
              if( dump /= 0.0_rp ) dump = 1.0_rp / dump
              !
              ! DDP(:)  = DUMP * DDP(:)
              ! DDP1(:) = DDP(:) . M(:)^{-1/2}
              !
              do kk = 1,nbnodes
                 ddp(kk)  = ddp(kk) * dump
                 ddp1(kk) = dt(kk)  * ddp(kk)
              end do
              !
              ! DX1 = A . DDP1
              ! DX  = ( A . DDP1 ) . M(:)^{-1/2}
              !
              if( solve_sol(1) % kfl_symme == 1 ) then        
                 call bsymax(1_ip,nbnodes,1_ip,an,c_sym,r_sym,ddp1,dx1)
              else
                 call bcsrax(1_ip,nbnodes,1_ip,an,c_dom,r_dom,ddp1,dx1)
              end if

              do kk = 1,nbnodes
                 dx(kk) = dt(kk) * dx1(kk)
              end do

              call prxyxz(1_ip,ddp,ax,dx,numer,denom)
              numer =  2.0_rp * numer
              denom = -alanda + denom

              TITA = 0.5_rp * abs(atan(-numer/denom))

              CTHETA = cos(TITA)
              STHETA = sin(TITA)

              do kk = 1,nbnodes
                 x(kk)  = CTHETA * x(kk)  + STHETA * ddp(kk)
                 Ax(kk) = CTHETA * Ax(kk) + STHETA * dx(kk)
              end do

           end do

           do kk = 1,nbnodes
              eigen(kk,im) = x(kk)
           end do
           eigva(im) = alanda

           if( INOTSLAVE ) WRITE(6,*)  'AUTOVALOR : ',im,ALANDA,ErrOR 
           if( INOTSLAVE ) WRITE(12,*) 'AUTOVALOR : ',npaso,im,ALANDA,ErrOR 

        else

           Napunte(im) = 1

           if( INOTSLAVE ) WRITE(6,*)  'convergido ',im,eigva(im),Napunte(im)
           if( INOTSLAVE ) WRITE(12,*) 'convergido ',npaso,im,eigva(im),Napunte(im)

        endif

        nsuma = 0
        do im2 = 1,nestates
           nsuma = nsuma + Napunte(im2)
        end do
        if( nsuma == nestates ) pasaje = 1

     end do
     !
     !  RAYLEIGGH-RITS ON SPAN[x] 
     !
     call ortogeig3(nbnodes,nestates,napunte,eigen)

     do im = 1,nestates  

        !       if(napunte(im).eq.0) then
        !              im2=im2+1
        !         nvalor(im2)=im
        do kk=1,nbnodes
           x1(kk) = dt(kk)  * eigen(kk,im)
           !     x1(kk) = dt(kk)  * eigen(kk,nvalor(im2))
        end do
        !
        !     CALL MATxVECSim(IA,JA,An,ad,x1,ax1, np)
        !
        if( solve_sol(1) % kfl_symme == 1 ) then        
           call bsymax(1_ip,nbnodes,1_ip,an,c_sym,r_sym,x1,ax1)
        else
           call bcsrax(1_ip,nbnodes,1_ip,an,c_dom,r_dom,x1,ax1)
        end if

        do kk = 1,nbnodes
           ax(kk) = dt(kk) * ax1(kk)
        end do

        call prodxy(1_ip,nbnodes,eigen(1,im),Ax,dump)

        H_esp(im+(im-1)*im/2) = dump  !  fundot(np,eigen(:,im),Ax)

        do kk=1,im-1
           npos = kk+(im-1_ip)*im/2_ip  
           !     do kk=1,im2-1
           !     npos = kk+(im2-1)*im2/2  
           call prodxy(1_ip,nbnodes,eigen(1,kk),Ax,dump)
           H_esp(npos) = dump    ! fundot(np,eigen(:,kk),Ax)
           !     H_esp(npos)=fundot(NP,eigen(:,nvalor(kk)),Ax)
        end do

     end do
     call eigensystem(nestates,H_esp,lan_esp,eig_esp) 

     call ortog_minimo(nestates,napunte,eig_esp)

     !   ESTO        

     do kk = 1,nestates
        eigva(kk) = lan_esp(kk)
        do jj = 1,nbnodes
           w_esp(jj,kk) = 0.0_rp
           do ii = 1,nestates
              w_esp(jj,kk) = w_esp(jj,kk) + eigen(jj,ii) * eig_esp(ii,kk)
           end do
        end do
     end do
 
  end do
  !
  ! Re-scaleo los autovectores
  !
  do im = 1,Nestates
     do kk = 1,nbnodes
        eigen(kk,im) = eigen(kk,im) * dt(kk)
     end do
  end do

  DEALLOCATE(x,Y,rr,prr,dd,dx,Ax,ddp,ad,dt,x1,Ax1,ddp1,dx1,H_esp,lan_esp,eig_esp,w_esp,napunte)

END subroutine eigit3

subroutine eigensystem(nmin,H_esp,lan_esp,eig_esp) 
  use def_kintyp, only : ip,rp
  imPLICIT NONE
  integer(ip)  :: nmin
  real(rp)     :: H_esp(*),lan_esp(nmin),eig_esp(nmin,nmin)
  CHARACTER(1) :: JOBZ, UPLO
  INTEGER(ip)  :: INFO
  real(rp)     :: WORK( 3*nmin )

  JOBZ = 'V'
  UPLO = 'U'
call runend('C EST HONTEUX')
  !call DSPEV( JOBZ, UPLO, Nmin, H_esp, lan_esp,eig_esp, nmin, WORK , INFO )

  !*  INFO    (output) INTEGER
  !*          = 0:  successful exit.
  !*          < 0:  if INFO = -i, the i-th argument had an illegal value.
  !*          > 0:  if INFO = i, the algorithm failed to converge; i
  !*                off-diagonal elements of an intermediate tridiagonal
  !*                form did not converge to zero.
  !*
  if(info.ne.0) then
     write(6,*) 'error  en dspev ',info
     stop ' '
  endif

end subroutine eigensystem

SUBROUTINE  ortog_minimo(MM,napunte,autovec) 
  use def_kintyp, only : ip,rp
  implicit none
  INTEGER(ip)             :: MM,ii,jj,kk
  INTEGER(ip)             :: NAPUNTE(mm)
  real(rp)                :: autovec(mm,mm),NORMA1,NORMA2,xmax
  real(rp),   ALLOCATABLE :: x(:)

  ALLOCATE(x(MM))

  do ii = 2,mm

     ! IF(NAPUNTE(ii).EQ.0) THEN
     do JJ = 1,MM
        x(JJ) = autovec(JJ,ii)
     end do

     do kk = 1,ii-1
        NORMA1 = 0.0_rp
        NORMA2 = 0.0_rp
        do JJ = 1,MM
           NORMA1 = NORMA1 + autovec(JJ,kk) * autovec(JJ,kk)
           NORMA2 = NORMA2 + autovec(JJ,kk) * autovec(JJ,ii)
        end do
        xMAx = 0.0_rp
        if( NORMA1 /= 0.0_rp ) NORMA1 = 1.0_rp / NORMA1
        do JJ = 1,MM
           x(JJ) = x(JJ) - NORMA2 * autovec(JJ,kk) * NORMA1
           IF( xMAx < ABS(x(JJ)) ) xMAx = ABS(x(JJ)) 
        end do
     end do

     do JJ = 1,MM
        autovec(JJ,ii) = x(JJ)
     end do

  end do

  DEALLOCATE(x)

  RETURN
END SUBROUTINE ortog_minimo
