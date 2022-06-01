subroutine gendiapre(amatr,ia,ja,dia,ND,V)
  implicit none
  integer*4 :: ia(*),ja(*),ND,V
  real*8 :: amatr(V,V,*),dia(*)
  integer*4 :: i,j,k
  !$OMP PARALLEL  DO                                     &
  !$OMP SCHEDULE ( STATIC )                              &
  !$OMP DEFAULT  ( NONE )                                &
  !$OMP SHARED   ( amatr, ia, ja, ND,V,dia )  &
  !$OMP PRIVATE  ( i, j,k)           
  do i=1,ND
     do j=ia(i),ia(i+1)-1
        if(ja(j)==i) then
           do k=1,V
               dia((i-1)*V + k) = amatr(k,k,j)
           end do
        end if
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine gendiapre

subroutine gendiapreupdate(amatr,ia,ja,dia,ND,V)
  implicit none
  integer*4 :: ia(*),ja(*),ND,V
  real*8 :: amatr(V,V,*),dia(*)
  integer*4 :: i,j,k
  real*8 :: rzero= 1e-12
  !$OMP PARALLEL  DO                                     &
  !$OMP SCHEDULE ( STATIC )                              &
  !$OMP DEFAULT  ( NONE )                                &
  !$OMP SHARED   ( amatr, ia, ja, ND,V,dia )  &
  !$OMP PRIVATE  ( i, j,k)           
  do i=1,ND
     do j=ia(i),ia(i+1)-1
        if(ja(j)==i) then
           do k=1,V
              if(abs( dia((i-1)*V + k))<rzero) then
                 dia((i-1)*V + k) = 0.0d0 
              else
                 dia((i-1)*V + k) = 1.0d0/dia((i-1)*V + k) 
              end if
           end do
        end if
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine gendiapreupdate



subroutine linearcombocpu(x,v,y,NV,res)
  implicit none
  real*8    :: x(*),y(*),v(NV,*)
  integer*4 :: NV,res,i,j

  do i =1,res
     !$OMP PARALLEL  DO                                     &
     !$OMP SCHEDULE ( STATIC )                              &
     !$OMP DEFAULT  ( NONE )                                &
     !$OMP SHARED   ( i,x,y,v,res,NV)  &
     !$OMP PRIVATE  ( j)           
     do j = 1,NV
        x(j) = x(j) + v(j,i) * y(i)
     end do
     !$OMP END PARALLEL DO
  end do
end subroutine linearcombocpu

subroutine cdivbyelem(offset,x,y,z,NV)
  implicit none
  integer*4 :: NV,offset
  real*8 :: x(*),y(*),z(*)
  integer*4 :: i
  !$OMP PARALLEL  DO                                     &
  !$OMP SCHEDULE ( STATIC )                              &
  !$OMP DEFAULT  ( NONE )                                &
  !$OMP SHARED   ( x,y,z,NV,offset)  &
  !$OMP PRIVATE  ( i)           
  do i=1+offset,NV
     z(i) = y(i)/x(i)
  end do
  !$OMP END PARALLEL DO
end subroutine cdivbyelem

subroutine cmulbyelem(offset,x,y,z,NV)
  implicit none
  integer*4 :: NV,offset
  real*8 :: x(*),y(*),z(*)
  integer*4 :: i
  !$OMP PARALLEL  DO                                     &
  !$OMP SCHEDULE ( STATIC )                              &
  !$OMP DEFAULT  ( NONE )                                &
  !$OMP SHARED   ( x,y,z,NV,offset)  &
  !$OMP PRIVATE  ( i)           
  do i=1+offset,NV
     z(i) = y(i)*x(i)
  end do
  !$OMP END PARALLEL DO
end subroutine cmulbyelem


subroutine cpudbsrmv(npoi1,nbnodes,nbvar,an,ia,ja,xx,yy) 

  use def_kintyp, only             :  ip,rp
  implicit none
  integer(ip), intent(in)          :: nbnodes,nbvar,npoi1
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(in)          :: xx(nbvar,*)
  real(rp),    intent(out), target :: yy(nbvar,*)
  integer(ip)                      :: ii,jj,kk,ll,col
  real(rp)                         :: raux,raux1,raux2,raux3


  if( nbvar == 1 ) then
     !
     ! NBVAR=1
     !
     !$OMP PARALLEL  DO                                     &
     !$OMP SCHEDULE ( STATIC )                              &
     !$OMP DEFAULT  ( NONE )                                &
     !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
     !$OMP PRIVATE  ( col, ii, jj, raux )            
     !
     do ii = npoi1+1,nbnodes
        yy(1,ii) = 0.0_rp
        do jj   = ia(ii),ia(ii+1)-1
           col  = ja(jj)
           raux = xx(1,col)
           yy(1,ii) = yy(1,ii) +an(1,1,jj) * raux
        end do
     end do
     !$OMP END PARALLEL DO
  else if( nbvar == 2 ) then
     !
     ! NBVAR=2
     !
     !$OMP PARALLEL  DO                                     &
     !$OMP SCHEDULE ( STATIC )                              &
     !$OMP DEFAULT  ( NONE )                                &
     !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
     !$OMP PRIVATE  ( col, ii, jj, raux1, raux2 )    
     !
     do ii = npoi1+1,nbnodes
        yy(1,ii) = 0.0_rp
        yy(2,ii) = 0.0_rp
        do jj       = ia(ii),ia(ii+1)-1
           col      = ja(jj)
           raux1    = xx(1,col)
           raux2    = xx(2,col)
           yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
           yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
           yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
           yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
        end do
     end do
     !$OMP END PARALLEL DO
  else if( nbvar == 3 ) then
     !
     ! NBVAR=3
     !
     !$OMP PARALLEL  DO                                       &
     !$OMP SCHEDULE ( STATIC )                                &
     !$OMP DEFAULT  ( NONE )                                  &
     !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )    &
     !$OMP PRIVATE  ( col, ii, jj, raux1, raux2, raux3 )    
     !
     do ii = npoi1+1,nbnodes
        yy(1,ii) = 0.0_rp
        yy(2,ii) = 0.0_rp
        yy(3,ii) = 0.0_rp
        do jj       = ia(ii),ia(ii+1)-1
           col      = ja(jj)
           raux1    = xx(1,col)
           raux2    = xx(2,col)
           raux3    = xx(3,col)
           yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
           yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
           yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux3
           yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
           yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
           yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux3
           yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1
           yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux2
           yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux3
        end do
        
     end do
     !$OMP END PARALLEL DO
  else if( nbvar == 4 ) then
     !
     ! NBVAR=4
     !
     !$OMP PARALLEL  DO                                     &
     !$OMP SCHEDULE ( STATIC )                              &
     !$OMP DEFAULT  ( NONE )                                &
     !$OMP SHARED   ( an, ia, ja, nbnodes, npoi1, xx, yy )  &
     !$OMP PRIVATE  ( col, ii, jj, raux )            
     !
     do ii = npoi1+1,nbnodes
        yy(1,ii) = 0.0_rp
        yy(2,ii) = 0.0_rp
        yy(3,ii) = 0.0_rp
        yy(4,ii) = 0.0_rp
        do jj       = ia(ii),ia(ii+1)-1
           col      = ja(jj)
           raux     = xx(1,col)
           yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
           yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
           yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
           yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
           raux     = xx(2,col)
           yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
           yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
           yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
           yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
           raux     = xx(3,col)
           yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
           yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
           yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
           yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
           raux     = xx(4,col)
           yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
           yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
           yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
           yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
        end do
        
     end do
     !$OMP END PARALLEL DO
  else
     !
     ! NBVAR = whatever
     !
     !$OMP PARALLEL  DO                                            &
     !$OMP SCHEDULE ( STATIC )                                     &
     !$OMP DEFAULT  ( NONE )                                       &
     !$OMP SHARED   ( an, ia, ja, nbnodes, nbvar, npoi1, xx, yy )  &
     !$OMP PRIVATE  ( col, ii, jj, kk, ll, raux )           
     !
     do ii = npoi1+1,nbnodes
        do kk = 1,nbvar
           yy(kk,ii) = 0.0_rp
        end do
        do jj  = ia(ii),ia(ii+1)-1
           col = ja(jj)
           do ll = 1,nbvar
              raux = xx(ll,col)
              do kk = 1,nbvar
                 yy(kk,ii) = yy(kk,ii) + an(ll,kk,jj) * raux
              end do
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
end subroutine cpudbsrmv

subroutine cpuddot(offset,npoi1,npoi2,npoi3,V,x,y,ans)
  implicit none
  real*8 :: x(*),y(*),ans
  real*8 :: ans1,ans2
  integer*4 :: npoi1,npoi2,npoi3,V,offset
  integer*4 :: i
  ans=0

  !$OMP PARALLEL DO REDUCTION(+:ans1) PRIVATE(i) DEFAULT(SHARED)
  do i=1+offset,npoi1*V
     ans1 = ans1 + (x(i)*y(i))
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO REDUCTION(+:ans2) PRIVATE(i) DEFAULT(SHARED)
  do i=npoi2*V,npoi3*V
     ans2 = ans2 + (x(i)*y(i))
  end do  
  !$OMP END PARALLEL DO

  ans = ans1 + ans2
  
end subroutine cpuddot

subroutine cpudaxpy(offset,N,x,y,alpha)
  implicit none
  real*8 :: x(*),y(*),alpha

  integer*4 :: N,offset
  integer*4 :: i

  
  !$OMP PARALLEL  DO                                            &
  !$OMP SCHEDULE ( STATIC )                                     &
  !$OMP DEFAULT  ( NONE )                                       &
  !$OMP SHARED   ( x, y,alpha,N,offset )  &
  !$OMP PRIVATE  ( i)           
  do i=1+offset,N
     y(i) = alpha*x(i) + y(i)
  end do
  !$OMP END PARALLEL DO
end subroutine cpudaxpy

subroutine cpudscal(offset,N,alpha,x)
  implicit none
  real*8 :: x(*),alpha

  integer*4 :: N,offset
  integer*4 :: i

  
  !$OMP PARALLEL  DO                                            &
  !$OMP SCHEDULE ( STATIC )                                     &
  !$OMP DEFAULT  ( NONE )                                       &
  !$OMP SHARED   ( x,alpha,N,offset )  &
  !$OMP PRIVATE  ( i)           
  do i=1+offset,N
     x(i) = alpha*x(i) 
  end do
  !$OMP END PARALLEL DO
end subroutine cpudscal

subroutine cpudcopy(offset,NV,x,y)
  implicit none
  real*8 :: x(*),y(*)
  integer*4 :: NV,offset
  integer*4 :: i
  !$OMP PARALLEL  DO                                            &
  !$OMP SCHEDULE ( STATIC )                                     &
  !$OMP DEFAULT  ( NONE )                                       &
  !$OMP SHARED   ( x, y,NV,offset)  &
  !$OMP PRIVATE  ( i)           
  do i=1+offset,NV
     y ( i ) = x ( i )
  end do
  !$OMP END PARALLEL DO
end subroutine cpudcopy
