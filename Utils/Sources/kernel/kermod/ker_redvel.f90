!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_redvel.f90
!> @author  Hadrien Calmet
!> @date    05/07/2012
!> @brief   Compute reduced velocity
!> @details Compute reduced velocity
!> @}
!-----------------------------------------------------------------------
subroutine ker_redvel(vered)

  use def_parame
  use def_master
  use def_domain
  use mod_gradie
  use mod_memory


  implicit none
  real(rp),   intent(in) :: vered(ndime,*)
  integer(ip)            :: ipoin
  real(rp),   pointer    :: gradu(:,:,:)

  if( INOTMASTER ) then


     nullify(gradu)
     call memory_alloca(mem_modul(1:2,modul),'GRADU','redvel',gradu,ndime,ndime,npoin)
     call gradie(veloc(:,:,1),gradu)

     points: do ipoin = 1,npoin

        call hadri3(gradu(1,1,ipoin),veloc(1,ipoin,1),veloc(1,ipoin,3),ndime,vered(1,ipoin))

     end do points

     call memory_deallo(mem_modul(1:2,modul),'GRADU','redvel',gradu)

  end if

end subroutine ker_redvel

subroutine hadri3(gvelo,velo,veli,ndime,vered)

  use def_kintyp, only : ip,rp,lg
  use def_kermod

  implicit none

  integer(ip),intent(in)  :: ndime
  real(rp),   intent(in)  :: gvelo(ndime,ndime)
  real(rp),   intent(in)  :: velo(ndime)
  real(rp),   intent(in)  :: veli(ndime)
  real(rp),   intent(out) :: vered(ndime)
  real(rp)                :: accel(ndime)
  integer(ip)             :: ierr,jdime,idime,i
  real(rp)                :: jaco(ndime,ndime)
  real(rp)                :: autovec(ndime,ndime)
  real(rp)                :: lambda_r(ndime)
  real(rp)                :: lambda_i(ndime)
  real(rp)                :: fv1(ndime)
  real(rp)                :: fv2(ndime),m,vr(ndime),n(ndime),velon
  logical(lg)             :: matzz

  i = 0
  do idime = 1,ndime
     lambda_r(idime) = 0.0_rp
     lambda_i(idime) = 0.0_rp
     fv1(idime)      = 0.0_rp
     fv2(idime)      = 0.0_rp
     vr(idime)       = 0.0_rp
     n(idime)        = 0.0_rp
     do jdime = 1,ndime
        jaco(idime,jdime)    = gvelo(idime,jdime)
        autovec(idime,jdime) = 0.0_rp
     end do
  end do


  matzz = .true.

  call hadri4(&
       ndime,ndime,jaco,lambda_r,lambda_i,matzz,&
       autovec,fv1,ierr)

  if(      lambda_i(1) == 0.0_rp .and. lambda_i(2) == -lambda_i(3) ) then
     i = 1
  else if( lambda_i(2) == 0.0_rp .and. lambda_i(1) == -lambda_i(3) ) then
     i = 2
  else if( lambda_i(3) == 0.0_rp .and. lambda_i(1) == -lambda_i(2) ) then
     i = 3
  end if

  if( i /= 0 ) then
     !
     ! calcul of Vr reduced velocity
     !
     ! normalized real eigenvector n
     !
     m = sqrt(&
          autovec(i,1)*autovec(i,1) + &
          autovec(i,2)*autovec(i,2) + &
          autovec(i,3)*autovec(i,3) )
     !
     !
     !
     n(1) = autovec(i,1) / m
     n(2) = autovec(i,2) / m
     n(3) = autovec(i,3) / m
     !
     ! calcul del scalar product velo*n
     !
     if( kfl_vortx_thres == 0 .or. kfl_vortx_thres == 1) then
        !
        ! Eigen method
        !
        velon = (velo(1)*n(1)+velo(2)*n(2)+velo(3)*n(3))
        !
        vered(1) = velo(1) - velon * n(1)
        vered(2) = velo(2) - velon * n(2)
        vered(3) = velo(3) - velon * n(3)

     elseif ( kfl_vortx_thres == 2) then
        !
        ! Parallel Method
        !
        do idime=1,ndime
           accel(idime) = (gvelo(1,idime)*velo(1)+gvelo(2,idime)*velo(2)+gvelo(3,idime)*velo(3))+(velo(idime)-veli(idime))
        enddo
        !
        velon = (accel(1)*n(1)+accel(2)*n(2)+accel(3)*n(3))
        !
        vered(1) = accel(1) - velon * n(1)
        vered(2) = accel(2) - velon * n(2)
        vered(3) = accel(3) - velon * n(3)
        !
     end if
  else
     vered(1) = 0.0_rp
     vered(2) = 0.0_rp
     vered(3) = 0.0_rp
  end if

  if (ierr/=0) then
     write(*,*)' todo mal '
     stop
  end if

end subroutine hadri3

subroutine hadri4(nm,n,a,wr,wi,matz,z,fv1,ierr)
  use      def_kintyp , only : ip,rp,lg
  implicit none

  integer(ip) :: n,nm,is1,is2,ierr
  logical(lg) :: matz
  real(rp)    :: a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
  integer(ip) :: iv1(n)

  if (n <= nm) go to 10
  ierr = 10 * n
  go to 50
  !
10 continue
  call balanc(nm,n,a,is1,is2,fv1)
  call ielmhes(nm,n,is1,is2,a,iv1)
  if( matz ) go to 20
  !
  ! Find eigenvalues only
  !
  call hqr(nm,n,is1,is2,a,wr,wi,ierr)
  go to 50
  !
  ! Find both eigenvalues and eigenvectors
  !
20 continue
  call eltran(nm,n,is1,is2,a,iv1,z)
  call hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
  if( ierr /= 0 ) go to 50
  call balbak(nm,n,is1,is2,fv1,n,z)
50 return

end subroutine hadri4

subroutine balanc(nm,n,a,low,igh,xscal)

  use def_kintyp, only : ip,rp,lg
  implicit none
  integer(ip)  :: i,j,k,l,m,n,jj,nm,igh,low,iexc
  real(rp)     :: a(nm,n),xscal(n)
  real(rp)     :: c,f,g,r,s,b2,radixr
  logical (lg) :: noconv

  radixr = 16.0_rp
  !
  b2 = radixr * radixr
  k = 1
  l = n
  go to 100
  !     .......... in-line procedure for row and
  !                column exchange ..........
20 xscal(m) = real(j,rp)
  if (j == m) go to 50
  !
  do i = 1, l
     f = a(i,j)
     a(i,j) = a(i,m)
     a(i,m) = f
  end do
!
  do i = k, n
     f = a(j,i)
     a(j,i) = a(m,i)
     a(m,i) = f
  end do
!
   50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 if (l == 1) go to 280
      l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
!
         do 110 i = 1, l
            if (i == j) go to 110
            if (a(j,i) .ne. 0.0_rp) go to 120
  110    continue
!
         m = l
         iexc = 1
         go to 20
  120 continue
!
      go to 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1
!
  140 do 170 j = k, l
!
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0_rp) go to 170
  150    continue
!
         m = k
         iexc = 2
         go to 20
  170 continue
!     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 xscal(i) = 1.0_rp
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!
      do 270 i = k, l
         c = 0.0_rp
         r = 0.0_rp
!
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  200    continue
!     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0_rp .or. r .eq. 0.0_rp) go to 270
         g = r / radixr
         f = 1.0_rp
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radixr
         c = c * b2
         go to 210
  220    g = r * radixr
  230    if (c .lt. g) go to 240
         f = f / radixr
         c = c / b2
         go to 230
!     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95_rp * s) go to 270
         g = 1.0_rp / f
         xscal(i) = xscal(i) * f
         noconv = .true.
!
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
!
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
!
  270 continue
!
      if (noconv) go to 190
!
  280 low = k
      igh = l
      return

end subroutine balanc


subroutine balbak(nm,n,low,igh,xscal,m,z)


  use def_kintyp, only : ip,rp
  implicit none
  integer(ip) :: i,j,k,m,n,ii,nm,igh,low
  real(rp) :: xscal(n),z(nm,m)
  real (rp):: s

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = xscal(i)
!    .......... left hand eigenvectors are back transformed
!               if the foregoing statement is replaced by
!               s=1.0d0/xscal(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
!
 110 continue
!    ......... for i=low-1 step -1 until 1,
!              igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = int(xscal(i),ip)
         if (k .eq. i) go to 140
!
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
!
  140 continue
!
  200 return

end subroutine balbak

subroutine ielmhes(nm,n,low,igh,a,inti)

  use def_kintyp, only : ip,rp
  implicit none
  integer(ip) :: i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
  real(rp) :: a(nm,n)
  real(rp) :: x,y
  integer(ip) :: inti(igh)


      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0_rp
         i = m
!
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
!
         inti(m) = i
         if (i .eq. m) go to 130
!     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
!
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
!     .......... end interchange ..........
  130    if (x .eq. 0.0_rp) go to 180
         mp1 = m + 1
!
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0_rp) go to 160
            y = y / x
            a(i,mm1) = y
!
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
!
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
!
  160    continue
!
  180 continue
!
  200 return
end subroutine ielmhes

subroutine eltran(nm,n,low,igh,a,inti,z)

  use def_kintyp, only : ip,rp
  implicit none
  integer(ip) :: i,j,n,kl,mm,mp,nm,igh,low,mp1
  real(rp) :: a(nm,igh),z(nm,n)
  integer(ip) :: inti(igh)

!
!     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
!
         do 60 i = 1, n
   60    z(i,j) = 0.0_rp
!
         z(j,j) = 1.0_rp
   80 continue
!
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
!
         i = inti(mp)
         if (i .eq. mp) go to 140
!
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0_rp
  130    continue
!
         z(i,mp) = 1.0_rp
  140 continue
!
  200 return

end subroutine eltran

subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
  !
  !  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
  !
  use def_kintyp, only : ip,rp,lg
  implicit none
  integer(ip) :: i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
  real(rp) :: h(nm,n),wr(n),wi(n)
  real(rp) :: p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
  logical(lg) :: notlas

!
  ierr = 0
  norm = 0.0_rp
  k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
  do 50 i = 1, n
!
         do 40 j = k, n
40          norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0_rp
50       continue
!
      en = igh
      t = 0.0_rp
      itn = 30*n
!     .......... search for next eigenvalues ..........
60    if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
70    do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0_rp) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
80       continue
!     .......... form shift ..........
100      x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
120      h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75_rp * s
      y = x
      w = -0.4375_rp * s * s
130   its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
140      continue
!
150      mp2 = m + 2
!
         do 160 i = mp2, en
         h(i,i-2) = 0.0_rp
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0_rp
160      continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
         do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0_rp
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0_rp) go to 260
         p = p / x
         q = q / x
         r = r / x
170      s = dsign(dsqrt(real(p*p+q*q+r*r, 8)),real(p,8))
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
180      if (l .ne. m) h(k,k-1) = -h(k,k-1)
190      p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
200         continue
!
         j = min(en,k+3)
!     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
210         continue
            go to 255
225      continue
!     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
230         continue
!
            j = min(en,k+3)
!     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
240         continue
255         continue
!
260         continue
!
            go to 70
!     .......... one root found ..........
270   wr(en) = x + t
      wi(en) = 0.0_rp
      en = na
      go to 60
!     .......... two roots found ..........
280   p = (y - x) / 2.0_rp
      q = p * p + w
      zz = sqrt(abs(q))
      x = x + t
      if (q .lt. 0.0_rp) go to 320
!     .......... real pair ..........
      zz = p + dsign(real(zz,8),real(p,8))
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0_rp) wr(en) = x - w / zz
      wi(na) = 0.0_rp
      wi(en) = 0.0_rp
      go to 330
      !     .......... complex pair ..........
320   wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
330   en = enm2
      go to 60
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
1000  ierr = en
1001  return

    end subroutine hqr

    subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)

      use def_kintyp, only : ip,rp,lg
      implicit none
      integer(ip) :: i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn, &
                     igh,itn,its,low,mp2,enm2,ierr
      real(rp) :: h(nm,n),wr(n),wi(n),z(nm,n)
      real(rp) :: p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical(lg) :: notlas

      ierr = 0
      norm = 0.0_rp
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0_rp
   50 continue
!
      en = igh
      t = 0.0_rp
      itn = 30*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0_rp) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75_rp * s
      y = x
      w = -0.4375_rp * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two conse!utive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0_rp
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0_rp
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0_rp
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0_rp) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(sqrt(real(p*p+q*q+r*r,8)),real(p,8))
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min(en,k+3)
!     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
!     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min(en,k+3)
!     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
!     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0_rp
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0_rp
      q = p * p + w
      zz = dsqrt(real(abs(q), 8))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0_rp) go to 320
!     .......... real pair ..........
      zz = p + dsign(real(zz,8),real(p,8))
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0_rp) wr(en) = x - w / zz
      wi(na) = 0.0_rp
      wi(en) = 0.0_rp
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(real(p*p+q*q,8))
      p = p / r
      q = q / r
!     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
!     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
!     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
!
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0_rp) go to 1001
!     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
!     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0_rp
         if (na .eq. 0) go to 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0_rp
!
            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)
!
            if (wi(i) .ge. 0.0_rp) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0_rp) go to 640
            t = w
            if (t .ne. 0.0_rp) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01_rp * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
!     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
!
!     .......... overflow control ..........
  680       t = abs(h(i,en))
            if (t .eq. 0.0_rp) go to 700
            tst1 = t
            tst2 = tst1 + 1.0_rp/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
!
  700    continue
!     .......... end real vector ..........
         go to 800
!     .......... complex vector ..........
  710    m = na
!     .......... last vector !omponent chosen imaginary so that
!                eigenvector matrix is triangular ..........
         if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0_rp,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0_rp
         h(en,en) = 1.0_rp
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0_rp
            sa = 0.0_rp
!
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. 0.0_rp) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0_rp) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
!     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0_rp * q
            if (vr .ne. 0.0_rp .or. vi .ne. 0.0_rp) go to 784
               tst1 = norm * (abs(w) + abs(q) + abs(x)  &
                         + abs(y) + abs(zz))
               vr = tst1
  783          vr = 0.01_rp * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi, &
                     h(i,na),h(i,en))
            if (abs(x) .le. abs(zz) + abs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,  &
                     h(i+1,na),h(i+1,en))
!
!     .......... overflow control ..........
  790       t = dmax1(abs(h(i,na)), abs(h(i,en)))
            if (t .eq. 0.0_rp) go to 795
            tst1 = t
            tst2 = tst1 + 1.0_rp/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
!
  795    continue
!     .......... end complex vector ..........
  800 continue
!     .......... end back substitution.
!                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
!
         do 820 j = i, n
  820    z(i,j) = h(i,j)
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min(j,igh)
!
         do 880 i = low, igh
            zz = 0.0_rp
!
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
!
            z(i,j) = zz
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
end subroutine hqr2

subroutine hadri5(vered)
  !
  !  linearly interpolated each components of the reduced velocity
  !  to find the center,we set vered in equationt to 0
  !  a1 + b1.r + c1.s + d1.t = 0  plan1
  !
  !  a2 + b2.r + c2.s + d2.t = 0  plan2
  !
  !  a3 + b3.r + c3.s + d3.t = 0  plan3
  !
  !
  use def_parame
  use def_master
  use def_domain
  use def_kintyp
  use mod_memory

  implicit none

  real(rp),    intent(in)  :: vered(ndime,npoin)
  integer(ip)              :: idime,ielem,counti,ncount,is,totnp,ii
  integer(ip)              :: pnode,pelty,ipoin1,ipoin2,ipoin3,ipoin4
  real(rp)                 :: a(3),b(3),c(3),d(3),r,s,t
  real(rp)                 :: vn1(3),vn2(3),vn3(3)
  real(rp)                 :: ps,mn1,mn2,para,det
  real(rp)                 :: ma(2,2),mai(2,2),ms(2),slop
  real(rp)                 :: coord_pt(3,3)
  character(8)             :: chtim
  integer(ip), target      :: dummp(1)
  integer(ip), pointer     :: icount(:),dplnp(:)
  real(rp),    pointer     :: coord_tmp(:,:)
  real(rp),    pointer     :: coord_jpt(:),npdat(:)
  real(rp)                 :: veloc_point(3),veloc_mod
  real(rp)                 :: qvorti_point


  ncount=0

  if( INOTSLAVE ) then ! Master allocates memory for indexes of displacements and number of points
     allocate (icount(npart+1),dplnp(npart+1))
     !
     ! file for vu , writing the header of .vu
     !
     if(ittim<10) then
        write(chtim,'(a,i1)') '0000000',ittim
     else if(ittim<100) then
        write(chtim,'(a,i2)') '000000',ittim
     else if(ittim<1000) then
        write(chtim,'(a,i3)') '00000',ittim
     else if(ittim<10000) then
        write(chtim,'(a,i4)') '0000',ittim
     else if(ittim<100000) then
        write(chtim,'(a,i5)') '000',ittim
     else if(ittim<1000000) then
        write(chtim,'(a,i6)') '00',ittim
     else if(ittim<10000000) then
        write(chtim,'(a,i7)') '0',ittim
     end if

     open(unit=51,file=trim(chtim)//'.res.vu',status='unknown')
     write(51,*)"CHAMP Coo( ) ="
     write(51,*)"{"
     write(51,*)"// Donnees x, y, z"

  end if

  if( INOTMASTER ) then
     !
     !initialise the qvorti array
     !
     call memory_alloca(mem_modul(1:2,modul),'COORD_TMP','hadri5',coord_tmp,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'VORTI'    ,'hadri5',vorti,ndime+1_ip,npoin)
     call vortic(1_ip)
     !
     !
     !
     elements:do ielem = 1,nelem
        pelty = ltype(ielem)
        pnode = nnode(pelty)

        if (pnode==4) then
           !
           !Shape function of tetra(shape3.f90 in kernel/mathru)
           !
           ipoin1 = lnods(1,ielem)
           if (lpoty(ipoin1)/= 0 ) then
              cycle
           end if
           ipoin2 = lnods(2,ielem)
           if (lpoty(ipoin2)/= 0 ) then
              cycle
           end if
           ipoin3 = lnods(3,ielem)
           if (lpoty(ipoin3)/= 0 ) then
              cycle
           end if
           ipoin4 = lnods(4,ielem)
           if (lpoty(ipoin4)/= 0 ) then
              cycle
           end if

           a(1) = vered(1,ipoin1)
           a(2) = vered(2,ipoin1)
           a(3) = vered(3,ipoin1)

           b(1) = - vered(1,ipoin1) + vered(1,ipoin2)
           b(2) = - vered(2,ipoin1) + vered(2,ipoin2)
           b(3) = - vered(3,ipoin1) + vered(3,ipoin2)

           c(1) = - vered(1,ipoin1) + vered(1,ipoin3)
           c(2) = - vered(2,ipoin1) + vered(2,ipoin3)
           c(3) = - vered(3,ipoin1) + vered(3,ipoin3)

           d(1) = - vered(1,ipoin1) + vered(1,ipoin4)
           d(2) = - vered(2,ipoin1) + vered(2,ipoin4)
           d(3) = - vered(3,ipoin1) + vered(3,ipoin4)
           !
           ! the line of intersection produced by the 2 plans  vn1^vn2 = vn3
           !
           vn1(1) = b(1)
           vn1(2) = c(1)
           vn1(3) = d(1)

           vn2(1) = b(2)
           vn2(2) = c(2)
           vn2(3) = d(2)
           !
           !  condition of non colinear abs(n1.n2)/||n1||*||n2||<1 if parallel =1
           !
           ps  = abs(  vn1(1) * vn2(1) + vn1(2) * vn2(2) + vn1(3) * vn2(3) )
           mn1 = sqrt( vn1(1) * vn1(1) + vn1(2) * vn1(2) + vn1(3) * vn1(3) )
           mn2 = sqrt( vn2(1) * vn2(1) + vn2(2) * vn2(2) + vn2(3) * vn2(3) )

           if( mn1 == 0.0_rp .or. mn2 == 0.0_rp ) then
              cycle
           end if

           para = ps / (mn1*mn2)

           if (para < 1) then
              call vecpro(vn1,vn2,vn3,ndime)
           else if (para ==1) then
              !write(*,*)"2 plans are parallel,next element"
              !write(*,*)"element=",ielem
              cycle
           end if
           !
           !  finding a point is belowed vn3 and the 2 plans
           !  solve a system of 2 eq 2 unknow
           !  ex r=0
           !  a1+c1.s+d1.t=0
           !  a2+c2.s+d2.t=0
           !
           ! |c1 d1| |s| = |-a1|  ===>  A . S = R ==> S=A-1 R
           ! |c2 d2| |t| = |-a2|
           !
           !  mai inverse matrix of A
           !
           ma(1,1) = c(1)
           ma(1,2) = d(1)
           ma(2,1) = c(2)
           ma(2,2) = d(2)

           call invmtx(ma,mai,det,ndime-1)

           ms(1) = mai(1,1)*(-a(1)) + mai(1,2)*(-a(2))
           ms(2) = mai(2,1)*(-a(1)) + mai(2,2)*(-a(2))
           !
           ! point P (0,s,t)
           !
           !r=0
           !s=ms(1)
           !t=ms(2)
           !
           !solve y=ax+b => y= k vn3 + P
           !
           counti=0
           !
           ! face1 t=0 => k=-p(3)/vn3(3)
           !
           !if (b(1)==c(1).and.c(1)==0) then
           ! plan parallel to rs so
           !
           !check if the line pass through the face 1
           !
           if( vn3(3) == 0.0_rp ) then
              cycle
           else
              slop = -ms(2)/vn3(3)
              !write(*,*)"slop face 1",slop
              r = 0.0_rp + slop*vn3(1)
              s =  ms(1) + slop*vn3(2)
              t = 0.0_rp
              if( r >= 0.0_rp .and. s >= 0.0_rp .and. (s+r) <= 1.0_rp ) then
                 !write(*,*)"intersection with face 1"
                 counti=counti+1
                 do idime=1,ndime
                    coord_pt(idime,counti)=&
                         (1-r-s-t)*coord(idime,ipoin1) &
                         +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                         +t*coord(idime,ipoin4)
                 end do
              end if
           end if
           !
           ! face2 s=0 => k=-p(2)/vn3(2)
           !
           if( vn3(2) == 0.0_rp ) then
              cycle
           else
              slop=-ms(1)/vn3(2)
              !write(*,*)"slop face 2",slop
              r= 0.0_rp + slop*vn3(1)
              s= 0.0_rp
              t= ms(2)  + slop*vn3(3)
              if( r >= 0.0_rp .and. t >= 0.0_rp .and. (r+t)<= 1.0_rp ) then
                 !write(*,*)"intersection with face 2"
                 counti=counti+1
                 do idime=1,ndime
                    coord_pt(idime,counti)=(1-r-s-t)*coord(idime,ipoin1) &
                         +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                         +t*coord(idime,ipoin4)
                 end do
              end if
           end if
           !
           ! face3 r=0 => k=-p(1)/vn3(1)
           !
           if( vn3(1) == 0.0_rp ) then
              cycle
           else
              slop=0.0_rp/vn3(2)
              !write(*,*)"slop face 3",slop
              r= 0.0_rp
              s= ms(1)  + slop*vn3(2)
              t= ms(2)  + slop*vn3(3)
              if( s >= 0.0_rp .and. t >= 0.0_rp .and. (s+t) <= 1.0_rp ) then
                 !write(*,*)"intersection with face 3"
                 counti=counti+1
                 do idime=1,ndime
                    coord_pt(idime,counti)=(1-r-s-t)*coord(idime,ipoin1) &
                         +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                         +t*coord(idime,ipoin4)
                 end do
                 if (counti==3) then
                    write(*,*) "the line pass through the 3 faces !"
                 end if
              end if
           end if
           !
           ! face4 r+s+t=1 => k=[1-p(1)-p(2)-p(3)]/[vn3(1)+vn(2)+vn(3)]
           !
           if( (vn3(1)+vn3(2)+vn3(3)) == 0.0_rp ) then
              cycle
           else
              slop=( 1 -s -t ) / (vn3(1)+vn3(2)+vn3(3))
              !write(*,*)"slop face 4",slop
              r= 0.0_rp + slop*vn3(1)
              s= ms(1)  + slop*vn3(2)
              t= ms(2)  + slop*vn3(3)
              if (r <= 1.0_rp.and.s <= 1.0_rp.and.  &
                   t <= 1.0_rp .and.                &
                   r >= 0.0_rp.and.s >= 0.0_rp.and. &
                   t >= 0.0_rp ) then
                 !write(*,*)"intersection with face 4"
                 counti=counti+1
                 do idime=1,ndime
                    coord_pt(idime,counti)=(1-r-s-t)*coord(idime,ipoin1) &
                         +r*coord(idime,ipoin2)+s*coord(idime,ipoin3)&
                         +t*coord(idime,ipoin4)
                 end do
                 if (counti==3) then
                    write(*,*) "the line pass through the 3 faces !"
                 else if (counti==4) then
                    write(*,*) "la concha ! algo raro pasa !"
                 end if
              end if
           end if

           !
           ! if counti>2 means that the line pass through the element
           !
           if (counti == 2) then
           !
           !  first threshold with velocity criterio ( < 0.1 )
           !
              do idime=1,ndime
                  veloc_point(idime)=(1-r-s-t)*veloc(idime,ipoin1,1) &
                         +r*veloc(idime,ipoin2,1)+s*veloc(idime,ipoin3,1)&
                         +t*veloc(idime,ipoin4,1)
              end do
              veloc_mod=sqrt(veloc_point(1)**2+veloc_point(2)**2+veloc_point(3)**2)
           !   write(*,*)'veloc_mod',veloc_mod
           !
           !  second threshold with qvorticity criterio ( >10 )
           !
                 qvorti_point=(1-r-s-t)*vorti(ndime+1,ipoin1) &
                      +r*vorti(ndime+1,ipoin2)+s*vorti(ndime+1,ipoin3)&
                      +t*vorti(ndime+1,ipoin4)
          !
          !
          !    if ((veloc_mod < 0.1_rp ).AND.(qvorti_point > 10.0_rp)) then
                 ncount=ncount+1
                 do idime=1,ndime
                    coord_tmp(idime,ncount)=(coord_pt(idime,1)+coord_pt(idime,2))/2.0_rp
                 end do
         !     end if
           end if
           !
           !
        else
           write(*,*)"error in redvel"
           stop
        end if
     end do elements
     !write(*,*)"kfl_paral, ncount",kfl_paral, ncount
     !
     !allocation the points array with the right size for each proc
     !and npdat equal but for the mpi communication
     !
     allocate( npdat(ncount*ndime) )

     do ii=1,ncount
        npdat(((ii-1)*3)+1)  =  coord_tmp(1,ii)
        npdat(((ii-1)*3)+2)  =  coord_tmp(2,ii)
        npdat(((ii-1)*3)+3)  =  coord_tmp(3,ii)
     end do
     call memory_deallo(mem_modul(1:2,modul),'VORTI'    ,'hadri5',vorti)
     call memory_deallo(mem_modul(1:2,modul),'COORD_TMP','hadri5',coord_tmp)
     !
     !Transfer number of point found for each proc
     !
     dummp(1) =  ncount
     paris    => dummp
     parig    => nul1i
     npasi    =  1
     !
     !master receive number of point from slaves
     !
  else if( IMASTER ) then

     call memgen(1_ip,npart+1,0_ip)
     dummp(1) =  0  !Master sends npts=0
     paris    => dummp
     parig    => gisca
     npasi    =  1
     npari    =  1
  end if
  !
  ! mpi_gather
  !
  call Parall(705_ip)
  !
  !Construct displacement list
  !
  if( IMASTER ) then
     icount = gisca*ndime
     dplnp(1) = 0
     dplnp(2) = 0
     do is=3,npart+1
        dplnp(is) = dplnp(is-1)+icount(is-1)
     enddo
     call memgen(3_ip,npart+1,0_ip)
     !
     !Total number of points detected by slaves
     !
     totnp = dplnp(npart+1)+icount(npart+1)
  else ! We are sequential
     totnp = ncount
  end if
  !
  ! Receive points index list from slaves
  !
  if( ISLAVE ) then

     parrs => npdat
     npasr =  ncount*ndime
     parre => nul1r
     parig => nul1i
     pari1 => nul1i

  else if( IMASTER ) then

     allocate( coord_jpt(totnp))
     parrs => nul1r
     npasr =  0          ! We know that master has nvoxl = 0
     parre => coord_jpt
     parig => icount     ! Receive count array
     pari1 => dplnp      ! Displacement list

  end if
  !
  ! MPI_GATHERV
  !
  call Parall(706_ip)

  if( IMASTER ) then
     do ii=1,totnp,3
        write(51,*)coord_jpt(ii),coord_jpt(ii+1),coord_jpt(ii+2)
     end do
     deallocate(coord_jpt)

     write(51,*)"};"
     write(51,*)
     write(51,*)"MAILLAGE MonMaillage( ) ="
     write(51,*)"{"
     write(51,*)"   ZONE Zone1( LagrPoint01, Coo%3,",totnp/3,");"
     write(51,*)"};"
     close(51)
  else if( ISLAVE ) then
     deallocate( npdat )
  end if

end subroutine hadri5

subroutine cdiv(ar,ai,br,bi,cr,ci)
  use def_kintyp, only : rp
  implicit none
  real(rp) :: ar,ai,br,bi,cr,ci
  !
  !     complex division, (cr,ci) = (ar,ai)/(br,bi)
  !
  real(rp) :: s,ars,ais,brs,bis
  s = abs(br) + abs(bi)
  ars = ar/s
  ais = ai/s
  brs = br/s
  bis = bi/s
  s = brs**2 + bis**2
  cr = (ars*brs + ais*bis)/s
  ci = (ais*brs - ars*bis)/s

end subroutine cdiv
