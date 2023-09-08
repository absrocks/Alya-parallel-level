module mod_intpol

  use def_kintyp, only       :  ip,rp
  implicit none

contains


  subroutine krigin(coord,uncoo,limit,lnode,shapl,incoo)
    !--------------------------------------------------------------------------
    !****f* ibm_fielib
    ! NAME
    !    ibm_fielib
    ! DESCRIPTION
    !    Find the best element formed by a given node and their neiborgs that 
    !    contains the point of intersection with the particle surface    
    !
    !    Returns the point of prejection, a list of elemental nodes 
    !    and the value of the shape functions
    ! USED BY
    !    nastin/nsi_embedd
    !--------------------------------------------------------------------------
    use def_domain,  only    :  ndime

    real(rp),    intent(in)             :: coord(ndime,*)
    real(rp),    intent(in)             :: uncoo(ndime)
    integer(ip), intent(in)             :: limit
    integer(ip), intent(in)             :: lnode(limit)
    real(rp),    intent(out)            :: shapl(limit+1)
    real(rp),    intent(in),  optional  :: incoo(ndime)


    integer(ip)                         :: idime,iinte,jinte
    integer(ip)                         :: dummi,infor,nukno,ndata
    integer(ip), pointer                :: linde(:)   => null()
    real(rp)                            :: coori(ndime),coorj(ndime)
    real(rp),    pointer                :: covma(:,:) => null()
    real(rp),    pointer                :: covve(:)   => null()

    ndata = limit
    if( present(incoo) ) then
       ndata = ndata + 1
    end if

    ! Use linear mean for kriging
    nukno = ndata+ndime+1
    ! Use constant mean for kriging
    if (ndata < ndime+1) nukno = ndata+1

    allocate(covma(nukno,nukno))
    allocate(covve(nukno))
    allocate(linde(nukno))

    do iinte = 1,nukno
       !do iinte = 1,limit+1
       covve(iinte) = 0.0_rp
       do jinte = 1,nukno
          !do jinte = 1,limit+1
          covma(iinte,jinte) = 0.0_rp
       end do
    end do
    !
    ! Assembly the covariance matrix
    !
    do iinte = 1,ndata
       do jinte = iinte,ndata
          do idime = 1,ndime
             if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
             if (iinte >  limit) coori(idime) = incoo(idime)
             if (jinte <= limit) coorj(idime) = coord(idime,lnode(jinte))
             if (jinte >  limit) coorj(idime) = incoo(idime)
          end do
          covma(iinte,jinte) = 0.0_rp
          do idime = 1,ndime
             covma(iinte,jinte) = covma(iinte,jinte) + (coori(idime) - coorj(idime))**2.0_rp           
          end do
          covma(iinte,jinte) = covma(iinte,jinte)**1.5_rp
       end do
    end do

    do iinte = 1,ndata-1
       do jinte = iinte+1,ndata
          covma(jinte,iinte) = covma(iinte,jinte)
       end do
    end do
    !
    ! Nugget effect
    !
    !do iinte = 1,limit
    !   covma(iinte,iinte) =  covma(iinte,iinte) + (1.0_rp)**2.0_rp
    !end do
    !
    ! Assembly the mean value basis
    !
    do iinte = 1,ndata
       covma(iinte,ndata+1) = 1.0_rp
       if (nukno == ndata+ndime+1) then
          do idime = 1,ndime
             if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
             if (iinte >  limit) coori(idime) = incoo(idime)
          end do
          do idime = 1,ndime
             covma(iinte,ndata+idime+1)       = coori(idime)
          end do
       end if
    end do
    do iinte = 1,ndata
       covma(ndata+1,iinte) = covma(iinte,ndata+1)
       if (nukno == ndata+ndime+1) then
          do idime = 1,ndime     
             covma(ndata+idime+1,iinte)       = covma(iinte,ndata+idime+1)
          end do
       end if
    end do
    !
    ! Assembly the right hand side
    !
    do iinte = 1,ndata
       do idime = 1,ndime
          if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
          if (iinte >  limit) coori(idime) = incoo(idime)
       end do
       covve(iinte) = 0.0_rp
       do idime = 1,ndime
          covve(iinte) =  covve(iinte) + ( uncoo(idime)  - coori(idime) )**2.0_rp
       end do
       covve(iinte) = covve(iinte)**1.5_rp
    end do

    covve(ndata+1) = 1.0_rp
    if (nukno == ndata+ndime+1) then
       do idime = 1,ndime
          covve(ndata+idime+1)       = uncoo(idime)
       end do
    end if
    !
    ! Obtain the krigging interpolation coefficinets
    !  
    call ludeco(covma,nukno,linde,dummi,infor)     
    call lusolv(covma,nukno,linde,covve)
    do iinte = 1,ndata
       shapl(iinte) = covve(iinte)
    end do

    deallocate(covma)
    deallocate(covve)
    deallocate(linde)

  end subroutine krigin
  !  ***************************************************************
  !  * Given an N x N matrix A, this routine replaces it by the LU *
  !  * decomposition of a rowwise permutation of itself. A and N   *
  !  * are input. INDX is an output vector which records the row   *
  !  * permutation effected by the partial pivoting; D is output   *
  !  * as -1 or 1, depending on whether the number of row inter-   *
  !  * changes was even or odd, respectively. This routine is used *
  !  * in combination with LUBKSB to solve linear equations or to  *
  !  * invert a matrix. Return code is 1, if matrix is singular.   *
  !  ***************************************************************
  subroutine ludeco(a,n,indx,d,code)
    use def_master, only :  zeror

    integer(ip), intent(in)    :: n
    integer(ip), intent(out)   :: indx(n),d,code
    real(rp),    intent(inout) :: a(n,n)
    integer(ip)                :: i,j,k,imax
    real(rp)                   :: amax,dum,rsum,rtiny,vv(100)

    rtiny = zeror
    !rtiny = 1.5e-10_rp
    d=1; code=0

    do i=1,n
       amax=0.0_rp
       do j=1,n
          if (abs(a(i,j)) > amax) amax=abs(a(i,j))
       end do ! j loop
       if(amax < rtiny) then
          call runend('KRIGING: THE COVARIANCE MATRIX IS SINGULAR')
       end if
       vv(i) = 1.0_rp / amax
    end do ! i loop

    do j=1,n
       do i=1,j-1
          rsum = a(i,j)
          do k=1,i-1
             rsum = rsum - a(i,k)*a(k,j) 
          end do ! k loop
          a(i,j) = rsum
       end do ! i loop
       amax = 0.0_rp
       do i=j,n
          rsum = a(i,j)
          do k=1,j-1
             rsum = rsum - a(i,k)*a(k,j) 
          end do ! k loop
          a(i,j) = rsum
          dum = vv(i)*abs(rsum)
          if(dum >= amax) then
             imax = i
             amax = dum
          end if
       end do ! i loop  

       if(j /= imax) then
          do k=1,n
             dum = a(imax,k)
             a(imax,k) = a(j,k)
             a(j,k) = dum
          end do ! k loop
          d = -d
          vv(imax) = vv(j)
       end if

       indx(j) = imax
       if(abs(a(j,j)) < rtiny) a(j,j) = rtiny

       if(j /= n) then
          dum = 1.0_rp / a(j,j)
          do i=j+1,n
             a(i,j) = a(i,j)*dum
          end do ! i loop
       end if
    end do ! j loop

    return
  end subroutine ludeco
  !  ******************************************************************
  !  * solves the set of n linear equations a . x = b.  here a is     *
  !  * input, not as the matrix a but rather as its lu decomposition, *
  !  * determined by the routine ludcmp. indx is input as the permuta-*
  !  * tion vector returned by ludcmp. b is input as the right-hand   *
  !  * side vector b, and returns with the solution vector x. a, n and*
  !  * indx are not modified by this routine and can be used for suc- *
  !  * cessive calls with different right-hand sides. this routine is *
  !  * also efficient for plain matrix inversion.                     *
  !  ******************************************************************
  subroutine lusolv(a,n,indx,bb)

    integer(ip), intent(in)    :: n
    integer(ip), intent(in)    :: indx(n)
    real(rp),    intent(in)    :: a(n,n)
    real(rp),    intent(inout) :: bb(n)      
    integer(ip)                :: ii,i,ll,j
    real(rp)                   :: rsum

    ii = 0

    do i=1,n
       ll = indx(i)
       rsum = bb(ll)
       bb(ll) = bb(i)
       if(ii.ne.0) then
          do j=ii,i-1
             rsum = rsum - a(i,j)*bb(j)
          end do ! j loop
       else if(rsum.ne.0.0_rp) then
          ii = i
       end if
       bb(i) = rsum
    end do ! i loop

    do i=n,1,-1
       rsum = bb(i)
       if(i < n) then
          do j=i+1,n
             rsum = rsum - a(i,j)*bb(j)
          end do ! j loop
       end if
       bb(i) = rsum / a(i,i)
    end do ! i loop

    return
  end subroutine lusolv


end module mod_intpol
