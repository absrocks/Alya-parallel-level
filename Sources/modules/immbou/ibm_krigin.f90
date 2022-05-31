subroutine ibm_krigin(uncoo,incoo,limit,lnode,shapl)
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
  use def_kintyp,  only    :  ip,rp
  use def_immbou
  use def_master
  use def_domain
  implicit none

  real(rp),    intent(in)  :: incoo(ndime),uncoo(ndime)
  integer(ip), intent(in)  :: limit
  integer(ip), intent(in)  :: lnode(limit)
  real(rp),    intent(out) :: shapl(limit+1)

  integer(ip)              :: idime,iinte,jinte
  integer(ip)              :: dummi,infor,nukno
  integer(ip), pointer     :: linde(:)
  real(rp)                 :: coori(ndime),coorj(ndime)
  real(rp),    pointer     :: covma(:,:),covve(:)


  ! Use linear mean for kriging
  nukno = limit+ndime+2
  ! Use constant mean for kriging
  if (limit+1 < ndime+1) nukno = limit+2

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
  do iinte = 1,limit+1
     do jinte = iinte,limit+1
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

  do iinte = 1,limit
     do jinte = iinte+1,limit+1
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
  do iinte = 1,limit+1
     covma(iinte,limit+2) = 1.0_rp
     if (nukno == limit+ndime+2) then
        do idime = 1,ndime
           if (iinte <= limit) coori(idime) = coord(idime,lnode(iinte)) 
           if (iinte >  limit) coori(idime) = incoo(idime)
        end do
        do idime = 1,ndime
           covma(iinte,limit+idime+2)       = coori(idime)
        end do
     end if
  end do
  do iinte = 1,limit+1
     covma(limit+2,iinte) = covma(iinte,limit+2)
     if (nukno == limit+ndime+2) then
        do idime = 1,ndime     
           covma(limit+idime+2,iinte)       = covma(iinte,limit+idime+2)
        end do
     end if
  end do
  !
  ! Assembly the right hand side
  !
  do iinte = 1,limit+1
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

  covve(limit+2) = 1.0_rp
  if (nukno == limit+ndime+2) then
     do idime = 1,ndime
        covve(limit+idime+2)       = uncoo(idime)
     end do
  end if
  !
  ! Obtain the krigging interpolation coefficinets
  !  
  call ibm_ludeco(covma,nukno,linde,dummi,infor)     
  call ibm_solver(covma,nukno,linde,covve)
  do iinte = 1,limit+1
     shapl(iinte) = covve(iinte)
  end do

  deallocate(covma)
  deallocate(covve)
  deallocate(linde)


end subroutine ibm_krigin

