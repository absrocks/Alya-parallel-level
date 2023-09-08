subroutine eigdir(eigen,unkno,amatr,cmass,npoin,ncomp_eig)

  !-----------------------------------------------------------------------
  !
  ! This routine solves generalized eigenvalues systems using a direct method.
  !  The method is taken (stolen) of eigenpack ans used the system complete.
  ! then, is more slow than the sex of the people of eighty yeards old!
  ! It is assumed that the mesh graph has been properly renumbered.      
  !
  !-----------------------------------------------------------------------
  use      def_solver
  implicit none

  integer(ip),intent(in)  :: npoin,ncomp_eig
  real(rp), intent(in)    :: amatr(npoin,npoin),cmass(npoin,npoin)
  real(rp), intent(out)   :: unkno(npoin),eigen(npoin*npoin)
  real(rp)                :: alphar(npoin),alphai(npoin),beta(npoin),&
       autoval(npoin,npoin),aux 
  real(rp)                :: cpu_refe1,cpu_refe2,epsil
  integer(ip), save       :: ipass=0,ierr=0,kpoin,niter,icont,jcont
  integer, allocatable    :: mini(:)
  logical                 :: matzz 

  !
  ! Initializations
  !
  call cputim(cpu_refe1)
  if(ipass==0) then
     smemo(1) = 0
     smemo(2) = 0
     ipass    = 1
  end if

  ! Call solver routine
  !      

  epsil=0.0D0
  matzz=.true.
  niter = eigen_sol(1)%diter

  call QZHES(npoin,amatr,cmass,matzz,autoval)
  call QZIT(npoin,amatr,cmass,epsil,matzz,autoval,ierr,niter)

  if(ierr/=0) then
     eigen_sol(1)%error = ierr
  endif

  call QZVAL(npoin,amatr,cmass,alphar,alphai,beta,matzz,autoval)

  do kpoin=1,npoin
     unkno(kpoin)= alphar(kpoin)/beta(kpoin)
  enddo

  call QZVEC(npoin,amatr,cmass,alphar,alphai,beta,autoval)


  ! ordered in the eigen-values

  do jcont=1,npoin-1
     do kpoin=jcont+1,npoin
        if(unkno(jcont) > unkno(kpoin)) then
           aux = unkno(jcont)
           unkno(jcont)=unkno(kpoin)
           unkno(kpoin) = aux
           do icont=1,npoin
              aux = autoval(icont,jcont)  
              autoval(icont,jcont) = autoval(icont,kpoin)
              autoval(icont,kpoin) = aux
           enddo
        endif
     enddo
  enddo

  !   search by the  eigen_sol(1)%neiva minimal values 

  allocate (mini(npoin))

  eigen_sol(1)%neiva=0
  do jcont=1,ncomp_eig
     eigen_sol(1)%neiva=eigen_sol(1)%neiva+1
     mini(eigen_sol(1)%neiva)=jcont
     unkno(eigen_sol(1)%neiva) = unkno(jcont) 
  enddo

  do jcont=1,eigen_sol(1)%neiva 
     do icont=1,npoin
        eigen(npoin*(jcont-1)+icont) = autoval(icont,mini(jcont))
     enddo
  enddo

  !
  ! Solver statistics
  !
  eigen_sol(1)%nsolv = eigen_sol(1)%nsolv + 1

  !
  ! Compute CPU time 
  !
  call cputim(cpu_refe2)
  cpu_solve = cpu_solve + cpu_refe2 - cpu_refe1

end subroutine eigdir
