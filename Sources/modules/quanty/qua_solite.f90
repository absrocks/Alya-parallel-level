subroutine qua_solite()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_solite
  ! NAME 
  !    qua_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the Schrodinger equations.
  ! USES
  !    qua_matrix
  !    qua_elmope
  !    Soldir
  !    Solite
  ! USED BY
  !    qua_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use def_solver
  use mod_outfor, only : outfor
  implicit none
  real(rp)    :: cpu_refe1,cpu_refe2,ainte,ai1,ai2,dummr
  integer(ip) :: kk,k,ielem,kpos,nk,nderiv

  call cputim(cpu_refe1)
  !
  ! Update inner iteration counter and write headings in the solver
  ! file.
  !
  itinn(modul) = itinn(modul) + 1

  call outfor(5_ip,lun_solve_qua,' ')
  !
  ! Construct the system matrix and right-hand-side
  !
  print*,'KS: matrix'
  if( INOTMASTER ) call qua_matrix()
  !
  ! Solve the algebraic system
  !       
  print*,'KS: solver'
  call soleig(amatr,eigva,eigen,bmatr,itinn(modul))

  ! normalizo los autovectores relevantes integrando !!! ojo esta en 1D!!
  kpos=-1
  if(kfl_spher==1) then
     DO KK=1,eigen_sol(1)%neiva
        AINTE = 0.0
        kpos=kpos+1
        DO ielem=1,nelem
           AI1 = (eigen(lnods(1,ielem)+kpos*npoin) * coord(1,lnods(1,ielem)) )**2
           AI2 = (eigen(lnods(3,ielem)+kpos*npoin) * coord(1,lnods(3,ielem)) )**2
           AINTE=AINTE+(AI1+AI2)*0.5*(coord(1,lnods(3,ielem))-coord(1,lnods(1,ielem)) )

           AI1 = (eigen(lnods(3,ielem)+kpos*npoin) * coord(1,lnods(3,ielem)) )**2
           AI2 = (eigen(lnods(2,ielem)+kpos*npoin) * coord(1,lnods(2,ielem)) )**2
           AINTE=AINTE+(AI1+AI2)*0.5*(coord(1,lnods(2,ielem))-coord(1,lnods(3,ielem)) )


        ENDDO
        DO K=1,npoin
           eigen(K + kpos*npoin) =eigen(K + kpos*npoin)/SQRT(AINTE)
        ENDDO

     ENDDO

  !ELSE
     !kpos = 0
     !DO KK = 1,eigen_sol(1)%neiva
     !   NK     =  0
     !   nderiv =  0     
     !   if( INOTMASTER ) then
     !      CALL qua_3dinte_R(NK,NDERIV,EIGEN(kpos+1),AINTE)
     !      print*,'EIGEN a= ',KK,SQRT(AINTE)  
     !      AINTE = 1.0_rp / SQRT(AINTE)
     !      DO K = 1,npoin
     !         eigen( K + kpos ) = eigen(K + kpos) * AINTE
     !      ENDDO
     !      print*,'EIGEN b= ',KK,SQRT(AINTE)  
     !      kpos  =  kpos + npoin
     !   else
     !      CALL qua_3dinte_R(NK,NDERIV,dummr,AINTE)
     !   end if
     !ENDDO

  endif

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine qua_solite

