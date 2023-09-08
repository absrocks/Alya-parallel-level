subroutine Solmum(rhsid,unkno,amatr)
  !-----------------------------------------------------------------------
  !
  ! This routine solves linear systems using MUMPS sparse direct method.     
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_solver, only    :  smemo,cpu_solve,memdi,solve_sol
  use def_domain, only    :  npoin,r_sol,c_sol,nzsol
  use mod_memchk
  implicit none
  real(rp),       target  :: rhsid(npoin*solve_sol(1)%ndofn),unkno(npoin*solve_sol(1)%ndofn)
  real(rp),       target  :: amatr(solve_sol(1)%nzmat)
  real(rp)                :: cpu_refe1,cpu_refe2
  integer(4),     save    :: ipass=0
  integer(4)              :: nnnn2,nzso2,ii,istat
  integer(4),     pointer :: j_col(:),i_row(:)

  INCLUDE 'dmumps_struc.h'
  INCLUDE 'mpif.h'  
  TYPE (DMUMPS_STRUC) mumps_par
  INTERFACE
     SUBROUTINE DMUMPS( id )
       !DEC$ ATTRIBUTES C,REFERENCE,NOMIXED_STR_LEN_ARG, ALIAS:'_DMUMPS'   :: DMUMPS
       INCLUDE 'dmumps_struc.h'
       TYPE (DMUMPS_STRUC) :: id
     END SUBROUTINE DMUMPS
  END INTERFACE
  !
  ! Initializations
  !
  call cputim(cpu_refe1)
  if(ipass==0) then
     smemo(1) = 0
     smemo(2) = 0
     ipass    = 1
  end if
  !
  ! obtain i_row & j_col from r_sol & c_sol
  !
  nnnn2=int(npoin,4)*int(solve_sol(1)%ndofn,4)
  nzso2=int(nzsol,4)*int(solve_sol(1)%ndofn,4)*int(solve_sol(1)%ndofn,4)

  allocate(j_col(nzso2),stat=istat)
  call memchk(0_ip,istat,memdi,'J_COL','Solmum',j_col)
  allocate(i_row(nzso2),stat=istat)
  call memchk(0_ip,istat,memdi,'I_ROW','Solmum',i_row)

  call matspr(i_row,j_col,nzso2)
  !
  ! Preliminaries  
  !
  ! Define a communicator for the package.
  mumps_par%COMM = PAR_COMM_MYCODE
  !
  ! Initialize an instance of the package
  ! for L U factorization (sym = 0, with working host)
  !
  mumps_par%JOB = -1_4
  mumps_par%SYM =  0_4
  mumps_par%PAR =  1_4
  CALL DMUMPS(mumps_par)
  !
  ! Define problem on the host (processor 0)
  !
  IF ( mumps_par%MYID==0 ) THEN
     ! mumps_par%N=npoin
     mumps_par%N=nnnn2
     ! mumps_par%NZ=nzsol
     mumps_par%NZ=nzso2
     ! ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )     !now i_row is used
     ! ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )     !now j_col is used
     ! ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )     !hhhrhs
  END IF
  !            
  mumps_par%IRN => i_row
  mumps_par%JCN => j_col
  if(solve_sol(1)%nrhss/=1) call runend('Solmum: nsist.ne.1, needs thinking')
  !
  ! Obtain mumps_par%RHS from rhsid
  !
#ifdef R4
  call runend('MUMPS DOES NOT WORK IN SIMPLE PRECISION')
#else
  mumps_par%RHS      => rhsid
  mumps_par%A        => amatr
  mumps_par%JOB      =  6_4
  mumps_par%ICNTL(1) =  0_4  ! output file 
  mumps_par%ICNTL(3) =  0_4  ! output file 
#endif
  !
  !  mumps_par%ICNTL(7)=5  ! chh reod metis  
  !  chh ojo con la version debug no usa metis????
  !  mumps_par%CNTL(1)=0.0 ! 0.0 no pivoting 
  !
  CALL DMUMPS(mumps_par)
  ! 
  ! Obtain unkno from mumps_par%RHS 
  !
  do ii=1,int(mumps_par%N,ip)
     unkno(ii)= mumps_par%RHS(ii)
  end do
  !
  ! Deallocate user data
  !
  IF ( mumps_par%MYID==0 )THEN
     !        DEALLOCATE( mumps_par%IRN )
     !        DEALLOCATE( mumps_par%JCN )
     nullify( mumps_par%IRN )
     nullify( mumps_par%JCN )

     call memchk(2_ip,istat,memdi,'I_ROW','Solmum',i_row)
     deallocate(i_row,stat=istat)
     if(istat/=0) call memerr(2_ip,'I_ROW','Solmum',0_ip)
     call memchk(2_ip,istat,memdi,'J_COL','Solmum',j_col)
     deallocate(j_col,stat=istat)
     if(istat/=0) call memerr(2_ip,'J_COL','Solmum',0_ip)

     nullify( mumps_par%A )
     nullify( mumps_par%RHS )
     !        DEALLOCATE( mumps_par%A   )
     !        DEALLOCATE( mumps_par%RHS )
  END IF
  !
  !  Destroy the instance (deallocate internal data structures)
  !
  mumps_par%JOB = -2_4
  !
  ! Solver statistics
  !
  solve_sol(1)%nsolv = solve_sol(1)%nsolv + 1
  CALL DMUMPS(mumps_par)
  !
  ! Compute CPU time 
  !
  call cputim(cpu_refe2)
  cpu_solve = cpu_solve + cpu_refe2 - cpu_refe1

end subroutine Solmum

! MODIFICACIONES QUE SE INTRODUCEN EN EL DSP  para compaq v fortran
!1)crear carpeta solmum a la misma altura de soldir y agregar arch
!
!2)seleccionar todos salvo solmu0 y dmumps_distribution.c
!y
!
!a) set/fort/ext_proc  /arg_pass_conv : cby ref
!		      /string_lenth: after all arg
!
!b) set/fort/prepos_symb   /use fpp
!			  /pred_prec_symb :metis	
!
!  NO LONGER NEEDED			  /include&usepath  : 
!			  C:\Herbert\cvs\ok\Zephyr\Sources\kernel\solmum
!			  
!3)sett/link  a)input agregar :  ./Libraries/libgoto_p4_256-r0.94.lib  ./Libraries/libmetis.a
!	     b)output  reserve: 0x3b9aca00
!    (la libgoto tiene una .lib y una dll las dll hay que copiarla en las 
!    carpetas debug y release pero en el dsp no se debe indicar nada
!       * si en lugar de usar la goto quiero usar la atlas pongo 
!       *  ./Libraries/libf77blas.a ./Libraries/libatlas.a
!       * y las agrego en libraries
!
!4) en la version release compilar mumps_cv.F con global optimizations
! ver comentario en readme de mumps donde dice no usar -03 para compaq
! por lo tanto uso global optimizations que de acuerdo al help entiendo
! que corresponde a -o2
!   
!5) ademas al compilar la version release tuve problemas me aparecerieron NAN
!   para solucionarlo tuve bajar la optimizacion de solmum.f a global optimizations
!   seguir analizando 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cosas a revisar
! 1) ver si tengo que usar npoin o n_eqs
! 2) j_col is a pointer and is received as an array en matspr 
!    can it introduce any problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Slight changes made to the original mumps subroutines
! 1) dmumps_c.c is not used . It was only called from dmumps_part8.F. In the places 
!    where it was called I introduced a stop in case it ever happens.  
! .......=======
!    can it introduce any problem
