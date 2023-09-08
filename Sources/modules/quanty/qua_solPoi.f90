subroutine qua_solPoi()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_solPoi
  ! NAME 
  !    qua_solPoi
  ! DESCRIPTION
  !    This routine solves the Poisson ecuation over the real domain 
  !    in order to obtain the nuclear portion of the Schrodinger equations.
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
  use mod_postpr
  use mod_outfor, only : outfor
  implicit none
  real(rp)    :: cpu_refe1,cpu_refe2,vhmax,vhmin
  integer(ip) :: inode
  character(5) :: wopos(2)

  call cputim(cpu_refe1)

  !  call outfor(5_ip,lun_solve_qua,' ')
  !
  ! Construct the system matrix and right-hand-side
  !
  call inisol()
  !
  ! calclula boundaries actualizadas
  !
  print*,'Poisson: boundary actualization=',kfl_paral
  call qua_poibou() 

  if( INOTMASTER ) then
     !
     ! Element assembly
     !
     print*,'Poisson: element assembly=',kfl_paral
     call qua_elmope(1_ip) 

  endif
  !
  ! Solve the algebraic system
  !
  print*,'Poisson: solver=',kfl_paral
  call solver(rhsid,unkno,amatr,pmatr)


  if( INOTMASTER ) then
     !wopos(1)='POISS'
     !wopos(2)='SCALA'
     !call postpr(unkno,wopos,ittim,cutim)
     ! ditribuyo la solucion en V_hartree
     !  FILE="data\\"//'PSPOT'//NUM(1:1)//'.3D'    
     !open(unit=111,file='data\\V_hartree.3D',STATUS='UNKNOWN')
     !write(111,*) ' x     y     z     V_H '
     if( kfl_spher == 0 ) then
        vhmax = unkno(1)
        vhmin = unkno(1)
        do inode = 1,npoin
           v_hartree(inode) = unkno(inode)
           if(v_hartree(inode) > vhmax) vhmax = v_hartree(inode)
           if(v_hartree(inode) < vhmin) vhmin = v_hartree(inode)
           !write(111,'(4e15.5)') coord(1,inode),coord(2,inode),coord(3,inode),unkno(inode)  
        end do
        !write(lun_witne_qua,*)  'v harttre:  ',vhmax,vhmin
        !else
        !   do inode=1,npoin
        !      v_hartree(inode) = unkno(inode)
        !      write(111,'(2e15.5)') coord(1,inode),unkno(inode)  
        !   enddo
     endif
     !close(111)
     !
  end if

  call cputim(cpu_refe2)
  cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1 

end subroutine qua_solPoi

