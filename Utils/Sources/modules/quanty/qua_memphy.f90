subroutine qua_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_memphy
  ! NAME
  !    qua_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES

  ! USED BY
  !    qua_reabcs
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: jj,ii
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     ! alloco estructuras de especies          
     allocate(especie_qua(atomo_qua(1)%nespec),stat=istat)          
     !read(111,'(i3)') nat 

     ! alloco atomo
     allocate(atomo_qua(1)%espe(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'NATOM','qua_reaphy',atomo_qua(1)%espe)
     allocate(atomo_qua(1)%tipoPP(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'tipPP','qua_reaphy',atomo_qua(1)%tipoPP)
     allocate(atomo_qua(1)%coorx(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'coorx','qua_reaphy',atomo_qua(1)%coorx)
     allocate(atomo_qua(1)%coory(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'coory','qua_reaphy',atomo_qua(1)%coory)
     allocate(atomo_qua(1)%coorz(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'coorz','qua_reaphy',atomo_qua(1)%coorz)
     allocate(atomo_qua(1)%spin(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'spin ','qua_reaphy',atomo_qua(1)%spin)
     allocate(atomo_qua(1)%Tiespe(atomo_qua(1)%natoms),stat=istat) 
     call memchk(zero,istat,mem_modul(1:2,modul),'Tiesp','qua_reaphy',atomo_qua(1)%Tiespe)

  case(2_ip)

     do ii = 1,atomo_qua(1)%nespec
        allocate(especie_qua(ii)%nocupa(100)) 
        do jj = 1,100
           especie_qua(ii)%nocupa(jj) = 0.0_rp
        end do
     end do

  end select

end subroutine qua_memphy
