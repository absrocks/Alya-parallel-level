!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_reaphy.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Allocate memory for the physical arrays
!> @details Allocate memory for the physical arrays
!> @} 
!------------------------------------------------------------------------
subroutine por_memphy(itask)
  use def_master
  use def_domain
  use def_porous
  use mod_memchk
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iwell,ntime,nposi

  select case ( itask )

  case ( 1_ip )
     !
     ! Properties
     ! 
  case ( 3_ip )
     !
     ! tkrel_por: Stores both k_rw & k_ro
     !
     call memory_alloca(mem_modul(1:2,modul),'TKREL_POR','por_memphy',tkrel_por,3_ip,mrows_por,nmate_por)    
     call memory_alloca(mem_modul(1:2,modul),'NROWS_POR','por_memphy',nrows_por,nmate_por)    

  case ( 4_ip )
     !
     ! well
     !
     call memory_alloca(mem_modul(1:2,modul),'WVALU_POR','por_memphy',wvalu_por,mheiw_por+1,mroww_por,nwell_por)    
     call memory_alloca(mem_modul(1:2,modul),'NROWW_POR','por_memphy',nroww_por,nwell_por)    
     call memory_alloca(mem_modul(1:2,modul),'RWELL_POR','por_memphy',rwell_por,nwell_por)    
     call memory_alloca(mem_modul(1:2,modul),'PBOHO_POR','por_memphy',pboho_por,nwell_por)    
     call memory_alloca(mem_modul(1:2,modul),'NHEIW_POR','por_memphy',nheiw_por,nwell_por)    
     call memory_alloca(mem_modul(1:2,modul),'KFL_WELLC_POR','por_memphy',kfl_wellc_por,nwell_por)    
     !
     ! New well
     !
     allocate( tywel_por(nwell_por) )
     do iwell = 1,nwell_por
        tywel_por(iwell) % itype = -1
        tywel_por(iwell) % ntime =  0
        tywel_por(iwell) % nposi =  0
        tywel_por(iwell) % radiu =  0.0_rp
        nullify( tywel_por(iwell) % q_table   )
        nullify( tywel_por(iwell) % pbh_table )
        nullify( tywel_por(iwell) % pbh_coord )
     end do

  case ( 5_ip )
     !
     ! Allocate memory for well iwell
     !
     iwell = igene
     ntime = tywel_por(iwell) % ntime
     nposi = tywel_por(iwell) % nposi
     call memory_alloca(mem_modul(1:2,modul),'Q_TABLE_POR'  ,'por_memphy',tywel_por(iwell) % q_table  ,ntime,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'PBH_TABLE_POR','por_memphy',tywel_por(iwell) % pbh_table,ntime,nposi+1_ip)
     call memory_alloca(mem_modul(1:2,modul),'PBH_COORD_POR','por_memphy',tywel_por(iwell) % pbh_coord,nposi)

  case ( -5_ip )
     !
     ! Deallocate memory for well iwell
     !
     do iwell = 1,nwell_por
        ntime = tywel_por(iwell) % ntime
        nposi = tywel_por(iwell) % nposi
        call memory_deallo(mem_modul(1:2,modul),'Q_TABLE_POR'  ,'por_memphy',tywel_por(iwell) % q_table  )
        call memory_deallo(mem_modul(1:2,modul),'PBH_TABLE_POR','por_memphy',tywel_por(iwell) % pbh_table)
        call memory_deallo(mem_modul(1:2,modul),'PBH_COORD_POR','por_memphy',tywel_por(iwell) % pbh_coord)
     end do
     deallocate( tywel_por )

  end select

end subroutine por_memphy
