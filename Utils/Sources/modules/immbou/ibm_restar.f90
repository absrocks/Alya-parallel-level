subroutine ibm_restar(itask)
  !------------------------------------------------------------------------
  !****f* immbou/ibm_restar
  ! NAME 
  !    ibm_restar
  ! DESCRIPTION
  !    Restart file for IB
  ! USES
  ! USED BY
  !    restar
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: kfl_ptask_tmp,kfl_gores
  integer(4)              :: istat

  if( nimbo == 0 ) return
  if( ISLAVE )     return

  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  !
  ! Open file
  !
  modul         = -1
  kfl_ptask_tmp = kfl_ptask

  if( itask == WRITE_RESTART_FILE ) kfl_ptask = 0
  if( itask == READ_RESTART_FILE )  kfl_ptask = 2

  if( nimbo > 0 ) then

     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0
        call ibm_defini( 1_ip )
        call ibm_defini( 7_ip )
        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'PARIN','ibm_restar',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'PARRE','ibm_restar',parre)
           if( itask == READ_RESTART_FILE ) then
              if( npari > 0 ) read(lun_rstib) parin(1:npari)
              if( nparr > 0 ) read(lun_rstib) parre(1:nparr)
           end if
        end if
     end do

     if( itask == WRITE_RESTART_FILE ) then
        if( npari > 0 ) write(lun_rstib) parin(1:npari)
        if( nparr > 0 ) write(lun_rstib) parre(1:nparr)
     end if

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'PARIN','ibm_restar',parin)
     deallocate(parin,stat=istat)
     if( istat /= 0 ) call memerr(two,'PARIN','ibm_restar',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'PARRE','ibm_restar',parre)
     deallocate(parre,stat=istat)
     if( istat /= 0 ) call memerr(two,'PARRE','ibm_restar',0_ip)

  end if
  !
  ! Close restart file
  !
  modul     = 0
  kfl_ptask = kfl_ptask_tmp
  call respre(3_ip,kfl_gores)

end subroutine ibm_restar
