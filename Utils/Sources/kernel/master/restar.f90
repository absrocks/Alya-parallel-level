subroutine restar(itask)
  !------------------------------------------------------------------------
  !****f* Master/restar
  ! NAME 
  !    restar
  ! DESCRIPTION
  !    Restart file
  ! USES
  ! USED BY
  !    Turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_iofile
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: dummi(10)
  real(rp)                :: dummr(10)
  integer(ip)             :: kfl_gores,npapa

  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  select case ( itask )

  case ( READ_RESTART_FILE )
     !
     ! Read continue restart file
     !
     if( kfl_rstar == 2 ) then

        if( INOTSLAVE ) then
           read(lun_rstar) ittim,cutim
           read(lun_rstar) npapa
           read(lun_rstar, end =100) prthe(2),prthe(3),prthe(4), dpthe  ! for lowmac flows
        else
           npapa = 0
        end if

100     if( IPARALL ) then

           dummi(1) = ittim
           dummi(2) = npapa

           dummr(1) = cutim
           dummr(2) = prthe(2)
           dummr(3) = prthe(3)
           dummr(4) = prthe(4)
           dummr(5) = dpthe

           call PAR_BROADCAST(2_ip,dummi)
           call PAR_BROADCAST(5_ip,dummr)

           ittim    = dummi(1) 
           npapa    = dummi(2)
           cutim    = dummr(1) 
           prthe(2) = dummr(2)  
           prthe(3) = dummr(3)  
           prthe(4) = dummr(4)  
           dpthe    = dummr(5)  

        end if

        if( cutim >= timef )     kfl_gotim = 0
        if( ittim >= mitim )     kfl_gotim = 0
        if( npapa /= npart ) &
             call runend('RESTAR: RESTART CAN ONLY BE DONE WITH THE SAME NUMBER OF SUBDOMAINS')
 
     end if

  case ( WRITE_RESTART_FILE )
     !
     ! Write restart file
     !
     if( INOTSLAVE ) then

        if( kfl_timei == 0 ) then
           write(lun_rstar) 0_ip,0.0_rp
        else
           write(lun_rstar) ittim,cutim
        end if
        write(lun_rstar) npart
        write(lun_rstar) prthe(1), prthe(3), prthe(4), dpthe
     end if

  end select
  !
  ! Close restart file
  !
  call respre(3_ip,kfl_gores)
  !
  ! Close restart file Hdf5
  !
  call Hdfpos(12_ip)
  
  
end subroutine restar
