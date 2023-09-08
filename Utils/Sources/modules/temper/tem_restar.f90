subroutine tem_restar(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_restar
  ! NAME 
  !    tem_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod, only : kfl_adj_prob
  use def_domain
  use def_temper
  use mod_postpr
  use mod_memchk
  use mod_ADR, only : AR_OSS
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores, ielem, pelty, pgaus, igaus, istat, ncomp
  integer(ip)             :: ipoin, nnqua
  integer(ip)             :: dummi(10), nreal
  real(rp)                :: dummr(10)
 
  !
  ! Check if restart file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = min(3_ip,ncomp_tem)  
  else
     icomp = 1 
  end if

  !----------------------------------------------------------------------
  !
  ! Temperature
  !
  !----------------------------------------------------------------------


  gesca => therm(:,icomp)
  call postpr(gesca,postp(1)%wopos(1:3,1),ittim,cutim)
  if (kfl_tisch_tem==2) then  ! BDF
     if (ncomp_tem > 3) then  ! BDF2 and bigger
        gesca => therm(:,4)
        call postpr(gesca,postp(1)%wopos(1:3,16),ittim,cutim)

        if (itask == READ_RESTART_FILE) then
           !
           ! Broadcast of file_opened
           !
           if( IPARALL ) then 
              if (IMASTER) then
                 if (file_opened) then 
                    dummi(1) = 1
                 else
                    dummi(1) = 0              
                 end if
              end if

              nreal = 1
              call parari('BCT',0_ip,nreal,dummi)
              if (dummi(1)==1) then
                 file_opened = .true.
              else if (dummi(1) == 0) then
                 file_opened = .false.
              end if

           end if
           kfl_rstar_two = file_opened ! last run was BDF2
        end if
     
        if (itask == READ_RESTART_FILE .and. .not. kfl_rstar_two ) then  
           neule_tem = ittim +1
           if (INOTMASTER) then
              do ipoin =1, npoin
                 therm(ipoin, 4) = therm(ipoin,icomp)
              end do
           end if
        end if
     end if
  end if  
  !----------------------------------------------------------------------
  !
  ! Projections for OSS stabilization
  !
  !----------------------------------------------------------------------

  if( kfl_ortho_tem > 0 )       call postpr(ADR_tem % proje1,postp(1)%wopos(1:3,10),ittim,cutim)
  if( kfl_ortho_tem == AR_OSS ) call postpr(ADR_tem % proje2,postp(1)%wopos(1:3,30),ittim,cutim)

  !----------------------------------------------------------------------
  !
  ! SGS for tracking
  !
  !----------------------------------------------------------------------

  if( kfl_sgsno_tem /= 0 .or. kfl_sgsti_tem /= 0 ) then
     !
     ! Temperature SGS
     ! 
     call postpr(ADR_tem % sgs,postp(1)%wopos(1:3,3),ittim,cutim,1_ip) !ADR_tem % sgs 1

     if( itask == READ_RESTART_FILE.and. INOTMASTER ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           do igaus = 1,pgaus
              ADR_tem % sgs(ielem) % a(1,igaus,2) = ADR_tem % sgs(ielem) % a(1,igaus,1)              
           end do
        end do
     end if

     !
     ! BDF 
     !
     if (kfl_tisch_tem==2) then        
        call postpr(ADR_tem % sgs,postp(1)%wopos(1:3,17),ittim,cutim,3_ip)  !ADR_tem % sgs 2   

        if(itask == READ_RESTART_FILE.and. IPARALL ) then 
           !
           ! Broadcast of file_opened
           !
           if (IMASTER) then
              if (file_opened) then 
                 dummi(1) = 1
              else
                 dummi(1) = 0              
              end if
           end if

           nreal = 1
           call parari('BCT',0_ip,nreal,dummi)
           if (dummi(1)==1) then
              file_opened = .true.
           else if (dummi(1) == 0) then
              file_opened = .false.
           end if

        end if
        if(.not.file_opened.and.itask == READ_RESTART_FILE.and.INOTMASTER) then
           do ielem = 1,nelem
              pelty = abs(ltype(ielem))
              pgaus = ngaus(pelty)
              do igaus = 1,pgaus
                 ADR_tem % sgs(ielem)%a(1,igaus,3) = ADR_tem % sgs(ielem)%a(1,igaus,2)              
              end do
           end do
        end if
     end if
     
  end if

  !----------------------------------------------------------------------
  !
  ! Restart file tem.rst
  !
  !----------------------------------------------------------------------
  if( itask == READ_RESTART_FILE ) then  
     if( INOTSLAVE ) then      
        read ( momod(modul) % lun_rstar) vinvt_tem(2)
        if (kfl_tisch_tem==2.and. kfl_rstar_two)   read ( momod(modul) % lun_rstar) vinvt_tem(3)
     end if
     if( IPARALL ) then 
        !
        ! Broadcast
        !
        dummr(1) = vinvt_tem(2)
        nreal = 1
        if (kfl_tisch_tem==2.and. kfl_rstar_two) then
           dummr(2) = vinvt_tem(3)
           nreal = nreal +1
        end if

        call pararr('BCT',0_ip,nreal,dummr)

        vinvt_tem(2) = dummr(1)
        if (kfl_tisch_tem==2.and. kfl_rstar_two) &
             vinvt_tem(3) = dummr(2) 

     end if
    

  else    
     !
     !  Writes restart
     !        
     if( INOTSLAVE ) then 
        write( momod(modul) % lun_rstar) vinvt_tem(2)
        if (kfl_tisch_tem==2)  write( momod(modul) % lun_rstar) vinvt_tem(3)
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! assign constant tempe forward values for adjoint
  !
  !----------------------------------------------------------------------  
   
  if( itask == READ_RESTART_FILE .and. kfl_adj_prob == 1 ) then
    do ipoin = 1, npoin
	tempe_forw(ipoin,1) = tempe(ipoin,icomp)
	tempe(ipoin,icomp) = 0.0_rp
    end do    
  endif
  
  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------  
  call respre(3_ip,kfl_gores)
  
  !----------------------------------------------------------------------
  !
  ! Updutes some scalars (needed in master)
  !
  !----------------------------------------------------------------------  
   if( itask == READ_RESTART_FILE ) then  
      prthe(1) = prthe(2)
      if ( .not. kfl_rstar_two ) prthe(3) = prthe(2)   
   end if

end subroutine tem_restar
 
