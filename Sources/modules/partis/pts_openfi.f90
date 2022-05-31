subroutine pts_openfi(itask)
  !------------------------------------------------------------------------
  !****f* partis/pts_openfi
  ! NAME 
  !    pts_openfi  
  ! DESCRIPTION
  !    itask = 1 -- open all files
  !    itask = 2 -- close previously open pts.res and open pts.000<timestep>.res 
  ! USES
  ! USED BY
  !    pts_turnon
  !***
  !------------------------------------------------------------------------
  use def_partis
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_iofile
  use mod_outfor, only : outfor
  use mod_opfpos, only: postpr_intto8
  implicit none
  integer(ip), intent(in) :: itask 
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit
  character(200)          :: messa

  if( INOTSLAVE .or. itask == 10 ) then

     if( kfl_rstar == 2 ) then
        statu = 'old'
        forma = 'formatted'
        posit = 'append'
     else
        statu = 'unknown'
        forma = 'formatted'
        posit = 'asis'
     end if

     select case ( itask )

     case (   1_ip )
        !
        ! Open result file
        !
        if ( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR1725',fil_resul_pts)
           call GET_ENVIRONMENT_VARIABLE('FOR1728',fil_oudep_pts)
           call GET_ENVIRONMENT_VARIABLE('FOR1729',fil_depsu_pts)
        else if ( kfl_naked == 1 ) then
           fil_resul_pts      = adjustl(trim(namda))//'.'           //exmod(modul)//'.res'
           fil_oudep_pts      = adjustl(trim(namda))//'-deposition.'//exmod(modul)//'.csv'
           fil_depsu_pts      = adjustl(trim(namda))//'-surface_depo.'//exmod(modul)//'.res'
   !!!        fil_num_partis_pts = adjustl(trim(namda))//'-particles.'//exmod(modul)//'.res' 
           !           
           !c etait mieux avant
           !
           !fil_oudep_pts = adjustl(trim(namda))//'-deposition.'//exmod(modul)//'.res' 
        end if


        if( kfl_rstar == 2 ) then 
           kfl_reawr = 1
        
           if (kfl_ptsres_binary==0) then
              call iofile(4_ip,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION','old','unformatted')
           else
              open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
           end if

           if( kfl_reawr == 1 ) then
              if (kfl_ptsres_binary==0) then
                 call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION','old','formatted','append')
              else
                 open(lun_resul_pts, file=fil_resul_pts, access='stream', position='append')
              end if
           else
              if (kfl_ptsres_binary==0) then
                 call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION')
                 call outfor(53_ip,lun_resul_pts,' ')
              else
                 open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
              end if
 
           end if
           kfl_reawr = 0
        else

           if (kfl_ptsres_binary==0) then
              call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION')
              !call outfor(53_ip,lun_resul_pts,' ')
           else
              open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
           end if

           call pts_output_header()

        end if
        !
        ! Convergence file
        !
        if( kfl_thermo_pts /= 0 ) then
           call outfor(55_ip,momod(modul) % lun_conve,' ')
        else
           call outfor(54_ip,momod(modul) % lun_conve,' ')
        endif

        !
        ! Deposition file
        !
        if( kfl_oudep_pts /= 0 ) then
           if( kfl_rstar == 2 ) then 
              kfl_reawr = 1
              call iofile(4_ip,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION','old','unformatted')
              if( kfl_reawr == 1 ) then
                 call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION','old','formatted','append')
              else
                 call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION')
                 igene = kfl_oudep_pts
                 call outfor(57_ip,lun_oudep_pts,' ')
              end if
              kfl_reawr = 0
           else
              call iofile(zero,lun_oudep_pts,fil_oudep_pts,'LAGRANGIAN PARTICLES POSITION')
              !igene = kfl_oudep_pts
              !call outfor(57_ip,lun_oudep_pts,' ')
              call pts_deposition_header()
              
           end if
        end if
        !
        ! Surface Deposition file
        !
        if( kfl_depos_surface_pts /= 0 ) then           
           call iofile(zero,lun_depsu_pts,fil_depsu_pts,'DEPOSITION SURFACE')
           call outfor(101_ip,lun_depsu_pts,' ')
        end if

     case (2_ip) 
        !
        ! this task is for splitted writing of pts.res
        ! open pts.res.000000<time_step_number> for writing
        ! ignore restart flag since there is no appending
        ! 
        fil_resul_pts = adjustl(trim(namda))//'.'//exmod(modul)//'.'//postpr_intto8(ittim)//'.res'
        !print *,'cutim: ',cutim

        close(lun_resul_pts)

        if (kfl_ptsres_binary==0) then
            call iofile(zero,lun_resul_pts,fil_resul_pts,'LAGRANGIAN PARTICLES POSITION')
        else
            open(lun_resul_pts, file=fil_resul_pts, access='stream',POSITION='REWIND',STATUS='REPLACE')
        end if

        call pts_output_header()

     end select

  end if

end subroutine pts_openfi

