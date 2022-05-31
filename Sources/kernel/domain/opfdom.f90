!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    opfdom.f90
!> @author  houzeaux
!> @date    2020-02-27
!> @brief   Open domain files
!> @details This subroutine gets ALL the file names to be used by Alya in two
!>          possible ways and them open them:
!>       
!>          1. Recalling them from the environment, when Alya is launched
!>          encapsulated in a shell script, or
!>       
!>          2. Composing the names out of the problem name which is given as argument
!>          when the binary file Alya is launched "naked". 
!> @} 
!-----------------------------------------------------------------------

subroutine opfdom(itask)
  
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_iofile
  use mod_output, only : output_file_names
  use def_parall, only : kfl_repart_par
  use def_parall, only : kfl_repart_post_par
  use def_parall, only : num_repart_par
  use def_AMR,    only : kfl_amr
  use def_AMR,    only : kfl_amr_post
  use def_AMR,    only : num_amr
  use mod_std
  
  implicit none
  
  integer(ip), intent(in) :: itask
  character(150)          :: fil_pdata_dom,fil_elsta_dom
  character(150)          :: fil_elmsh_dom,fil_elres_dom
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     select case ( itask )

     case ( 1_ip )
        !
        ! Get file names
        !    
        if( kfl_naked == 0 ) then
           
           call GET_ENVIRONMENT_VARIABLE('FOR021',fil_pdata_dom) 
           call GET_ENVIRONMENT_VARIABLE('FOR025',fil_elsta_dom)
           call GET_ENVIRONMENT_VARIABLE('FOR026',fil_elmsh_dom)
           call GET_ENVIRONMENT_VARIABLE('FOR027',fil_elres_dom)
           call GET_ENVIRONMENT_VARIABLE('FOR030',fil_pos00) 
           call GET_ENVIRONMENT_VARIABLE('FOR031',fil_pos01)
           call GET_ENVIRONMENT_VARIABLE('FOR032',fil_pos02)           

        else if( kfl_naked == 1 ) then
           
           fil_pdata_dom = adjustl(trim(namda))//'.dom.dat'
           fil_elsta_dom = adjustl(trim(namda))//'.els.log'
           fil_elmsh_dom = adjustl(trim(namda))//'.els.post.msh'
           fil_elres_dom = adjustl(trim(namda))//'.els.post.res'
           fil_pos00     = adjustl(trim(namda))//'.post.alyapar'
           fil_pos01     = adjustl(trim(namda))//'.post.alyafil'
           fil_pos02     = adjustl(trim(namda))//'.post.alyalog'
           
        end if
        !
        ! Open main file: *.DOM.DAT
        !
        call iofile(zero,lun_pdata_dom,fil_pdata_dom,'DOMAIN','old')

     case ( 4_ip )
        !
        ! Define unit opening option if this is a restart run for alyafil
        !
        fil_pos00_save = fil_pos00
        fil_pos01_save = fil_pos01
        fil_pos02_save = fil_pos02

        call output_file_names(REPARTITIONING=.true.,AMR=.true.)
        if( kfl_repart_par /= 0 .and. kfl_repart_post_par == 2 ) then 
           call execute_command_line ('mkdir -p repartition'//trim(intost(num_repart_par)))
        end if
        if( kfl_amr /= 0 .and. kfl_amr_post == 2 ) then 
           call execute_command_line ('mkdir -p amr'//trim(intost(num_amr)))
        end if
        !
        ! Open *.POST.ALYAFIL
        !
        if( kfl_rstar == 2 .and. iofile_file_exists(fil_pos01) ) then
           call iofile_restart_run(statu,forma,posit)
        else
           call iofile_normal_run(statu,forma,posit)
        end if
        call iofile(zero,lun_pos01,fil_pos01,'POSTPROCESS FILE NAMES',statu,forma,posit)
        !
        ! Open *.POST.ALYALOG
        !
        if( kfl_rstar == 2 .and. iofile_file_exists(fil_pos02) ) then
           call iofile_restart_run(statu,forma,posit)
        else
           call iofile_normal_run(statu,forma,posit)
        end if
        call iofile(zero,lun_pos02,fil_pos02,'POSTPROCESS LOG',statu,forma,posit)

     case ( 3_ip )
        !
        ! Close domain data file
        !
        call iofile(two,lun_pdata_dom,'','DOMAIN')
        lispa = 0
        lisda = lun_pdata      ! Recover data file
        lisre = lun_outpu      ! Recover results file

     end select

  end if

end subroutine opfdom
